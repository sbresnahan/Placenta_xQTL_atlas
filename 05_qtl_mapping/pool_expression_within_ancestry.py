#!/usr/bin/env python3
"""
pool_expression_within_ancestry.py — Pool pre-normalization expression BEDs
across cohorts within each ancestry stratum.

Reads the unnorm gene-level expression BED files produced by the PANTRY
pipeline (step 10_aggregate_expression.sh writes output/unnorm/expression.bed)
and pools them within ancestry strata, taking the gene union across cohorts
by default (NaN-filled where a cohort lacks a feature; the pooled detection
filters at the normalization stage decide the final feature set).
`--pool-mode intersection` restores the legacy behavior. The output is the
input for QN + within-cohort INT + HCP estimation.

Usage:
  python3 pool_expression_within_ancestry.py \
      --ancestry-map pooled_sample_ancestry_RNAseq.tsv \
      --config config.yml \
      --output-dir pooled/ \
      [--modality expression]

The ancestry map is a 3-column TSV: sample_id, assigned_ancestry, cohort.

Output: <ancestry>_pooled_<modality>.bed per ancestry stratum (e.g.
expression -> <ancestry>_pooled_expression.bed; expression_counts ->
<ancestry>_pooled_expression_counts.bed, the TMM->VST input written per
cohort by tximport_counts.R).
"""

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd


BED_META_COLS = ["#chr", "start", "end", "phenotype_id"]


def load_ancestry_map(path):
    """Load the ancestry map TSV.

    Returns a DataFrame with columns: sample_id, assigned_ancestry, cohort.
    """
    df = pd.read_csv(path, sep="\t")
    expected = {"sample_id", "assigned_ancestry", "cohort"}
    missing = expected - set(df.columns)
    if missing:
        sys.stderr.write(
            f"ERROR: ancestry map missing columns: {missing}\n"
            f"  Expected: {expected}\n  Got: {list(df.columns)}\n")
        sys.exit(1)
    return df


def load_cohort_bed(cohort_dir, modality):
    """Load an unnorm expression BED file for a cohort.

    Returns a DataFrame with the 4 metadata columns + sample columns.
    """
    bed_path = os.path.join(cohort_dir, "output", "unnorm", f"{modality}.bed")
    if not os.path.exists(bed_path):
        sys.stderr.write(f"  WARN: {bed_path} not found\n")
        return None
    df = pd.read_csv(bed_path, sep="\t",
                     dtype={"#chr": str, "start": int, "end": int, "phenotype_id": str})
    return df


def load_config_cohorts(config_path):
    """Extract cohort names and output_base from config.yml via config_get.py."""
    scripts_dir = os.path.dirname(os.path.abspath(__file__))
    config_get = os.path.join(scripts_dir, "config_get.py")
    if not os.path.exists(config_get):
        # repo layout: config_get.py lives in ../03_phenotyping
        config_get = os.path.join(scripts_dir, "..", "03_phenotyping", "config_get.py")
    if not os.path.exists(config_get):
        config_get = "config_get.py"  # last resort: cwd/PATH (flat deploy)

    result = subprocess.run(
        ["python3", config_get, config_path],
        capture_output=True, text=True, check=True
    )
    config_vars = {}
    for line in result.stdout.splitlines():
        if line.startswith("export "):
            m = re.match(r"export (\w+)=(.*)", line)
            if m:
                val = m.group(2)
                if val.startswith("'") and val.endswith("'"):
                    val = val[1:-1].replace("'\"'\"'", "'")
                config_vars[m.group(1)] = val

    output_base = config_vars.get("OUTPUT_BASE", "")

    # Extract cohort names by parsing config.yml directly for the cohorts block
    with open(config_path) as f:
        config_text = f.read()

    # Simple parse: find lines under "cohorts:" that are 2-space indented keys
    cohorts = []
    in_cohorts = False
    for line in config_text.splitlines():
        stripped = line.split("#")[0].rstrip()
        if not stripped.strip():
            continue
        if stripped.startswith("cohorts:"):
            in_cohorts = True
            continue
        if in_cohorts:
            # Cohort keys are indented exactly 2 spaces and end with ':'
            if re.match(r"^  \w", stripped) and stripped.strip().endswith(":"):
                cohort_name = stripped.strip().rstrip(":")
                cohorts.append(cohort_name)
            elif not stripped.startswith(" "):
                # Left-aligned line means we've left the cohorts block
                in_cohorts = False
                break

    return output_base, cohorts


def pool_ancestry_stratum(ancestry, ancestry_samples, cohort_dirs, modality,
                          output_dir, pool_mode="union"):
    """Pool expression BEDs for one ancestry stratum.

    Args:
        ancestry: ancestry label (e.g. "EUR")
        ancestry_samples: DataFrame of samples in this stratum (sample_id, cohort)
        cohort_dirs: dict cohort -> cohort_dir
        modality: which unnorm BED to pool
        output_dir: where to write the pooled BED
        pool_mode: "union" (default; devBrain-style — pooled detection
            filters at normalization decide the final feature set) or
            "intersection" (legacy behavior)

    Returns:
        dict with pooling statistics
    """
    # Group samples by cohort within this ancestry
    samples_by_cohort = defaultdict(list)
    for _, row in ancestry_samples.iterrows():
        samples_by_cohort[row["cohort"]].append(row["sample_id"])

    print(f"\n{'='*60}")
    print(f"Ancestry stratum: {ancestry} (N={len(ancestry_samples)})")
    print(f"  Cohorts: {dict(samples_by_cohort)}")

    # Load each cohort's BED and extract relevant samples
    cohort_dfs = {}
    gene_sets = {}
    samples_found = defaultdict(list)
    samples_missing = defaultdict(list)

    for cohort, sample_list in samples_by_cohort.items():
        cohort_dir = cohort_dirs.get(cohort)
        if cohort_dir is None:
            sys.stderr.write(f"  WARN: no directory for cohort {cohort}, skipping\n")
            continue

        bed = load_cohort_bed(cohort_dir, modality)
        if bed is None:
            sys.stderr.write(f"  WARN: no {modality} BED for cohort {cohort}, skipping\n")
            continue

        # Check which samples are present in the BED
        bed_samples = [c for c in bed.columns if c not in BED_META_COLS]
        for s in sample_list:
            if s in bed_samples:
                samples_found[cohort].append(s)
            else:
                samples_missing[cohort].append(s)

        # Extract metadata + found samples
        cols = BED_META_COLS + samples_found[cohort]
        cohort_dfs[cohort] = bed[cols].copy()
        gene_sets[cohort] = set(bed["phenotype_id"])

    if not cohort_dfs:
        sys.stderr.write(f"  ERROR: no cohort BEDs loaded for ancestry {ancestry}\n")
        return None

    # Report missing samples
    for cohort, missing in samples_missing.items():
        if missing:
            print(f"  Missing samples in {cohort}: {len(missing)}")
            for s in missing:
                print(f"    {s}")

    # Feature set: union (default) or legacy intersection across cohorts.
    # Union is the devBrain-style schema: per-cohort prefilters no longer
    # decide the pooled feature set; the pooled detection filters at the
    # normalization stage (TPM>0.1 in >=25% / >=40% non-NA + no-variance +
    # per-cohort >=2 non-NA guards) decide. Features absent from a cohort
    # get NaN for that cohort's samples (ratio modalities only — expression
    # BEDs share the GTF feature space, so union == intersection there once
    # per-cohort prefilters are skipped).
    gene_intersection = set.intersection(*gene_sets.values())
    gene_union = set.union(*gene_sets.values())
    n_genes_per_cohort = {c: len(g) for c, g in gene_sets.items()}
    print(f"  Genes per cohort: {n_genes_per_cohort}")
    print(f"  Gene intersection: {len(gene_intersection)}")
    print(f"  Gene union: {len(gene_union)} "
          f"(+{len(gene_union) - len(gene_intersection)} vs intersection)")

    feature_set = gene_union if pool_mode == "union" else gene_intersection
    union_sorted = sorted(feature_set)

    # Merge: reindex every cohort to the pooled feature set; metadata
    # (#chr/start/end) come from the first cohort in which the feature
    # appears (combine_first).
    all_samples = []
    merged_meta = None
    merged_data = []
    n_na_cells = 0

    for cohort, df in cohort_dfs.items():
        df_idx = df.set_index("phenotype_id")
        df_re = df_idx.reindex(union_sorted)

        meta_c = df_re[["#chr", "start", "end"]]
        if merged_meta is None:
            merged_meta = meta_c.copy()
        else:
            merged_meta = merged_meta.combine_first(meta_c)

        sample_cols = [c for c in df.columns if c not in BED_META_COLS]
        n_na_cells += int(df_re[sample_cols].isna().sum().sum())
        merged_data.append(df_re[sample_cols])
        all_samples.extend(sample_cols)

    # All frames are indexed by phenotype_id here; concat aligns on it.
    pooled = pd.concat([merged_meta] + merged_data, axis=1)
    pooled = pooled.reset_index()  # phenotype_id back as a column
    pooled = pooled[BED_META_COLS + all_samples]
    if n_na_cells:
        print(f"  NaN cells introduced by union pooling: {n_na_cells} "
              f"(features absent from a cohort; imputed/filtered at "
              f"normalization)")

    # Verify column count
    assert len(pooled.columns) == 4 + len(all_samples), \
        f"Column count mismatch: {len(pooled.columns)} vs {4 + len(all_samples)}"

    # Union reindexing pads absent features with NaN, leaving start/end as
    # float64; float_format="%g" would then write large coordinates in
    # scientific notation (e.g. 1.0064e+08), which breaks downstream readers
    # expecting integer BED coordinates. combine_first guarantees every
    # pooled feature has metadata from >=1 cohort, so cast back to int.
    for coord in ("start", "end"):
        if pooled[coord].isna().any():
            bad = pooled.loc[pooled[coord].isna(), "phenotype_id"].head(5).tolist()
            raise ValueError(
                f"Missing {coord} metadata for {int(pooled[coord].isna().sum())} "
                f"features (e.g. {bad}); cannot write integer BED coordinates")
        pooled[coord] = pooled[coord].astype(int)

    # Write output
    out_path = os.path.join(output_dir, f"{ancestry}_pooled_{modality}.bed")
    pooled.to_csv(out_path, sep="\t", index=False, float_format="%g")
    print(f"  Output: {out_path}")
    print(f"  Shape: {pooled.shape[0]} genes x {len(all_samples)} samples")

    return {
        "ancestry": ancestry,
        "n_samples": len(all_samples),
        "n_genes": len(feature_set),
        "n_genes_intersection": len(gene_intersection),
        "n_genes_union": len(gene_union),
        "n_gained_vs_intersection": len(gene_union) - len(gene_intersection),
        "n_na_cells_union": n_na_cells,
        "pool_mode": pool_mode,
        "n_cohorts": len(cohort_dfs),
        "genes_per_cohort": n_genes_per_cohort,
        "samples_per_cohort": {c: len(s) for c, s in samples_found.items()},
        "missing_samples": {c: len(s) for c, s in samples_missing.items() if s},
        "output": out_path,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Pool unnorm expression BEDs across cohorts within ancestry strata")
    parser.add_argument("--ancestry-map", required=True,
                        help="TSV: sample_id, assigned_ancestry, cohort")
    parser.add_argument("--config", help="config.yml (for cohort paths)")
    parser.add_argument("--output-base", help="Override output_base from config")
    parser.add_argument("--cohort-dirs", nargs="*",
                        help="Explicit cohort dirs (overrides config). "
                             "Format: cohort_name=/path/to/cohort_dir")
    parser.add_argument("--modality", default="expression",
                        help="Which unnorm BED to pool (default: expression)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for pooled matrices")
    parser.add_argument("--ancestries", nargs="*",
                        help="Subset of ancestries to process (default: all)")
    parser.add_argument("--pool-mode", choices=["union", "intersection"],
                        default="union",
                        help="Feature pooling across cohorts: union (default, "
                             "devBrain-style; pooled detection filters decide) "
                             "or intersection (legacy)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load ancestry map
    ancestry_map = load_ancestry_map(args.ancestry_map)
    print(f"Ancestry map: {len(ancestry_map)} samples")
    print(f"  Ancestries: {sorted(ancestry_map['assigned_ancestry'].unique())}")
    print(f"  Cohorts: {sorted(ancestry_map['cohort'].unique())}")

    # Resolve cohort directories
    cohort_dirs = {}
    if args.cohort_dirs:
        for item in args.cohort_dirs:
            if "=" in item:
                name, path = item.split("=", 1)
                cohort_dirs[name] = path
            else:
                sys.stderr.write(f"  WARN: ignoring malformed --cohort-dirs entry: {item}\n")
    elif args.config:
        output_base, cohort_names = load_config_cohorts(args.config)
        if args.output_base:
            output_base = args.output_base
        for c in cohort_names:
            cohort_dirs[c] = os.path.join(output_base, c)
        print(f"  Config cohorts: {cohort_names}")
        print(f"  Output base: {output_base}")
    else:
        sys.stderr.write("ERROR: need --config or --cohort-dirs to locate cohort BEDs\n")
        sys.exit(1)

    # Filter to requested ancestries
    ancestries = args.ancestries or sorted(ancestry_map["assigned_ancestry"].unique())

    # Pool each ancestry stratum
    all_stats = []
    for ancestry in ancestries:
        stratum_samples = ancestry_map[ancestry_map["assigned_ancestry"] == ancestry]
        if len(stratum_samples) == 0:
            print(f"\nNo samples for ancestry {ancestry}, skipping")
            continue
        stats = pool_ancestry_stratum(
            ancestry, stratum_samples, cohort_dirs, args.modality, args.output_dir,
            pool_mode=args.pool_mode)
        if stats:
            all_stats.append(stats)

    # Summary report
    print(f"\n{'='*60}")
    print("POOLING SUMMARY")
    print(f"{'='*60}")
    for s in all_stats:
        print(f"  {s['ancestry']}: {s['n_samples']} samples, "
              f"{s['n_genes']} genes, {s['n_cohorts']} cohorts")
        print(f"    Samples/cohort: {s['samples_per_cohort']}")
        if s["missing_samples"]:
            print(f"    Missing samples: {s['missing_samples']}")

    # Write summary TSV
    summary_path = os.path.join(args.output_dir, "pooling_summary.tsv")
    summary_rows = []
    for s in all_stats:
        for cohort, n in s["samples_per_cohort"].items():
            summary_rows.append({
                "ancestry": s["ancestry"],
                "cohort": cohort,
                "n_samples": n,
                "n_genes_pooled": s["n_genes"],
                "n_genes_intersection": s["n_genes_intersection"],
                "n_genes_union": s["n_genes_union"],
                "n_gained_vs_intersection": s["n_gained_vs_intersection"],
                "n_na_cells_union": s["n_na_cells_union"],
                "pool_mode": s["pool_mode"],
                "n_missing": s["missing_samples"].get(cohort, 0),
            })
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    print(f"\nSummary: {summary_path}")


if __name__ == "__main__":
    main()
