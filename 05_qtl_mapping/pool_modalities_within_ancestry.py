#!/usr/bin/env python3
"""
pool_modalities_within_ancestry.py — Pool pre-normalization RNA phenotype BEDs
across cohorts within each ancestry stratum, for the non-expression modalities.

This is the unified pooler for roadmap step 4 (cross-cohort ComBat + pooling).
Gene-level expression is handled separately by the HCP pipeline
(pool_expression_within_ancestry.py + combat_normalize_hcp.R); this script
covers the remaining modalities:

    isoforms, alt_TSS, alt_polyA, splicing, intron_retention,
    RNA_editing, stability

Two input conventions are supported:

  1. Per-cohort unnorm BEDs (default):
     <cohort_dir>/output/unnorm/<modality>.bed
     Sample columns use ORIGINAL sample IDs. Features are pooled by
     phenotype_id INTERSECTION across cohorts (all cohorts share the HPLRv2
     reference, so IDs are directly comparable). Used for:
       isoforms, alt_TSS, alt_polyA, RNA_editing, stability

  2. Pre-pooled BEDs (splicing, intron_retention):
     The harmonize_within_ancestry.py step already produced a single pooled
     unnorm BED per ancestry with NAMESPACED sample IDs ({cohort}_{sample}),
     because splicing cluster numbers and IR event IDs are cohort-specific and
     require coordinate harmonization before pooling. For these modalities,
     pass --pre-pooled-dir pointing at the harmonize output root; this script
     locates <pre_pooled_dir>/<ancestry>/<modality>/unnorm/<modality>.bed and
     passes it through (no re-pooling), parsing cohort labels from the
     {cohort}_ prefix.

Cohort labels (needed downstream for ComBat batch) are resolved by:
  - namespaced IDs  -> split on first underscore
  - original IDs    -> join to the ancestry map on sample_id

Usage:
  # Per-cohort pooling (alt_TSS, alt_polyA, RNA_editing, stability, isoforms)
  python3 pool_modalities_within_ancestry.py \
      --ancestry-map pooled_sample_ancestry_RNAseq.tsv \
      --config config.yml \
      --modality alt_TSS \
      --output-dir pooled/

  # Pre-pooled pass-through (splicing, intron_retention)
  python3 pool_modalities_within_ancestry.py \
      --ancestry-map pooled_sample_ancestry_RNAseq.tsv \
      --config config.yml \
      --modality splicing \
      --pre-pooled-dir /path/to/harmonize_output \
      --output-dir pooled/

The ancestry map is a 3-column TSV: sample_id, assigned_ancestry, cohort.

Output: <ancestry>_<modality>_pooled.bed per ancestry stratum, plus a
pooling_summary.tsv.

Dependency: config_get.py (shipped with the PANTRY pipeline). It is resolved
relative to THIS file's location first, then falls back to PATH, so the pooler
works whether it lives alongside the main PANTRY scripts or in its own
subdirectory.
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

# Modalities whose unnorm BEDs are produced per-cohort and pooled here by
# phenotype_id intersection (original sample IDs).
INTERSECTION_MODALITIES = {
    "isoforms", "alt_TSS", "alt_polyA", "RNA_editing", "stability",
}

# Modalities already pooled by harmonize_within_ancestry.py (namespaced IDs).
PREPOOLED_MODALITIES = {"splicing", "intron_retention"}

ALL_MODALITIES = INTERSECTION_MODALITIES | PREPOOLED_MODALITIES


# ---------------------------------------------------------------------------
# Config / ancestry map loading
# ---------------------------------------------------------------------------

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


def resolve_config_get(scripts_dir=None):
    """Resolve the path to config_get.py.

    Order of precedence:
      1. Explicit scripts_dir argument
      2. Directory of THIS file (so the pooler works from a subdirectory)
      3. Repo layout: ../03_phenotyping relative to THIS file
      4. PATH fallback ("config_get.py")
    """
    if scripts_dir is None:
        scripts_dir = os.path.dirname(os.path.abspath(__file__))
    candidate = os.path.join(scripts_dir, "config_get.py")
    if os.path.exists(candidate):
        return candidate
    # repo layout: config_get.py lives in ../03_phenotyping
    candidate = os.path.join(scripts_dir, "..", "03_phenotyping", "config_get.py")
    if os.path.exists(candidate):
        return candidate
    return "config_get.py"  # last resort: cwd/PATH


def load_config(config_path, scripts_dir=None):
    """Load config.yml via config_get.py and parse exported vars.

    Returns a dict of config vars. config_get.py is resolved relative to THIS
    file's location first (so the pooler works from any directory), then falls
    back to PATH.
    """
    config_get = resolve_config_get(scripts_dir)

    result = subprocess.run(
        ["python3", config_get, config_path],
        capture_output=True, text=True, check=True
    )
    cfg = {}
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line.startswith("export "):
            continue
        line = line[len("export "):]
        key, _, value = line.partition("=")
        if value.startswith("'") and value.endswith("'"):
            value = value[1:-1].replace("'\"'\"'", "'")
        cfg[key] = value
    return cfg


def load_cohorts(config_path):
    """Extract cohort names from config.yml (the 'cohorts:' block)."""
    with open(config_path) as f:
        text = f.read()
    m = re.search(r"^cohorts:\s*$", text, re.MULTILINE)
    if not m:
        return []
    start = m.end()
    cohorts = []
    for line in text[start:].splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if not line.startswith(" ") and not line.startswith("\t"):
            break
        if stripped.endswith(":") and ":" not in stripped[:-1]:
            cohorts.append(stripped[:-1])
    return cohorts


# ---------------------------------------------------------------------------
# Per-cohort pooling (phenotype_id intersection)
# ---------------------------------------------------------------------------

def load_cohort_bed(cohort_dir, modality):
    """Load an unnorm BED file for a cohort.

    Returns a DataFrame with the 4 metadata columns + sample columns, or None
    if the file is absent.
    """
    bed_path = os.path.join(cohort_dir, "output", "unnorm", f"{modality}.bed")
    if not os.path.exists(bed_path):
        sys.stderr.write(f"  WARN: {bed_path} not found\n")
        return None
    df = pd.read_csv(bed_path, sep="\t",
                     dtype={"#chr": str, "start": int, "end": int,
                            "phenotype_id": str})
    return df


def pool_intersection_stratum(ancestry, ancestry_samples, cohort_dirs,
                              modality, output_dir):
    """Pool per-cohort BEDs for one ancestry stratum by phenotype_id intersection.

    Used for isoforms, alt_TSS, alt_polyA, RNA_editing, stability.
    """
    samples_by_cohort = defaultdict(list)
    for _, row in ancestry_samples.iterrows():
        samples_by_cohort[row["cohort"]].append(row["sample_id"])

    print(f"\n{'='*60}")
    print(f"Ancestry stratum: {ancestry} (N={len(ancestry_samples)})")
    print(f"  Modality: {modality} (intersection pooling)")
    print(f"  Cohorts: {dict(samples_by_cohort)}")

    cohort_dfs = {}
    feature_sets = {}
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

        bed_samples = [c for c in bed.columns if c not in BED_META_COLS]
        for s in sample_list:
            if s in bed_samples:
                samples_found[cohort].append(s)
            else:
                samples_missing[cohort].append(s)

        cols = BED_META_COLS + samples_found[cohort]
        cohort_dfs[cohort] = bed[cols].copy()
        feature_sets[cohort] = set(bed["phenotype_id"])

    if not cohort_dfs:
        sys.stderr.write(f"  ERROR: no cohort BEDs loaded for ancestry {ancestry}\n")
        return None

    for cohort, missing in samples_missing.items():
        if missing:
            print(f"  Missing samples in {cohort}: {len(missing)}")

    feature_intersection = set.intersection(*feature_sets.values())
    n_per_cohort = {c: len(g) for c, g in feature_sets.items()}
    print(f"  Features per cohort: {n_per_cohort}")
    print(f"  Feature intersection: {len(feature_intersection)}")

    all_samples = []
    merged_meta = None
    merged_data = []

    for cohort, df in cohort_dfs.items():
        df_f = df[df["phenotype_id"].isin(feature_intersection)].copy()
        df_f = df_f.sort_values("phenotype_id").reset_index(drop=True)

        if merged_meta is None:
            merged_meta = df_f[BED_META_COLS].copy()
        else:
            assert merged_meta["phenotype_id"].equals(df_f["phenotype_id"]), \
                f"phenotype_id mismatch when merging cohort {cohort}"

        sample_cols = [c for c in df_f.columns if c not in BED_META_COLS]
        merged_data.append(df_f[sample_cols])
        all_samples.extend(sample_cols)

    pooled_data = pd.concat(merged_data, axis=1)
    pooled = pd.concat([merged_meta, pooled_data], axis=1)

    assert len(pooled.columns) == 4 + len(all_samples)

    out_path = os.path.join(output_dir, f"{ancestry}_{modality}_pooled.bed")
    pooled.to_csv(out_path, sep="\t", index=False, float_format="%g")
    print(f"  Output: {out_path}")
    print(f"  Shape: {pooled.shape[0]} features x {len(all_samples)} samples")

    return {
        "ancestry": ancestry,
        "modality": modality,
        "pooling": "intersection",
        "n_samples": len(all_samples),
        "n_features": len(feature_intersection),
        "n_cohorts": len(cohort_dfs),
        "features_per_cohort": n_per_cohort,
        "samples_per_cohort": {c: len(s) for c, s in samples_found.items()},
        "missing_samples": {c: len(s) for c, s in samples_missing.items() if s},
        "output": out_path,
    }


# ---------------------------------------------------------------------------
# Pre-pooled pass-through (splicing, intron_retention)
# ---------------------------------------------------------------------------

def find_prepooled_bed(pre_pooled_dir, ancestry, modality):
    """Locate the harmonize-step pooled BED for a modality/ancestry.

    Expected layout (from harmonize_within_ancestry.py):
        <pre_pooled_dir>/<ancestry>/<modality>/unnorm/<modality>.bed
    """
    candidates = [
        os.path.join(pre_pooled_dir, ancestry, modality, "unnorm", f"{modality}.bed"),
        os.path.join(pre_pooled_dir, ancestry, modality, f"{modality}.bed"),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return None


def parse_cohort_from_namespaced(sample_ids):
    """Parse cohort labels from namespaced sample IDs ({cohort}_{sample}).

    Returns a dict sample_id -> cohort. Splits on the FIRST underscore; sample
    IDs are assumed not to contain underscores in the cohort prefix. If a
    sample ID has no underscore, it is assigned cohort "unknown".
    """
    labels = {}
    for s in sample_ids:
        if "_" in s:
            cohort, _, rest = s.partition("_")
            labels[s] = cohort
        else:
            labels[s] = "unknown"
    return labels


def pool_prepooled_stratum(ancestry, ancestry_samples, pre_pooled_dir,
                           modality, output_dir):
    """Pass through a pre-pooled (harmonized) BED for splicing/IR.

    No re-pooling is performed; the BED is copied to the output dir with the
    canonical naming convention. Cohort labels are parsed from the namespaced
    sample IDs and written to a sidecar TSV for the ComBat step.
    """
    print(f"\n{'='*60}")
    print(f"Ancestry stratum: {ancestry}")
    print(f"  Modality: {modality} (pre-pooled pass-through)")

    bed_path = find_prepooled_bed(pre_pooled_dir, ancestry, modality)
    if bed_path is None:
        sys.stderr.write(
            f"  ERROR: pre-pooled BED not found for {ancestry}/{modality} "
            f"under {pre_pooled_dir}\n")
        return None

    print(f"  Input: {bed_path}")
    df = pd.read_csv(bed_path, sep="\t",
                     dtype={"#chr": str, "start": int, "end": int,
                            "phenotype_id": str})
    sample_cols = [c for c in df.columns if c not in BED_META_COLS]
    print(f"  Shape: {df.shape[0]} features x {len(sample_cols)} samples")

    # Parse cohort labels from namespaced IDs
    cohort_labels = parse_cohort_from_namespaced(sample_cols)
    n_unknown = sum(1 for v in cohort_labels.values() if v == "unknown")
    if n_unknown:
        sys.stderr.write(
            f"  WARN: {n_unknown} sample IDs have no cohort prefix (underscore); "
            f"assigned cohort='unknown'\n")

    cohort_counts = pd.Series(cohort_labels).value_counts().to_dict()
    print(f"  Cohorts (from namespaced IDs): {cohort_counts}")

    # Write pooled BED (canonical name)
    out_path = os.path.join(output_dir, f"{ancestry}_{modality}_pooled.bed")
    df.to_csv(out_path, sep="\t", index=False, float_format="%g")
    print(f"  Output: {out_path}")

    # Write cohort-label sidecar (consumed by the ComBat R script)
    label_path = os.path.join(output_dir, f"{ancestry}_{modality}_cohort_labels.tsv")
    pd.DataFrame([
        {"sample_id": s, "cohort": cohort_labels[s]}
        for s in sample_cols
    ]).to_csv(label_path, sep="\t", index=False)
    print(f"  Cohort labels: {label_path}")

    return {
        "ancestry": ancestry,
        "modality": modality,
        "pooling": "prepooled",
        "n_samples": len(sample_cols),
        "n_features": df.shape[0],
        "n_cohorts": len(set(cohort_labels.values())),
        "samples_per_cohort": cohort_counts,
        "missing_samples": {},
        "output": out_path,
        "cohort_labels": label_path,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Pool unnorm RNA phenotype BEDs across cohorts within "
                    "ancestry strata (non-expression modalities)")
    parser.add_argument("--ancestry-map", required=True,
                        help="TSV: sample_id, assigned_ancestry, cohort")
    parser.add_argument("--config", help="config.yml (for cohort paths)")
    parser.add_argument("--output-base", help="Override output_base from config")
    parser.add_argument("--cohort-dirs", nargs="*",
                        help="Explicit cohort dirs (overrides config). "
                             "Format: cohort_name=/path/to/cohort_dir")
    parser.add_argument("--modality", required=True, choices=sorted(ALL_MODALITIES),
                        help="Which modality to pool")
    parser.add_argument("--pre-pooled-dir",
                        help="Root dir of harmonize_within_ancestry.py outputs "
                             "(required for splicing/intron_retention)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for pooled matrices")
    parser.add_argument("--ancestries", nargs="*",
                        help="Subset of ancestries to process (default: all)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Validate modality / mode compatibility
    if args.modality in PREPOOLED_MODALITIES:
        if not args.pre_pooled_dir:
            sys.stderr.write(
                f"ERROR: --pre-pooled-dir required for modality '{args.modality}' "
                f"(splicing/intron_retention are pre-pooled by harmonize_within_ancestry.py)\n")
            sys.exit(1)
    else:
        if args.pre_pooled_dir:
            sys.stderr.write(
                f"ERROR: --pre-pooled-dir not applicable for modality '{args.modality}' "
                f"(only splicing/intron_retention are pre-pooled)\n")
            sys.exit(1)

    # Load ancestry map
    ancestry_map = load_ancestry_map(args.ancestry_map)
    print(f"Ancestry map: {len(ancestry_map)} samples")
    print(f"  Ancestries: {sorted(ancestry_map['assigned_ancestry'].unique())}")
    print(f"  Cohorts: {sorted(ancestry_map['cohort'].unique())}")

    # Resolve cohort directories (only needed for intersection modalities)
    cohort_dirs = {}
    if args.modality in INTERSECTION_MODALITIES:
        if args.cohort_dirs:
            for item in args.cohort_dirs:
                if "=" in item:
                    name, path = item.split("=", 1)
                    cohort_dirs[name] = path
                else:
                    sys.stderr.write(f"  WARN: ignoring malformed --cohort-dirs entry: {item}\n")
        elif args.config:
            cfg = load_config(args.config)
            output_base = args.output_base or cfg.get("OUTPUT_BASE", "")
            cohort_names = load_cohorts(args.config)
            for c in cohort_names:
                cohort_dirs[c] = os.path.join(output_base, c)
            print(f"  Config cohorts: {cohort_names}")
            print(f"  Output base: {output_base}")
        else:
            sys.stderr.write("ERROR: need --config or --cohort-dirs for intersection modalities\n")
            sys.exit(1)

    ancestries = args.ancestries or sorted(ancestry_map["assigned_ancestry"].unique())

    all_stats = []
    for ancestry in ancestries:
        stratum_samples = ancestry_map[ancestry_map["assigned_ancestry"] == ancestry]
        if len(stratum_samples) == 0:
            print(f"\nNo samples for ancestry {ancestry}, skipping")
            continue

        if args.modality in PREPOOLED_MODALITIES:
            stats = pool_prepooled_stratum(
                ancestry, stratum_samples, args.pre_pooled_dir,
                args.modality, args.output_dir)
        else:
            stats = pool_intersection_stratum(
                ancestry, stratum_samples, cohort_dirs,
                args.modality, args.output_dir)
        if stats:
            all_stats.append(stats)

    # Summary
    print(f"\n{'='*60}")
    print("POOLING SUMMARY")
    print(f"{'='*60}")
    for s in all_stats:
        print(f"  {s['ancestry']} / {s['modality']}: {s['n_samples']} samples, "
              f"{s['n_features']} features, {s['n_cohorts']} cohorts "
              f"({s['pooling']})")
        print(f"    Samples/cohort: {s['samples_per_cohort']}")
        if s["missing_samples"]:
            print(f"    Missing samples: {s['missing_samples']}")

    summary_path = os.path.join(args.output_dir, "pooling_summary.tsv")
    summary_rows = []
    for s in all_stats:
        for cohort, n in s["samples_per_cohort"].items():
            summary_rows.append({
                "ancestry": s["ancestry"],
                "modality": s["modality"],
                "cohort": cohort,
                "n_samples": n,
                "n_features": s["n_features"],
                "pooling": s["pooling"],
                "n_missing": s["missing_samples"].get(cohort, 0),
            })
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    print(f"\nSummary: {summary_path}")


if __name__ == "__main__":
    main()
