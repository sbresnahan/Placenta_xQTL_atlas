#!/usr/bin/env python3
"""
harmonize_within_ancestry.py — Cross-cohort coordinate harmonization for
splicing (leafCutter) and intron retention (MAJIQ) within an ancestry stratum.

Pools per-cohort intermediates, harmonizes features by stable genomic
coordinates, and writes pooled intermediates in the tool's native format so
that assemble_bed.py can produce a unified unnorm/ BED ready for downstream
ComBat + normalization.

Splicing strategy:
    leafCutter cluster numbers are cohort-specific. Junctions are harmonized
    by (chrom, start, end, strand). Per-cohort clusters sharing >=1 junction
    are linked in a graph; connected components become meta-clusters.
    Proportions are recomputed within meta-clusters from pooled raw counts.

IR strategy:
    MAJIQ IR event IDs are build-specific. Events are matched by exact
    (seqid, start, end, strand, gene_id_base). All cohorts use the same MAJIQ
    GFF3 from Tier 1 shared reference prep, so coordinates are expected to be
    identical across builds. Events present in only one cohort are included
    with NaN PSI for other cohorts' samples.

Usage:
    harmonize_within_ancestry.py splicing \\
        --config config.yml \\
        --ancestry-map samples_ancestry.tsv \\
        --ancestry EUR \\
        --scripts-dir /path/to/scripts \\
        --output-dir /path/to/output/EUR/splicing

    harmonize_within_ancestry.py intron-retention \\
        --config config.yml \\
        --ancestry-map samples_ancestry.tsv \\
        --ancestry EUR \\
        --scripts-dir /path/to/scripts \\
        --output-dir /path/to/output/EUR/intron_retention
"""

import argparse
import gzip
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


# ---------------------------------------------------------------------------
# Config loading (reuses config_get.py)
# ---------------------------------------------------------------------------

def load_config(config_path: str, scripts_dir: str) -> dict:
    """Load config.yml via config_get.py and parse the exported vars."""
    result = subprocess.run(
        ["python3", os.path.join(scripts_dir, "config_get.py"), config_path],
        capture_output=True, text=True, check=True,
    )
    cfg = {}
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line.startswith("export "):
            continue
        line = line[len("export "):]
        key, _, value = line.partition("=")
        # Strip surrounding single quotes (config_get.py shell-quotes values)
        if value.startswith("'") and value.endswith("'"):
            value = value[1:-1].replace("'\"'\"'", "'")
        cfg[key] = value
    return cfg


def load_cohorts(config_path: str) -> list:
    """Extract cohort names from config.yml (the 'cohorts:' block)."""
    import re
    with open(config_path) as f:
        text = f.read()
    # Find the cohorts: block and parse 2-space-indented cohort names
    m = re.search(r"^cohorts:\s*$", text, re.MULTILINE)
    if not m:
        return []
    start = m.end()
    # Cohort names are indented under 'cohorts:' and end with ':'
    cohorts = []
    for line in text[start:].splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        # Stop if we hit a top-level (non-indented) key
        if not line.startswith(" ") and not line.startswith("\t"):
            break
        # Cohort lines look like "  cohort1:" or "  GUSTO:"
        if stripped.endswith(":") and ":" not in stripped[:-1]:
            cohorts.append(stripped[:-1])
    return cohorts


# ---------------------------------------------------------------------------
# Ancestry / sample loading
# ---------------------------------------------------------------------------

def load_ancestry_map(path: str, sample_col: str = "sample_id",
                      ancestry_col: str = "ancestry") -> pd.DataFrame:
    """Load the sample -> ancestry TSV."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if sample_col not in df.columns or ancestry_col not in df.columns:
        raise ValueError(
            f"Ancestry map {path} must contain columns '{sample_col}' and "
            f"'{ancestry_col}'. Found: {list(df.columns)}"
        )
    return df[[sample_col, ancestry_col]]


def get_cohort_samples_in_ancestry(cohort_samples_file: str,
                                   ancestry_df: pd.DataFrame,
                                   ancestry: str,
                                   sample_col: str = "sample_id",
                                   ancestry_col: str = "ancestry") -> list:
    """Return sample IDs in a cohort that belong to the target ancestry."""
    with open(cohort_samples_file) as f:
        cohort_samples = [line.strip() for line in f if line.strip()]
    ancestry_samples = set(
        ancestry_df.loc[ancestry_df[ancestry_col] == ancestry, sample_col]
    )
    return [s for s in cohort_samples if s in ancestry_samples]


# ---------------------------------------------------------------------------
# Splicing harmonization
# ---------------------------------------------------------------------------

def parse_leafcutter_counts(path: str) -> tuple:
    """Parse leafcutter_perind_numers.counts.gz.

    Returns:
        sample_ids: list of sample column names
        rows: list of (junction_key, cohort_cluster_id, counts_array)
              where junction_key = (chrom, start, end, strand)
              and cohort_cluster_id = int (the clu_N from the row index)
    """
    with gzip.open(path, "rt") as f:
        header = f.readline().strip().split(" ")
        # header[0] is "chrom", rest are sample IDs
        sample_ids = header[1:]
        rows = []
        for line in f:
            parts = line.strip().split(" ")
            if len(parts) < 2:
                continue
            row_id = parts[0]
            counts = np.array([float(x) for x in parts[1:]], dtype=np.float64)
            # row_id format: chrom:start:end:clu_N_strand
            # Chromosomes may contain colons in some builds, but leafCutter
            # uses the format chrom:start:end:clu_N_strand where the last
            # three colon-fields are start, end, and clu_N_strand.
            fields = row_id.split(":")
            if len(fields) < 4:
                continue
            cluster_field = fields[-1]  # clu_N_strand
            end = int(fields[-2])
            start = int(fields[-3])
            chrom = ":".join(fields[:-3])
            # cluster_field = clu_N_strand -> split by _
            sub = cluster_field.split("_")
            # Expected: ["clu", N, strand]
            if len(sub) < 3:
                continue
            strand = sub[-1]
            clu_n = int(sub[-2])
            junction_key = (chrom, start, end, strand)
            rows.append((junction_key, clu_n, counts))
    return sample_ids, rows


def harmonize_splicing(config: dict, cohorts: list, ancestry_samples: dict,
                       output_dir: Path, scripts_dir: str,
                       ref_anno: str) -> None:
    """Harmonize splicing across cohorts within an ancestry stratum."""
    output_base = config["OUTPUT_BASE"]
    pantry_scripts = config["PANTRY_SCRIPTS"]

    intermediate_dir = output_dir / "intermediate"
    unnorm_dir = output_dir / "unnorm"
    intermediate_dir.mkdir(parents=True, exist_ok=True)
    unnorm_dir.mkdir(parents=True, exist_ok=True)

    # --- Load per-cohort leafCutter counts ---
    # cohort_data[cohort] = {
    #   'sample_ids': [...],  # namespaced
    #   'rows': [(junction_key, clu_n, counts_subset), ...]
    # }
    cohort_data = {}
    # junction -> set of (cohort, clu_n) that contain it
    junction_to_clusters = {}

    for cohort in cohorts:
        counts_path = Path(output_base) / cohort / "intermediate" / "splicing" / "leafcutter_perind_numers.counts.gz"
        if not counts_path.exists():
            print(f"[WARN] {cohort}: splicing counts not found at {counts_path}, skipping", file=sys.stderr)
            continue
        cohort_samples = ancestry_samples.get(cohort, [])
        if not cohort_samples:
            print(f"[WARN] {cohort}: 0 samples in ancestry, skipping splicing", file=sys.stderr)
            continue

        sample_ids, rows = parse_leafcutter_counts(str(counts_path))
        # Subset columns to ancestry samples
        col_idx = []
        kept_samples = []
        for s in cohort_samples:
            if s in sample_ids:
                col_idx.append(sample_ids.index(s))
                kept_samples.append(s)
        if not col_idx:
            print(f"[WARN] {cohort}: no overlap between cohort samples and leafCutter columns, skipping", file=sys.stderr)
            continue

        col_idx = np.array(col_idx)
        namespaced_samples = [f"{cohort}_{s}" for s in kept_samples]
        subset_rows = []
        for junction_key, clu_n, counts in rows:
            counts_subset = counts[col_idx]
            subset_rows.append((junction_key, clu_n, counts_subset))
            junction_to_clusters.setdefault(junction_key, set()).add((cohort, clu_n))

        cohort_data[cohort] = {
            "sample_ids": namespaced_samples,
            "rows": subset_rows,
        }
        print(f"[INFO] {cohort}: loaded {len(subset_rows)} junctions, {len(kept_samples)} samples")

    if not cohort_data:
        raise RuntimeError("No cohort splicing data available for this ancestry")

    # --- Build connected components of (cohort, clu_n) by shared junctions ---
    # Index all (cohort, clu_n) pairs
    all_clusters = set()
    for cohort, data in cohort_data.items():
        for _, clu_n, _ in data["rows"]:
            all_clusters.add((cohort, clu_n))
    all_clusters = sorted(all_clusters)
    cluster_index = {c: i for i, c in enumerate(all_clusters)}
    n_clusters = len(all_clusters)

    # Build adjacency via shared junctions
    if n_clusters == 0:
        raise RuntimeError("No clusters to harmonize")

    edges_src = []
    edges_dst = []
    for junction_key, cluster_set in junction_to_clusters.items():
        cluster_list = sorted(cluster_set)
        for i in range(len(cluster_list)):
            for j in range(i + 1, len(cluster_list)):
                edges_src.append(cluster_index[cluster_list[i]])
                edges_dst.append(cluster_index[cluster_list[j]])

    if HAVE_SCIPY and edges_src:
        adj = coo_matrix(
            (np.ones(len(edges_src)), (edges_src, edges_dst)),
            shape=(n_clusters, n_clusters),
        )
        n_components, labels = connected_components(adj, directed=False)
    else:
        # Fallback: union-find
        n_components, labels = _union_find(n_clusters, edges_src, edges_dst)

    print(f"[INFO] {n_clusters} per-cohort clusters -> {n_components} meta-clusters")

    # Map (cohort, clu_n) -> meta_cluster_id
    cluster_to_meta = {}
    for (cohort, clu_n), idx in cluster_index.items():
        cluster_to_meta[(cohort, clu_n)] = int(labels[idx])

    # Renumber meta-clusters as clu_1, clu_2, ...
    unique_meta = sorted(set(cluster_to_meta.values()))
    meta_renumber = {old: new + 1 for new, old in enumerate(unique_meta)}

    # --- Build pooled counts matrix ---
    # Collect all junctions across all cohorts, keyed by (junction_key, meta_cluster)
    # meta_cluster determines the proportion group
    pooled = {}  # (junction_key, meta_cluster) -> {sample: count}
    all_samples = []
    for cohort, data in cohort_data.items():
        all_samples.extend(data["sample_ids"])

    for cohort, data in cohort_data.items():
        for junction_key, clu_n, counts_subset in data["rows"]:
            meta = meta_renumber[cluster_to_meta[(cohort, clu_n)]]
            key = (junction_key, meta)
            if key not in pooled:
                pooled[key] = {s: 0.0 for s in all_samples}
            for i, s in enumerate(data["sample_ids"]):
                pooled[key][s] += counts_subset[i]

    # Write pooled counts in leafCutter format
    pooled_counts_path = intermediate_dir / "leafcutter_harmonized_perind_numers.counts.gz"
    with gzip.open(pooled_counts_path, "wt") as f:
        # Header must be sample names ONLY (no "chrom" token). assemble_bed.py
        # reads this with pd.read_csv(sep=' '), and pandas auto-uses column 1
        # as the string row index only when the header has one fewer field
        # than the data rows -- matching per-cohort leafcutter numers format.
        f.write(" ".join(all_samples) + "\n")
        for (junction_key, meta), sample_counts in pooled.items():
            chrom, start, end, strand = junction_key
            row_id = f"{chrom}:{start}:{end}:clu_{meta}_{strand}"
            counts_str = " ".join(str(int(sample_counts[s])) for s in all_samples)
            f.write(row_id + " " + counts_str + "\n")

    print(f"[INFO] Wrote pooled counts: {pooled_counts_path} ({len(pooled)} junctions)")

    # Write harmonization map
    harm_map_path = intermediate_dir / "harmonization_map.tsv"
    with open(harm_map_path, "w") as f:
        f.write("cohort\tcohort_cluster\tmeta_cluster\tn_junctions\tn_samples\n")
        # Count junctions per (cohort, clu_n)
        cohort_cluster_junctions = {}
        for cohort, data in cohort_data.items():
            for _, clu_n, _ in data["rows"]:
                key = (cohort, clu_n)
                cohort_cluster_junctions[key] = cohort_cluster_junctions.get(key, 0) + 1
        for (cohort, clu_n), idx in cluster_index.items():
            meta = meta_renumber[cluster_to_meta[(cohort, clu_n)]]
            n_junc = cohort_cluster_junctions.get((cohort, clu_n), 0)
            n_samples = len(cohort_data[cohort]["sample_ids"])
            f.write(f"{cohort}\t{clu_n}\t{meta}\t{n_junc}\t{n_samples}\n")

    # --- Call assemble_bed.py splicing ---
    splicing_bed = unnorm_dir / "splicing.bed"
    cmd = [
        "python3", os.path.join(pantry_scripts, "assemble_bed.py"), "splicing",
        "--input", str(pooled_counts_path),
        "--ref-anno", ref_anno,
        "--output", str(splicing_bed),
    ]
    print(f"[INFO] Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # Write phenotype_groups.txt
    groups_path = unnorm_dir / "splicing.phenotype_groups.txt"
    _write_phenotype_groups(splicing_bed, groups_path)

    print(f"[INFO] Done: splicing harmonized -> {splicing_bed}")


# ---------------------------------------------------------------------------
# Intron retention harmonization
# ---------------------------------------------------------------------------

def harmonize_intron_retention(config: dict, cohorts: list,
                               ancestry_samples: dict, output_dir: Path,
                               scripts_dir: str, ref_anno: str) -> None:
    """Harmonize intron retention across cohorts within an ancestry stratum."""
    output_base = config["OUTPUT_BASE"]
    pantry_scripts = config["PANTRY_SCRIPTS"]

    intermediate_dir = output_dir / "intermediate"
    unnorm_dir = output_dir / "unnorm"
    intermediate_dir.mkdir(parents=True, exist_ok=True)
    unnorm_dir.mkdir(parents=True, exist_ok=True)

    # Metadata columns in retained_intron_psi.tsv.gz (from extract_ir_psi.py)
    metadata_cols = [
        "majiq_ir_id", "ec_idx", "seqid", "start", "end", "strand",
        "gene_id", "gene_id_base", "gene_name", "event_type",
        "is_denovo", "event_denovo", "ref_exon_start", "ref_exon_end",
        "other_exon_start", "other_exon_end",
    ]

    # --- Load per-cohort IR PSI tables ---
    # event_data[event_key] = {
    #   'metadata': {...},          # from first cohort that has it (for output)
    #   'cohort_meta': {cohort: {'majiq_ir_id': ..., ...}, ...},  # per-cohort original IDs
    #   'cohort_psi': {cohort: {namespaced_sample: psi, ...}, ...}
    # }
    event_data = {}
    all_samples = []

    for cohort in cohorts:
        psi_path = Path(output_base) / cohort / "intermediate" / "intron_retention" / "retained_intron_psi.tsv.gz"
        if not psi_path.exists():
            print(f"[WARN] {cohort}: IR PSI not found at {psi_path}, skipping", file=sys.stderr)
            continue
        cohort_samples = ancestry_samples.get(cohort, [])
        if not cohort_samples:
            print(f"[WARN] {cohort}: 0 samples in ancestry, skipping IR", file=sys.stderr)
            continue

        df = pd.read_csv(psi_path, sep="\t", dtype={"seqid": str})
        # Identify sample columns (everything not in metadata_cols)
        sample_cols = [c for c in df.columns if c not in metadata_cols]
        # Subset to ancestry samples
        kept_cols = [c for c in cohort_samples if c in sample_cols]
        if not kept_cols:
            print(f"[WARN] {cohort}: no overlap between cohort samples and IR columns, skipping", file=sys.stderr)
            continue

        namespaced = [f"{cohort}_{s}" for s in kept_cols]
        all_samples.extend(namespaced)

        for _, row in df.iterrows():
            gene_id_base = str(row.get("gene_id_base", "")).split(".")[0]
            event_key = (
                str(row["seqid"]),
                int(row["start"]),
                int(row["end"]),
                str(row["strand"]),
                gene_id_base,
            )
            if event_key not in event_data:
                # Store metadata from first cohort that has this event
                meta = {}
                for mc in metadata_cols:
                    if mc in row.index:
                        meta[mc] = row[mc]
                event_data[event_key] = {
                    "metadata": meta,
                    "cohort_meta": {},
                    "cohort_psi": {},
                }
            # Store this cohort's original metadata (for the harmonization map)
            cohort_meta = {}
            for mc in metadata_cols:
                if mc in row.index:
                    cohort_meta[mc] = row[mc]
            event_data[event_key]["cohort_meta"][cohort] = cohort_meta
            # Store this cohort's PSI values
            psi_vals = {}
            for s, ns in zip(kept_cols, namespaced):
                v = row[s]
                psi_vals[ns] = v if pd.notna(v) else np.nan
            event_data[event_key]["cohort_psi"][cohort] = psi_vals

        print(f"[INFO] {cohort}: loaded {len(df)} IR events, {len(kept_cols)} samples")

    if not event_data:
        raise RuntimeError("No cohort IR data available for this ancestry")

    # --- Build unified PSI table ---
    # Sort event keys for deterministic ordering
    sorted_events = sorted(event_data.keys(), key=lambda k: (k[0], k[1], k[2], k[3], k[4]))

    # Assign unified majiq_ir_id
    unified_rows = []
    harm_map_rows = []
    for unified_idx, event_key in enumerate(sorted_events):
        ed = event_data[event_key]
        seqid, start, end, strand, gene_id_base = event_key
        unified_ir_id = (
            f"IR:{gene_id_base}:{seqid}:{start}:{end}:{strand}:{unified_idx}"
        )
        meta = ed["metadata"].copy()
        meta["majiq_ir_id"] = unified_ir_id
        meta["ec_idx"] = unified_idx

        row = {}
        for mc in metadata_cols:
            row[mc] = meta.get(mc, "")
        # Fill in PSI values from all cohorts
        for s in all_samples:
            row[s] = np.nan
        for cohort, psi_vals in ed["cohort_psi"].items():
            for ns, v in psi_vals.items():
                row[ns] = v
        unified_rows.append(row)

        # Harmonization map: every cohort's original event ID -> unified event
        for cohort in ed["cohort_psi"]:
            cohort_original_id = ed["cohort_meta"][cohort].get("majiq_ir_id", "")
            harm_map_rows.append({
                "cohort": cohort,
                "cohort_majiq_ir_id": cohort_original_id,
                "unified_majiq_ir_id": unified_ir_id,
                "event_key": f"{seqid}:{start}:{end}:{strand}:{gene_id_base}",
            })

    unified_df = pd.DataFrame(unified_rows, columns=metadata_cols + all_samples)

    # Write harmonized PSI TSV
    harmonized_psi_path = intermediate_dir / "retained_intron_psi_harmonized.tsv.gz"
    unified_df.to_csv(harmonized_psi_path, sep="\t", index=False,
                      compression="gzip", float_format="%g")
    print(f"[INFO] Wrote harmonized PSI: {harmonized_psi_path} ({len(unified_df)} events)")

    # Write harmonization map
    harm_map_path = intermediate_dir / "harmonization_map.tsv"
    harm_df = pd.DataFrame(harm_map_rows)
    harm_df.to_csv(harm_map_path, sep="\t", index=False)
    print(f"[INFO] Wrote harmonization map: {harm_map_path} ({len(harm_map_rows)} rows)")

    # --- Call assemble_bed.py intron-retention ---
    ir_bed = unnorm_dir / "intron_retention.bed"
    cmd = [
        "python3", os.path.join(pantry_scripts, "assemble_bed.py"), "intron-retention",
        "--input", str(harmonized_psi_path),
        "--ref-anno", ref_anno,
        "--output", str(ir_bed),
    ]
    print(f"[INFO] Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # Write phenotype_groups.txt
    groups_path = unnorm_dir / "intron_retention.phenotype_groups.txt"
    _write_phenotype_groups(ir_bed, groups_path)

    print(f"[INFO] Done: intron retention harmonized -> {ir_bed}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _union_find(n: int, edges_src: list, edges_dst: list) -> tuple:
    """Fallback connected components via union-find (no scipy)."""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for s, d in zip(edges_src, edges_dst):
        union(s, d)

    # Compress and renumber
    roots = {}
    labels = np.zeros(n, dtype=int)
    for i in range(n):
        r = find(i)
        if r not in roots:
            roots[r] = len(roots)
        labels[i] = roots[r]
    return len(roots), labels


def _write_phenotype_groups(bed_path: Path, groups_path: Path) -> None:
    """Write phenotype_groups.txt from a BED file (gene grouping for tensorQTL).

    phenotype_id format: {gene_id}__{rest} -> group by gene_id.
    """
    import subprocess
    # Use the same awk logic as scripts 12/13
    cmd = (
        f"cat {bed_path} | tail -n +2 | cut -f4 | "
        f"awk '{{ g=$1; sub(/__.*$/, \"\", g); print $1 \"\\t\" g }}' > {groups_path}"
    )
    subprocess.run(cmd, shell=True, check=True)


def write_sample_namespace_map(cohort_samples: dict, output_path: Path) -> None:
    """Write {cohort}_{sample} -> sample mapping for traceability."""
    rows = []
    for cohort, samples in cohort_samples.items():
        for s in samples:
            rows.append({"namespaced_id": f"{cohort}_{s}", "sample_id": s, "cohort": cohort})
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Cross-cohort coordinate harmonization for splicing & IR within ancestry"
    )
    sub = parser.add_subparsers(dest="modality", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--config", required=True, help="config.yml path")
    common.add_argument("--ancestry-map", required=True,
                        help="TSV mapping sample_id -> ancestry")
    common.add_argument("--ancestry", required=True,
                        help="Ancestry group to process (e.g. EUR, EAS, AFR, HIS, SAS)")
    common.add_argument("--scripts-dir", required=True,
                        help="Directory containing config_get.py and stage scripts")
    common.add_argument("--output-dir", required=True,
                        help="Output directory for this ancestry/modality")
    common.add_argument("--sample-col", default="sample_id",
                        help="Column name for sample IDs in ancestry map (default: sample_id)")
    common.add_argument("--ancestry-col", default="ancestry",
                        help="Column name for ancestry in ancestry map (default: ancestry)")

    sub.add_parser("splicing", parents=[common], help="Harmonize splicing")
    sub.add_parser("intron-retention", parents=[common], help="Harmonize intron retention")

    args = parser.parse_args()

    # Load config
    config = load_config(args.config, args.scripts_dir)
    cohorts = load_cohorts(args.config)
    if not cohorts:
        raise RuntimeError("No cohorts found in config.yml")

    ref_anno = config.get("NORMALIZED_GTF")
    if not ref_anno:
        raise RuntimeError("NORMALIZED_GTF not found in config")

    # Load ancestry map
    ancestry_df = load_ancestry_map(args.ancestry_map, args.sample_col, args.ancestry_col)

    # Get samples per cohort within this ancestry
    output_base = config["OUTPUT_BASE"]
    ancestry_samples = {}
    all_namespaced = {}
    for cohort in cohorts:
        samples_file = Path(output_base) / cohort / "samples.txt"
        if not samples_file.exists():
            print(f"[WARN] {cohort}: samples.txt not found, skipping", file=sys.stderr)
            continue
        cohort_samples = get_cohort_samples_in_ancestry(
            str(samples_file), ancestry_df, args.ancestry,
            args.sample_col, args.ancestry_col,
        )
        if cohort_samples:
            ancestry_samples[cohort] = cohort_samples
            all_namespaced[cohort] = cohort_samples
            print(f"[INFO] {cohort}: {len(cohort_samples)} samples in ancestry {args.ancestry}")

    if not ancestry_samples:
        raise RuntimeError(f"No samples found for ancestry '{args.ancestry}' across cohorts")

    # Write sample namespace map
    output_dir = Path(args.output_dir)
    namespace_map_path = output_dir.parent / "sample_namespace_map.tsv"
    namespace_map_path.parent.mkdir(parents=True, exist_ok=True)
    write_sample_namespace_map(all_namespaced, namespace_map_path)
    print(f"[INFO] Wrote sample namespace map: {namespace_map_path}")

    # Dispatch
    if args.modality == "splicing":
        harmonize_splicing(config, cohorts, ancestry_samples, output_dir,
                           args.scripts_dir, ref_anno)
    elif args.modality == "intron-retention":
        harmonize_intron_retention(config, cohorts, ancestry_samples, output_dir,
                                   args.scripts_dir, ref_anno)


if __name__ == "__main__":
    main()
