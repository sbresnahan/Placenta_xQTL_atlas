#!/usr/bin/env python3
"""Script 5: Pool per-cohort deconvolution outputs and generate QC report.

Reads the per-cohort output TSVs from run_music_deconvolution.R and produces:
  - all_cohorts_cell_proportions_full.tsv       (pooled samples x 27 types)
  - all_cohorts_cell_proportions_collapsed.tsv  (pooled samples x 8 types)
  - all_cohorts_maternal_flag.tsv               (pooled maternal flags)
  - deconvolution_qc_report.md                  (QC summary report)

Usage:
  python3 pool_deconvolution_outputs.py \
    --input-dir /mnt/results/deconvolution \
    --output-dir /mnt/results/deconvolution \
    --ancestry-map /path/to/sample_ancestry.tsv  # optional: sample -> ancestry

The ancestry map is an optional TSV with columns: sample_id, ancestry.
If provided, the QC report includes ancestry-stratified summaries.
"""

import argparse
import glob
import os
from pathlib import Path

import pandas as pd
import numpy as np


def find_cohort_files(input_dir):
    """Find all per-cohort output files."""
    full_files = sorted(glob.glob(os.path.join(input_dir, "*_cell_proportions_full.tsv")))
    collapsed_files = sorted(glob.glob(os.path.join(input_dir, "*_cell_proportions_collapsed.tsv")))
    flag_files = sorted(glob.glob(os.path.join(input_dir, "*_maternal_flag.tsv")))

    # Extract cohort names from filenames (strip suffix)
    cohorts = []
    for f in full_files:
        base = os.path.basename(f)
        cohort = base.replace("_cell_proportions_full.tsv", "")
        if cohort != "all_cohorts":  # skip previously pooled files
            cohorts.append(cohort)

    return cohorts, full_files, collapsed_files, flag_files


def pool_files(files, label):
    """Concatenate a list of TSV files, adding a 'cohort' column."""
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        base = os.path.basename(f)
        # Extract cohort from filename
        if "_cell_proportions_full" in base:
            cohort = base.replace("_cell_proportions_full.tsv", "")
        elif "_cell_proportions_collapsed" in base:
            cohort = base.replace("_cell_proportions_collapsed.tsv", "")
        elif "_maternal_flag" in base:
            cohort = base.replace("_maternal_flag.tsv", "")
        else:
            cohort = "unknown"
        df.insert(1, 'cohort', cohort)
        dfs.append(df)
    if not dfs:
        raise ValueError(f"No {label} files found")
    pooled = pd.concat(dfs, ignore_index=True)
    return pooled, [os.path.basename(f) for f in files]


def generate_qc_report(full, collapsed, flags, cohorts, ancestry_map=None):
    """Generate a markdown QC report."""
    lines = []
    lines.append("# Cell-Type Deconvolution QC Report\n")
    lines.append(f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n")
    lines.append(f"**Method:** MuSiC (Wang et al. 2019) with Campbell et al. 2023 reference\n")
    lines.append(f"**Reference:** GSE182381 — 16,003 genes, 40,494 cells, 27 cell types (19 fetal + 8 maternal)\n")

    # --- Overview ---
    lines.append("## Overview\n")
    lines.append(f"| Metric | Value |")
    lines.append(f"|---|---|")
    lines.append(f"| Cohorts | {len(cohorts)} |")
    lines.append(f"| Total samples | {len(full)} |")
    lines.append(f"| Full cell types | {full.shape[1] - 2} |")  # minus sample + cohort cols
    lines.append(f"| Collapsed cell types | {collapsed.shape[1] - 2} |")
    n_flagged = flags['flagged'].sum() if 'flagged' in flags.columns else 0
    lines.append(f"| Samples flagged (>10% maternal) | {n_flagged} ({100*n_flagged/len(flags):.1f}%) |")
    lines.append("")

    # --- Per-cohort sample counts ---
    lines.append("## Per-cohort sample counts\n")
    counts = full.groupby('cohort').size().reset_index(name='n_samples')
    counts['n_flagged'] = counts['cohort'].map(
        flags.groupby('cohort')['flagged'].sum() if 'flagged' in flags.columns else 0
    )
    lines.append("| Cohort | N samples | N flagged |")
    lines.append("|---|---|---|")
    for _, row in counts.iterrows():
        lines.append(f"| {row['cohort']} | {row['n_samples']} | {row['n_flagged']} |")
    lines.append("")

    # --- Collapsed proportion summary (median per cohort) ---
    lines.append("## Collapsed cell-type proportions (median per cohort)\n")
    prop_cols = [c for c in collapsed.columns if c not in ['sample', 'cohort']]
    summary = collapsed.groupby('cohort')[prop_cols].median().round(4)
    lines.append("| Cohort | " + " | ".join(prop_cols) + " |")
    lines.append("|---|" + "---|" * len(prop_cols))
    for cohort, row in summary.iterrows():
        lines.append(f"| {cohort} | " + " | ".join(f"{row[c]:.4f}" for c in prop_cols) + " |")
    lines.append("")

    # --- Overall collapsed proportion distribution ---
    lines.append("## Overall collapsed proportion distribution\n")
    lines.append("| Cell type | Min | Median | Mean | Max |")
    lines.append("|---|---|---|---|---|")
    for col in prop_cols:
        vals = collapsed[col]
        lines.append(f"| {col} | {vals.min():.4f} | {vals.median():.4f} | {vals.mean():.4f} | {vals.max():.4f} |")
    lines.append("")

    # --- Maternal fraction distribution ---
    lines.append("## Maternal fraction distribution\n")
    if 'maternal_fraction' in flags.columns:
        mf = flags['maternal_fraction']
        lines.append(f"| Statistic | Value |")
        lines.append(f"|---|---|")
        lines.append(f"| Min | {mf.min():.4f} |")
        lines.append(f"| Median | {mf.median():.4f} |")
        lines.append(f"| Mean | {mf.mean():.4f} |")
        lines.append(f"| Max | {mf.max():.4f} |")
        lines.append(f"| Samples >5% | {(mf > 0.05).sum()} |")
        lines.append(f"| Samples >10% | {(mf > 0.10).sum()} |")
        lines.append(f"| Samples >20% | {(mf > 0.20).sum()} |")
    lines.append("")

    # --- Ancestry-stratified summary (if map provided) ---
    if ancestry_map is not None:
        lines.append("## Ancestry-stratified summary\n")
        merged = collapsed.merge(ancestry_map, on='sample', how='left')
        if 'ancestry' in merged.columns:
            anc_summary = merged.groupby('ancestry')[prop_cols].median().round(4)
            lines.append("| Ancestry | " + " | ".join(prop_cols) + " |")
            lines.append("|---|" + "---|" * len(prop_cols))
            for ancestry, row in anc_summary.iterrows():
                lines.append(f"| {ancestry} | " + " | ".join(f"{row[c]:.4f}" for c in prop_cols) + " |")
            lines.append("")

    # --- Proportion sum check ---
    lines.append("## Proportion sum validation\n")
    row_sums = collapsed[prop_cols].sum(axis=1)
    lines.append(f"| Check | Result |")
    lines.append(f"|---|---|")
    lines.append(f"| All row sums = 1.0 (±0.01) | {((row_sums > 0.99) & (row_sums < 1.01)).all()} |")
    lines.append(f"| Min row sum | {row_sums.min():.6f} |")
    lines.append(f"| Max row sum | {row_sums.max():.6f} |")
    lines.append("")

    # --- Assumptions and limitations ---
    lines.append("## Assumptions and limitations\n")
    lines.append("- **Pseudo-subject variance:** MuSiC's cross-subject weighting uses 5 pseudo-subjects "
                 "(random cell splits within each cell type) rather than true biological replicates. "
                 "The integrated Seurat object with biorep metadata was not deposited publicly (GEO GSE182381 "
                 "only has a KC-only object with 2 bioreps and no trophoblasts). Gene weights reflect "
                 "within-reference technical variation, not inter-individual biological variation.\n")
    lines.append("- **Gene symbol overlap:** The reference uses GENCODE v32 symbols; bulk uses HPLRv2 "
                 "gene names (SQANTI3 associated_gene + GENCODE v45 override). Novel HPLRv2 genes without "
                 "GENCODE matches are absent from deconvolution. See per-cohort gene overlap in deconvolution logs.\n")
    lines.append("- **Platform mismatch:** Reference is scRNA-seq counts (10x Chromium); bulk is Salmon TPM "
                 "(short-read RNA-seq). MuSiC handles cross-platform deconvolution but normalization differences "
                 "may affect accuracy.\n")
    lines.append("- **No validation:** In-silico and DNA-based VerifyBamID validation deferred per user decision. "
                 "Proportion estimates should be interpreted with caution until validated.\n")
    lines.append("- **Healthy reference, diseased bulk:** Reference is from healthy term placentas; bulk includes "
                 "GDM/obesity/exposure samples. MuSiC does not correct for condition-specific DE genes (unlike MuSiC2).\n")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Pool per-cohort deconvolution outputs and generate QC report"
    )
    parser.add_argument('--input-dir', type=str, required=True,
                        help='Directory containing per-cohort output TSVs')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='Directory for pooled outputs and QC report')
    parser.add_argument('--ancestry-map', type=str, default=None,
                        help='Optional TSV: sample_id, ancestry')
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== Pooling deconvolution outputs ===")
    print(f"Input dir:  {input_dir}")
    print(f"Output dir: {output_dir}\n")

    # --- Find cohort files ---
    cohorts, full_files, collapsed_files, flag_files = find_cohort_files(input_dir)
    print(f"Cohorts found: {len(cohorts)}")
    for c in cohorts:
        print(f"  - {c}")
    print()

    if not cohorts:
        print("ERROR: No per-cohort output files found.")
        return

    # --- Pool outputs ---
    print("Pooling full proportions...")
    full_pooled, full_names = pool_files(full_files, "full")
    full_path = output_dir / "all_cohorts_cell_proportions_full.tsv"
    full_pooled.to_csv(full_path, sep='\t', index=False)
    print(f"  {full_pooled.shape[0]} samples x {full_pooled.shape[1]-2} cell types -> {full_path}")

    print("Pooling collapsed proportions...")
    collapsed_pooled, collapsed_names = pool_files(collapsed_files, "collapsed")
    collapsed_path = output_dir / "all_cohorts_cell_proportions_collapsed.tsv"
    collapsed_pooled.to_csv(collapsed_path, sep='\t', index=False)
    print(f"  {collapsed_pooled.shape[0]} samples x {collapsed_pooled.shape[1]-2} cell types -> {collapsed_path}")

    print("Pooling maternal flags...")
    flags_pooled, flag_names = pool_files(flag_files, "flag")
    flags_path = output_dir / "all_cohorts_maternal_flag.tsv"
    flags_pooled.to_csv(flags_path, sep='\t', index=False)
    print(f"  {flags_pooled.shape[0]} samples -> {flags_path}")

    # --- Load ancestry map if provided ---
    ancestry_map = None
    if args.ancestry_map and os.path.exists(args.ancestry_map):
        ancestry_map = pd.read_csv(args.ancestry_map, sep='\t')
        print(f"\nLoaded ancestry map: {len(ancestry_map)} samples")

    # --- Generate QC report ---
    print("\nGenerating QC report...")
    report = generate_qc_report(full_pooled, collapsed_pooled, flags_pooled, cohorts, ancestry_map)
    report_path = output_dir / "deconvolution_qc_report.md"
    with open(report_path, 'w') as f:
        f.write(report)
    print(f"  -> {report_path}")

    print(f"\n=== Pooling complete ===")
    print(f"Total samples: {len(full_pooled)}")
    if 'flagged' in flags_pooled.columns:
        print(f"Flagged samples: {flags_pooled['flagged'].sum()}")


if __name__ == '__main__':
    main()
