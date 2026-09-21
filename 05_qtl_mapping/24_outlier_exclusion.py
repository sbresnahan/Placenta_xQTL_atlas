#!/usr/bin/env python3
"""
24_outlier_exclusion.py — Select genotype PCs and remove PC outliers

REVISED (GTEx-conventions round): PC selection is now simply the FIRST 5
genotype PCs (PC1-PC5), matching GTEx/PANTRY convention. The previous
within-stratum variance selection (which picked high-numbered PCs from a
1KG-projected PC-AiR solution) is removed — PCs now come from cohort-only
PCA (22_genotype_pca.sh), so eigenvalue order is meaningful.

For each ancestry:
  1. Loads genotype PCs (from 22_genotype_pca.sh).
  2. Selects PC1..PC5 (--n-pcs, default 5).
  3. Writes {ANC}_selected_pcs.txt (consumed by 25_build_covariates.py) and
     {ANC}_pc_variance.tsv (full variance table, diagnostic only).
  4. Flags samples where |PC_i - mean| > --sd-threshold * SD for any
     SELECTED PC.
  5. Writes outlier list {ANC}_outliers.tsv.
  6. Removes outliers from: genotype pgen, expression BED, HCP factors,
     deconvolution proportions, metadata. The pgen rewrite also applies
     --mac (minor-allele-count floor, default 5) so low-carrier variants
     (the round-1 artifact class) never enter testing.

Usage:
  python3 24_outlier_exclusion.py \
      --qtl-dir <qtl_inputs dir> \
      --pcair-dir <genotype_pcs dir> \
      --ancestries EAS EUR \
      --sd-threshold 6 \
      --n-pcs 5 \
      --mac 5
"""

import argparse
import os
import subprocess
import sys
import pandas as pd
import numpy as np


def filter_bed_samples(bed_path, out_path, keep_samples):
    """Keep only columns in keep_samples (plus BED meta columns)."""
    with open(bed_path, 'r') as f:
        header = f.readline()
        data = f.read()

    cols = header.strip().split('\t')
    meta = {'#chr', 'start', 'end', 'phenotype_id', 'chr'}
    keep_indices = []
    for i, c in enumerate(cols):
        if c in meta or c in keep_samples:
            keep_indices.append(i)

    new_header = '\t'.join(cols[i] for i in keep_indices) + '\n'
    new_data_lines = []
    for line in data.strip().split('\n'):
        fields = line.split('\t')
        new_data_lines.append('\t'.join(fields[i] for i in keep_indices))
    new_data = '\n'.join(new_data_lines) + '\n'

    with open(out_path, 'w') as f:
        f.write(new_header)
        f.write(new_data)

    dropped = len(cols) - len(keep_indices)
    # Dropped columns are all samples (meta columns are always kept)
    meta_count = sum(1 for c in cols if c in meta)
    return len(keep_indices) - meta_count, dropped


def select_first_pcs(pc_cols_all, n_pcs):
    """GTEx convention: select the first n_pcs genotype PCs by eigenvalue
    order (PC1..PCn). pc_cols_all must be sorted by PC number."""
    return pc_cols_all[:n_pcs]


def pc_variance_table(pcs, pc_cols_all, selected):
    """Diagnostic variance table across all computed PCs (not used for
    selection — retained for the scree-style diagnostic only)."""
    variances = pcs[pc_cols_all].var()
    total_var = variances.sum()
    pct_var = 100.0 * variances / total_var
    return pd.DataFrame({
        'pc': pc_cols_all,
        'variance': variances[pc_cols_all].values,
        'pct_variance': pct_var[pc_cols_all].values,
        'selected': [c in selected for c in pc_cols_all],
    })


def main():
    parser = argparse.ArgumentParser(description="PC selection + outlier exclusion on genotype PCs")
    parser.add_argument("--qtl-dir", required=True, help="QTL inputs directory (from Step 3)")
    parser.add_argument("--pcair-dir", required=True, help="Genotype PCs directory (from Step 2)")
    parser.add_argument("--ancestries", default="EAS EUR", help="Space-separated ancestry labels")
    parser.add_argument("--sd-threshold", type=float, default=6.0,
                        help="SD threshold for outlier flagging (default: 6)")
    parser.add_argument("--n-pcs", type=int, default=5,
                        help="Number of top genotype PCs to select (GTEx convention: 5)")
    parser.add_argument("--mac", type=int, default=5,
                        help="Minor-allele-count floor applied to the post-outlier pgen "
                             "(default: 5; 0 disables). Prevents low-carrier artifacts.")
    args = parser.parse_args()

    ancestries = args.ancestries.split()

    for anc in ancestries:
        print(f"\n{'='*60}")
        print(f"Ancestry: {anc}")
        print(f"{'='*60}")

        # ---- Load genotype PCs ----
        pcs_path = os.path.join(args.pcair_dir, f"{anc}_genotype_pcs.tsv")
        if not os.path.exists(pcs_path):
            print(f"  ERROR: genotype PCs not found: {pcs_path}")
            continue
        pcs = pd.read_csv(pcs_path, sep='\t')
        print(f"  Loaded {len(pcs)} samples, {len(pcs.columns)-1} PCs")

        # ---- Sanity-check PC sample IDs against the qtl psam ----
        qtl_psam_path = os.path.join(args.qtl_dir, f"{anc}_qtl.psam")
        if not os.path.exists(qtl_psam_path):
            print(f"  ERROR: qtl psam not found: {qtl_psam_path} (run Step 3 first)")
            continue
        with open(qtl_psam_path) as f:
            psam_header = f.readline().strip()
        psam = pd.read_csv(qtl_psam_path, sep='\t', comment='#', header=None)
        if psam_header.startswith('#'):
            cols = psam_header.lstrip('#').split('\t')
            if len(cols) == psam.shape[1]:
                psam.columns = cols
        if 'IID' not in psam.columns:
            psam.columns = ['FID', 'IID'] + [f'extra_{i}' for i in range(psam.shape[1] - 2)]
        if 'FID' not in psam.columns:
            # plink2 uses '0' as FID when the .psam has no FID column
            psam['FID'] = '0'

        pc_ids = set(pcs['sample_id'])
        psam_iids = set(psam['IID'])
        overlap = pc_ids & psam_iids
        print(f"  PC sample IDs vs psam IIDs: {len(overlap)} of {len(pc_ids)} PC samples match")
        if len(overlap) == 0:
            print(f"  ERROR: no overlap between PC sample IDs and psam IIDs.")
            print(f"    PC examples:    {sorted(pc_ids)[:3]}")
            print(f"    psam examples:  {sorted(psam_iids)[:3]}")
            print(f"    ID formats differ — report these to adjust the harmonization.")
            continue
        # Restrict PCs to samples actually in the qtl pgen
        pcs = pcs[pcs['sample_id'].isin(psam_iids)].reset_index(drop=True)

        # ---- Select the first N PCs (GTEx convention) ----
        pc_cols_all = sorted([c for c in pcs.columns if c.startswith('PC')],
                             key=lambda c: int(c[2:]))
        selected = select_first_pcs(pc_cols_all, args.n_pcs)
        var_table = pc_variance_table(pcs, pc_cols_all, selected)

        print(f"\n  Within-stratum PC variance (n={len(pcs)} samples; diagnostic only):")
        for _, row in var_table.iterrows():
            mark = ' *' if row['selected'] else ''
            print(f"    {row['pc']:>5}: {row['pct_variance']:6.2f}%{mark}")
        print(f"\n  Selected first {len(selected)} PCs (GTEx convention): "
              f"{', '.join(selected)}")

        sel_path = os.path.join(args.qtl_dir, f"{anc}_selected_pcs.txt")
        with open(sel_path, 'w') as f:
            f.write('\n'.join(selected) + '\n')
        print(f"  Written: {sel_path}")

        var_path = os.path.join(args.qtl_dir, f"{anc}_pc_variance.tsv")
        var_table.to_csv(var_path, sep='\t', index=False)
        print(f"  Written: {var_path}")

        # ---- Flag outliers on SELECTED PCs ----
        outlier_mask = pd.Series(False, index=pcs.index)
        outlier_records = []

        for pc in selected:
            mean = pcs[pc].mean()
            sd = pcs[pc].std()
            if sd == 0:
                print(f"  {pc}: zero variance within stratum, skipping outlier check")
                continue
            sd_dist = (pcs[pc] - mean).abs() / sd
            is_outlier = sd_dist > args.sd_threshold
            n_outliers = is_outlier.sum()
            if n_outliers > 0:
                print(f"  {pc}: mean={mean:.4f}, sd={sd:.4f}, outliers (> {args.sd_threshold} SD): {n_outliers}")
                for idx in pcs.index[is_outlier]:
                    outlier_records.append({
                        'sample_id': pcs.loc[idx, 'sample_id'],
                        'pc': pc,
                        'value': pcs.loc[idx, pc],
                        'mean': mean,
                        'sd': sd,
                        'sd_distance': sd_dist[idx]
                    })
            outlier_mask = outlier_mask | is_outlier

        outlier_samples = set(pcs.loc[outlier_mask, 'sample_id'])
        clean_samples = set(pcs.loc[~outlier_mask, 'sample_id'])

        print(f"\n  Total outliers: {len(outlier_samples)}")
        print(f"  Clean samples: {len(clean_samples)}")

        # ---- Write outlier list ----
        outlier_path = os.path.join(args.qtl_dir, f"{anc}_outliers.tsv")
        if outlier_records:
            outlier_df = pd.DataFrame(outlier_records)
            outlier_df.to_csv(outlier_path, sep='\t', index=False)
            print(f"  Written: {outlier_path}")
        else:
            print(f"  No outliers found. Writing empty file.")
            pd.DataFrame(columns=['sample_id', 'pc', 'value', 'mean', 'sd', 'sd_distance']).to_csv(
                outlier_path, sep='\t', index=False)

        if not outlier_samples:
            print(f"  No outliers to remove. All files unchanged."
                  + (f" (pgen retains the MAC >= {args.mac} floor from 23_prepare_intersection.py)"
                     if args.mac > 0 else ""))
            continue

        # ---- Remove outliers from genotype pgen ----
        keep_psam = os.path.join(args.qtl_dir, f"{anc}_keep_clean.txt")
        # Write actual FID/IID pairs from the qtl psam (plink2 --keep matches both)
        keep_df = psam[psam['IID'].isin(clean_samples)][['FID', 'IID']]
        keep_df.to_csv(keep_psam, sep='\t', index=False, header=False)

        qtl_prefix = os.path.join(args.qtl_dir, f"{anc}_qtl")
        tmp_prefix = os.path.join(args.qtl_dir, f"{anc}_qtl_clean")
        mac_flag = f"--mac {args.mac} " if args.mac > 0 else ""
        cmd = (
            f"plink2 --pgen {qtl_prefix}.pgen "
            f"--pvar {qtl_prefix}.pvar "
            f"--psam {qtl_prefix}.psam "
            f"--keep {keep_psam} "
            f"{mac_flag}"
            f"--make-pgen --threads 4 "
            f"--out {tmp_prefix}"
        )
        print(f"\n  Filtering genotype pgen...")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ERROR: plink2 failed: {result.stderr}")
            continue
        # Replace original with clean version
        for ext in ['.pgen', '.pvar', '.psam']:
            os.replace(f"{tmp_prefix}{ext}", f"{qtl_prefix}{ext}")
        print(f"    Genotype: {len(clean_samples)} samples (removed {len(outlier_samples)})"
              + (f"; MAC >= {args.mac} floor re-applied" if args.mac > 0 else ""))

        # ---- Remove outliers from expression BED ----
        expr_path = os.path.join(args.qtl_dir, f"{anc}_expression_harmonized.bed")
        if os.path.exists(expr_path):
            kept, dropped = filter_bed_samples(expr_path, expr_path, clean_samples)
            print(f"    Expression: {kept} samples (removed {dropped})")

        # ---- Remove outliers from HCP factors ----
        hcp_path = os.path.join(args.qtl_dir, f"{anc}_hcp_factors_harmonized.tsv")
        if os.path.exists(hcp_path):
            hcp_df = pd.read_csv(hcp_path, sep='\t', index_col=0)
            drop_cols = [c for c in hcp_df.columns if c not in clean_samples]
            hcp_df = hcp_df.drop(columns=drop_cols)
            hcp_df.to_csv(hcp_path, sep='\t')
            print(f"    HCP factors: {len(hcp_df.columns)} samples (removed {len(drop_cols)})")

        # ---- Remove outliers from deconvolution ----
        deconv_path = os.path.join(args.qtl_dir, f"{anc}_deconvolution_harmonized.tsv")
        if os.path.exists(deconv_path):
            deconv_df = pd.read_csv(deconv_path, sep='\t')
            deconv_df = deconv_df[deconv_df['sample_id'].isin(clean_samples)]
            deconv_df.to_csv(deconv_path, sep='\t', index=False)
            print(f"    Deconvolution: {len(deconv_df)} samples")

        # ---- Remove outliers from metadata ----
        meta_path = os.path.join(args.qtl_dir, f"{anc}_metadata.tsv")
        if os.path.exists(meta_path):
            meta_df = pd.read_csv(meta_path, sep='\t')
            meta_df = meta_df[meta_df['array_id'].isin(clean_samples)]
            meta_df.to_csv(meta_path, sep='\t', index=False)
            print(f"    Metadata: {len(meta_df)} samples")

        print(f"\n  Done: {anc} — {len(clean_samples)} samples after outlier exclusion")


if __name__ == "__main__":
    main()
