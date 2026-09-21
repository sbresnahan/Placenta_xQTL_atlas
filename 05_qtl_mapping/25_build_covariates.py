#!/usr/bin/env python3
"""
25_build_covariates.py — Assemble and optimize tensorQTL covariate table

REVISED (GTEx-conventions round): ppBMI removed per study decision; new
hard cap of --max-covariates (default 25) applied after correlation pruning.

Combines:
  - HCP factors (HCP_1 … HCP_15) from harmonized HCP file
  - Genotype PCs selected by 24_outlier_exclusion.py ({ANC}_selected_pcs.txt;
    now the first 5 PCs, GTEx convention)
  - Sex (M→0, F→1) from metadata
  - Gestational age (mean-centered) from metadata
  - Cell-type proportions (collapsed types, arcsinh-transformed,
    mean-centered, DOMINANT type dropped as compositional reference)

Optimization (in order):
  1. Pre-filter: drop near-zero-variance covariates (sd < --min-sd).
     Catches degenerate HCP factors from singular HCP runs.
  2. Correlation pruning: iteratively drop the lower-priority member of the
     max-|r| pair while any |r| > --cor-threshold. Keep priority:
     sex/GA > genotype PCs > cell types > HCP factors; within a block,
     lower index kept.
  3. Cap: while more than --max-covariates remain, drop the lowest-priority
     covariate (highest tier, then highest within-block index — i.e.
     HCP_15 first). sex/GA and genotype PCs are effectively never
     cap-dropped at the default cap of 25.

Diagnostics per ancestry:
  - {ANC}_covariate_correlation.png — before/after correlation heatmaps
  - {ANC}_covariate_pruning.tsv — every dropped covariate with reason
  - Console: per-block counts, final count, condition number

Output: {ANC}_covariates.tsv in tensorQTL format
  (rows = covariates, columns = samples, tab-delimited, first column = name)

Usage:
  python3 25_build_covariates.py \
      --qtl-dir <qtl_inputs dir> \
      --pcair-dir <genotype_pcs dir> \
      --ancestries EAS EUR \
      --cor-threshold 0.9 \
      --min-sd 1e-8 \
      --max-covariates 25
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']

# Priority tiers for correlation pruning (lower = higher priority)
TIER_DEMOGRAPHIC = 0   # sex, GA
TIER_GENO_PC = 1       # genotype PCs
TIER_CELLTYPE = 2      # cell-type proportions
TIER_HCP = 3           # HCP factors


def plot_cor_heatmaps(cor_before, cor_after, out_path, anc):
    """Side-by-side correlation heatmaps before/after pruning."""
    n_before = cor_before.shape[0]
    n_after = cor_after.shape[0]
    panel_w = max(5.5, 0.28 * n_before + 2)
    fig, axes = plt.subplots(1, 2, figsize=(2 * panel_w, panel_w))
    im = None
    for ax, cor, title in zip(
            axes, [cor_before, cor_after],
            [f'Before pruning ({n_before} covariates)',
             f'After pruning ({n_after} covariates)']):
        im = ax.imshow(cor.values, cmap='RdBu_r', vmin=-1, vmax=1,
                       aspect='auto')
        ax.set_xticks(range(cor.shape[1]))
        ax.set_xticklabels(cor.columns, rotation=90, fontsize=6)
        ax.set_yticks(range(cor.shape[0]))
        ax.set_yticklabels(cor.index, fontsize=6)
        ax.set_title(f'{anc}: {title}', fontsize=10)
    # Dedicated colorbar axis to the right of both panels (no overlap)
    fig.subplots_adjust(left=0.08, right=0.90, wspace=0.35,
                        top=0.92, bottom=0.18)
    cax = fig.add_axes([0.92, 0.30, 0.015, 0.45])
    fig.colorbar(im, cax=cax, label='Pearson r')
    plt.savefig(out_path, dpi=150)
    plt.close()


def prune_correlated(covariates, priority, threshold):
    """Iteratively drop the lower-priority member of the max-|r| pair until
    no |r| > threshold remains. covariates: rows = covariates, cols = samples.
    Returns (pruned_df, dropped_records)."""
    dropped = []
    current = covariates.copy()
    while current.shape[0] > 1:
        cor = current.T.corr().abs()
        cor_vals = np.array(cor.values, copy=True)  # pandas .values may be read-only
        np.fill_diagonal(cor_vals, 0.0)
        max_r = cor_vals.max()
        if not np.isfinite(max_r) or max_r <= threshold:
            break
        i, j = np.unravel_index(cor_vals.argmax(), cor_vals.shape)
        c1, c2 = cor.index[i], cor.index[j]
        r_val = cor.loc[c1, c2]
        # Drop lower priority: higher tier; within tier, higher index
        if priority[c1] <= priority[c2]:
            drop, keep = c2, c1
        else:
            drop, keep = c1, c2
        dropped.append({
            'covariate': drop,
            'reason': 'correlated',
            'correlated_with': keep,
            'r_value': round(float(r_val), 4),
        })
        current = current.drop(index=drop)
    return current, dropped


def condition_number(covariates):
    """Condition number of the design matrix (intercept + z-scored
    covariates). Scale-free measure of collinearity."""
    X = covariates.T.values.astype(float)  # samples x covariates
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    Xs = (X - X.mean(axis=0)) / sd
    Xd = np.column_stack([np.ones(Xs.shape[0]), Xs])
    return np.linalg.cond(Xd)


def main():
    parser = argparse.ArgumentParser(description="Build + optimize tensorQTL covariate table")
    parser.add_argument("--qtl-dir", required=True, help="QTL inputs directory (from Steps 3-4)")
    parser.add_argument("--pcair-dir", required=True, help="Genotype PCs directory (from Step 2)")
    parser.add_argument("--ancestries", default="EAS EUR", help="Space-separated ancestry labels")
    parser.add_argument("--cor-threshold", type=float, default=0.9,
                        help="Prune covariate pairs with |r| above this (default: 0.9)")
    parser.add_argument("--min-sd", type=float, default=1e-8,
                        help="Drop covariates with sd below this (default: 1e-8)")
    parser.add_argument("--max-covariates", type=int, default=25,
                        help="Hard cap on final covariate count, applied after "
                             "correlation pruning (default: 25). Lowest-priority "
                             "covariates (highest-index HCPs first) are dropped.")
    args = parser.parse_args()

    ancestries = args.ancestries.split()

    for anc in ancestries:
        print(f"\n{'='*60}")
        print(f"Ancestry: {anc}")
        print(f"{'='*60}")

        covariate_blocks = []
        priority = {}  # covariate name -> (tier, within-block index)
        dropped_records = []

        # ---- 1. HCP factors ----
        hcp_path = os.path.join(args.qtl_dir, f"{anc}_hcp_factors_harmonized.tsv")
        if not os.path.exists(hcp_path):
            print(f"  ERROR: HCP factors not found: {hcp_path}")
            continue
        hcp_df = pd.read_csv(hcp_path, sep='\t', index_col=0)
        # HCP factors: rows = HCP_1..15, columns = samples
        print(f"  HCP factors: {hcp_df.shape[0]} factors, {hcp_df.shape[1]} samples")
        for name in hcp_df.index:
            try:
                idx = int(str(name).split('_')[1])
            except (IndexError, ValueError):
                idx = 999
            priority[name] = (TIER_HCP, idx)
        covariate_blocks.append(hcp_df)

        # ---- 2. Genotype PCs (selected by script 24) ----
        sel_path = os.path.join(args.qtl_dir, f"{anc}_selected_pcs.txt")
        if not os.path.exists(sel_path):
            print(f"  ERROR: selected PC list not found: {sel_path}")
            print(f"         Run 24_outlier_exclusion.py first (it writes this file).")
            continue
        with open(sel_path) as f:
            selected_pcs = [line.strip() for line in f if line.strip()]
        pcs_path = os.path.join(args.pcair_dir, f"{anc}_genotype_pcs.tsv")
        if not os.path.exists(pcs_path):
            print(f"  ERROR: genotype PCs not found: {pcs_path}")
            continue
        pcs_df = pd.read_csv(pcs_path, sep='\t')
        # Filter to samples in HCP (post-outlier-exclusion)
        hcp_samples = set(hcp_df.columns)
        pcs_df = pcs_df[pcs_df['sample_id'].isin(hcp_samples)]
        pcs_df = pcs_df.set_index('sample_id')
        pc_cols = [c for c in selected_pcs if c in pcs_df.columns]
        missing_pcs = [c for c in selected_pcs if c not in pcs_df.columns]
        if missing_pcs:
            print(f"  WARN: selected PCs missing from PC file: {missing_pcs}")
        pcs_df = pcs_df[pc_cols].T  # rows = PCs, columns = samples
        print(f"  Genotype PCs: {len(pc_cols)} selected PCs "
              f"({', '.join(pc_cols)}), {pcs_df.shape[1]} samples")
        for name in pc_cols:
            priority[name] = (TIER_GENO_PC, int(name[2:]))
        covariate_blocks.append(pcs_df)

        # ---- 3. Sex and GA from metadata ----
        meta_path = os.path.join(args.qtl_dir, f"{anc}_metadata.tsv")
        if not os.path.exists(meta_path):
            print(f"  ERROR: metadata not found: {meta_path}")
            continue
        meta_df = pd.read_csv(meta_path, sep='\t')
        meta_df = meta_df[meta_df['array_id'].isin(hcp_samples)]
        meta_df = meta_df.set_index('array_id')

        # Sex: M→0, F→1
        sex_map = {'M': 0, 'F': 1, 'm': 0, 'f': 1, 'Male': 0, 'Female': 1}
        sex_row = meta_df['sex'].map(sex_map)
        sex_df = pd.DataFrame({'sex': sex_row}).T
        print(f"  Sex: {int(sex_row.sum())} F, {int((sex_row == 0).sum())} M, "
              f"{int(sex_row.isna().sum())} unknown")
        if sex_row.isna().any():
            mode = sex_row.mode()[0]
            sex_df = sex_df.fillna(mode)
            print(f"    Imputed {int(sex_row.isna().sum())} missing with mode ({mode})")
        priority['sex'] = (TIER_DEMOGRAPHIC, 0)
        covariate_blocks.append(sex_df)

        # GA: mean-centered
        if 'GA' in meta_df.columns:
            ga = meta_df['GA'].astype(float)
            ga_centered = ga - ga.mean()
            ga_df = pd.DataFrame({'GA': ga_centered}).T
            print(f"  GA: mean={ga.mean():.2f} weeks, range=[{ga.min():.1f}, {ga.max():.1f}]")
            if ga_centered.isna().any():
                ga_df = ga_df.fillna(0)
                print(f"    Imputed {int(ga_centered.isna().sum())} missing with mean (0 after centering)")
            priority['GA'] = (TIER_DEMOGRAPHIC, 1)
            covariate_blocks.append(ga_df)
        else:
            print(f"  WARN: GA column not found in metadata, skipping")

        # ppBMI: REMOVED per study decision (GTEx-conventions round) — not
        # included in the covariate pool.

        # ---- 4. Cell-type proportions (collapsed types) ----
        deconv_path = os.path.join(args.qtl_dir, f"{anc}_deconvolution_harmonized.tsv")
        if os.path.exists(deconv_path):
            deconv_df = pd.read_csv(deconv_path, sep='\t')
            deconv_df = deconv_df[deconv_df['sample_id'].isin(hcp_samples)]
            deconv_df = deconv_df.set_index('sample_id')
            prop_cols = [c for c in deconv_df.columns if c not in ('cohort',)]
            deconv_df = deconv_df[prop_cols]
            # Drop the dominant (most abundant) cell type as the
            # compositional reference: proportions sum to ~1, so keeping all
            # types makes the block rank-deficient with the intercept.
            dominant = deconv_df.mean().idxmax()
            dropped_records.append({
                'covariate': dominant,
                'reason': 'reference_cell_type',
                'correlated_with': '',
                'r_value': '',
            })
            print(f"  Dropping dominant cell type as reference: {dominant} "
                  f"(mean proportion {deconv_df.mean().max():.3f})")
            deconv_df = deconv_df.drop(columns=[dominant])
            # Arcsinh transform + mean-center
            deconv_transformed = np.arcsinh(deconv_df)
            deconv_transformed = deconv_transformed - deconv_transformed.mean()
            deconv_df_t = deconv_transformed.T
            # Prefix cell type names with "ct_" to avoid collisions with
            # HCP/PC/demographic covariate names
            deconv_df_t.index = ['ct_' + str(c) for c in deconv_df_t.index]
            print(f"  Cell-type proportions: {deconv_df_t.shape[0]} types "
                  f"(after dropping reference), {deconv_df_t.shape[1]} samples")
            for i, name in enumerate(deconv_df_t.index):
                priority[name] = (TIER_CELLTYPE, i)
            covariate_blocks.append(deconv_df_t)
        else:
            print(f"  WARN: deconvolution not found: {deconv_path}, skipping cell-type covariates")

        # ---- Combine all covariate blocks ----
        all_samples = set(hcp_df.columns)
        for block in covariate_blocks[1:]:
            all_samples = all_samples & set(block.columns)

        print(f"\n  Final sample count (intersection of all covariate blocks): {len(all_samples)}")

        aligned_blocks = []
        for bi, block in enumerate(covariate_blocks):
            # Deduplicate columns (sample IDs) — keep first occurrence
            if block.columns.duplicated().any():
                dup_cols = block.columns[block.columns.duplicated(keep=False)].unique().tolist()
                print(f"  WARN: block {bi} has duplicate sample columns: {dup_cols[:5]}, keeping first")
                block = block.loc[:, ~block.columns.duplicated(keep='first')]
            # Also check for duplicate row indices (covariate names)
            if block.index.duplicated().any():
                dup_idx = block.index[block.index.duplicated(keep=False)].unique().tolist()
                print(f"  WARN: block {bi} has duplicate covariate names: {dup_idx[:5]}, keeping first")
                block = block[~block.index.duplicated(keep='first')]
            aligned_blocks.append(block[list(all_samples)])
        # Check for duplicate covariate names (row indices) before concatenating
        all_names = []
        for b in aligned_blocks:
            all_names.extend(b.index.tolist())
        dupes = [n for n in set(all_names) if all_names.count(n) > 1]
        if dupes:
            sys.exit(f"  ERROR: Duplicate covariate names across blocks: {dupes}. "
                     f"Cannot concatenate. Check for name collisions.")
        try:
            covariates = pd.concat(aligned_blocks, axis=0)
        except Exception as e:
            print(f"\n  ERROR during pd.concat: {e}")
            for bi, b in enumerate(aligned_blocks):
                print(f"    block {bi}: shape={b.shape}, "
                      f"index_dups={b.index.duplicated().sum()}, "
                      f"col_dups={b.columns.duplicated().sum()}, "
                      f"index[:3]={list(b.index[:3])}")
            sys.exit(1)

        if covariates.isna().any().any():
            n_na = int(covariates.isna().sum().sum())
            print(f"  WARN: {n_na} NaN values in covariate table, filling with 0")
            covariates = covariates.fillna(0)

        n_assembled = covariates.shape[0]
        block_counts = {
            'HCP': sum(1 for c in covariates.index if priority.get(c, (None,))[0] == TIER_HCP),
            'geno_PC': sum(1 for c in covariates.index if priority.get(c, (None,))[0] == TIER_GENO_PC),
            'demographic': sum(1 for c in covariates.index if priority.get(c, (None,))[0] == TIER_DEMOGRAPHIC),
            'cell_type': sum(1 for c in covariates.index if priority.get(c, (None,))[0] == TIER_CELLTYPE),
        }
        print(f"\n  Assembled {n_assembled} covariates: "
              + ', '.join(f"{k}={v}" for k, v in block_counts.items()))

        # ---- Optimization 1: near-zero-variance pre-filter ----
        sds = covariates.std(axis=1)
        zero_var = sds[sds < args.min_sd].index.tolist()
        for name in zero_var:
            dropped_records.append({
                'covariate': name,
                'reason': 'zero_variance',
                'correlated_with': '',
                'r_value': '',
            })
            print(f"  Pre-filter: dropping {name} (sd={sds[name]:.2e} < {args.min_sd})")
        covariates = covariates.drop(index=zero_var)

        # Correlation matrix before correlation pruning
        cor_before = covariates.T.corr()

        # ---- Optimization 2: correlation pruning ----
        covariates, cor_dropped = prune_correlated(
            covariates, priority, args.cor_threshold)
        dropped_records.extend(cor_dropped)
        for rec in cor_dropped:
            print(f"  Pruning: dropping {rec['covariate']} "
                  f"(|r|={rec['r_value']} with {rec['correlated_with']})")

        # ---- Optimization 3: hard cap at --max-covariates ----
        # Drop lowest-priority covariates first: highest tier, then highest
        # within-block index (HCP_15 before HCP_14, ...). At the default cap
        # of 25 this only ever reaches HCPs/cell types — sex/GA and genotype
        # PCs are lower tiers and are effectively protected.
        while covariates.shape[0] > args.max_covariates:
            drop = max(covariates.index,
                       key=lambda c: priority.get(c, (99, 999)))
            dropped_records.append({
                'covariate': drop,
                'reason': 'cap_25' if args.max_covariates == 25 else f'cap_{args.max_covariates}',
                'correlated_with': '',
                'r_value': '',
            })
            print(f"  Cap ({args.max_covariates}): dropping {drop} "
                  f"(lowest priority)")
            covariates = covariates.drop(index=drop)

        cor_after = covariates.T.corr()

        # ---- Diagnostics ----
        cond = condition_number(covariates)
        final_counts = {
            'HCP': sum(1 for c in covariates.index if priority.get(c, (None,))[0] == TIER_HCP),
            'geno_PC': sum(1 for c in covariates.index if priority.get(c, (None,))[0] == TIER_GENO_PC),
            'demographic': sum(1 for c in covariates.index if priority.get(c, (None,))[0] == TIER_DEMOGRAPHIC),
            'cell_type': sum(1 for c in covariates.index if priority.get(c, (None,))[0] == TIER_CELLTYPE),
        }
        print(f"\n  Final: {covariates.shape[0]} covariates "
              f"(dropped {n_assembled - covariates.shape[0]}): "
              + ', '.join(f"{k}={v}" for k, v in final_counts.items()))
        print(f"  Condition number of final design matrix: {cond:.1f}")
        if cond > 1000:
            print(f"  WARN: condition number > 1000 — residual collinearity remains")

        png_path = os.path.join(args.qtl_dir, f"{anc}_covariate_correlation.png")
        plot_cor_heatmaps(cor_before, cor_after, png_path, anc)
        print(f"  Written: {png_path}")

        pruning_path = os.path.join(args.qtl_dir, f"{anc}_covariate_pruning.tsv")
        pd.DataFrame(dropped_records,
                     columns=['covariate', 'reason', 'correlated_with', 'r_value']
                     ).to_csv(pruning_path, sep='\t', index=False)
        print(f"  Written: {pruning_path}")

        # ---- Write covariate table (tensorQTL format) ----
        out_path = os.path.join(args.qtl_dir, f"{anc}_covariates.tsv")
        covariates.index.name = 'covariate_id'
        covariates.to_csv(out_path, sep='\t')

        print(f"\n  Written: {out_path}")
        print(f"  Shape: {covariates.shape[0]} covariates x {covariates.shape[1]} samples")
        print(f"  Covariates: {list(covariates.index)}")


if __name__ == "__main__":
    main()
