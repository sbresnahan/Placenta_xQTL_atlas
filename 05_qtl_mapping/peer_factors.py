#!/usr/bin/env python3
"""
peer_factors.py — Extract latent factors from ComBat+INT expression using PEER

Replaces the HCP step in combat_normalize_hcp.R when Rhcpp fails (singular
matrices, zero-variance factors). PEER is the standard latent factor method
for eQTL studies (used by GTEx). Falls back to SVD if the peer package is
not installed.

Reads the ComBat+INT expression BED written by combat_normalize_hcp.R,
extracts K latent factors, and writes them in the same tensorQTL covariate
format as the HCP script ({ANC}_hcp_factors.tsv).

Usage:
  python3 peer_factors.py \
      --expression EUR_combat_int_expression.bed \
      --output-dir hcp_factors/ \
      --ancestry EUR \
      --k 15
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']


def read_expression_bed(bed_path):
    """Read expression BED → (samples, genes) numpy matrix + sample IDs."""
    with open(bed_path, 'r') as f:
        header = f.readline().strip().split('\t')
    meta = {'#chr', 'start', 'end', 'phenotype_id', 'chr'}
    sample_ids = [c for c in header if c not in meta]
    df = pd.read_csv(bed_path, sep='\t', usecols=sample_ids)
    # samples x genes
    mat = df.values.T.astype(float)
    return mat, sample_ids


def run_peer(mat, k):
    """Run PEER factor extraction. Returns (factors, method_name)."""
    try:
        import peer
    except ImportError:
        return None, 'peer_not_installed'

    model = peer.PEER()
    model.setNk(k)
    # PEER expects samples x genes
    model.setPheno(mat)
    model.update()
    factors = model.getX()  # samples x k
    return factors, 'peer'


def run_svd(mat, k):
    """SVD-based factor extraction (PCA). Returns factors (samples x k).

    This is a deterministic fallback when PEER is not available. The factors
    are the top-k left singular vectors scaled by singular values.
    """
    # Center each gene (column) to mean 0
    mat_centered = mat - mat.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(mat_centered, full_matrices=False)
    k = min(k, U.shape[1])
    factors = U[:, :k] * S[:k]
    return factors, 'svd'


def main():
    parser = argparse.ArgumentParser(description="PEER latent factor extraction")
    parser.add_argument("--expression", required=True,
                        help="ComBat+INT expression BED file")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory (same as HCP output)")
    parser.add_argument("--ancestry", required=True,
                        help="Ancestry label (e.g. EAS, EUR)")
    parser.add_argument("--k", type=int, default=15,
                        help="Number of latent factors (default: 15)")
    args = parser.parse_args()

    print(f"\n{'='*60}")
    print(f"PEER factor extraction: {args.ancestry}")
    print(f"{'='*60}")

    # ---- Load expression ----
    print(f"  Loading expression: {args.expression}")
    mat, sample_ids = read_expression_bed(args.expression)
    print(f"  Expression: {mat.shape[0]} samples x {mat.shape[1]} genes")

    # ---- Run PEER (or SVD fallback) ----
    factors, method = run_peer(mat, args.k)
    if factors is None:
        print(f"  PEER package not available, falling back to SVD (PCA).")
        print(f"  To use PEER: pip install peer (requires Cython)")
        factors, method = run_svd(mat, args.k)

    k_actual = factors.shape[1]
    print(f"  Method: {method}")
    print(f"  Factors: {factors.shape[0]} samples x {k_actual} factors")

    # Check for zero-variance factors
    factor_vars = factors.var(axis=0)
    n_zero = int((factor_vars < 1e-10).sum())
    if n_zero > 0:
        print(f"  WARN: {n_zero} of {k_actual} factors have ~zero variance")
    else:
        print(f"  All factors have non-zero variance (OK)")

    # ---- Write factors as tensorQTL covariate table ----
    # Format: rows = factors, columns = samples, first col = covariate name
    factor_names = [f"HCP_{i+1}" for i in range(k_actual)]
    factors_df = pd.DataFrame(factors, index=sample_ids, columns=factor_names).T
    factors_df.index.name = 'covariate'

    out_path = os.path.join(args.output_dir, f"{args.ancestry}_hcp_factors.tsv")
    factors_df.to_csv(out_path, sep='\t')
    print(f"  Written: {out_path} ({k_actual} factors x {len(sample_ids)} samples)")

    # ---- Scree plot ----
    var_explained = factor_vars / factor_vars.sum() * 100
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(range(1, k_actual + 1), var_explained, color='#0279EE')
    ax.set_xlabel('Factor')
    ax.set_ylabel('Variance explained (%)')
    ax.set_title(f'{args.ancestry}: {method.upper()} factor scree')
    plt.tight_layout()
    scree_path = os.path.join(args.output_dir, f"{args.ancestry}_peer_scree.png")
    plt.savefig(scree_path, dpi=150)
    plt.close()
    print(f"  Written: {scree_path}")

    print(f"\n  Done: {args.ancestry} — {k_actual} factors via {method}")


if __name__ == "__main__":
    main()
