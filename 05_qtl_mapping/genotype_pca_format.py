#!/usr/bin/env python3
"""
genotype_pca_format.py — Format plink2 PCA output for the QTL pipeline

Reads plink2 --pca output (.eigenvec + .eigenval) and writes:
  1. {ANC}_genotype_pcs.tsv — sample_id, PC1..PCn (the schema
     24_outlier_exclusion.py already reads; replaces the PC-AiR output).
  2. {ANC}_genotype_pcs_scree.png — scree plot (per-PC % variance +
     cumulative %) from the eigenvalues.

Usage:
  python3 genotype_pca_format.py \
      --eigenvec <work>/EAS_pca.eigenvec \
      --eigenval <work>/EAS_pca.eigenval \
      --ancestry EAS --n-pcs 20 --output-dir <genotype_pcs dir>
"""

import argparse
import os
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]


def main():
    ap = argparse.ArgumentParser(description="Format plink2 PCA output + scree plot")
    ap.add_argument("--eigenvec", required=True, help="plink2 .eigenvec file")
    ap.add_argument("--eigenval", required=True, help="plink2 .eigenval file")
    ap.add_argument("--ancestry", required=True, help="Ancestry label (e.g. EAS)")
    ap.add_argument("--n-pcs", type=int, default=20, help="Number of PCs (default: 20)")
    ap.add_argument("--output-dir", default=".", help="Output directory")
    args = ap.parse_args()

    # ---- Load eigenvectors ----
    # plink2 writes a whitespace-delimited file with header '#FID IID PC1 ...'
    ev = pd.read_csv(args.eigenvec, sep=r"\s+")
    ev.columns = [c.lstrip("#") for c in ev.columns]
    if "IID" not in ev.columns:
        sys.exit(f"ERROR: no IID column in {args.eigenvec}; columns: {list(ev.columns)}")
    pc_cols = [c for c in ev.columns if c.startswith("PC")]
    if len(pc_cols) < args.n_pcs:
        print(f"  WARN: only {len(pc_cols)} PCs in eigenvec (requested {args.n_pcs})")
    pc_cols = sorted(pc_cols, key=lambda c: int(c[2:]))[: args.n_pcs]

    out = ev[["IID"] + pc_cols].rename(columns={"IID": "sample_id"})
    os.makedirs(args.output_dir, exist_ok=True)
    tsv_path = os.path.join(args.output_dir, f"{args.ancestry}_genotype_pcs.tsv")
    out.to_csv(tsv_path, sep="\t", index=False)
    print(f"  Written: {tsv_path} ({out.shape[0]} samples x {len(pc_cols)} PCs)")

    # ---- Scree plot from eigenvalues ----
    eigenvalues = pd.read_csv(args.eigenval, header=None)[0].values
    n = min(len(eigenvalues), len(pc_cols))
    if n == 0:
        print("  WARN: no eigenvalues found, skipping scree plot")
        return
    pct = 100.0 * eigenvalues[:n] / eigenvalues.sum()
    cum = pct.cumsum()

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(range(1, n + 1), pct, color="#0279EE", alpha=0.7, label="Per-PC variance")
    ax.plot(range(1, n + 1), cum, color="#FF9400", marker="o", markersize=3,
            linewidth=0.8, label="Cumulative")
    ax.set_xticks(range(1, n + 1))
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance explained (%)")
    ax.set_title(f"Genotype PC scree plot — {args.ancestry} (cohort-only PCA)")
    ax.legend(frameon=False)
    fig.tight_layout()
    png_path = os.path.join(args.output_dir, f"{args.ancestry}_genotype_pcs_scree.png")
    fig.savefig(png_path, dpi=150)
    print(f"  Written: {png_path}")

    print("  Eigenvalue summary:")
    for i in range(min(n, 10)):
        print(f"    PC{i+1}: {eigenvalues[i]:.4f} ({pct[i]:.2f}%, cum {cum[i]:.2f}%)")


if __name__ == "__main__":
    main()
