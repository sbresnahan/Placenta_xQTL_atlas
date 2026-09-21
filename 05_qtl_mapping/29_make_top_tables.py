#!/usr/bin/env python3
"""
29_make_top_tables.py — Regenerate *_cisqtl_top.tsv from existing parquets

The top tables written by earlier versions of 27_run_tensorqtl.py were
unsorted head(1000) subsets with the phenotype/group IDs dropped (a 'pval'
column check that never matches tensorQTL's pval_nominal). The parquet
files are complete and correct, so this script rebuilds the top tables
from them directly — no tensorQTL reruns needed.

For each {ANC}_{modality}_cisqtl.parquet in the results directory:
  - restores the index (phenotype/group IDs) as a column
  - sorts by pval_nominal
  - writes {ANC}_{modality}_cisqtl_top.tsv
  - prints association counts at common thresholds

Usage:
  python3 29_make_top_tables.py --results-dir /path/to/qtl_results
"""

import argparse
import glob
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Regenerate sorted cis-xQTL top tables from parquets")
    parser.add_argument("--results-dir", required=True,
                        help="Directory containing *_cisqtl.parquet files")
    args = parser.parse_args()

    parquets = sorted(glob.glob(os.path.join(args.results_dir, "*_cisqtl.parquet")))
    if not parquets:
        raise SystemExit(f"No *_cisqtl.parquet files found in {args.results_dir}")

    for pq in parquets:
        base = pq.replace("_cisqtl.parquet", "")
        df = pd.read_parquet(pq)

        # Restore IDs from the index (map_cis puts phenotype/group IDs there)
        if not isinstance(df.index, pd.RangeIndex):
            df = df.reset_index()
            if 'index' in df.columns and 'phenotype_id' not in df.columns:
                df = df.rename(columns={'index': 'phenotype_id'})

        pcol = next((c for c in ['pval_nominal', 'pval'] if c in df.columns), None)
        if pcol is not None:
            df = df.sort_values(pcol)

        out_path = f"{base}_cisqtl_top.tsv"
        df.to_csv(out_path, sep='\t', index=False)

        name = os.path.basename(base)
        msg = f"  {name}: {len(df)} associations"
        if pcol is not None:
            msg += (f" | {pcol} < 5e-8: {(df[pcol] < 5e-8).sum()}"
                    f" | < 1e-5: {(df[pcol] < 1e-5).sum()}")
        if 'group_size' in df.columns:
            msg += f" | multi-phenotype groups: {(df['group_size'] > 1).sum()}"
        print(msg + f"\n    Written: {out_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
