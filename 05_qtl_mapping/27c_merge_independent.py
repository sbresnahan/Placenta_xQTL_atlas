#!/usr/bin/env python3
"""
27c_merge_independent.py — merge chunked cis.map_independent outputs into the
final per-ancestry × modality independent-xQTL tables.

Reads {output_dir}/independent_chunks/{ANC}_{label}/chunk_KKKK.parquet (written
by 27_run_tensorqtl.py --independent-only --chunk-index, submitted as an LSF
array by 28b_submit_independent.sh), verifies completeness against the number
of FDR-significant rows in {ANC}_{label}_cisqtl.parquet, and writes the same
final outputs the monolithic --independent mode produces:

  {ANC}_{label}_cisqtl_independent.parquet
  {ANC}_{label}_cisqtl_independent_top.tsv

Each chunk carries a chunk_KKKK.json sidecar recording the chunk_size and
independent_fdr it was computed with; the merge refuses to mix chunks whose
parameters differ from its own (stale chunks from a run with different
chunking parameters would otherwise be silently misassembled).

Usage:
  python3 27c_merge_independent.py \
      --output-dir <qtl_results dir> \
      --ancestry EAS \
      --modality splicing \
      [--no-groups] [--chunk-size 100] [--independent-fdr 0.05]

Exit 1 (printing the missing chunk indices and a resubmission command) if any
expected chunk is missing or parameter-mismatched. With zero significant
phenotypes/groups, empty outputs are written (exit 0), mirroring the
monolithic --independent behavior.
"""

import argparse
import glob
import json
import math
import os
import sys

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Merge chunked cis.map_independent outputs")
    parser.add_argument("--output-dir", required=True,
                        help="qtl_results directory")
    parser.add_argument("--ancestry", required=True, help="Ancestry label")
    parser.add_argument("--modality", required=True, help="Modality label")
    parser.add_argument("--no-groups", action="store_true",
                        help="Merge the ungrouped ({mod}_ungrouped) outputs")
    parser.add_argument("--chunk-size", type=int, default=100,
                        help="Chunk size the array was submitted with "
                             "(default: 100)")
    parser.add_argument("--independent-fdr", type=float, default=0.05,
                        help="FDR threshold the array was submitted with "
                             "(default: 0.05)")
    args = parser.parse_args()

    anc = args.ancestry
    label = f"{args.modality}_ungrouped" if args.no_groups else args.modality
    output_dir = args.output_dir

    cis_path = os.path.join(output_dir, f"{anc}_{label}_cisqtl.parquet")
    ind_path = os.path.join(
        output_dir, f"{anc}_{label}_cisqtl_independent.parquet")
    ind_top_path = os.path.join(
        output_dir, f"{anc}_{label}_cisqtl_independent_top.tsv")
    chunk_dir = os.path.join(output_dir, "independent_chunks",
                             f"{anc}_{label}")

    print(f"[{anc} / {label}] merge chunked independent-scan results")

    if not os.path.exists(cis_path):
        sys.exit(f"ERROR: map_cis results not found: {cis_path}")

    qval = pd.to_numeric(
        pd.read_parquet(cis_path, columns=["qval"])["qval"],
        errors="coerce")
    n_sig = int((qval <= args.independent_fdr).sum())
    n_chunks = int(np.ceil(n_sig / args.chunk_size)) if n_sig else 0
    print(f"  {n_sig} significant phenotypes/groups at FDR <= "
          f"{args.independent_fdr} -> {n_chunks} expected chunks "
          f"(chunk size {args.chunk_size})")

    if n_chunks == 0:
        pd.DataFrame().to_parquet(ind_path)
        pd.DataFrame().to_csv(ind_top_path, sep="\t", index=False)
        print(f"  Nothing significant — wrote empty outputs:")
        print(f"    {ind_path}")
        print(f"    {ind_top_path}")
        return

    # ---- Verify all expected chunks exist with matching parameters ----
    missing, mismatched = [], []
    chunk_paths = []
    for k in range(1, n_chunks + 1):
        p = os.path.join(chunk_dir, f"chunk_{k:04d}.parquet")
        chunk_paths.append(p)
        if not os.path.exists(p):
            missing.append(k)
            continue
        sidecar = p.replace(".parquet", ".json")
        if not os.path.exists(sidecar):
            mismatched.append((k, "missing sidecar .json"))
            continue
        with open(sidecar) as f:
            meta = json.load(f)
        if meta.get("chunk_size") != args.chunk_size:
            mismatched.append(
                (k, f"chunk_size {meta.get('chunk_size')} != "
                    f"{args.chunk_size}"))
        elif not math.isclose(float(meta.get("independent_fdr", -1)),
                              args.independent_fdr):
            mismatched.append(
                (k, f"independent_fdr {meta.get('independent_fdr')} != "
                    f"{args.independent_fdr}"))

    if missing or mismatched:
        if missing:
            idx = ",".join(map(str, missing))
            print(f"  ERROR: {len(missing)} of {n_chunks} chunk files missing "
                  f"from {chunk_dir}", file=sys.stderr)
            print(f"    missing chunk indices: {idx}", file=sys.stderr)
            print(f"  Resubmit just those chunks, then the merge, with:",
                  file=sys.stderr)
            print(f"    CHUNK_INDICES=\"{idx}\" ANCESTRIES={anc} "
                  f"MODALITIES={args.modality} bash 28b_submit_independent.sh",
                  file=sys.stderr)
        if mismatched:
            print(f"  ERROR: {len(mismatched)} chunk(s) were computed with "
                  f"different parameters (stale chunks?):", file=sys.stderr)
            for k, why in mismatched[:20]:
                print(f"    chunk {k}: {why}", file=sys.stderr)
            print(f"  Delete the stale chunk directory and resubmit the full "
                  f"array:", file=sys.stderr)
            print(f"    rm -r {chunk_dir}", file=sys.stderr)
            print(f"    ANCESTRIES={anc} MODALITIES={args.modality} "
                  f"bash 28b_submit_independent.sh", file=sys.stderr)
        sys.exit(1)

    # ---- Concatenate chunks (chunk order == significant-row order == the
    # order the monolithic run processes phenotypes/groups) ----
    # ignore_index=True: each chunk's map_independent output carries its own
    # 0-based RangeIndex; resetting reproduces the monolithic output's clean
    # RangeIndex exactly.
    dfs = [pd.read_parquet(p) for p in chunk_paths]
    ind_result = pd.concat(dfs, ignore_index=True)
    n_rows = len(ind_result)

    ind_result.to_parquet(ind_path)
    print(f"  Written: {ind_path} ({n_rows} independent associations)")

    ind_top = ind_result.copy()
    if len(ind_top) and not isinstance(ind_top.index, pd.RangeIndex):
        ind_top = ind_top.reset_index()
    if len(ind_top) and 'pval_nominal' in ind_top.columns:
        ind_top = ind_top.sort_values('pval_nominal')
    ind_top.to_csv(ind_top_path, sep='\t', index=False)
    print(f"  Written: {ind_top_path} ({len(ind_top)} associations)")

    if len(ind_result) and 'rank' in ind_result.columns:
        if 'phenotype_id' in ind_result.columns:
            gid = ind_result['phenotype_id']
        elif not isinstance(ind_result.index, pd.RangeIndex):
            gid = pd.Series(ind_result.index)
        else:
            gid = None
        if gid is not None:
            max_rank = ind_result['rank'].groupby(gid.values).max()
            print(f"  Genes/groups with >1 conditionally independent xQTL: "
                  f"{(max_rank > 1).sum()}")

    n_chunk_files = len(glob.glob(os.path.join(chunk_dir, "chunk_*.parquet")))
    if n_chunk_files > n_chunks:
        print(f"  NOTE: {n_chunk_files - n_chunks} extra chunk file(s) in "
              f"{chunk_dir} beyond the {n_chunks} expected — ignored. "
              f"(Stale files from an earlier submission; safe to delete.)")


if __name__ == "__main__":
    main()
