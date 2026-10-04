#!/usr/bin/env python3
"""
picard_persample.py — Scan/merge helper for sample-parallel Picard QC.

Companion to 19a_picard_sharded.sh (sample-parallel version). Two subcommands:

  scan   — build the per-cohort todo list: samples.txt minus samples already
           present in legacy chunk outputs (${COHORT}.chunk*.qcmetrics.tsv)
           or in valid per-sample outputs (per_sample/${COHORT}/<sample>.qcmetrics.tsv).
           Writes the todo file and prints the todo count on stdout (last line);
           all diagnostics go to stderr so the count can be captured with $(...).

  merge  — combine legacy chunk rows (kept verbatim — they carry the 19b
           fixes) with per-sample rows for samples not in any chunk file
           (chunk wins on conflict). Missing cells in per-sample-derived rows
           are imputed with cohort-level column medians (computed over the
           full merged cohort). Writes ${COHORT}_qc_metrics.tsv per cohort —
           the same files 19_hcp_factors.sh Stage 1/2 expects.

Backwards compatibility: chunk files are never modified. Per-sample outputs
live in the per_sample/ subdirectory so the Stage-2 glob
(${QC_DIR}/*_qc_metrics.tsv) never picks them up.

Usage:
  python3 picard_persample.py scan  --samples samples.txt --qc-dir QC_DIR \
      --cohort cohort1 --todo-out TODO_FILE
  python3 picard_persample.py merge --qc-dir QC_DIR [--output-base DIR] \
      [--cohorts cohort1 cohort3]
"""

import argparse
import glob
import os
import sys

PER_SAMPLE_DIRNAME = "per_sample"
PER_SAMPLE_SUFFIX = ".qcmetrics.tsv"


# =============================================================================
# Scan (no pandas required — runs on the login node at submit time)
# =============================================================================

def _read_chunk_samples(qc_dir, cohort):
    """Sample IDs (first column) from all legacy chunk files for a cohort."""
    done = set()
    pattern = os.path.join(qc_dir, f"{cohort}.chunk*.qcmetrics.tsv")
    for path in sorted(glob.glob(pattern)):
        with open(path) as f:
            header = f.readline()
            if not header.split("\t")[0].strip() == "sample":
                sys.stderr.write(f"  WARN: unexpected header in {path}, skipping\n")
                continue
            for line in f:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                done.add(line.split("\t")[0].strip())
    return done


def _valid_per_sample_file(path):
    """A per-sample output counts as done only if it has the expected header
    and at least one data row (guards against partial/corrupt writes)."""
    try:
        with open(path) as f:
            header = f.readline()
            if not header:
                return False
            first_col = header.split("\t")[0].strip()
            if first_col != "sample":
                return False
            for line in f:
                if line.strip():
                    return True
        return False
    except OSError:
        return False


def _read_per_sample_done(qc_dir, cohort):
    """Sample IDs with a valid per-sample output file."""
    done = set()
    ps_dir = os.path.join(qc_dir, PER_SAMPLE_DIRNAME, cohort)
    for path in sorted(glob.glob(os.path.join(ps_dir, f"*{PER_SAMPLE_SUFFIX}"))):
        if _valid_per_sample_file(path):
            sample = os.path.basename(path)[: -len(PER_SAMPLE_SUFFIX)]
            done.add(sample)
        else:
            sys.stderr.write(f"  WARN: invalid/partial per-sample file ignored: {path}\n")
    return done


def done_samples(qc_dir, cohort):
    """Union of samples completed by legacy chunks and per-sample jobs."""
    return _read_chunk_samples(qc_dir, cohort) | _read_per_sample_done(qc_dir, cohort)


def write_todo(samples_file, qc_dir, cohort, todo_out):
    """Write the todo list (samples.txt minus done-set, deduped, original
    order preserved). Returns (todo, n_done)."""
    done = done_samples(qc_dir, cohort)
    todo = []
    seen = set()
    with open(samples_file) as f:
        for line in f:
            sample = line.strip()
            if not sample or sample in seen:
                continue
            seen.add(sample)
            if sample not in done:
                todo.append(sample)
    with open(todo_out, "w") as out:
        for sample in todo:
            out.write(sample + "\n")
    return todo, len(done & seen)


# =============================================================================
# Merge (pandas required)
# =============================================================================

def merge_cohort(qc_dir, cohort, samples_file=None):
    """Merge legacy chunk rows + new per-sample rows for one cohort.

    Returns (merged_df, stats dict). Writes nothing — caller persists.
    """
    import numpy as np
    import pandas as pd

    stats = {"cohort": cohort, "n_chunk_rows": 0, "n_new_rows": 0,
             "n_chunk_wins_skipped": 0, "n_imputed_cells": 0,
             "n_missing_samples": 0}

    # ---- Legacy chunk rows (authoritative; kept verbatim) ----
    chunk_files = sorted(glob.glob(os.path.join(
        qc_dir, f"{cohort}.chunk*.qcmetrics.tsv")))
    chunk_dfs = [pd.read_csv(f, sep="\t", index_col=0) for f in chunk_files]
    chunk_df = (pd.concat(chunk_dfs) if chunk_dfs
                else pd.DataFrame())
    if not chunk_df.empty:
        chunk_df = chunk_df[~chunk_df.index.duplicated(keep="first")]
    stats["n_chunk_files"] = len(chunk_files)
    stats["n_chunk_rows"] = len(chunk_df)

    # ---- Per-sample rows for samples not in any chunk (chunk wins) ----
    ps_dir = os.path.join(qc_dir, PER_SAMPLE_DIRNAME, cohort)
    new_rows = []
    for path in sorted(glob.glob(os.path.join(ps_dir, f"*{PER_SAMPLE_SUFFIX}"))):
        sample = os.path.basename(path)[: -len(PER_SAMPLE_SUFFIX)]
        if not _valid_per_sample_file(path):
            sys.stderr.write(f"  WARN: skipping invalid per-sample file: {path}\n")
            continue
        df = pd.read_csv(path, sep="\t", index_col=0)
        if len(df) != 1:
            sys.stderr.write(
                f"  WARN: {path} has {len(df)} rows (expected 1), skipping\n")
            continue
        if sample in chunk_df.index or df.index[0] in chunk_df.index:
            stats["n_chunk_wins_skipped"] += 1
            continue
        new_rows.append(df)
    new_df = pd.concat(new_rows) if new_rows else pd.DataFrame()
    if not new_df.empty:
        new_df = new_df[~new_df.index.duplicated(keep="first")]
    stats["n_new_rows"] = len(new_df)

    if chunk_df.empty and new_df.empty:
        return pd.DataFrame(), stats

    merged = pd.concat([chunk_df, new_df])
    merged.index.name = "sample"

    # ---- Cohort-level median imputation, new rows only ----
    if not new_df.empty:
        num_cols = merged.select_dtypes(include=[np.number]).columns
        medians = merged[num_cols].median().fillna(0.0)
        block = merged.loc[new_df.index, num_cols]
        n_missing = int(block.isna().sum().sum())
        if n_missing:
            merged.loc[new_df.index, num_cols] = block.fillna(medians)
            stats["n_imputed_cells"] = n_missing

    # Safety net: schema drift (a column absent from all chunk files) would
    # leave NaN in chunk rows. Should not happen (same extractor throughout);
    # fill loudly rather than write NaN into the HCP prior matrix.
    num_cols = merged.select_dtypes(include=[np.number]).columns
    leftover = int(merged[num_cols].isna().sum().sum())
    if leftover:
        medians = merged[num_cols].median().fillna(0.0)
        merged[num_cols] = merged[num_cols].fillna(medians)
        sys.stderr.write(
            f"  WARN: {cohort}: {leftover} NaN cells remained outside new rows "
            f"(schema drift?); filled with cohort medians\n")

    # ---- Missing-sample report (vs samples.txt, if available) ----
    if samples_file and os.path.exists(samples_file):
        with open(samples_file) as f:
            expected = {line.strip() for line in f if line.strip()}
        missing = sorted(expected - set(merged.index))
        stats["n_missing_samples"] = len(missing)
        if missing:
            miss_path = os.path.join(qc_dir, f"{cohort}.missing_samples.txt")
            with open(miss_path, "w") as out:
                out.write("\n".join(missing) + "\n")
            stats["missing_file"] = miss_path

    return merged, stats


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_scan = sub.add_parser("scan", help="write todo list; print todo count")
    p_scan.add_argument("--samples", required=True, help="cohort samples.txt")
    p_scan.add_argument("--qc-dir", required=True)
    p_scan.add_argument("--cohort", required=True)
    p_scan.add_argument("--todo-out", required=True)

    p_merge = sub.add_parser("merge", help="merge chunks + per-sample outputs")
    p_merge.add_argument("--qc-dir", required=True)
    p_merge.add_argument("--output-base",
                         help="cohort parent dir (for samples.txt missing report)")
    p_merge.add_argument("--cohorts", nargs="*",
                         help="cohorts to merge (default: all discovered)")

    args = parser.parse_args()

    if args.cmd == "scan":
        todo, n_done = write_todo(args.samples, args.qc_dir, args.cohort,
                                  args.todo_out)
        sys.stderr.write(
            f"  {args.cohort}: {n_done} done, {len(todo)} to do -> {args.todo_out}\n")
        print(len(todo))  # stdout: machine-readable count (last/only line)
        return

    # ---- merge ----
    import pandas as pd

    qc_dir = args.qc_dir
    if args.cohorts:
        cohorts = args.cohorts
    else:
        cohorts = set()
        for f in glob.glob(os.path.join(qc_dir, "*.chunk*.qcmetrics.tsv")):
            cohorts.add(os.path.basename(f).split(".chunk")[0])
        ps_root = os.path.join(qc_dir, PER_SAMPLE_DIRNAME)
        if os.path.isdir(ps_root):
            for d in os.listdir(ps_root):
                if os.path.isdir(os.path.join(ps_root, d)):
                    cohorts.add(d)
        cohorts = sorted(cohorts)
    if not cohorts:
        raise SystemExit(f"ERROR: no chunk files or per-sample dirs in {qc_dir}")

    for co in cohorts:
        samples_file = None
        if args.output_base:
            candidate = os.path.join(args.output_base, co, "samples.txt")
            if os.path.exists(candidate):
                samples_file = candidate
        merged, stats = merge_cohort(qc_dir, co, samples_file=samples_file)
        if merged.empty:
            sys.stderr.write(f"  {co}: nothing to merge, skipping\n")
            continue
        out = os.path.join(qc_dir, f"{co}_qc_metrics.tsv")
        merged.to_csv(out, sep="\t")
        print(f"{co}: {stats['n_chunk_rows']} legacy chunk rows "
              f"({stats['n_chunk_files']} files) + {stats['n_new_rows']} new "
              f"per-sample rows -> {merged.shape[0]} samples x "
              f"{merged.shape[1]} metrics")
        if stats["n_chunk_wins_skipped"]:
            print(f"  {stats['n_chunk_wins_skipped']} per-sample files ignored "
                  f"(sample already in a chunk file; chunk row kept)")
        if stats["n_imputed_cells"]:
            print(f"  {stats['n_imputed_cells']} missing cells in new rows "
                  f"imputed with cohort-level medians")
        if stats["n_missing_samples"]:
            print(f"  {stats['n_missing_samples']} samples still missing "
                  f"(see {stats['missing_file']}) — resubmit to fill them")
        print(f"  Written: {out}")
    print("Merge complete. Now run 19_hcp_factors.sh (it will skip Stage 1).")


if __name__ == "__main__":
    main()
