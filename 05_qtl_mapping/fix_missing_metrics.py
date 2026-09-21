#!/usr/bin/env python3
"""
fix_missing_metrics.py — Re-run only the Picard tools that failed or had
incorrect column names in sharded QC output from 19a_picard_sharded.sh.

Re-runs per sample:
  1. CollectAlignmentSummaryMetrics  — fixes empty READS_ALIGNED /
     PCT_READS_ALIGNED (Picard 2.27 renamed to PF_READS_ALIGNED /
     PCT_PF_READS_ALIGNED; old code looked for the wrong key).
  2. CollectRnaSeqMetrics             — was failing entirely (now with
     VALIDATION_STRINGENCY=LENIENT; full stderr captured for diagnosis).
  3. CollectGcBiasMetrics             — was failing (needed R=<fasta>;
     now parses the summary file for AT_DROPOUT / GC_DROPOUT / MEAN_COVERAGE).

Existing InsertSize and DupMetrics columns are preserved unchanged.
The output file is written in-place (updated TSV with all columns).

Usage:
  python3 fix_missing_metrics.py \
    --chunk-file  /path/to/cohort1.chunk1.qcmetrics.tsv \
    --bam-dir     /path/to/cohort1/intermediate/star_out \
    --refflat     /path/to/HPLRv2.refFlat \
    --fasta       /path/to/GRCh38.fa \
    --picard-cmd  picard
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Picard metrics parsing (same logic as picard_qc.py)
# ---------------------------------------------------------------------------

def parse_picard_metrics(metrics_file):
    """Parse a PicardTools metrics file into {column: value}."""
    if not os.path.exists(metrics_file):
        return {}
    with open(metrics_file) as f:
        lines = f.readlines()
    start = None
    for i, line in enumerate(lines):
        if line.startswith("## METRICS CLASS"):
            start = i
            break
    if start is None or start + 2 >= len(lines):
        return {}
    header = lines[start + 1].strip().split("\t")
    data = lines[start + 2].strip().split("\t")
    if len(data) < len(header):
        data.extend([""] * (len(header) - len(data)))
    return dict(zip(header, data))


def coerce_numeric(val):
    if val is None or val == "" or val == "?":
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def run_picard(cmd, label, sample):
    """Run a Picard command. Returns True on success, False on failure.
    Prints FULL stderr on failure (not just last 3 lines)."""
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600, check=False
        )
        if result.returncode != 0:
            sys.stderr.write(
                f"  WARN [{label}] {sample}: picard exited {result.returncode}\n")
            if result.stderr:
                for line in result.stderr.strip().split("\n"):
                    sys.stderr.write(f"        {line}\n")
            return False
        return True
    except subprocess.TimeoutExpired:
        sys.stderr.write(f"  WARN [{label}] {sample}: timed out (3600s)\n")
        return False
    except Exception as e:
        sys.stderr.write(f"  WARN [{label}] {sample}: {e}\n")
        return False


# ---------------------------------------------------------------------------
# BAM discovery (same logic as picard_qc.py)
# ---------------------------------------------------------------------------

def find_bam(sample, star_out_dir, bam_dir=None):
    pre_shrink = os.path.join(star_out_dir, f"{sample}.Aligned.sortedByCoord.out.bam")
    if os.path.exists(pre_shrink):
        return pre_shrink, "pre-shrink"
    if bam_dir:
        shrunk = os.path.join(bam_dir, f"{sample}.bam")
        if os.path.exists(shrunk):
            return shrunk, "shrunk"
    return None, None


# ---------------------------------------------------------------------------
# Picard tool runners (fixed versions)
# ---------------------------------------------------------------------------

def collect_alignment_summary(picard_cmd, bam, tmpdir, sample):
    out = os.path.join(tmpdir, "alignment_metrics.txt")
    cmd = [picard_cmd, "CollectAlignmentSummaryMetrics",
           f"I={bam}", f"O={out}", "VALIDATION_STRINGENCY=LENIENT"]
    if run_picard(cmd, "AlignSummary", sample):
        m = parse_picard_metrics(out)
        return {
            "AlignMetrics.TOTAL_READS": coerce_numeric(m.get("TOTAL_READS")),
            "AlignMetrics.READS_ALIGNED": coerce_numeric(
                m.get("PF_READS_ALIGNED", m.get("READS_ALIGNED"))),
            "AlignMetrics.MEAN_READ_LENGTH": coerce_numeric(m.get("MEAN_READ_LENGTH")),
            "AlignMetrics.STRAND_BALANCE": coerce_numeric(m.get("STRAND_BALANCE")),
            "AlignMetrics.PCT_READS_ALIGNED": coerce_numeric(
                m.get("PCT_PF_READS_ALIGNED", m.get("PCT_READS_ALIGNED"))),
        }
    return {}


def collect_rna_seq_metrics(picard_cmd, bam, tmpdir, sample, refflat):
    out = os.path.join(tmpdir, "rna_metrics.txt")
    cmd = [picard_cmd, "CollectRnaSeqMetrics",
           f"I={bam}", f"O={out}", f"REF_FLAT={refflat}",
           "STRAND=NONE", "VALIDATION_STRINGENCY=LENIENT"]
    if run_picard(cmd, "RnaSeqMetrics", sample):
        m = parse_picard_metrics(out)
        return {
            "RnaMetrics.PCT_RIBOSOMAL_BASES": coerce_numeric(m.get("PCT_RIBOSOMAL_BASES")),
            "RnaMetrics.PCT_CODING_BASES": coerce_numeric(m.get("PCT_CODING_BASES")),
            "RnaMetrics.PCT_UTR_BASES": coerce_numeric(m.get("PCT_UTR_BASES")),
            "RnaMetrics.PCT_INTRONIC_BASES": coerce_numeric(m.get("PCT_INTRONIC_BASES")),
            "RnaMetrics.PCT_INTERGENIC_BASES": coerce_numeric(m.get("PCT_INTERGENIC_BASES")),
            "RnaMetrics.MEDIAN_5PRIME_BIAS": coerce_numeric(m.get("MEDIAN_5PRIME_BIAS")),
            "RnaMetrics.MEDIAN_3PRIME_BIAS": coerce_numeric(m.get("MEDIAN_3PRIME_BIAS")),
            "RnaMetrics.MEDIAN_5PRIME_TO_3PRIME_BIAS": coerce_numeric(
                m.get("MEDIAN_5PRIME_TO_3PRIME_BIAS")),
        }
    return {}


def collect_gc_bias(picard_cmd, bam, tmpdir, sample, fasta):
    out = os.path.join(tmpdir, "gc_bias_metrics.txt")
    summary = os.path.join(tmpdir, "gc_bias_summary.txt")
    chart = os.path.join(tmpdir, "gc_bias_chart.pdf")
    cmd = [picard_cmd, "CollectGcBiasMetrics",
           f"I={bam}", f"O={out}", f"S={summary}", f"CHART={chart}",
           f"R={fasta}", "VALIDATION_STRINGENCY=LENIENT"]
    if run_picard(cmd, "GcBias", sample):
        m = parse_picard_metrics(summary)
        return {
            "GcBias.AT_DROPOUT": coerce_numeric(m.get("AT_DROPOUT")),
            "GcBias.GC_DROPOUT": coerce_numeric(m.get("GC_DROPOUT")),
            "GcBias.MEAN_COVERAGE": coerce_numeric(m.get("MEAN_COVERAGE")),
        }
    return {}


# ---------------------------------------------------------------------------
# Diagnostic: compare chromosome names between BAM and refFlat
# ---------------------------------------------------------------------------

def diagnose_chr_mismatch(bam, refflat):
    """Print first few chromosome names from BAM header and refFlat."""
    try:
        import pysam
        bam_chroms = set()
        with pysam.AlignmentFile(bam, "rb") as bf:
            bam_chroms = set(bf.references[:5])
        ref_chroms = set()
        with open(refflat) as f:
            for i, line in enumerate(f):
                if i >= 5:
                    break
                ref_chroms.add(line.split("\t")[1])
        print(f"  BAM chroms (first 5):  {sorted(bam_chroms)}")
        print(f"  refFlat chroms (first 5): {sorted(ref_chroms)}")
        if bam_chroms and ref_chroms:
            if not any(c in ref_chroms for c in bam_chroms):
                print("  WARNING: no overlap between BAM and refFlat chromosome names!")
    except ImportError:
        print("  (pysam not available, skipping chr diagnostic)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Re-run failed Picard tools and merge into existing chunk TSV")
    parser.add_argument("--chunk-file", required=True,
                        help="Existing chunk TSV from 19a")
    parser.add_argument("--bam-dir", required=True,
                        help="star_out directory containing pre-shrink BAMs")
    parser.add_argument("--refflat", required=True)
    parser.add_argument("--fasta", required=True,
                        help="Reference FASTA (required for CollectGcBiasMetrics)")
    parser.add_argument("--picard-cmd", default="picard")
    parser.add_argument("--output", default=None,
                        help="Output path (default: overwrite --chunk-file in place)")
    args = parser.parse_args()

    out_path = args.output or args.chunk_file

    # ---- Read existing chunk TSV ----
    df = pd.read_csv(args.chunk_file, sep="\t", index_col=0, dtype=str)
    # Drop duplicate "sample" column if present (from earlier bug)
    if "sample" in df.columns:
        df = df.drop(columns=["sample"])
    df.index.name = "sample"
    samples = list(df.index)
    print(f"Loaded {len(samples)} samples from {args.chunk_file}")
    print(f"  Existing columns: {list(df.columns)}")

    # ---- Diagnostic on first sample ----
    if samples:
        first_bam, bam_type = find_bam(samples[0], args.bam_dir)
        if first_bam:
            print(f"\nDiagnostic (first sample {samples[0]}):")
            diagnose_chr_mismatch(first_bam, args.refflat)

    # ---- Re-run 3 Picard tools per sample ----
    new_rows = []
    for i, sample in enumerate(samples):
        bam, bam_type = find_bam(sample, args.bam_dir)
        if bam is None:
            sys.stderr.write(f"  [{i+1}/{len(samples)}] {sample}: BAM not found, skipping\n")
            new_rows.append({})
            continue

        print(f"  [{i+1}/{len(samples)}] {sample} ({bam_type} BAM)")
        with tempfile.TemporaryDirectory() as tmpdir:
            metrics = {}
            metrics.update(
                collect_alignment_summary(args.picard_cmd, bam, tmpdir, sample))
            metrics.update(
                collect_rna_seq_metrics(args.picard_cmd, bam, tmpdir, sample, args.refflat))
            metrics.update(
                collect_gc_bias(args.picard_cmd, bam, tmpdir, sample, args.fasta))
        new_rows.append(metrics)

    # ---- Merge new columns into existing DataFrame ----
    new_df = pd.DataFrame(new_rows, index=df.index)
    new_df.index.name = "sample"

    # Update existing columns (AlignMetrics fixes) and add new ones
    for col in new_df.columns:
        df[col] = new_df[col]

    # ---- Report ----
    metric_groups = ["AlignMetrics", "InsertSize", "RnaMetrics",
                     "GcBias", "DupMetrics", "SubjectBias"]
    print(f"\nFinal columns ({len(df.columns)}):")
    for g in metric_groups:
        cols = [c for c in df.columns if c.startswith(g + ".")]
        status = "OK" if cols else "MISSING"
        print(f"  {g}: {len(cols)} cols [{status}]")

    # ---- Write ----
    df.to_csv(out_path, sep="\t", index=True)
    print(f"\nWrote {len(df)} samples x {len(df.columns)} metrics to {out_path}")


if __name__ == "__main__":
    main()
