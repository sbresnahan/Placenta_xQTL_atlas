#!/usr/bin/env python3
"""
picard_qc.py — Collect PicardTools sequencing QC metrics per sample.

Standalone script that finds pre-shrink STAR BAMs and runs a curated subset of
PicardTools Collect* tools to produce the HCP prior matrix F (sequencing QC
metrics used as priors for Hidden Covariates with Prior estimation).

Metric panel: devBrain-style full field set — every numeric field of each
tool's summary table is emitted as `Tool.FIELD`, dropping only
identifier/categorical columns (QC_EXCLUDE_FIELDS, mirroring devBrain's
picard_merge_metrics.R). Multi-row tables (AlignmentSummaryMetrics) emit the
library-level PAIR row unsuffixed and other categories suffixed
(.FIRST_OF_PAIR / .SECOND_OF_PAIR). The original 24 column names are
preserved (FIELD_RENAME maps Picard >=2.27 PF_* names back), but
AlignMetrics.* values are now library-level (PAIR row) rather than
FIRST_OF_PAIR. Downstream zero-variance + |r|>0.9 pruning in
hcp_from_matrix.R handles redundancy.

Raw Picard outputs persist under <output dir>/raw/<sample>/ so the panel can
be re-extracted later with --parse-only (no BAM access needed).

Also computes two subject-specific covariates from the original HCP paper
(Mostafavi et al. 2013):
  - Subject-specific GC bias: Pearson corr(log2(TPM+1), gene GC content) per sample
  - Subject-specific gene-length bias: Pearson corr(log2(TPM+1), gene length) per sample

These require the gene-level TPM matrix (from Salmon quant.sf) and a gene
GC/length annotation table.

Usage:
  # Collect QC metrics for one cohort
  python3 picard_qc.py \
      --config config.yml \
      --cohort cohort1 \
      --output qc_metrics.tsv \
      --refflat HPLRv2.refFlat

  # Generate refFlat from the HPLRv2 PANTRY-normalized GTF
  python3 picard_qc.py --generate-refflat \
      --gtf HPLRv2.0.annotated.PANTRY.gtf \
      --output HPLRv2.refFlat

  # Compute gene GC/length annotation table (for subject-specific bias)
  python3 picard_qc.py --generate-gene-annot \
      --gtf HPLRv2.0.annotated.PANTRY.gtf \
      --fasta genome.fasta \
      --output gene_gc_length.tsv

Requires PicardTools on seadragon:
  module load picard
  eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
  conda activate picard-2.27.4
"""

import argparse
import os
import re
import subprocess
import sys
# tempfile no longer needed: raw Picard outputs persist under --raw-dir
from collections import defaultdict

import numpy as np
import pandas as pd


# =============================================================================
# PicardTools metric parsing
# =============================================================================

def parse_picard_metrics(metrics_file):
    """Parse a PicardTools metrics file into a dict of {column: value}.

    PicardTools metrics files have a header section, then a METRICS CLASS line,
    then a column header line, then one or more data rows. We extract the first
    data row (most tools produce one row per BAM).
    """
    if not os.path.exists(metrics_file):
        return {}

    with open(metrics_file) as f:
        lines = f.readlines()

    metrics_start = None
    for i, line in enumerate(lines):
        if line.startswith("## METRICS CLASS"):
            metrics_start = i
            break

    if metrics_start is None:
        return {}

    if metrics_start + 2 >= len(lines):
        return {}

    header = lines[metrics_start + 1].strip().split("\t")
    data = lines[metrics_start + 2].strip().split("\t")

    if len(data) < len(header):
        data.extend([""] * (len(header) - len(data)))

    return dict(zip(header, data))


def coerce_numeric(val):
    """Try to convert a string to float; return None if not numeric."""
    if val is None or val == "" or val == "?":
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


# =============================================================================
# Full-panel extraction (devBrain-style)
# =============================================================================
# The original pipeline extracted a 24-metric whitelist. The devBrain schema
# (picard_merge_metrics.R) instead keeps the full MultiQC-style field set,
# dropping only identifier/categorical columns and letting downstream
# zero-variance + |r|>0.9 pruning handle redundancy. We mirror that here:
# every numeric field of each tool's summary table is emitted as
# `Tool.FIELD`, with the exclusions below.

QC_EXCLUDE_FIELDS = {
    "SAMPLE", "LIBRARY", "READ_GROUP", "CATEGORY",
    "ACCUMULATION_LEVEL", "READS_USED", "PAIR_ORIENTATION",
}

# Picard >= 2.27 renamed these fields; map the new raw names back to the
# column names the pipeline has always emitted (applied to the primary row
# only, so the original 24 column names are preserved).
FIELD_RENAME = {
    "PF_READS_ALIGNED": "READS_ALIGNED",
    "PCT_PF_READS_ALIGNED": "PCT_READS_ALIGNED",
}


def parse_picard_metrics_all(metrics_file):
    """Parse a PicardTools metrics file into a list of {column: value} dicts,
    one per data row. Multi-row tables (e.g. AlignmentSummaryMetrics) carry
    one row per CATEGORY. Parsing stops at the end of the metrics table
    (blank line or the next ## section, e.g. a HISTOGRAM_CLASS block)."""
    if not os.path.exists(metrics_file):
        return []

    with open(metrics_file) as f:
        lines = f.readlines()

    metrics_start = None
    for i, line in enumerate(lines):
        if line.startswith("## METRICS CLASS"):
            metrics_start = i
            break

    if metrics_start is None or metrics_start + 2 >= len(lines):
        return []

    header = lines[metrics_start + 1].strip().split("\t")
    rows = []
    for line in lines[metrics_start + 2:]:
        if not line.strip() or line.startswith("#"):
            break
        data = line.strip().split("\t")
        if len(data) < len(header):
            data.extend([""] * (len(header) - len(data)))
        rows.append(dict(zip(header, data)))
    return rows


def extract_full_panel(rows, prefix):
    """Flatten ALL numeric fields of a Picard metrics table into
    `prefix.FIELD` columns (devBrain-style full panel).

    - ID/categorical columns in QC_EXCLUDE_FIELDS are dropped.
    - The primary row (CATEGORY=PAIR when present, else UNPAIRED, else the
      first row) is emitted unsuffixed; additional rows (e.g.
      FIRST_OF_PAIR / SECOND_OF_PAIR, or extra LIBRARY rows) are suffixed
      with their category value.
    - FIELD_RENAME maps Picard >= 2.27 names back to the pipeline's
      historical column names on the primary row.

    NOTE: for multi-category tables the primary row is the library-level
    PAIR row. The legacy whitelist read the FIRST data row (FIRST_OF_PAIR
    for paired-end BAMs), so AlignMetrics.* values are library-level now —
    column names are unchanged.
    """
    if not rows:
        return {}
    cats = [r.get("CATEGORY", "") for r in rows]
    primary_idx = 0
    for want in ("PAIR", "UNPAIRED"):
        if want in cats:
            primary_idx = cats.index(want)
            break
    multi = len(rows) > 1
    out = {}
    for i, row in enumerate(rows):
        suffix = ""
        if multi and i != primary_idx:
            cat = row.get("CATEGORY") or row.get("LIBRARY") or f"ROW{i + 1}"
            suffix = f".{cat}"
        for field, val in row.items():
            if field in QC_EXCLUDE_FIELDS:
                continue
            num = coerce_numeric(val)
            if num is None:
                continue
            if not suffix:
                field = FIELD_RENAME.get(field, field)
            out[f"{prefix}.{field}{suffix}"] = num
    return out


# =============================================================================
# PicardTools Collect* tool runners
# =============================================================================

def run_picard(cmd, label, sample):
    """Run a PicardTools command, capturing stderr. Returns True on success."""
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600, check=False
        )
        if result.returncode != 0:
            sys.stderr.write(
                f"  WARN [{label}] {sample}: picard exited {result.returncode}\n"
            )
            if result.stderr:
                err_lines = result.stderr.strip().split("\n")
                for line in err_lines:
                    sys.stderr.write(f"        {line}\n")
            return False
        return True
    except subprocess.TimeoutExpired:
        sys.stderr.write(f"  WARN [{label}] {sample}: timed out (3600s)\n")
        return False
    except Exception as e:
        sys.stderr.write(f"  WARN [{label}] {sample}: {e}\n")
        return False


def collect_alignment_summary(picard_cmd, bam, out_dir, sample, parse_only=False):
    """CollectAlignmentSummaryMetrics -> alignment_metrics.txt (full panel)"""
    out_file = os.path.join(out_dir, "alignment_metrics.txt")
    cmd = [picard_cmd, "CollectAlignmentSummaryMetrics", f"I={bam}", f"O={out_file}"]
    if parse_only or run_picard(cmd, "AlignSummary", sample):
        return extract_full_panel(parse_picard_metrics_all(out_file), "AlignMetrics")
    return {}


def collect_insert_size(picard_cmd, bam, out_dir, sample, parse_only=False):
    """CollectInsertSizeMetrics -> insert_size_metrics.txt (full panel)"""
    out_file = os.path.join(out_dir, "insert_size_metrics.txt")
    hist_file = os.path.join(out_dir, "insert_size_histogram.pdf")
    cmd = [picard_cmd, "CollectInsertSizeMetrics",
           f"I={bam}", f"O={out_file}", f"H={hist_file}"]
    if parse_only or run_picard(cmd, "InsertSize", sample):
        return extract_full_panel(parse_picard_metrics_all(out_file), "InsertSize")
    return {}


def collect_rna_seq_metrics(picard_cmd, bam, out_dir, sample, refflat,
                            parse_only=False):
    """CollectRnaSeqMetrics -> rna_metrics.txt (full panel)"""
    out_file = os.path.join(out_dir, "rna_metrics.txt")
    cmd = [picard_cmd, "CollectRnaSeqMetrics",
           f"I={bam}", f"O={out_file}", f"REF_FLAT={refflat}", "STRAND=NONE"]
    if parse_only or run_picard(cmd, "RnaSeqMetrics", sample):
        return extract_full_panel(parse_picard_metrics_all(out_file), "RnaMetrics")
    return {}


def collect_gc_bias(picard_cmd, bam, out_dir, sample, fasta=None,
                    parse_only=False):
    """CollectGcBiasMetrics -> gc_bias_metrics.txt (+ summary)

    R= (reference FASTA) is a REQUIRED argument for CollectGcBiasMetrics;
    without it the tool errors out and no GcBias columns are produced.
    The dropout/coverage summaries live in the S= summary metrics file,
    not the per-GC-bin detail table (O=).
    """
    out_file = os.path.join(out_dir, "gc_bias_metrics.txt")
    summary_file = os.path.join(out_dir, "gc_bias_summary.txt")
    chart_file = os.path.join(out_dir, "gc_bias_chart.pdf")
    cmd = [picard_cmd, "CollectGcBiasMetrics",
           f"I={bam}", f"O={out_file}", f"S={summary_file}", f"CHART={chart_file}"]
    if fasta:
        cmd.append(f"R={fasta}")
    if parse_only or run_picard(cmd, "GcBias", sample):
        return extract_full_panel(parse_picard_metrics_all(summary_file), "GcBias")
    return {}


def mark_duplicates(picard_cmd, bam, out_dir, sample, parse_only=False):
    """MarkDuplicates -> dup_metrics.txt (metrics only; temp BAM discarded)"""
    out_file = os.path.join(out_dir, "dup_metrics.txt")
    tmp_bam = os.path.join(out_dir, "tmp_dedup.bam")
    cmd = [picard_cmd, "MarkDuplicates",
           f"I={bam}", f"O={tmp_bam}", f"M={out_file}",
           "ASSUME_SORTED=true", "VALIDATION_STRINGENCY=LENIENT"]
    if parse_only or run_picard(cmd, "MarkDup", sample):
        result = extract_full_panel(parse_picard_metrics_all(out_file), "DupMetrics")
        if os.path.exists(tmp_bam):
            os.remove(tmp_bam)
        return result
    if os.path.exists(tmp_bam):
        os.remove(tmp_bam)
    return {}


# =============================================================================
# Subject-specific GC and gene-length bias (from HCP paper)
# =============================================================================

def compute_subject_specific_bias(tpm_matrix, gene_annot):
    """Compute subject-specific GC bias and gene-length bias.

    Following Mostafavi et al. 2013:
      - GC bias: Pearson corr(log2(TPM+1), gene GC content) per sample
      - Gene-length bias: Pearson corr(log2(TPM+1), gene length) per sample
    """
    common_genes = tpm_matrix.index.intersection(gene_annot.index)
    if len(common_genes) == 0:
        sys.stderr.write(
            "  WARN: no common genes between TPM matrix and gene annotation\n")
        # Keep 'sample' as the index (not a column) so the downstream join
        # aligns and no duplicate 'sample' column appears in the output
        return pd.DataFrame(
            columns=["SubjectBias.GC", "SubjectBias.GENE_LENGTH"]
        ).rename_axis("sample")

    tpm_sub = tpm_matrix.loc[common_genes]
    gc = gene_annot.loc[common_genes, "gc_content"].values
    length = gene_annot.loc[common_genes, "gene_length"].values
    log_tpm = np.log2(tpm_sub.values + 1)

    results = []
    for i, sample in enumerate(tpm_sub.columns):
        expr = log_tpm[:, i]
        gc_bias = np.corrcoef(expr, gc)[0, 1] if np.std(expr) > 0 and np.std(gc) > 0 else np.nan
        len_bias = np.corrcoef(expr, length)[0, 1] if np.std(expr) > 0 and np.std(length) > 0 else np.nan
        results.append({"sample": sample, "SubjectBias.GC": gc_bias,
                        "SubjectBias.GENE_LENGTH": len_bias})
    return pd.DataFrame(results).set_index("sample")


def load_tx_to_gene(gtf_path):
    """Parse transcript_id -> gene_id mapping from GTF transcript features.

    Needed because HPLRv2/SQANTI3 transcript IDs do not share a parseable
    prefix with gene IDs, so heuristic string-splitting of Salmon Names
    produces gene IDs that never match the GTF gene_id space.
    """
    tx_to_gene = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attr_dict = {}
            for m in re.finditer(r'(\w+)\s+"([^"]*)"', fields[8]):
                attr_dict[m.group(1)] = m.group(2)
            tid = attr_dict.get("transcript_id")
            gid = attr_dict.get("gene_id")
            if tid and gid:
                tx_to_gene[tid] = gid
    return tx_to_gene


def load_tpm_matrix(salmon_dir, samples, tx_to_gene=None):
    """Load gene-level TPM matrix from per-sample Salmon quant.sf files.

    Aggregates transcript-level TPM to gene-level by summing transcripts per gene.
    If tx_to_gene (from load_tx_to_gene) is given, transcript Names are mapped
    exactly to GTF gene_ids; otherwise falls back to a heuristic
    (Name prefix before the first '.' or '__').
    """
    tpm_dfs = []
    for sample in samples:
        quant_file = os.path.join(salmon_dir, sample, "quant.sf")
        if not os.path.exists(quant_file):
            sys.stderr.write(f"  WARN: quant.sf not found for {sample}: {quant_file}\n")
            continue
        df = pd.read_csv(quant_file, sep="\t")
        if tx_to_gene:
            df["gene_id"] = df["Name"].map(tx_to_gene)
            df = df.dropna(subset=["gene_id"])
        else:
            df["gene_id"] = df["Name"].str.split(".").str[0]
            df["gene_id"] = df["gene_id"].str.split("__").str[0]
        gene_tpm = df.groupby("gene_id")["TPM"].sum()
        gene_tpm.name = sample
        tpm_dfs.append(gene_tpm)

    if not tpm_dfs:
        return pd.DataFrame()

    tpm_matrix = pd.concat(tpm_dfs, axis=1)
    tpm_matrix.index.name = "gene_id"
    return tpm_matrix


# =============================================================================
# refFlat generation from GTF
# =============================================================================

def generate_refflat(gtf_path, output_path):
    """Convert a GTF file to refFlat format for PicardTools CollectRnaSeqMetrics.

    refFlat format (tab-delimited, one line per transcript, 11 fields):
      geneName name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds

    Coordinates are 0-based. exonStarts and exonEnds are comma-separated lists
    ending with a comma.
    """
    genes = {}  # gene_name -> {transcript_id -> tx dict}
    tx_to_gene = {}  # transcript_id -> gene_name

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            start = int(start) - 1  # GTF 1-based -> refFlat 0-based
            end = int(end)

            attr_dict = {}
            for m in re.finditer(r'(\w+)\s+"([^"]*)"', attrs):
                attr_dict[m.group(1)] = m.group(2)

            gene_name = attr_dict.get("gene_name") or attr_dict.get("gene_id", "")
            transcript_id = attr_dict.get("transcript_id", "")

            if feature == "transcript" and transcript_id:
                if gene_name not in genes:
                    genes[gene_name] = {}
                genes[gene_name][transcript_id] = {
                    "chrom": chrom, "strand": strand,
                    "txStart": start, "txEnd": end,
                    "cdsStart": start, "cdsEnd": end, "exons": [],
                }
                tx_to_gene[transcript_id] = gene_name

            elif feature == "exon" and transcript_id:
                gname = tx_to_gene.get(transcript_id)
                if gname and transcript_id in genes.get(gname, {}):
                    genes[gname][transcript_id]["exons"].append((start, end))

            elif feature == "CDS" and transcript_id:
                gname = tx_to_gene.get(transcript_id)
                if gname and transcript_id in genes.get(gname, {}):
                    tx = genes[gname][transcript_id]
                    if tx["cdsStart"] == tx["txStart"]:
                        tx["cdsStart"] = start
                    tx["cdsEnd"] = end

    with open(output_path, "w") as out:
        for gene_name in sorted(genes):
            txs = genes[gene_name]
            for tx_id in sorted(txs):
                tx = txs[tx_id]
                if not tx["exons"]:
                    continue
                exons = sorted(tx["exons"], key=lambda x: x[0])
                exon_count = len(exons)
                exon_starts = ",".join(str(e[0]) for e in exons) + ","
                exon_ends = ",".join(str(e[1]) for e in exons) + ","
                out.write("\t".join([
                    gene_name, tx_id, tx["chrom"], tx["strand"],
                    str(tx["txStart"]), str(tx["txEnd"]),
                    str(tx["cdsStart"]), str(tx["cdsEnd"]),
                    str(exon_count), exon_starts, exon_ends,
                ]) + "\n")

    n_genes = len(genes)
    n_tx = sum(len(txs) for txs in genes.values())
    print(f"refFlat written: {output_path} ({n_genes} genes, {n_tx} transcripts)")


# =============================================================================
# Gene GC/length annotation table generation
# =============================================================================

def generate_gene_gc_length(gtf_path, fasta_path, output_path):
    """Compute per-gene GC content and gene length from GTF + genome FASTA.

    Gene length = sum of exon lengths (union of merged exons per gene).
    GC content = fraction of G/C bases in the gene's exonic sequence.
    """
    try:
        import pysam
    except ImportError:
        sys.stderr.write(
            "ERROR: pysam required for --generate-gene-annot. "
            "Install with: uv pip install pysam\n")
        sys.exit(1)

    gene_exons = defaultdict(list)  # gene_id -> [(chrom, start, end)]

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            if feature != "exon":
                continue
            start = int(start) - 1
            end = int(end)

            attr_dict = {}
            for m in re.finditer(r'(\w+)\s+"([^"]*)"', attrs):
                attr_dict[m.group(1)] = m.group(2)
            gene_id = attr_dict.get("gene_id", "")
            if not gene_id:
                continue
            gene_exons[gene_id].append((chrom, start, end))

    print(f"Loaded {len(gene_exons)} genes from GTF")

    fasta = pysam.FastaFile(fasta_path)
    fasta_chroms = set(fasta.references)

    results = []
    n_done = 0
    for gene_id in sorted(gene_exons):
        exons = gene_exons[gene_id]
        by_chrom = defaultdict(list)
        for chrom, start, end in exons:
            by_chrom[chrom].append((start, end))

        total_length = 0
        total_gc = 0

        for chrom, intervals in by_chrom.items():
            if chrom not in fasta_chroms:
                continue
            intervals.sort()
            merged = []
            for start, end in intervals:
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            for start, end in merged:
                seq = fasta.fetch(chrom, start, end).upper()
                total_length += len(seq)
                total_gc += seq.count("G") + seq.count("C")

        gc_content = total_gc / total_length if total_length > 0 else np.nan
        results.append({"gene_id": gene_id, "gc_content": gc_content,
                        "gene_length": total_length})

        n_done += 1
        if n_done % 2000 == 0:
            print(f"  Processed {n_done}/{len(gene_exons)} genes")

    df = pd.DataFrame(results).set_index("gene_id")
    df.to_csv(output_path, sep="\t")
    print(f"Gene GC/length annotation written: {output_path} ({len(df)} genes)")


# =============================================================================
# Main: collect QC metrics for a cohort
# =============================================================================

def find_bam(sample, star_out_dir, bam_dir):
    """Find the pre-shrink STAR BAM, falling back to shrunk BAM."""
    pre_shrink = os.path.join(star_out_dir, f"{sample}.Aligned.sortedByCoord.out.bam")
    if os.path.exists(pre_shrink):
        return pre_shrink, "pre-shrink"
    shrunk = os.path.join(bam_dir, f"{sample}.bam")
    if os.path.exists(shrunk):
        return shrunk, "shrunk"
    return None, "missing"


def load_config_vars(config_path, cohort):
    """Load config.yml via config_get.py and return a dict of shell vars."""
    scripts_dir = os.path.dirname(os.path.abspath(__file__))
    config_get = os.path.join(scripts_dir, "config_get.py")
    if not os.path.exists(config_get):
        # repo layout: config_get.py lives in ../03_phenotyping
        config_get = os.path.join(scripts_dir, "..", "03_phenotyping", "config_get.py")
    if not os.path.exists(config_get):
        config_get = "config_get.py"  # last resort: cwd/PATH (flat deploy)

    result = subprocess.run(
        ["python3", config_get, config_path, "--cohort", cohort],
        capture_output=True, text=True, check=True
    )
    config_vars = {}
    for line in result.stdout.splitlines():
        if line.startswith("export "):
            m = re.match(r"export (\w+)=(.*)", line)
            if m:
                val = m.group(2)
                if val.startswith("'") and val.endswith("'"):
                    val = val[1:-1].replace("'\"'\"'", "'")
                config_vars[m.group(1)] = val
    return config_vars


def main():
    parser = argparse.ArgumentParser(
        description="Collect curated PicardTools QC metrics for HCP priors")
    parser.add_argument("--config", help="config.yml path (for cohort paths)")
    parser.add_argument("--cohort", help="Cohort name (key in config.yml)")
    parser.add_argument("--samples", help="Sample list file (one ID per line)")
    parser.add_argument("--output", required=True, help="Output QC metrics TSV")
    parser.add_argument("--bam-dir", help="Override star_out directory")
    parser.add_argument("--refflat", help="refFlat file for CollectRnaSeqMetrics")
    parser.add_argument("--picard-cmd", default="picard",
                        help="PicardTools command (default: picard)")
    parser.add_argument("--salmon-dir", help="Salmon output dir for subject-specific bias")
    parser.add_argument("--gene-annot", help="Gene GC/length annotation TSV")
    parser.add_argument("--generate-refflat", action="store_true",
                        help="Generate refFlat from GTF and exit")
    parser.add_argument("--generate-gene-annot", action="store_true",
                        help="Generate gene GC/length table and exit")
    parser.add_argument("--gtf", help="GTF file (for --generate-refflat/--generate-gene-annot)")
    parser.add_argument("--fasta", help="Genome FASTA (for --generate-gene-annot)")
    parser.add_argument("--raw-dir",
                        help="Directory for persisted per-sample raw Picard outputs "
                             "(default: <output dir>/raw). Raw outputs are kept so the "
                             "metric panel can be re-extracted later with --parse-only "
                             "without re-running Picard on the BAMs.")
    parser.add_argument("--parse-only", action="store_true",
                        help="Do not run Picard; re-extract metrics from the persisted "
                             "raw outputs under --raw-dir (BAMs not required).")
    args = parser.parse_args()

    # ---- Subcommand: generate refFlat ----
    if args.generate_refflat:
        if not args.gtf:
            sys.stderr.write("--gtf required for --generate-refflat\n")
            sys.exit(1)
        generate_refflat(args.gtf, args.output)
        return

    # ---- Subcommand: generate gene GC/length annotation ----
    if args.generate_gene_annot:
        if not args.gtf or not args.fasta:
            sys.stderr.write("--gtf and --fasta required for --generate-gene-annot\n")
            sys.exit(1)
        generate_gene_gc_length(args.gtf, args.fasta, args.output)
        return

    # ---- Main: collect QC metrics ----
    if not args.config or not args.cohort:
        sys.stderr.write("--config and --cohort required for QC collection\n")
        sys.exit(1)
    if not args.refflat:
        sys.stderr.write("--refflat required for QC collection (CollectRnaSeqMetrics)\n")
        sys.exit(1)

    config_vars = load_config_vars(args.config, args.cohort)
    output_base = config_vars.get("OUTPUT_BASE", "")
    cohort_dir = os.path.join(output_base, args.cohort)
    interm_dir = os.path.join(cohort_dir, "intermediate")
    star_out_dir = args.bam_dir or os.path.join(interm_dir, "star_out")
    bam_dir = os.path.join(interm_dir, "bam")

    # Load samples
    if args.samples:
        with open(args.samples) as f:
            samples = [line.strip() for line in f if line.strip()]
    else:
        samples_file = os.path.join(cohort_dir, "samples.txt")
        if os.path.exists(samples_file):
            with open(samples_file) as f:
                samples = [line.strip() for line in f if line.strip()]
        else:
            sys.stderr.write("ERROR: no samples file found. Use --samples.\n")
            sys.exit(1)

    print(f"Collecting PicardTools QC metrics for cohort {args.cohort}")
    print(f"  Samples: {len(samples)}")
    print(f"  star_out: {star_out_dir}")
    print(f"  refflat: {args.refflat}")
    print(f"  picard: {args.picard_cmd}")

    all_metrics = []
    missing_bams = []
    shrunk_bams = []

    # Raw Picard outputs persist per sample so the metric panel can be
    # re-extracted later (--parse-only) without re-running Picard on BAMs.
    raw_root = args.raw_dir or os.path.join(
        os.path.dirname(os.path.abspath(args.output)), "raw")
    if not args.parse_only:
        os.makedirs(raw_root, exist_ok=True)

    for i, sample in enumerate(samples):
        sample_dir = os.path.join(raw_root, sample)

        if args.parse_only:
            if not os.path.isdir(sample_dir):
                sys.stderr.write(
                    f"  [{i+1}/{len(samples)}] {sample}: no raw metrics at "
                    f"{sample_dir}, skipping\n")
                missing_bams.append(sample)
                continue
            print(f"  [{i+1}/{len(samples)}] {sample} (parse-only)")
            sample_metrics = {"sample": sample}
            sample_metrics.update(
                collect_alignment_summary(args.picard_cmd, None, sample_dir, sample,
                                          parse_only=True))
            sample_metrics.update(
                collect_insert_size(args.picard_cmd, None, sample_dir, sample,
                                    parse_only=True))
            sample_metrics.update(
                collect_rna_seq_metrics(args.picard_cmd, None, sample_dir, sample,
                                        args.refflat, parse_only=True))
            sample_metrics.update(
                collect_gc_bias(args.picard_cmd, None, sample_dir, sample,
                                parse_only=True))
            sample_metrics.update(
                mark_duplicates(args.picard_cmd, None, sample_dir, sample,
                                parse_only=True))
            all_metrics.append(sample_metrics)
            continue

        bam_path, bam_type = find_bam(sample, star_out_dir, bam_dir)

        if bam_path is None:
            sys.stderr.write(f"  [{i+1}/{len(samples)}] {sample}: BAM not found, skipping\n")
            missing_bams.append(sample)
            continue

        if bam_type == "shrunk":
            shrunk_bams.append(sample)

        print(f"  [{i+1}/{len(samples)}] {sample} ({bam_type} BAM)")

        os.makedirs(sample_dir, exist_ok=True)
        sample_metrics = {"sample": sample}
        sample_metrics.update(
            collect_alignment_summary(args.picard_cmd, bam_path, sample_dir, sample))
        sample_metrics.update(
            collect_insert_size(args.picard_cmd, bam_path, sample_dir, sample))
        sample_metrics.update(
            collect_rna_seq_metrics(args.picard_cmd, bam_path, sample_dir, sample, args.refflat))
        sample_metrics.update(
            collect_gc_bias(args.picard_cmd, bam_path, sample_dir, sample,
                            fasta=args.fasta))
        sample_metrics.update(
            mark_duplicates(args.picard_cmd, bam_path, sample_dir, sample))
        all_metrics.append(sample_metrics)

    if not all_metrics:
        sys.stderr.write("ERROR: no samples had collectible metrics\n")
        sys.exit(1)

    qc_df = pd.DataFrame(all_metrics).set_index("sample")

    # ---- Subject-specific GC and gene-length bias ----
    if args.salmon_dir and args.gene_annot:
        print("\nComputing subject-specific GC and gene-length bias...")
        gene_annot = pd.read_csv(args.gene_annot, sep="\t", index_col="gene_id")
        tx_to_gene = None
        if args.gtf and os.path.exists(args.gtf):
            tx_to_gene = load_tx_to_gene(args.gtf)
            print(f"  Loaded transcript->gene map: {len(tx_to_gene)} transcripts")
        else:
            sys.stderr.write(
                "  WARN: no --gtf given; falling back to heuristic gene-ID "
                "parsing of Salmon Names (may not match GTF gene_id space)\n")
        tpm_matrix = load_tpm_matrix(args.salmon_dir, samples, tx_to_gene=tx_to_gene)
        if not tpm_matrix.empty:
            bias_df = compute_subject_specific_bias(tpm_matrix, gene_annot)
            qc_df = qc_df.join(bias_df, how="left")
        else:
            sys.stderr.write("  WARN: could not load TPM matrix, skipping subject bias\n")
    else:
        print("\nSkipping subject-specific bias (need --salmon-dir and --gene-annot)")

    # ---- Impute missing values (column median) ----
    numeric_cols = qc_df.select_dtypes(include=[np.number]).columns
    for col in numeric_cols:
        n_missing = qc_df[col].isna().sum()
        if n_missing > 0:
            median_val = qc_df[col].median()
            if pd.isna(median_val):
                median_val = 0.0
            qc_df[col] = qc_df[col].fillna(median_val)
            if n_missing > len(qc_df) * 0.5:
                sys.stderr.write(
                    f"  WARN: {col} had {n_missing}/{len(qc_df)} missing values "
                    f"(>50%), imputed with median\n")

    # ---- Report ----
    print(f"\nQC metrics collected: {len(qc_df)} samples x {len(qc_df.columns)} metrics")
    if missing_bams:
        print(f"  Missing BAMs (skipped): {len(missing_bams)}")
        for s in missing_bams:
            print(f"    {s}")
    if shrunk_bams:
        print(f"  Shrunk BAMs (some metrics may be unreliable): {len(shrunk_bams)}")
    print(f"  Metrics: {', '.join(qc_df.columns)}")

    # ---- Write output ----
    qc_df.to_csv(args.output, sep="\t")
    print(f"\nOutput: {args.output}")


if __name__ == "__main__":
    main()
