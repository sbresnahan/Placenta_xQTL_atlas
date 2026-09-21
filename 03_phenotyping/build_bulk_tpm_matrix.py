#!/usr/bin/env python3
"""Script 2: Build a gene-symbol x sample TPM matrix from Salmon quant.sf files.

Aggregates per-sample Salmon transcript TPMs to gene-level TPMs using the
gene_name attribute from the normalized HPLRv2 GTF (produced by prep_PANTRY_gtf.R).
This matrix is the bulk mixture input for MuSiC deconvolution.

Using gene symbols (not Ensembl gene_id) maximizes overlap with the Campbell
et al. 2023 scRNA-seq reference (GSE182381), which uses GENCODE gene symbols.

Inputs:
  - Salmon quant.sf files: <salmon_dir>/<sample>/quant.sf
    Columns: Name, Length, EffectiveLength, TPM, NumReads
  - Normalized HPLRv2 GTF (has gene_name + transcript_id attributes on
    transcript feature lines, added by prep_PANTRY_gtf.R)
  - Samples file: one sample ID per line

Output:
  - TSV: genes (symbols) x samples, TPM values, non-log linear space

Usage:
  python3 build_bulk_tpm_matrix.py \
    --gtf        /path/to/HPLRv2.0.annotated.PANTRY.gtf \
    --salmon-dir /path/to/cohort/salmon_expression \
    --samples    /path/to/samples.txt \
    --output     /path/to/cohort_bulk_tpm_symbols.tsv

  # Or via config.yml (seadragon):
  python3 build_bulk_tpm_matrix.py \
    --config     config.yml \
    --cohort     cohort1 \
    --output     /path/to/cohort_bulk_tpm_symbols.tsv
"""

import argparse
import sys
import os
from pathlib import Path

import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# GTF parsing — extract transcript_id -> gene_name map
# ---------------------------------------------------------------------------

def build_tx_to_gene_name_map(gtf_path: Path) -> pd.DataFrame:
    """Parse GTF and build a transcript_id -> gene_name mapping table.

    Uses the 'transcript' feature lines, which carry both transcript_id and
    (after prep_PANTRY_gtf.R normalization) gene_name attributes.

    For transcripts missing gene_name, falls back to gene_id with the version
    suffix stripped (e.g. ENSG00000123456.1 -> ENSG00000123456).

    Returns a DataFrame with columns: transcript_id, gene_name
    """
    records = []
    with open(gtf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'transcript':
                continue
            attrs = fields[8]
            # Parse the GTF attribute column (key "value"; pairs separated by ';')
            attr_dict = {}
            for pair in attrs.split(';'):
                pair = pair.strip()
                if not pair:
                    continue
                # Handle both: key "value"  and  key "value";
                parts = pair.split(' ', 1)
                if len(parts) != 2:
                    continue
                key = parts[0]
                val = parts[1].strip().strip('"')
                attr_dict[key] = val
            tx_id = attr_dict.get('transcript_id')
            gene_name = attr_dict.get('gene_name')
            gene_id = attr_dict.get('gene_id')
            if tx_id is None:
                continue
            if gene_name is None or gene_name == '':
                # Fallback: strip version suffix from gene_id
                if gene_id is not None:
                    gene_name = gene_id.split('.')[0]
                else:
                    continue
            records.append({'transcript_id': tx_id, 'gene_name': gene_name})

    df = pd.DataFrame(records)
    # Drop duplicate transcript_ids (keep first)
    df = df.drop_duplicates(subset='transcript_id', keep='first')
    print(f"  GTF: {len(df)} transcript_id -> gene_name mappings")
    if df['gene_name'].duplicated().any():
        n_dup = df['gene_name'].duplicated().sum()
        print(f"  GTF: {n_dup} transcripts share gene_name with another transcript (expected)")
    return df


# ---------------------------------------------------------------------------
# Salmon quant.sf aggregation
# ---------------------------------------------------------------------------

def aggregate_sample_tpm(quant_sf: Path, tx_map: pd.DataFrame) -> pd.Series:
    """Read one sample's quant.sf and aggregate transcript TPMs to gene-level.

    Sums TPMs across all transcripts mapping to the same gene_name.
    Returns a Series indexed by gene_name with TPM values.
    """
    df = pd.read_csv(quant_sf, sep='\t')
    # Salmon quant.sf: Name, Length, EffectiveLength, TPM, NumReads
    # 'Name' is the transcript ID
    df = df[['Name', 'TPM']].rename(columns={'Name': 'transcript_id', 'TPM': 'tpm'})
    # Merge with transcript->gene map
    df = df.merge(tx_map, on='transcript_id', how='inner')
    # Aggregate: sum TPM across transcripts per gene
    gene_tpm = df.groupby('gene_name')['tpm'].sum()
    return gene_tpm


def build_bulk_matrix(salmon_dir: Path, samples: list, tx_map: pd.DataFrame) -> pd.DataFrame:
    """Build genes x samples TPM matrix from per-sample Salmon quant.sf files.

    Args:
        salmon_dir: Directory containing <sample>/quant.sf subdirectories.
        samples: List of sample IDs to process.
        tx_map: DataFrame with transcript_id, gene_name columns.

    Returns:
        DataFrame: rows = gene_name, columns = sample IDs, values = TPM.
    """
    matrices = []
    missing = []
    for i, sample in enumerate(samples):
        quant_path = salmon_dir / sample / 'quant.sf'
        if not quant_path.exists():
            missing.append(sample)
            continue
        gene_tpm = aggregate_sample_tpm(quant_path, tx_map)
        gene_tpm.name = sample
        matrices.append(gene_tpm)
        if (i + 1) % 50 == 0 or (i + 1) == len(samples):
            print(f"  Processed {i+1}/{len(samples)} samples")

    if missing:
        print(f"  WARNING: {len(missing)} samples missing quant.sf: {missing[:5]}{'...' if len(missing)>5 else ''}")

    if not matrices:
        raise ValueError(f"No quant.sf files found in {salmon_dir}")

    # Outer join to keep all genes (fill missing with 0)
    bulk = pd.concat(matrices, axis=1)
    bulk = bulk.fillna(0.0)

    # Handle duplicate gene names (sum them — standard for gene-level aggregation)
    if bulk.index.duplicated().any():
        n_before = len(bulk)
        bulk = bulk.groupby(level=0).sum()
        print(f"  Summed duplicate gene names: {n_before} -> {len(bulk)} unique genes")

    return bulk


# ---------------------------------------------------------------------------
# Config.yml integration (optional, for seadragon)
# ---------------------------------------------------------------------------

def read_config(config_path: str, cohort: str = None):
    """Minimal config.yml reader (mirrors config_get.py logic).

    Returns a nested dict. If cohort is given, returns cohort-specific keys
    merged with top-level keys.
    """
    import yaml
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    if cohort and cohort in cfg.get('cohorts', {}):
        cohort_cfg = cfg['cohorts'][cohort]
        # Merge: cohort-specific overrides top-level where keys overlap
        merged = {**cfg, **cohort_cfg}
        return merged
    return cfg


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build gene-symbol x sample TPM matrix from Salmon quant.sf for MuSiC deconvolution"
    )
    parser.add_argument('--gtf', type=str, default=None,
                        help='Path to normalized HPLRv2 GTF (with gene_name attribute)')
    parser.add_argument('--salmon-dir', type=str, default=None,
                        help='Directory containing <sample>/quant.sf subdirectories')
    parser.add_argument('--samples', type=str, default=None,
                        help='File with one sample ID per line')
    parser.add_argument('--output', type=str, required=True,
                        help='Output TSV path (genes x samples, TPM values)')
    # Config-based alternative (seadragon)
    parser.add_argument('--config', type=str, default=None,
                        help='config.yml path (alternative to --gtf/--salmon-dir/--samples)')
    parser.add_argument('--cohort', type=str, default=None,
                        help='Cohort name (used with --config)')
    args = parser.parse_args()

    # Resolve paths from config if --config given
    if args.config:
        cfg = read_config(args.config, args.cohort)
        gtf_path = Path(cfg['normalized_gtf'])
        salmon_dir = Path(cfg.get('salmon_expression_dir',
                                  os.path.join(cfg.get('output_base', ''), args.cohort,
                                               'salmon_expression')))
        samples_file = Path(cfg.get('samples_file',
                                    os.path.join(cfg.get('output_base', ''), args.cohort,
                                                 'samples.txt')))
    else:
        if not (args.gtf and args.salmon_dir and args.samples):
            parser.error("Either --config/--cohort or --gtf/--salmon-dir/--samples must be provided")
        gtf_path = Path(args.gtf)
        salmon_dir = Path(args.salmon_dir)
        samples_file = Path(args.samples)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"=== Building bulk TPM matrix ===")
    print(f"GTF:        {gtf_path}")
    print(f"Salmon dir: {salmon_dir}")
    print(f"Samples:    {samples_file}")
    print(f"Output:     {output_path}")
    print()

    # --- Load samples ---
    with open(samples_file) as fh:
        samples = [line.strip() for line in fh if line.strip()]
    print(f"Samples to process: {len(samples)}")

    # --- Build transcript -> gene_name map ---
    print(f"\nParsing GTF for transcript_id -> gene_name map...")
    tx_map = build_tx_to_gene_name_map(gtf_path)

    # --- Build bulk matrix ---
    print(f"\nAggregating Salmon TPMs to gene level...")
    bulk = build_bulk_matrix(salmon_dir, samples, tx_map)
    print(f"\nBulk matrix: {bulk.shape[0]} genes x {bulk.shape[1]} samples")

    # --- Report gene overlap with reference ---
    # The reference (GSE182381) uses GENCODE gene symbols. We report overlap
    # if a reference gene list is available alongside.
    ref_genes_path = output_path.parent / 'reference_gene_symbols.txt'
    if ref_genes_path.exists():
        ref_genes = set(open(ref_genes_path).read().splitlines())
        bulk_genes = set(bulk.index)
        shared = bulk_genes & ref_genes
        print(f"\nGene overlap with reference:")
        print(f"  Bulk genes:    {len(bulk_genes)}")
        print(f"  Reference genes: {len(ref_genes)}")
        print(f"  Shared:        {len(shared)} ({100*len(shared)/len(ref_genes):.1f}% of reference)")
    else:
        print(f"\n(No reference_gene_symbols.txt found — skipping overlap report)")
        print(f"  To enable: extract gene symbols from GSE182381_reference_sample.txt.gz")
        print(f"  and save as {ref_genes_path}")

    # --- Write output ---
    print(f"\nWriting {output_path}...")
    # Write TSV: first column = gene symbol, then sample columns
    bulk_out = bulk.copy()
    bulk_out.index.name = 'gene_symbol'
    bulk_out.to_csv(output_path, sep='\t', float_format='%.6g')
    print(f"Done. Size: {output_path.stat().st_size / 1e6:.1f} MB")
    print(f"\n=== Bulk TPM matrix complete ===")


if __name__ == '__main__':
    main()
