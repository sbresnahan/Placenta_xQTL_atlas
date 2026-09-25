#!/usr/bin/env python3
"""
23_prepare_intersection.py — Sample intersection + ID harmonization for tensorQTL

For each ancestry stratum:
  1. Identifies samples present in ALL data sources (expression BED, pooled
     pgen, HCP factors, deconvolution proportions, metadata).
  2. Filters the genotype pgen to intersection samples (via plink2 --keep),
     applying a minor-allele-count floor (--mac, default 5) so low-carrier
     variants (the round-1 artifact class) never enter QTL testing.
  3. Renames sample columns in phenotype/covariate files from rnaseq_id to
     array_id (matching the pgen).
  4. Reports intersection counts.

Usage:
  python3 23_prepare_intersection.py \
      --config <config.yml> \
      --metadata <placenta_QTL_cohort_metadata.tsv> \
      --ancestries EAS EUR \
      --output-dir <qtl_inputs dir>

Or via LSF wrapper (23_prepare_intersection.sh).
"""

import argparse
import os
import subprocess
import sys
import pandas as pd
import numpy as np


def load_config(config_path, scripts_dir):
    """Load config.yml via config_get.py."""
    config_get = os.path.join(scripts_dir, "config_get.py")
    if not os.path.exists(config_get):
        # repo layout: config_get.py lives in ../03_phenotyping
        config_get = os.path.join(scripts_dir, "..", "03_phenotyping", "config_get.py")
    result = subprocess.run(
        ["python3", config_get, config_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        sys.exit(f"ERROR: config_get.py failed: {result.stderr}")
    config = {}
    for line in result.stdout.strip().split('\n'):
        if '=' in line and not line.strip().startswith('#'):
            key, val = line.split('=', 1)
            key = key.strip()
            # config_get.py outputs "export KEY='value'" — strip the prefix
            if key.startswith('export '):
                key = key[len('export '):]
            config[key] = val.strip().strip("'\"")
    return config


def read_bed_samples(bed_path):
    """Read sample IDs from the header line of a BED file."""
    with open(bed_path, 'r') as f:
        header = f.readline().strip()
    cols = header.split('\t')
    # BED meta columns: #chr, start, end, phenotype_id
    meta = {'#chr', 'start', 'end', 'phenotype_id', 'chr'}
    return [c for c in cols if c not in meta]


def read_psam(psam_path):
    """Read a plink2 .psam file. Returns DataFrame with at least FID and IID.

    The psam header line starts with '#' (e.g. '#FID\\tIID\\t...'), which
    pandas comment='#' would drop — causing the first data row to be
    misread as the header. Parse the header explicitly instead.
    """
    with open(psam_path) as f:
        header_line = f.readline().strip()
    df = pd.read_csv(psam_path, sep='\t', comment='#', header=None)
    if header_line.startswith('#'):
        cols = header_line.lstrip('#').split('\t')
        if len(cols) == df.shape[1]:
            df.columns = cols
    if 'IID' not in df.columns:
        # Fallback: assume FID, IID, ... column order
        df.columns = ['FID', 'IID'] + [f'extra_{i}' for i in range(df.shape[1] - 2)]
    if 'FID' not in df.columns:
        # plink2 uses '0' as FID when the .psam has no FID column
        df['FID'] = '0'
    return df


def read_tsv_samples(tsv_path, id_col='sample_id'):
    """Read sample IDs from a TSV file."""
    df = pd.read_csv(tsv_path, sep='\t')
    if id_col in df.columns:
        return df[id_col].tolist()
    # Try first column
    return df.iloc[:, 0].tolist()


def rename_bed_columns(bed_path, out_path, id_map):
    """Rename sample columns in a BED file using rnaseq_id -> array_id map.
    
    id_map: dict {rnaseq_id: array_id}
    """
    with open(bed_path, 'r') as f:
        header = f.readline()
        data = f.read()
    
    cols = header.strip().split('\t')
    meta = {'#chr', 'start', 'end', 'phenotype_id', 'chr'}
    new_cols = []
    rename_count = 0
    drop_indices = []
    
    for i, c in enumerate(cols):
        if c in meta:
            new_cols.append(c)
        elif c in id_map:
            new_cols.append(id_map[c])
            rename_count += 1
        else:
            # Sample not in metadata — mark for dropping
            new_cols.append(c)
            drop_indices.append(i)
    
    if drop_indices:
        # Drop columns not in the intersection
        keep_indices = [i for i in range(len(cols)) if i not in drop_indices]
        new_cols = [new_cols[i] for i in keep_indices]
        # Filter data lines
        new_data_lines = []
        for line in data.strip().split('\n'):
            fields = line.split('\t')
            new_data_lines.append('\t'.join(fields[i] for i in keep_indices))
        new_data = '\n'.join(new_data_lines) + '\n'
    else:
        new_data = data
    
    new_header = '\t'.join(new_cols) + '\n'
    
    with open(out_path, 'w') as f:
        f.write(new_header)
        f.write(new_data)
    
    return rename_count, len(drop_indices)


def rename_tsv_samples(tsv_path, out_path, id_col, id_map, keep_samples=None):
    """Rename sample IDs in a TSV file (e.g., HCP factors, deconvolution).
    
    If the file has samples as rows (with an ID column), rename the ID column.
    If samples are columns, rename the columns.
    
    id_map: dict {rnaseq_id: array_id}
    keep_samples: set of array_ids to keep (if None, keep all mapped)
    """
    df = pd.read_csv(tsv_path, sep='\t')
    
    if id_col in df.columns:
        # Samples are rows — rename the ID column values
        df[id_col] = df[id_col].map(id_map)
        unmapped = df[id_col].isna().sum()
        df = df.dropna(subset=[id_col])
        if keep_samples is not None:
            df = df[df[id_col].isin(keep_samples)]
        df.to_csv(out_path, sep='\t', index=False)
        return len(df), unmapped
    else:
        # Samples might be columns — try renaming columns
        new_cols = []
        rename_count = 0
        drop_count = 0
        for c in df.columns:
            if c in id_map:
                new_cols.append(id_map[c])
                rename_count += 1
            elif c in keep_samples if keep_samples else True:
                new_cols.append(c)
            else:
                drop_count += 1
        df.columns = new_cols
        if keep_samples is not None:
            sample_cols = [c for c in new_cols if c in keep_samples]
            non_sample = [c for c in new_cols if c not in keep_samples]
            df = df[non_sample + sample_cols]
        df.to_csv(out_path, sep='\t', index=False)
        return rename_count, drop_count


def main():
    parser = argparse.ArgumentParser(description="Sample intersection + ID harmonization")
    parser.add_argument("--config", required=True, help="Path to config.yml")
    parser.add_argument("--metadata", required=True,
                        help="Path to placenta_QTL_cohort_metadata.tsv")
    parser.add_argument("--ancestries", default="EAS EUR",
                        help="Space-separated ancestry labels")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory (default: ${OUTPUT_BASE}/qtl_inputs)")
    parser.add_argument("--scripts-dir", default=None,
                        help="Scripts directory (for config_get.py)")
    parser.add_argument("--mac", type=int, default=5,
                        help="Minor-allele-count floor for the intersection pgen "
                             "(default: 5; 0 disables). Prevents low-carrier artifacts.")
    args = parser.parse_args()
    
    scripts_dir = args.scripts_dir or os.path.dirname(os.path.abspath(__file__))
    config = load_config(args.config, scripts_dir)
    output_base = config.get('OUTPUT_BASE', '')
    
    if not output_base:
        sys.exit("ERROR: OUTPUT_BASE not found in config")
    
    output_dir = args.output_dir or os.path.join(output_base, "qtl_inputs")
    os.makedirs(output_dir, exist_ok=True)
    
    geno_dir = "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/genotypes"
    hcp_dir = os.path.join(output_base, "hcp", "hcp_factors")
    deconv_dir = os.path.join(output_base, "deconvolution")
    
    # ---- Load metadata ----
    meta = pd.read_csv(args.metadata, sep='\t')
    print(f"Metadata: {len(meta)} rows")
    print(f"  Columns: {list(meta.columns)}")
    
    # Build rnaseq_id -> array_id map
    id_map = dict(zip(meta['rnaseq_id'], meta['array_id']))
    print(f"  ID map: {len(id_map)} entries")
    
    ancestries = args.ancestries.split()
    
    for anc in ancestries:
        print(f"\n{'='*60}")
        print(f"Ancestry: {anc}")
        print(f"{'='*60}")
        
        # ---- Load all sample sets ----
        # 1. Expression BED
        expr_bed = os.path.join(hcp_dir, f"{anc}_combat_int_expression.bed")
        if not os.path.exists(expr_bed):
            print(f"  ERROR: expression BED not found: {expr_bed}")
            continue
        expr_samples = read_bed_samples(expr_bed)
        print(f"  Expression BED: {len(expr_samples)} samples")
        
        # 2. Pooled pgen
        psam_path = os.path.join(geno_dir, f"{anc}_pooled.psam")
        if not os.path.exists(psam_path):
            print(f"  ERROR: psam not found: {psam_path}")
            continue
        psam = read_psam(psam_path)
        geno_samples = psam['IID'].tolist()
        print(f"  Genotype psam: {len(geno_samples)} samples (IID), FID example: {psam['FID'].iloc[0]}")
        
        # 3. HCP factors
        hcp_path = os.path.join(hcp_dir, f"{anc}_hcp_factors.tsv")
        if not os.path.exists(hcp_path):
            print(f"  ERROR: HCP factors not found: {hcp_path}")
            continue
        # HCP factors: rows = HCP_1..15, columns = samples (rnaseq_id)
        hcp_df = pd.read_csv(hcp_path, sep='\t', index_col=0)
        hcp_samples = list(hcp_df.columns)
        print(f"  HCP factors: {len(hcp_samples)} samples")
        
        # 4. Deconvolution
        deconv_path = os.path.join(deconv_dir, "all_cohorts_cell_proportions_collapsed.tsv")
        if not os.path.exists(deconv_path):
            print(f"  WARN: deconvolution not found: {deconv_path}")
            deconv_samples = []
        else:
            deconv_df = pd.read_csv(deconv_path, sep='\t')
            # Filter to this ancestry via metadata
            anc_meta = meta[meta['ancestry'] == anc]
            anc_rnaseq = set(anc_meta['rnaseq_id'])
            deconv_df = deconv_df[deconv_df['sample'].isin(anc_rnaseq)]
            deconv_samples = deconv_df['sample'].tolist()
            print(f"  Deconvolution: {len(deconv_samples)} samples (ancestry-filtered)")
        
        # 5. Metadata for this ancestry
        anc_meta = meta[meta['ancestry'] == anc]
        meta_rnaseq = set(anc_meta['rnaseq_id'])
        meta_array = set(anc_meta['array_id'])
        print(f"  Metadata ({anc}): {len(anc_meta)} samples")
        
        # ---- Compute intersection ----
        # RNA side: samples in expression, HCP, deconvolution, and metadata
        rna_intersection = set(expr_samples) & set(hcp_samples) & meta_rnaseq
        if deconv_samples:
            rna_intersection = rna_intersection & set(deconv_samples)
        
        # Map RNA intersection to array_id
        rna_to_array = {r: id_map[r] for r in rna_intersection if r in id_map}
        array_intersection = set(rna_to_array.values())
        
        # DNA side: intersection with genotype
        final_array = array_intersection & set(geno_samples)
        
        # Back-map to rnaseq_id for filtering
        final_rnaseq = {v: k for k, v in rna_to_array.items()}
        final_rnaseq = {final_rnaseq[a] for a in final_array if a in final_rnaseq}
        
        print(f"\n  Intersection:")
        print(f"    RNA samples (expr ∩ HCP ∩ deconv ∩ meta): {len(rna_intersection)}")
        print(f"    Mapped to array_id: {len(array_intersection)}")
        print(f"    DNA samples (genotype): {len(geno_samples)}")
        print(f"    FINAL intersection (RNA ∩ DNA): {len(final_array)}")
        
        # Report drops
        rna_only = array_intersection - set(geno_samples)
        dna_only = set(geno_samples) - array_intersection
        unmapped_rna = rna_intersection - set(rna_to_array.keys())
        print(f"    Dropped (RNA but no genotype): {len(rna_only)}")
        print(f"    Dropped (genotype but no RNA): {len(dna_only)}")
        print(f"    Dropped (RNA but no array_id mapping): {len(unmapped_rna)}")
        
        if len(final_array) < 10:
            print(f"  ERROR: Too few samples ({len(final_array)}) for QTL mapping. Skipping.")
            continue
        
        # ---- Filter genotype pgen to intersection ----
        qtl_prefix = os.path.join(output_dir, f"{anc}_qtl")
        keep_psam = os.path.join(output_dir, f"{anc}_keep.txt")
        # Write actual FID/IID pairs from the psam (plink2 --keep matches both)
        keep_df = psam[psam['IID'].isin(final_array)][['FID', 'IID']]
        keep_df.to_csv(keep_psam, sep='\t', index=False, header=False)
        if len(keep_df) != len(final_array):
            print(f"  WARN: keep file has {len(keep_df)} rows but intersection has {len(final_array)}")
        
        print(f"\n  Filtering genotype pgen to {len(final_array)} samples...")
        mac_flag = f"--mac {args.mac} " if args.mac > 0 else ""
        cmd = (
            f"plink2 --pgen {geno_dir}/{anc}_pooled.pgen "
            f"--pvar {geno_dir}/{anc}_pooled.pvar "
            f"--psam {geno_dir}/{anc}_pooled.psam "
            f"--keep {keep_psam} "
            f"{mac_flag}"
            f"--make-pgen --threads 4 "
            f"--out {qtl_prefix}"
        )
        print(f"    {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ERROR: plink2 failed: {result.stderr}")
            continue
        n_var = sum(1 for line in open(f"{qtl_prefix}.pvar") if not line.startswith('#'))
        print(f"    Written: {qtl_prefix}.pgen/pvar/psam ({n_var} variants"
              + (f", MAC >= {args.mac}" if args.mac > 0 else "") + ")")
        
        # ---- Harmonize expression BED ----
        # Build filtered id_map (only final intersection)
        final_id_map = {r: id_map[r] for r in final_rnaseq if r in id_map}
        final_array_set = final_array
        
        expr_out = os.path.join(output_dir, f"{anc}_expression_harmonized.bed")
        print(f"\n  Harmonizing expression BED...")
        renamed, dropped = rename_bed_columns(expr_bed, expr_out, final_id_map)
        print(f"    Renamed {renamed} columns, dropped {dropped} non-intersection samples")
        print(f"    Written: {expr_out}")
        
        # ---- Harmonize HCP factors ----
        hcp_out = os.path.join(output_dir, f"{anc}_hcp_factors_harmonized.tsv")
        print(f"  Harmonizing HCP factors...")
        # HCP factors: rows = covariates, columns = samples
        hcp_df = pd.read_csv(hcp_path, sep='\t', index_col=0)
        # Rename columns
        new_cols = {}
        drop_cols = []
        for c in hcp_df.columns:
            if c in final_id_map:
                new_cols[c] = final_id_map[c]
            else:
                drop_cols.append(c)
        hcp_df = hcp_df.drop(columns=drop_cols)
        hcp_df = hcp_df.rename(columns=new_cols)
        hcp_df.to_csv(hcp_out, sep='\t')
        print(f"    Renamed {len(new_cols)} columns, dropped {len(drop_cols)} samples")
        print(f"    Written: {hcp_out}")
        
        # ---- Harmonize deconvolution ----
        if os.path.exists(deconv_path):
            deconv_out = os.path.join(output_dir, f"{anc}_deconvolution_harmonized.tsv")
            print(f"  Harmonizing deconvolution proportions...")
            deconv_df = pd.read_csv(deconv_path, sep='\t')
            # Filter to ancestry + intersection
            deconv_df = deconv_df[deconv_df['sample'].isin(final_id_map.keys())]
            # Rename sample column
            deconv_df['sample'] = deconv_df['sample'].map(final_id_map)
            deconv_df = deconv_df.rename(columns={'sample': 'sample_id'})
            deconv_df.to_csv(deconv_out, sep='\t', index=False)
            print(f"    {len(deconv_df)} samples, written: {deconv_out}")
        
        # ---- Save metadata subset ----
        meta_out = os.path.join(output_dir, f"{anc}_metadata.tsv")
        anc_meta_final = anc_meta[anc_meta['array_id'].isin(final_array)]
        anc_meta_final.to_csv(meta_out, sep='\t', index=False)
        print(f"  Metadata subset: {len(anc_meta_final)} samples, written: {meta_out}")
        
        print(f"\n  Done: {anc} — {len(final_array)} samples ready for QTL mapping")


if __name__ == "__main__":
    main()
