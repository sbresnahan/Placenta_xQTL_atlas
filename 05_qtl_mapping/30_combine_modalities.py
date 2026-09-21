#!/usr/bin/env python3
"""
30_combine_modalities.py — build the PANTRY-style combined cross-modality BED

Per ancestry, concatenates ALL modalities' harmonized phenotype BEDs
(expression + the 7 non-expression modalities) into one combined BED, with:

  1. phenotype_id namespaced by modality:  {modality}__{original_id}
     (required: expression and stability both use bare ENSG IDs and would
     collide otherwise; namespacing also makes every phenotype's modality
     self-identifying in downstream results).
  2. Sample columns intersected across all modalities (a sample missing
     from ANY modality is dropped from the combined BED).
  3. A phenotype->gene groups file for tensorQTL grouped mapping (group_s):
     every phenotype of a gene, across all modalities, forms ONE group
     (PANTRY convention). Gene ID = original phenotype ID up to the first
     '__' (bare ENSG IDs map to themselves).
  4. A modality sidecar TSV for downstream reporting:
     namespaced phenotype_id -> modality, gene_id, original_id.

Inputs (in --qtl-dir), per ancestry ANC:
  {ANC}_expression_harmonized.bed
  {ANC}_{modality}_harmonized.bed   for the 7 non-expression modalities

Outputs (in --qtl-dir), per ancestry:
  {ANC}_combined_harmonized.bed            (input to 27_run_tensorqtl.sh)
  {ANC}_combined.phenotype_groups.txt      (no header: phenotype_id, gene_id)
  {ANC}_combined.phenotype_modalities.tsv  (header: phenotype_id, modality,
                                            gene_id, original_id)

Usage:
  python3 30_combine_modalities.py --qtl-dir <qtl_inputs dir> \
      --ancestries EAS EUR
"""

import argparse
import os
import sys

import pandas as pd

DEFAULT_MODALITIES = ["expression", "alt_polyA", "alt_TSS", "intron_retention",
                      "isoforms", "RNA_editing", "splicing", "stability"]


def load_bed(path, modality):
    """Read one harmonized BED; return (df, meta_cols, sample_cols).
    Requires the 4th column to be phenotype_id; preserves the file's own
    meta-column names (e.g. '#chr' vs 'chr')."""
    df = pd.read_csv(path, sep='\t')
    df.iloc[:, 0] = df.iloc[:, 0].astype(str)  # chromosome column as string
    meta_cols = list(df.columns[:4])
    if meta_cols[3] != 'phenotype_id':
        sys.exit(f"ERROR: {path}: 4th column must be 'phenotype_id', "
                 f"found '{meta_cols[3]}'")
    sample_cols = list(df.columns[4:])
    if not sample_cols:
        sys.exit(f"ERROR: {path}: no sample columns found")
    # Phenotype values must be numeric
    non_numeric = [c for c in sample_cols
                   if not pd.api.types.is_numeric_dtype(df[c])]
    if non_numeric:
        sys.exit(f"ERROR: {path}: non-numeric sample columns: "
                 f"{non_numeric[:5]}")
    # Within-modality duplicate phenotype IDs are never valid
    dups = df['phenotype_id'][df['phenotype_id'].duplicated()].unique()
    if len(dups):
        sys.exit(f"ERROR: {path}: {len(dups)} duplicate phenotype_id values "
                 f"(e.g. {list(dups[:3])})")
    return df, meta_cols, sample_cols


def main():
    parser = argparse.ArgumentParser(
        description="Build combined cross-modality phenotype BED (PANTRY-style)")
    parser.add_argument("--qtl-dir", required=True, help="QTL inputs directory")
    parser.add_argument("--ancestries", nargs='+', default=["EAS", "EUR"],
                        help="Ancestry labels (default: EAS EUR)")
    parser.add_argument("--modalities", nargs='+', default=DEFAULT_MODALITIES,
                        help="Modalities to combine (default: expression + all 7)")
    parser.add_argument("--label", default="combined",
                        help="Output label (default: combined)")
    args = parser.parse_args()

    qtl_dir = args.qtl_dir

    for anc in args.ancestries:
        print("=" * 60)
        print(f"[{anc}] Combining {len(args.modalities)} modalities")
        print("=" * 60)

        beds = {}
        meta_ref = None
        sample_sets = {}
        for mod in args.modalities:
            path = os.path.join(qtl_dir, f"{anc}_{mod}_harmonized.bed")
            if not os.path.exists(path):
                sys.exit(f"ERROR: required BED not found: {path}\n"
                         f"  (run 26_harmonize_modalities.py first; the "
                         f"combined run needs ALL modalities)")
            df, meta_cols, sample_cols = load_bed(path, mod)
            if meta_ref is None:
                meta_ref = meta_cols
            elif meta_cols != meta_ref:
                sys.exit(f"ERROR: {path}: meta columns {meta_cols} differ "
                         f"from {meta_ref} — harmonize the BED headers first")
            beds[mod] = df
            sample_sets[mod] = set(sample_cols)
            print(f"  {mod:18s} {df.shape[0]:>7,} phenotypes x "
                  f"{len(sample_cols)} samples")

        # ---- Sample intersection across ALL modalities ----
        common = set.intersection(*sample_sets.values())
        # Keep the column order of the first modality, filtered to common
        first_mod = args.modalities[0]
        common_ordered = [s for s in beds[first_mod].columns[4:]
                          if s in common]
        print(f"\n  Samples per modality: "
              f"{ {m: len(s) for m, s in sample_sets.items()} }")
        print(f"  Common across all modalities: {len(common_ordered)}")
        if len(common_ordered) < 10:
            sys.exit(f"ERROR: too few common samples "
                     f"({len(common_ordered)}) for QTL mapping")

        # ---- Namespace phenotype IDs + concatenate ----
        combined = []
        sidecar = []
        for mod in args.modalities:
            df = beds[mod]
            original_ids = df['phenotype_id'].astype(str)
            df = df[meta_ref + common_ordered].copy()
            df['phenotype_id'] = mod + '__' + original_ids.values
            combined.append(df)
            gene_ids = original_ids.str.split('__').str[0]
            sidecar.append(pd.DataFrame({
                'phenotype_id': df['phenotype_id'].values,
                'modality': mod,
                'gene_id': gene_ids.values,
                'original_id': original_ids.values,
            }))

        combined_df = pd.concat(combined, ignore_index=True)
        sidecar_df = pd.concat(sidecar, ignore_index=True)

        # ---- Global duplicate check (post-namespacing) ----
        dups = combined_df['phenotype_id'][
            combined_df['phenotype_id'].duplicated()].unique()
        if len(dups):
            sys.exit(f"ERROR: {len(dups)} duplicate phenotype_ids remain "
                     f"after namespacing (e.g. {list(dups[:3])}) — this "
                     f"should not happen; inspect the input BEDs")

        # ---- Write outputs ----
        bed_path = os.path.join(qtl_dir, f"{anc}_{args.label}_harmonized.bed")
        combined_df.to_csv(bed_path, sep='\t', index=False)
        print(f"\n  Written: {bed_path}")
        print(f"    {combined_df.shape[0]:,} phenotypes x "
              f"{len(common_ordered)} samples")

        groups_path = os.path.join(
            qtl_dir, f"{anc}_{args.label}.phenotype_groups.txt")
        sidecar_df[['phenotype_id', 'gene_id']].to_csv(
            groups_path, sep='\t', header=False, index=False)
        n_genes = sidecar_df['gene_id'].nunique()
        print(f"  Written: {groups_path}")
        print(f"    {len(sidecar_df):,} phenotypes -> {n_genes:,} genes "
              f"(cross-modality groups)")

        sidecar_path = os.path.join(
            qtl_dir, f"{anc}_{args.label}.phenotype_modalities.tsv")
        sidecar_df.to_csv(sidecar_path, sep='\t', index=False)
        print(f"  Written: {sidecar_path}")

        # Per-modality phenotype summary
        print(f"\n  Phenotypes per modality in combined BED:")
        for mod, n in sidecar_df['modality'].value_counts().items():
            print(f"    {mod:18s} {n:>7,}")
        print(f"\n  Done: {anc}\n")


if __name__ == "__main__":
    main()
