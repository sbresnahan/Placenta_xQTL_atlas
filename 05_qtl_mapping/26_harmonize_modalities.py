#!/usr/bin/env python3
"""
26_harmonize_modalities.py — Harmonize modality BEDs to array_id sample space

Brings the 7 non-expression modality BEDs from script 20
(combat_modalities/combat_int/{ANC}_{modality}_combat_int.bed) into the final
QTL sample set (array_id), using the post-intersection, post-outlier metadata
written by scripts 23+24 ({ANC}_metadata.tsv in qtl_inputs). Because the
metadata is already post-outlier, this single step handles both intersection
and outlier exclusion for the modality BEDs — scripts 23/24 do NOT need
rerunning.

Sample-ID conventions handled:
  - Intersection modalities (isoforms, alt_TSS, alt_polyA, RNA_editing,
    stability): sample columns are original rnaseq_ids — matched directly.
  - Pre-pooled modalities (splicing, intron_retention): sample columns are
    NAMESPACED as {cohort}_{rnaseq_id} by harmonize_within_ancestry.py
    (cohort-specific cluster/event IDs required coordinate harmonization
    before pooling). For any column that is not a known rnaseq_id, the
    {cohort}_ prefix is stripped (trying every underscore split point, so
    multi-underscore cohort names like NIEHS_RICHS_... also work) and the
    remainder is matched against the metadata.

Replicate handling: some individuals (array_id) have multiple RNA-seq
samples (e.g. NIGMS: two placental quadrant samples per individual).
QTL mapping requires one observation per genome, so when multiple
QC-passing rnaseq_ids map to the same array_id, their normalized
phenotype values are averaged per phenotype (the continuous-data
equivalent of sumTechReps used by the source NIGMS study). Merges are
logged to {ANC}_replicate_merges.tsv in the QTL directory.

Also copies each modality's phenotype_groups.txt into qtl_inputs.

Usage:
  python3 26_harmonize_modalities.py \
      --qtl-dir <qtl_inputs dir> \
      --modality-dir <combat_modalities/combat_int dir> \
      --ancestries "EAS EUR" \
      --modalities "alt_polyA alt_TSS intron_retention isoforms RNA_editing splicing stability"
"""

import argparse
import os
import shutil
import sys
import pandas as pd

BED_META_COLS = {'#chr', 'start', 'end', 'phenotype_id', 'chr'}


def strip_namespace(col, id_set):
    """Try stripping a leading '{prefix}_' namespace from a sample column;
    return the matching base sample ID, or None. Tries every underscore
    split point (shortest prefix first) so multi-underscore prefixes such
    as 'NIEHS_RICHS_SAMPLE1' are handled."""
    parts = col.split('_')
    for i in range(1, len(parts)):
        cand = '_'.join(parts[i:])
        if cand in id_set:
            return cand
    return None


def harmonize_bed(bed_path, out_path, id_map, keep_samples):
    """Rename sample columns to array_id (handling namespaced IDs), keep
    final samples, and average replicates (>1 rnaseq_id -> same array_id)
    per phenotype.

    Returns (n_kept_individuals, n_dropped_samples, kept_ids, merges,
             n_namespaced, ns_examples).
    """
    with open(bed_path, 'r') as f:
        header = f.readline().strip().split('\t')
    sample_cols = [c for c in header if c not in BED_META_COLS]

    id_set = set(id_map.keys())
    rename = {}
    dropped = []
    n_namespaced = 0
    ns_examples = []
    for c in sample_cols:
        if c in id_map:
            # Original ID: keep only if in the final post-outlier set
            if id_map[c] in keep_samples:
                rename[c] = id_map[c]
            else:
                dropped.append(c)
        else:
            # Namespaced ID ({cohort}_{rnaseq_id}): strip prefix and retry
            stripped = strip_namespace(c, id_set)
            if stripped is not None and id_map[stripped] in keep_samples:
                rename[c] = id_map[stripped]
                n_namespaced += 1
                if len(ns_examples) < 3:
                    ns_examples.append((c, stripped))
            else:
                dropped.append(c)

    usecols = [c for c in header if c in BED_META_COLS or c in rename]
    df = pd.read_csv(bed_path, sep='\t', usecols=usecols)
    df = df.rename(columns=rename)

    # Identify replicate merges (array_id <- multiple sample columns)
    inv = {}
    for orig, arr_id in rename.items():
        inv.setdefault(arr_id, []).append(orig)
    merges = {a: sorted(r) for a, r in inv.items() if len(r) > 1}

    meta_cols = [c for c in dict.fromkeys(df.columns) if c in BED_META_COLS]
    if df.columns.duplicated().any():
        # Average replicate sample columns per phenotype (row-wise mean)
        samp = df.drop(columns=meta_cols)
        samp = samp.T.groupby(level=0).mean().T
        # Restore first-appearance column order
        order = list(dict.fromkeys(rename.values()))
        samp = samp[order]
        df = pd.concat([df[meta_cols].reset_index(drop=True),
                        samp.reset_index(drop=True)], axis=1)

    # Safety: columns must be unique now
    dup = df.columns[df.columns.duplicated()].tolist()
    if dup:
        sys.exit(f"  ERROR: duplicate sample columns remain after averaging "
                 f"in {bed_path}: {dup[:5]}")

    df.to_csv(out_path, sep='\t', index=False)
    kept_ids = [c for c in df.columns if c not in BED_META_COLS]
    return len(kept_ids), len(dropped), kept_ids, merges, n_namespaced, ns_examples


def main():
    parser = argparse.ArgumentParser(description="Harmonize modality BEDs to array_id space")
    parser.add_argument("--qtl-dir", required=True,
                        help="QTL inputs directory (contains {ANC}_metadata.tsv from scripts 23+24)")
    parser.add_argument("--modality-dir", required=True,
                        help="Directory with {ANC}_{modality}_combat_int.bed and "
                             "{ANC}_{modality}.phenotype_groups.txt (script 20 output)")
    parser.add_argument("--ancestries", default="EAS EUR",
                        help="Space-separated ancestry labels")
    parser.add_argument("--modalities",
                        default="alt_polyA alt_TSS intron_retention isoforms RNA_editing splicing stability",
                        help="Space-separated modality labels")
    args = parser.parse_args()

    ancestries = args.ancestries.split()
    modalities = args.modalities.split()

    for anc in ancestries:
        print(f"\n{'='*60}")
        print(f"Ancestry: {anc}")
        print(f"{'='*60}")

        # ---- Final sample set: post-intersection, post-outlier metadata ----
        meta_path = os.path.join(args.qtl_dir, f"{anc}_metadata.tsv")
        if not os.path.exists(meta_path):
            print(f"  ERROR: metadata not found: {meta_path} (run scripts 23+24 first)")
            continue
        meta = pd.read_csv(meta_path, sep='\t')
        id_map = dict(zip(meta['rnaseq_id'], meta['array_id']))
        keep_samples = set(meta['array_id'])
        n_indiv = len(keep_samples)
        n_rows = len(meta)
        print(f"  Final samples (post-outlier): {n_indiv} unique individuals "
              f"from {n_rows} RNA-seq samples")
        if n_rows > n_indiv:
            print(f"  NOTE: {n_rows - n_indiv} replicate RNA-seq sample(s) will be "
                  f"averaged per individual where present in a modality BED")

        all_merges = []
        for mod in modalities:
            bed_path = os.path.join(args.modality_dir, f"{anc}_{mod}_combat_int.bed")
            groups_path = os.path.join(args.modality_dir, f"{anc}_{mod}.phenotype_groups.txt")
            out_bed = os.path.join(args.qtl_dir, f"{anc}_{mod}_harmonized.bed")
            out_groups = os.path.join(args.qtl_dir, f"{anc}_{mod}.phenotype_groups.txt")

            if not os.path.exists(bed_path):
                print(f"  WARN: {mod}: BED not found: {bed_path}, skipping")
                continue

            (n_kept, n_dropped, kept_ids, merges,
             n_namespaced, ns_examples) = harmonize_bed(
                bed_path, out_bed, id_map, keep_samples)
            n_pheno = sum(1 for _ in open(out_bed)) - 1
            print(f"  {mod}: {n_kept} individuals kept "
                  f"({len(merges)} averaged from replicates, "
                  f"{n_namespaced} matched via namespace stripping), "
                  f"{n_dropped} samples dropped (not in final set), "
                  f"{n_pheno} phenotypes")
            if ns_examples:
                print(f"    Namespace examples: " +
                      "; ".join(f"{a} -> {b}" for a, b in ns_examples))
            print(f"    Written: {out_bed}")

            for arr_id, orig_ids in merges.items():
                all_merges.append({'modality': mod, 'array_id': arr_id,
                                   'n_replicates': len(orig_ids),
                                   'sample_columns': ','.join(orig_ids)})

            missing_vs_final = len(keep_samples - set(kept_ids))
            if missing_vs_final > 0:
                print(f"    WARN: {missing_vs_final} final individuals have no {mod} data "
                      f"(tensorQTL will use the common subset)")

            if os.path.exists(groups_path):
                shutil.copy(groups_path, out_groups)
                print(f"    Written: {out_groups}")
            else:
                print(f"    WARN: groups file not found: {groups_path} "
                      f"(script 27 will run {mod} ungrouped)")

        if all_merges:
            merges_df = pd.DataFrame(all_merges)
            merges_path = os.path.join(args.qtl_dir, f"{anc}_replicate_merges.tsv")
            merges_df.to_csv(merges_path, sep='\t', index=False)
            n_indiv_merged = merges_df['array_id'].nunique()
            print(f"\n  Replicate merges: {n_indiv_merged} individuals averaged "
                  f"(logged to {merges_path})")

    print(f"\nDone.")


if __name__ == "__main__":
    main()
