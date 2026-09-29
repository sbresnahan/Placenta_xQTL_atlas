#!/usr/bin/env python3
"""34_pip_annotation_enrichment.py — functional-annotation enrichment of
high-PIP fine-mapped variants (Objective 1.5, second half).

Annotates the fine-mapped variant universe (all variants in the aggregated
PIP table) with:
  - VEP consequence classes (from a VEP run on the variant universe; see
    runbook — VEP input is generated with --make-vep-input),
  - ENCODE SCREEN cCRE classes (PLS, pELS, dELS, CTCF-bound, CA-TF BEDs),
  - placenta-specific open chromatin (ENCODE placenta DNase/ATAC peaks BED).

A variant appearing in multiple loci is counted once, carrying its maximum
PIP and minimum distance to a phenotype start across loci.

Enrichment per annotation class (high-PIP = pip_all >= --pip-threshold,
default 0.9; background = all other fine-mapped variants):
  - Fisher exact odds ratio + p-value (primary),
  - logistic regression high_pip ~ annotation + log10(dist_to_pheno_start+1)
    (sensitivity; MAF is not carried because the intersected pgens already
    enforce a uniform MAC>=5 floor),
  - PIP-weighted enrichment: (sum of PIPs in class / sum of all PIPs) /
    (fraction of variants in class).
BH-FDR is applied across classes for each test family.

Usage:
  python3 34_pip_annotation_enrichment.py \
      --pips finemap/aggregated/finemap_pips.tsv.gz \
      --vep vep_output.txt \
      --ccre GRCh38-cCREs.PLS.bed:PLS GRCh38-cCREs.pELS.bed:pELS ... \
      --placenta-ocr encode_placenta_dnase_peaks.bed \
      --out pip_annotation_enrichment.tsv

  python3 34_pip_annotation_enrichment.py --make-vep-input \
      --pips finemap/aggregated/finemap_pips.tsv.gz --out vep_input.tsv
"""

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# VEP Sequence Ontology consequence -> grouped class, in priority order
# (first match wins; a variant's class = highest-priority term across all
# its transcript rows).
VEP_PRIORITY = [
    ("splice_site", {"splice_acceptor_variant", "splice_donor_variant"}),
    ("coding_lof", {"stop_gained", "stop_lost", "start_lost",
                    "frameshift_variant"}),
    ("coding_nonsynonymous", {"missense_variant", "inframe_deletion",
                              "inframe_insertion"}),
    ("synonymous", {"synonymous_variant"}),
    ("utr", {"5_prime_UTR_variant", "3_prime_UTR_variant"}),
    ("intron", {"intron_variant"}),
    ("regulatory_region", {"regulatory_region_variant"}),
    ("flanking", {"upstream_gene_variant", "downstream_gene_variant"}),
    ("intergenic", {"intergenic_variant"}),
]
TERM_TO_CLASS = {t: c for c, terms in VEP_PRIORITY for t in terms}
CLASS_PRIORITY = {c: i for i, (c, _) in enumerate(VEP_PRIORITY)}


def make_vep_input(pips_path, out_path):
    """Write VEP region-format input: chrom start end ref/alt strand id."""
    pips = pd.read_csv(pips_path, sep="\t", usecols=["snp"])
    snps = pips["snp"].drop_duplicates()
    parts = snps.str.split(":", expand=True)
    vep = pd.DataFrame({
        "chrom": parts[0], "start": parts[1].astype(int),
        "end": parts[1].astype(int),
        "allele": parts[2] + "/" + parts[3], "strand": "+", "id": snps})
    vep = vep.sort_values(["chrom", "start"])
    vep.to_csv(out_path, sep="\t", header=False, index=False)
    print(f"wrote {len(vep)} variants -> {out_path}")
    print("Run VEP with default output, e.g.:")
    print(f"  vep -i {out_path} --cache --offline --assembly GRCh38 "
          f"--output_file vep_output.txt --force_overwrite")


def parse_vep(vep_path):
    """Parse VEP default output -> {variant_id: grouped_class}."""
    best = {}
    with open(vep_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            var_id, consequence = fields[0], fields[6]
            cls = None
            for term in consequence.split(","):
                term = term.strip()
                if term in TERM_TO_CLASS:
                    c = TERM_TO_CLASS[term]
                    if cls is None or CLASS_PRIORITY[c] < CLASS_PRIORITY[cls]:
                        cls = c
            if cls is None:
                cls = "other"
            if var_id not in best or CLASS_PRIORITY[cls] < CLASS_PRIORITY[best[var_id]]:
                best[var_id] = cls
    return best


def load_bed_intervals(bed_path):
    """BED -> {chrom: (starts, ends)} sorted numpy arrays (0-based)."""
    bed = pd.read_csv(bed_path, sep="\t", header=None, comment="#",
                      usecols=[0, 1, 2], names=["chrom", "start", "end"])
    bed["chrom"] = bed["chrom"].astype(str).str.replace("^chr", "", regex=True)
    out = {}
    for chrom, grp in bed.groupby("chrom"):
        g = grp.sort_values("start")
        out[str(chrom)] = (g["start"].to_numpy(), g["end"].to_numpy())
    return out


def overlap_flags(chroms, positions, intervals):
    """Boolean array: does each (chrom, pos) point fall in any interval?"""
    flags = np.zeros(len(positions), dtype=bool)
    df = pd.DataFrame({"chrom": chroms, "pos": positions})
    for chrom, idx in df.groupby("chrom").groups.items():
        if str(chrom) not in intervals:
            continue
        starts, ends = intervals[str(chrom)]
        pos = df.loc[idx, "pos"].to_numpy()
        # rightmost interval with start <= pos; overlap if its end > pos
        j = np.searchsorted(starts, pos, side="right") - 1
        hit = (j >= 0) & (ends[np.clip(j, 0, None)] > pos)
        flags[df.index.get_indexer(idx)] = hit
    return flags


def enrichment_table(df, anno_cols, pip_threshold):
    """Fisher + logistic + PIP-weighted enrichment per annotation class."""
    from scipy.stats import fisher_exact
    try:
        import statsmodels.api as sm
        have_sm = True
    except ImportError:
        have_sm = False

    high = (df["pip_all"] >= pip_threshold).to_numpy()
    log_dist = np.log10(df["min_dist_pheno"].to_numpy() + 1.0)
    rows = []
    for col in anno_cols:
        anno = df[col].to_numpy()
        n_anno = int(anno.sum())
        if n_anno == 0:
            continue
        a = int((high & anno).sum())       # high-PIP in class
        b = int((high & ~anno).sum())      # high-PIP not in class
        c = int((~high & anno).sum())      # background in class
        d = int((~high & ~anno).sum())     # background not in class
        odds, p_fisher = fisher_exact([[a, b], [c, d]])
        pip_w = (df.loc[anno, "pip_all"].sum() / df["pip_all"].sum()) / \
                (n_anno / len(df))
        row = {"annotation": col, "n_variants": n_anno,
               "n_high_pip": a,
               "cell_a_high_in": a, "cell_b_high_out": b,
               "cell_c_bg_in": c, "cell_d_bg_out": d,
               "frac_high_in_class": a / n_anno,
               "frac_high_outside": b / max(int((~anno).sum()), 1),
               "odds_ratio": odds, "fisher_p": p_fisher,
               "pip_weighted_enrichment": pip_w}
        if have_sm and 0 < a < n_anno and b > 0:
            try:
                X = pd.DataFrame({"anno": anno.astype(float),
                                  "log_dist": log_dist})
                X = sm.add_constant(X)
                fit = sm.Logit(high.astype(float), X).fit(disp=0)
                row["logistic_or"] = float(np.exp(fit.params["anno"]))
                row["logistic_p"] = float(fit.pvalues["anno"])
            except Exception:
                row["logistic_or"] = np.nan
                row["logistic_p"] = np.nan
        rows.append(row)
    res = pd.DataFrame(rows)
    for col in ("fisher_p", "logistic_p"):
        if col in res and res[col].notna().any():
            p = res[col].fillna(1.0).to_numpy()
            order = np.argsort(p)
            ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
            fdr = np.minimum.accumulate(ranked[::-1])[::-1]
            res[col.replace("_p", "_fdr")] = np.clip(
                fdr[np.argsort(order)], 0, 1)
    return res.sort_values("fisher_p")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pips", required=True,
                   help="Aggregated finemap_pips.tsv.gz from script 33")
    p.add_argument("--vep", default=None,
                   help="VEP default-format output for the variant universe")
    p.add_argument("--ccre", nargs="*", default=[],
                   help="BED:label pairs, e.g. GRCh38-cCREs.PLS.bed:PLS")
    p.add_argument("--placenta-ocr", default=None,
                   help="BED of placenta open-chromatin peaks")
    p.add_argument("--pip-threshold", type=float, default=0.9)
    p.add_argument("--make-vep-input", action="store_true",
                   help="Only write the VEP input file and exit")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    if args.make_vep_input:
        make_vep_input(args.pips, args.out)
        return

    print(f"Loading PIPs from {args.pips} ...")
    pips = pd.read_csv(args.pips, sep="\t")
    # one row per unique variant: max PIP, min distance to a phenotype start
    pips["dist"] = (pips["pos"] - pips["pheno_start"]).abs()
    uni = (pips.groupby(["chrom", "pos", "snp"], as_index=False)
           .agg(pip_all=("pip_all", "max"),
                min_dist_pheno=("dist", "min")))
    print(f"  {len(uni)} unique variants "
          f"({int((uni['pip_all'] >= args.pip_threshold).sum())} with "
          f"PIP >= {args.pip_threshold})")

    anno_cols = []

    if args.vep:
        print(f"Parsing VEP consequences from {args.vep} ...")
        vep_cls = parse_vep(args.vep)
        uni["vep_class"] = uni["snp"].map(vep_cls).fillna("not_annotated")
        for cls in ["coding_lof", "coding_nonsynonymous", "synonymous",
                    "splice_site", "utr", "intron", "regulatory_region",
                    "flanking", "intergenic"]:
            uni[f"vep_{cls}"] = uni["vep_class"] == cls
            anno_cols.append(f"vep_{cls}")
        uni["vep_coding_any"] = (uni["vep_coding_lof"]
                                 | uni["vep_coding_nonsynonymous"])
        anno_cols.append("vep_coding_any")

    for spec in args.ccre:
        bed_path, label = spec.rsplit(":", 1)
        col = f"ccre_{label}"
        print(f"Intersecting cCRE {label} ({bed_path}) ...")
        uni[col] = overlap_flags(uni["chrom"].to_numpy(),
                                 uni["pos"].to_numpy(),
                                 load_bed_intervals(bed_path))
        anno_cols.append(col)

    if args.placenta_ocr:
        print(f"Intersecting placenta OCR ({args.placenta_ocr}) ...")
        uni["placenta_ocr"] = overlap_flags(
            uni["chrom"].to_numpy(), uni["pos"].to_numpy(),
            load_bed_intervals(args.placenta_ocr))
        anno_cols.append("placenta_ocr")

    if not anno_cols:
        raise SystemExit("No annotations supplied (--vep, --ccre, "
                         "--placenta-ocr)")

    res = enrichment_table(uni, anno_cols, args.pip_threshold)
    res.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {len(res)} annotation tests -> {args.out}")
    show = res[["annotation", "n_variants", "n_high_pip", "odds_ratio",
                "fisher_p"]].head(20)
    print(show.to_string(index=False))


if __name__ == "__main__":
    main()
