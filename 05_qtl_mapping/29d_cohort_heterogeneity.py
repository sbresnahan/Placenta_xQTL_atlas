#!/usr/bin/env python3
"""
29d_cohort_heterogeneity.py — Cohort-heterogeneity battery for pooled
cis-QTL hits (module 05; part of the 29b-29f validation gate).

Pair set A (primary): pooled-significant pairs — by default the qval <=
--fdr leads from {results_dir}/{ANC}_{modality}_cisqtl.parquet, or an
explicit --pairs TSV with phenotype_id, variant_id columns.

  1. Per-cohort nominal tensorQTL scans (default: the two largest cohorts)
     on cohort subsets of the harmonized BED and covariates, with cohort
     indicator covariates dropped (constant within a cohort) and MAF
     re-filtered within each cohort by map_nominal.
  2. Inverse-variance fixed-effect meta-analysis across the per-cohort
     effect estimates -> Cochran Q + I2 per pair (BH-FDR across pairs).
  3. Pooled genotype x cohort interaction scan (tensorQTL map_nominal with
     interaction_df = second-cohort indicator on the two-cohort pooled
     sample set; the interaction main effect absorbs the cohort shift, so
     cohort dummies are dropped from the covariates here too) -> per-pair
     interaction p-value (BH-FDR across pairs).

  Output: {ANC}_{modality}_cohort_heterogeneity.tsv.gz +
          {ANC}_{modality}_cohort_heterogeneity_summary.tsv
          (n pairs, median I2, FDR-significant Cochran-Q and GxC counts —
          the metrics tracked by the validation gate).

Pair set B (optional, --extra-pairs): e.g. Choi-defined pairs from 29c.
Per-cohort nominal stats only (no meta-analysis), written to
{ANC}_{modality}_extra_pairs_withincohort.tsv.gz; 29e derives the
within-cohort nominal-P<0.05 fractions from it.

Runs in the tensorqtl conda env. Example:
  python3 29d_cohort_heterogeneity.py --qtl-dir $QTL_DIR \
      --results-dir $RESULTS_DIR --ancestry EAS --output-dir $DIAG_DIR/EAS \
      [--cohorts GUSTO,SNUH] [--extra-pairs $DIAG_DIR/EAS/Choi_vs_pooled_EAS.choi_defined_pairs.tsv]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

BED_META = ["#chr", "start", "end", "phenotype_id"]


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return out
    order = np.argsort(pv)
    q = pv[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    tmp = np.empty(n)
    tmp[order] = np.minimum(q, 1.0)
    out[ok] = tmp
    return out


def load_pairs_from_parquet(parquet, fdr):
    df = pd.read_parquet(parquet)
    if not isinstance(df.index, pd.RangeIndex):
        df = df.reset_index()
        if "index" in df.columns and "phenotype_id" not in df.columns:
            df = df.rename(columns={"index": "phenotype_id"})
    if "qval" not in df.columns or not df["qval"].notna().any():
        sys.exit(f"ERROR: no qval column in {parquet}; pass --pairs instead")
    sig = df[pd.to_numeric(df["qval"], errors="coerce") <= fdr]
    pairs = sig[["phenotype_id", "variant_id"]].drop_duplicates()
    print(f"  {len(pairs)} pooled-significant pairs "
          f"({sig['phenotype_id'].nunique()} phenotypes) at q <= {fdr}")
    return pairs


def load_pairs_tsv(path):
    d = pd.read_csv(path, sep="\t")
    if not {"phenotype_id", "variant_id"} <= set(d.columns):
        sys.exit(f"ERROR: {path} must have phenotype_id and variant_id "
                 f"columns (has: {list(d.columns)})")
    return d[["phenotype_id", "variant_id"]].drop_duplicates()


def prep_covariates(cov_path, samples, drop_cohort=True):
    """samples x covariates, aligned to samples; cohort dummies and
    zero-variance columns dropped."""
    cov = pd.read_csv(cov_path, sep="\t", index_col=0).T
    cov = cov.apply(pd.to_numeric, errors="coerce")
    cov = cov.reindex(samples)
    drop = []
    for c in cov.columns:
        if drop_cohort and str(c).startswith("cohort_"):
            drop.append(c)
        elif cov[c].isna().all() or cov[c].std(ddof=0) == 0:
            drop.append(c)
    if drop:
        print(f"    dropped covariates (cohort dummies / zero variance): "
              f"{drop}")
    cov = cov.drop(columns=drop)
    # any remaining NaN -> column mean (defensive; covariates are complete
    # in the canonical pipeline)
    cov = cov.fillna(cov.mean())
    return cov


def reconcile_chr(pheno_pos, variant_df):
    pheno_has = pheno_pos["chr"].astype(str).str.startswith("chr").any()
    var_has = variant_df["chrom"].astype(str).str.startswith("chr").any()
    if pheno_has and not var_has:
        pheno_pos = pheno_pos.copy()
        pheno_pos["chr"] = pheno_pos["chr"].str.replace("^chr", "", regex=True)
    elif var_has and not pheno_has:
        variant_df = variant_df.copy()
        variant_df["chrom"] = "chr" + variant_df["chrom"].astype(str)
    return pheno_pos, variant_df


def run_nominal(cis, genotype_df, variant_df, pheno, pheno_pos, cov,
                maf_threshold, interaction_df=None):
    """map_nominal with write_output=False; returns the full pair table."""
    tag = " + interaction" if interaction_df is not None else ""
    print(f"    map_nominal: {pheno.shape[0]} phenotypes x "
          f"{pheno.shape[1]} samples{tag}")
    res = cis.map_nominal(
        genotype_df=genotype_df,
        variant_df=variant_df,
        phenotype_df=pheno,
        phenotype_pos_df=pheno_pos,
        covariates_df=cov,
        interaction_df=interaction_df,
        maf_threshold=maf_threshold,
        write_output=False,
        verbose=False,
    )
    return res


def extract_pairs(nominal_df, pairs, prefix):
    """Pull the requested (phenotype_id, variant_id) rows out of a
    map_nominal result; rename stat columns with the given prefix."""
    nominal_df = nominal_df.set_index(
        pd.MultiIndex.from_frame(nominal_df[["phenotype_id", "variant_id"]]))
    want = pd.MultiIndex.from_frame(pairs)
    hit = nominal_df.reindex(want)
    n_hit = int(hit["pval_nominal"].notna().sum()) \
        if "pval_nominal" in hit.columns else 0
    out = pd.DataFrame(index=want)
    for src, dst in (("slope", f"beta_{prefix}"), ("slope_se", f"se_{prefix}"),
                     ("pval_nominal", f"p_{prefix}"), ("af", f"af_{prefix}"),
                     ("ma_samples", f"n_{prefix}")):
        if src in hit.columns:
            out[dst] = hit[src].values
    return out.reset_index(), n_hit


def find_interaction_pcol(cols):
    for c in cols:
        if c.lower() in ("pval_gi", "pval_gxi", "pval_gxc"):
            return c
    return None


def main():
    ap = argparse.ArgumentParser(
        description="Cohort-heterogeneity battery for pooled cis-QTL hits")
    ap.add_argument("--qtl-dir", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--ancestry", required=True)
    ap.add_argument("--modality", default="expression")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--pairs", default=None,
                    help="TSV with phenotype_id, variant_id (default: "
                         "qval <= --fdr leads from the results parquet)")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--cohorts", default=None,
                    help="Comma-separated pair of cohorts (default: the two "
                         "largest in {ANC}_metadata.tsv)")
    ap.add_argument("--extra-pairs", default=None,
                    help="Optional second pair set (e.g. Choi-defined pairs "
                         "from 29c): per-cohort nominal stats only")
    ap.add_argument("--maf-threshold", type=float, default=0.01)
    ap.add_argument("--cis-window", type=int, default=1000000)
    args = ap.parse_args()

    anc, mod = args.ancestry, args.modality
    qtl_dir = args.qtl_dir
    os.makedirs(args.output_dir, exist_ok=True)

    # ---- inputs ----
    bed_path = os.path.join(qtl_dir, f"{anc}_{mod}.bed.gz")
    if not os.path.exists(bed_path):
        bed_path = os.path.join(qtl_dir, f"{anc}_{mod}_harmonized.bed")
    meta_path = os.path.join(qtl_dir, f"{anc}_metadata.tsv")
    cov_path = os.path.join(qtl_dir, f"{anc}_covariates_{mod}.tsv")
    if not os.path.exists(cov_path):
        cov_path = os.path.join(qtl_dir, f"{anc}_covariates.tsv")
    plink_prefix = os.path.join(qtl_dir, f"{anc}_qtl")
    for f in (bed_path, meta_path, cov_path, f"{plink_prefix}.pgen"):
        if not os.path.exists(f):
            sys.exit(f"ERROR: required input not found: {f}")

    if args.pairs:
        pairs = load_pairs_tsv(args.pairs)
    else:
        pairs = load_pairs_from_parquet(
            os.path.join(args.results_dir,
                         f"{anc}_{mod}_cisqtl.parquet"), args.fdr)
    if len(pairs) == 0:
        sys.exit("ERROR: empty pair set")

    meta = pd.read_csv(meta_path, sep="\t")
    cohort_of = dict(zip(meta["array_id"], meta["cohort"]))

    print(f"[{anc} / {mod}] 29d cohort heterogeneity")
    print(f"  pairs: {len(pairs)}")

    import tensorqtl
    from tensorqtl import genotypeio, cis

    print("  Loading phenotypes + genotypes ...")
    pheno_all, pheno_pos_all = tensorqtl.read_phenotype_bed(bed_path)
    genotype_df, variant_df = genotypeio.load_genotypes(plink_prefix)
    pheno_pos_all, variant_df = reconcile_chr(pheno_pos_all, variant_df)

    # cohort assignments for BED samples
    bed_samples = list(pheno_all.columns)
    cohorts_series = pd.Series({s: cohort_of.get(s, "NA")
                                for s in bed_samples})
    cohorts_series = cohorts_series[cohorts_series != "NA"]
    if args.cohorts:
        coh_a, coh_b = [c.strip() for c in args.cohorts.split(",")]
    else:
        top2 = cohorts_series.value_counts().index[:2]
        coh_a, coh_b = sorted(top2,
                              key=lambda c: -int((cohorts_series == c).sum()))
    n_a = int((cohorts_series == coh_a).sum())
    n_b = int((cohorts_series == coh_b).sum())
    print(f"  cohorts: {coh_a} (n={n_a}) vs {coh_b} (n={n_b})")
    if n_a < 10 or n_b < 10:
        print("  WARN: a comparison cohort has < 10 samples; per-cohort "
              "estimates will be noisy")

    # subset phenotypes to the pair-set phenotypes
    ph_keep = pairs["phenotype_id"].unique()
    pheno = pheno_all.loc[pheno_all.index.intersection(ph_keep)]
    pheno_pos = pheno_pos_all.loc[pheno.index]
    n_missing = len(set(ph_keep) - set(pheno.index))
    if n_missing:
        print(f"  NOTE: {n_missing} pair phenotypes absent from the BED "
              f"(skipped)")
    pairs = pairs[pairs["phenotype_id"].isin(pheno.index)].reset_index(
        drop=True)
    print(f"  scanning {pheno.shape[0]} phenotypes / {len(pairs)} pairs")

    samples_a = [s for s in bed_samples if cohort_of.get(s) == coh_a]
    samples_b = [s for s in bed_samples if cohort_of.get(s) == coh_b]

    # ---- per-cohort nominal scans ----
    per_cohort = {}
    for tag, samples in ((coh_a, samples_a), (coh_b, samples_b)):
        print(f"  -- cohort {tag} (n={len(samples)}) --")
        cov_sub = prep_covariates(cov_path, samples)
        res = run_nominal(cis, genotype_df[samples], variant_df,
                          pheno[samples], pheno_pos, cov_sub,
                          args.maf_threshold)
        ext, n_hit = extract_pairs(res, pairs, tag)
        per_cohort[tag] = ext.set_index(["phenotype_id", "variant_id"])
        print(f"    exact pairs recovered: {n_hit} / {len(pairs)}")

    out = pairs.copy()
    for tag in (coh_a, coh_b):
        out = out.merge(per_cohort[tag].reset_index(),
                        on=["phenotype_id", "variant_id"], how="left")

    # ---- inverse-variance meta-analysis + Cochran Q + I2 ----
    ba, sa = out[f"beta_{coh_a}"].values, out[f"se_{coh_a}"].values
    bb, sb = out[f"beta_{coh_b}"].values, out[f"se_{coh_b}"].values
    ok = np.isfinite(ba) & np.isfinite(bb) & (sa > 0) & (sb > 0)
    wa = np.where(ok, 1.0 / sa**2, 0.0)
    wb = np.where(ok, 1.0 / sb**2, 0.0)
    wsum = wa + wb
    with np.errstate(divide="ignore", invalid="ignore"):
        b_meta = np.where(ok, (wa * ba + wb * bb) / wsum, np.nan)
        se_meta = np.where(ok, np.sqrt(1.0 / wsum), np.nan)
        q = np.where(ok, wa * (ba - b_meta)**2 + wb * (bb - b_meta)**2,
                     np.nan)
    from scipy.stats import chi2, norm
    p_q = chi2.sf(q, 1)
    i2 = np.where(np.isfinite(q) & (q > 0),
                  np.maximum(0.0, (q - 1.0) / q),
                  np.where(np.isfinite(q), 0.0, np.nan))
    z_meta = b_meta / se_meta
    out["beta_meta"] = b_meta
    out["se_meta"] = se_meta
    out["p_meta"] = 2 * norm.sf(np.abs(z_meta))
    out["cochran_q"] = q
    out["p_cochran_q"] = p_q
    out["i2"] = i2
    out["q_cochran_q"] = bh_fdr(p_q)

    # ---- pooled genotype x cohort interaction scan ----
    print(f"  -- GxC interaction scan ({coh_a}+{coh_b} pooled; "
          f"interaction = {coh_b} indicator) --")
    samples_ab = samples_a + samples_b
    cov_ab = prep_covariates(cov_path, samples_ab)
    interaction = pd.DataFrame(
        {f"cohort_{coh_b}": [0.0] * len(samples_a) + [1.0] * len(samples_b)},
        index=samples_ab)
    res_gxc = run_nominal(cis, genotype_df[samples_ab], variant_df,
                          pheno[samples_ab], pheno_pos, cov_ab,
                          args.maf_threshold, interaction_df=interaction)
    pcol_i = find_interaction_pcol(res_gxc.columns)
    if pcol_i is None:
        print(f"  WARN: no interaction p-value column found in map_nominal "
              f"output (columns: {list(res_gxc.columns)}); GxC skipped")
        out["p_gxc"] = np.nan
        out["q_gxc"] = np.nan
    else:
        bcol_i = pcol_i.replace("pval", "b")
        scol_i = pcol_i.replace("pval", "b") + "_se"
        gxc = res_gxc.set_index(["phenotype_id", "variant_id"])
        want = pd.MultiIndex.from_frame(pairs)
        hit = gxc.reindex(want)
        out["beta_gxc"] = hit[bcol_i].values if bcol_i in gxc.columns else np.nan
        out["se_gxc"] = hit[scol_i].values if scol_i in gxc.columns else np.nan
        out["p_gxc"] = hit[pcol_i].values
        out["q_gxc"] = bh_fdr(out["p_gxc"].values)
        print(f"    GxC pairs recovered: "
              f"{int(np.isfinite(out['p_gxc']).sum())} / {len(pairs)}")

    out_path = os.path.join(args.output_dir,
                            f"{anc}_{mod}_cohort_heterogeneity.tsv.gz")
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"  Written: {out_path}")

    # ---- run summary (validation-gate metrics) ----
    n_q = int((out["q_cochran_q"] <= args.fdr).sum())
    n_gxc = int((out["q_gxc"] <= args.fdr).sum()) \
        if "q_gxc" in out.columns else 0
    summary = pd.DataFrame({
        "metric": ["n_pairs", "n_pairs_both_cohorts", "median_i2",
                   "n_cochran_q_fdr", "n_gxc_fdr", "cohort_a", "cohort_b",
                   "n_cohort_a", "n_cohort_b"],
        "value": [len(pairs), int(ok.sum()),
                  float(np.nanmedian(out["i2"])), n_q, n_gxc,
                  coh_a, coh_b, n_a, n_b],
    })
    sum_path = os.path.join(
        args.output_dir, f"{anc}_{mod}_cohort_heterogeneity_summary.tsv")
    summary.to_csv(sum_path, sep="\t", index=False)
    print(f"  Written: {sum_path}")
    print(f"  median I2 = {summary.loc[2, 'value']:.3f}; "
          f"Cochran-Q FDR hits = {n_q}; GxC FDR hits = {n_gxc}")

    # ---- optional pair set B: per-cohort nominal stats only ----
    if args.extra_pairs:
        print("  -- extra pair set (per-cohort nominal only) --")
        extra = load_pairs_tsv(args.extra_pairs)
        ph_b = extra["phenotype_id"].unique()
        pheno_b = pheno_all.loc[pheno_all.index.intersection(ph_b)]
        pos_b = pheno_pos_all.loc[pheno_b.index]
        extra = extra[extra["phenotype_id"].isin(pheno_b.index)]
        print(f"    {len(extra)} pairs over {pheno_b.shape[0]} phenotypes "
              f"present in the BED")
        out_b = extra.reset_index(drop=True).copy()
        for tag, samples in ((coh_a, samples_a), (coh_b, samples_b)):
            cov_sub = prep_covariates(cov_path, samples)
            res = run_nominal(cis, genotype_df[samples], variant_df,
                              pheno_b[samples], pos_b, cov_sub,
                              args.maf_threshold)
            ext, n_hit = extract_pairs(res, extra.reset_index(drop=True),
                                       tag)
            out_b = out_b.merge(ext,
                                on=["phenotype_id", "variant_id"], how="left")
            print(f"    {tag}: exact pairs recovered: {n_hit} / {len(extra)}")
        b_path = os.path.join(
            args.output_dir, f"{anc}_{mod}_extra_pairs_withincohort.tsv.gz")
        out_b.to_csv(b_path, sep="\t", index=False, compression="gzip")
        print(f"  Written: {b_path}")

    print(f"\n[{anc} / {mod}] 29d done.")


if __name__ == "__main__":
    main()
