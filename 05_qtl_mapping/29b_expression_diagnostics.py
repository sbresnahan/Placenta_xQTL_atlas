#!/usr/bin/env python3
"""
29b_expression_diagnostics.py — Expression diagnostic inputs for the
within-cohort-INT validation gate (module 05; runs after 29, gates 31).

Generates, per ancestry stratum, the three inputs consumed by
29e_plot_diagnostics.R (panel set 1):

  1. {ANC}_expression_sample_PCs.tsv
     PCA of the standardized harmonized expression matrix (sample_id, PC1-10,
     cohort). After pooled QN + within-cohort INT, no cohort structure should
     remain on the major PCs.

  2. {ANC}_gene_residual_variance_by_cohort.tsv.gz
     Per-gene residuals from the pooled mapping-covariate fit (the same
     covariate file used for the expression scan), summarized per cohort as
     residual variance (var_<cohort>), the max:min ratio across cohorts, and
     log2(max:min). Residuals come from the pooled fit so the tiny cohorts
     (n=4/7) remain measurable. Within-cohort INT equilibrates per-cohort
     scale by construction, so post-fix ratios should center on 1.

  3. {ANC}_expression_hcp15_vs_hcp45_exact_pairs.tsv.gz
     Union of chr1 lead SNP-gene pairs from the 25a k=15 and k=45 parquets,
     refit at BOTH k by OLS on covariate-residualized genotype/phenotype
     (the tensorQTL nominal model: intercept + the per-k covariate file),
     with flags marking lead status at each k. The 25a parquets are
     lead-per-gene only, so cross-k evaluation requires this targeted refit;
     it needs the intersected pgen and the per-k covariates from the 25a
     staging dir and therefore must run in the tensorqtl conda env
     (genotypeio). Skip with --skip-hcp-comparison.

Inputs (per ancestry; all post-23/24/25 pipeline outputs):
  {qtl_dir}/{ANC}_expression_harmonized.bed
  {qtl_dir}/{ANC}_covariates_expression.tsv  (fallback: {ANC}_covariates.tsv)
  {qtl_dir}/{ANC}_metadata.tsv               (array_id -> cohort)
  {qtl_dir}/{ANC}_qtl.{pgen,pvar,psam}       (k-comparison only)
  25a staging, auto-detected for K = --k-low/--k-high:
    isolated per-k: {qtl_dir}/hcp_optimization_per_k/{ANC}/k{K}/work/{ANC}/
    legacy serial:  {hcp_opt_dir}/{ANC}/
    legacy flat:    {hcp_opt_dir}/

Usage:
  python3 29b_expression_diagnostics.py --qtl-dir $QTL_DIR --ancestry EAS \
      --output-dir $DIAG_DIR/EAS

  python3 29b_expression_diagnostics.py --qtl-dir $QTL_DIR --ancestry EAS \
      --output-dir $DIAG_DIR/EAS --hcp-opt-dir $QTL_DIR/hcp_optimization \
      --k-low 15 --k-high 45
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

BED_META = ["#chr", "start", "end", "phenotype_id"]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def load_expression_bed(path):
    bed = pd.read_csv(path, sep="\t")
    missing = [c for c in BED_META if c not in bed.columns]
    if missing:
        sys.exit(f"ERROR: {path} missing BED columns {missing}")
    return bed


def load_metadata_cohorts(meta_path):
    meta = pd.read_csv(meta_path, sep="\t")
    if "array_id" not in meta.columns or "cohort" not in meta.columns:
        sys.exit(f"ERROR: {meta_path} must carry array_id and cohort columns "
                 f"(has: {list(meta.columns)})")
    return dict(zip(meta["array_id"], meta["cohort"]))


def load_covariates(cov_path):
    cov = pd.read_csv(cov_path, sep="\t", index_col=0).T  # samples x covariates
    cov = cov.apply(pd.to_numeric, errors="coerce")
    return cov


def residualize(mat, cov):
    """Residualize rows of mat (features x samples) on [1, cov] via QR.

    cov: samples x p DataFrame aligned to mat.columns. Returns residuals
    (features x samples, same order) and the residual degrees of freedom.
    """
    C = np.column_stack([np.ones(cov.shape[0]), cov.values])
    Q, _ = np.linalg.qr(C)
    X = mat.values.astype(float)  # features x samples
    resid = X - X @ Q @ Q.T
    df = C.shape[0] - np.linalg.matrix_rank(C)
    return resid, df


def ols_refit(genotype_df, pheno_vec, cov, variant_ids):
    """Refit specific variants against one phenotype vector.

    genotype_df: variants x samples (dosages), pheno_vec: pd.Series indexed
    by sample, cov: samples x p covariates. Returns per-variant
    (slope, se, t, af, n). Mirrors the tensorQTL nominal model: intercept +
    covariates residualized out of both genotype and phenotype.
    """
    samples = list(pheno_vec.index)
    C = np.column_stack([np.ones(len(samples)), cov.loc[samples].values])
    Q, _ = np.linalg.qr(C)
    rank = np.linalg.matrix_rank(C)
    df = len(samples) - rank - 1  # -1 for the genotype itself
    y = pheno_vec.values.astype(float)
    y_r = y - Q @ (Q.T @ y)

    out = {}
    for v in variant_ids:
        if v not in genotype_df.index:
            out[v] = (np.nan, np.nan, np.nan, np.nan, 0)
            continue
        x = genotype_df.loc[v, samples].values.astype(float)
        n = int(np.sum(~np.isnan(x)))
        af = np.nanmean(x) / 2.0
        af = min(af, 1 - af)
        x = np.nan_to_num(x, nan=np.nanmean(x))
        x_r = x - Q @ (Q.T @ x)
        xx = float(x_r @ x_r)
        if xx <= 0 or df <= 0:
            out[v] = (np.nan, np.nan, np.nan, af, n)
            continue
        slope = float(x_r @ y_r) / xx
        resid = y_r - slope * x_r
        sigma2 = float(resid @ resid) / df
        se = np.sqrt(sigma2 / xx)
        out[v] = (slope, se, slope / se if se > 0 else np.nan, af, n)
    return out


# ---------------------------------------------------------------------------
# component 1: sample PCs
# ---------------------------------------------------------------------------

def write_sample_pcs(bed, cohort_of, out_path, n_pcs=10):
    samples = [c for c in bed.columns if c not in BED_META]
    X = bed[samples].values.astype(float)  # genes x samples
    # z-score per gene across the pooled stratum
    mu = np.nanmean(X, axis=1, keepdims=True)
    sd = np.nanstd(X, axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = 1.0
    Z = np.nan_to_num((X - mu) / sd)  # genes x samples
    # SVD of samples x genes
    U, S, _ = np.linalg.svd(Z.T, full_matrices=False)
    pcs = U[:, :n_pcs] * S[:n_pcs]
    k = min(n_pcs, pcs.shape[1])
    out = pd.DataFrame({"sample_id": samples})
    for i in range(k):
        out[f"PC{i+1}"] = pcs[:, i]
    out["cohort"] = [cohort_of.get(s, "NA") for s in samples]
    out.to_csv(out_path, sep="\t", index=False)
    print(f"  Written: {out_path} ({len(out)} samples, PC1-{k})")
    return out


# ---------------------------------------------------------------------------
# component 2: per-cohort residual variance
# ---------------------------------------------------------------------------

def write_residual_variance(bed, cov, cohort_of, out_path):
    samples = [c for c in bed.columns if c not in BED_META]
    common = [s for s in samples if s in cov.index]
    if len(common) < len(samples):
        print(f"  NOTE: {len(samples) - len(common)} BED samples lack "
              f"covariates and are dropped from the residual fit")
    bed = bed[["phenotype_id"] + common]
    cov = cov.loc[common]

    resid, df = residualize(bed.set_index("phenotype_id")[common], cov)
    genes = bed["phenotype_id"].values
    cohorts = pd.Series([cohort_of.get(s, "NA") for s in common])
    keep = cohorts != "NA"
    print(f"  Residual fit: {len(genes)} genes x {int(keep.sum())} samples "
          f"(df={df}); cohorts: "
          + ", ".join(f"{c}={int(n)}" for c, n in cohorts[keep].value_counts().items()))

    out = pd.DataFrame({"phenotype_id": genes})
    R = pd.DataFrame(resid, columns=common)
    var_cols = []
    for co in sorted(cohorts[keep].unique()):
        idx = np.where((cohorts == co).values)[0]
        v = np.var(R.values[:, idx], axis=1, ddof=1) if len(idx) > 1 else np.full(len(genes), np.nan)
        col = f"var_{co}"
        out[col] = v
        var_cols.append(col)
    vmat = out[var_cols].values
    with np.errstate(divide="ignore", invalid="ignore"):
        vmax = np.nanmax(vmat, axis=1)
        vmin = np.nanmin(np.where(vmat > 0, vmat, np.nan), axis=1)
        ratio = vmax / vmin
    out["max_min_variance_ratio"] = ratio
    out["log2_variance_ratio"] = np.log2(ratio)
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"  Written: {out_path} ({len(out)} genes, "
          f"{len(var_cols)} cohort variance columns)")
    return out


# ---------------------------------------------------------------------------
# component 3: k-low vs k-high exact-pair comparison (25a parquets + refit)
# ---------------------------------------------------------------------------

def write_hcp_k_comparison(anc, qtl_dir, hcp_opt_dir, k_low, k_high,
                           out_path):
    # Lazy import: requires the tensorqtl conda env (pgen reading).
    try:
        from tensorqtl import genotypeio
    except ImportError:
        sys.exit("ERROR: the k-comparison needs tensorqtl (genotypeio); "
                 "run in the tensorqtl conda env or use --skip-hcp-comparison")

    pq = {}
    cov = {}
    for k in (k_low, k_high):
        # 25a has had three staging layouts. Prefer the current isolated
        # per-k sandbox, then the legacy serial ancestry subdir, then the
        # oldest flat layout. Require BOTH the parquet and matching per-k
        # covariate table from the same layout.
        isolated_stage = os.path.join(
            qtl_dir, "hcp_optimization_per_k", anc, f"k{k}", "work", anc)
        serial_stage = os.path.join(hcp_opt_dir, anc)
        flat_stage = hcp_opt_dir

        candidates = [
            ("isolated-per-k", isolated_stage),
            ("legacy-serial", serial_stage),
            ("legacy-flat", flat_stage),
        ]

        p = c = layout = None
        attempted = []
        for label, stage in candidates:
            p_try = os.path.join(
                stage, f"results_k{k}", f"{anc}_expression_cisqtl.parquet")
            c_try = os.path.join(stage, f"{anc}_covariates_k{k}.tsv")
            attempted.append((label, p_try, c_try))
            if os.path.exists(p_try) and os.path.exists(c_try):
                p, c, layout = p_try, c_try, label
                break

        if p is None:
            msg = [f"ERROR: required 25a staging files not found for {anc} k={k}.",
                   "Tried:"]
            for label, p_try, c_try in attempted:
                msg.append(f"  [{label}] parquet: {p_try}")
                msg.append(f"  [{label}] covariates: {c_try}")
            sys.exit("\n".join(msg))

        print(f"  k={k}: using {layout} 25a staging")
        print(f"    parquet:   {p}")
        print(f"    covariates:{c}")

        d = pd.read_parquet(p)
        if not isinstance(d.index, pd.RangeIndex):
            d = d.reset_index()
            if "index" in d.columns and "phenotype_id" not in d.columns:
                d = d.rename(columns={"index": "phenotype_id"})
        pq[k] = d
        cov[k] = load_covariates(c)

    lead_low = set(zip(pq[k_low]["phenotype_id"], pq[k_low]["variant_id"]))
    lead_high = set(zip(pq[k_high]["phenotype_id"], pq[k_high]["variant_id"]))
    pair_union = sorted(lead_low | lead_high)
    print(f"  Lead pairs: {len(lead_low)} at k={k_low}, {len(lead_high)} at "
          f"k={k_high}, union {len(pair_union)}")

    # Genotypes + phenotypes (chr1 harmonized BED)
    plink_prefix = os.path.join(qtl_dir, f"{anc}_qtl")
    genotype_df, variant_df = genotypeio.load_genotypes(plink_prefix)
    bed = load_expression_bed(os.path.join(
        qtl_dir, f"{anc}_expression_harmonized.bed"))
    bed = bed[bed["#chr"].astype(str).str.replace("^chr", "", regex=True) == "1"]
    expr = bed.set_index("phenotype_id")

    rows = []
    by_gene = {}
    for ph, var in pair_union:
        by_gene.setdefault(ph, []).append(var)
    n_missing_pheno = 0
    for ph, variants in by_gene.items():
        if ph not in expr.index:
            n_missing_pheno += 1
            continue
        for k in (k_low, k_high):
            samples = [s for s in cov[k].index if s in expr.columns
                       and s in genotype_df.columns]
            y = expr.loc[ph, samples]
            res = ols_refit(genotype_df, y, cov[k], variants)
            for v in variants:
                slope, se, t, af, n = res[v]
                rows.append({
                    "phenotype_id": ph, "variant_id": v, "k": k,
                    "slope": slope, "slope_se": se, "t": t, "af": af, "n": n,
                })
    if n_missing_pheno:
        print(f"  NOTE: {n_missing_pheno} lead phenotypes absent from the "
              f"chr1 harmonized BED (skipped)")

    long = pd.DataFrame(rows)
    wide = long.pivot_table(index=["phenotype_id", "variant_id"], columns="k",
                            values="t", aggfunc="first").reset_index()
    wide.columns = ["phenotype_id", "variant_id",
                    f"t_k{k_low}", f"t_k{k_high}"]
    wide[f"abs_t_k{k_low}"] = wide[f"t_k{k_low}"].abs()
    wide[f"abs_t_k{k_high}"] = wide[f"t_k{k_high}"].abs()
    wide[f"lead_at_k{k_low}"] = [
        (p, v) in lead_low for p, v in zip(wide["phenotype_id"], wide["variant_id"])]
    wide[f"lead_at_k{k_high}"] = [
        (p, v) in lead_high for p, v in zip(wide["phenotype_id"], wide["variant_id"])]
    wide = wide[["phenotype_id", "variant_id",
                 f"abs_t_k{k_low}", f"abs_t_k{k_high}",
                 f"lead_at_k{k_low}", f"lead_at_k{k_high}"]]
    wide.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"  Written: {out_path} ({len(wide)} exact pairs)")
    return wide


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Expression diagnostic inputs (sample PCs, per-cohort "
                    "residual variance, HCP k-low vs k-high exact pairs)")
    ap.add_argument("--qtl-dir", required=True,
                    help="QTL inputs dir (post-23/24/25 outputs)")
    ap.add_argument("--ancestry", required=True, help="Ancestry label")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--hcp-opt-dir", default=None,
                    help="25a staging dir (default: {qtl-dir}/hcp_optimization)")
    ap.add_argument("--k-low", type=int, default=15)
    ap.add_argument("--k-high", type=int, default=45)
    ap.add_argument("--skip-hcp-comparison", action="store_true",
                    help="Skip the k-low vs k-high refit (no tensorqtl needed)")
    args = ap.parse_args()

    anc = args.ancestry
    qtl_dir = args.qtl_dir
    hcp_opt_dir = args.hcp_opt_dir or os.path.join(qtl_dir, "hcp_optimization")
    os.makedirs(args.output_dir, exist_ok=True)

    bed_path = os.path.join(qtl_dir, f"{anc}_expression_harmonized.bed")
    meta_path = os.path.join(qtl_dir, f"{anc}_metadata.tsv")
    cov_path = os.path.join(qtl_dir, f"{anc}_covariates_expression.tsv")
    if not os.path.exists(cov_path):
        cov_path = os.path.join(qtl_dir, f"{anc}_covariates.tsv")
    for f in (bed_path, meta_path, cov_path):
        if not os.path.exists(f):
            sys.exit(f"ERROR: required input not found: {f}")

    print(f"[{anc}] 29b expression diagnostics")
    print(f"  Expression BED: {bed_path}")
    print(f"  Covariates:     {cov_path}")
    print(f"  Metadata:       {meta_path}")

    bed = load_expression_bed(bed_path)
    cohort_of = load_metadata_cohorts(meta_path)
    cov = load_covariates(cov_path)
    print(f"  Loaded: {bed.shape[0]} genes x "
          f"{bed.shape[1] - len(BED_META)} samples; "
          f"{cov.shape[1]} covariates")

    print("\n-- sample PCs --")
    write_sample_pcs(bed, cohort_of,
                     os.path.join(args.output_dir,
                                  f"{anc}_expression_sample_PCs.tsv"))

    print("\n-- per-cohort residual variance --")
    write_residual_variance(
        bed, cov, cohort_of,
        os.path.join(args.output_dir,
                     f"{anc}_gene_residual_variance_by_cohort.tsv.gz"))

    if args.skip_hcp_comparison:
        print("\n-- k-comparison skipped (--skip-hcp-comparison) --")
    else:
        print(f"\n-- HCP k={args.k_low} vs k={args.k_high} exact pairs --")
        write_hcp_k_comparison(
            anc, qtl_dir, hcp_opt_dir, args.k_low, args.k_high,
            os.path.join(
                args.output_dir,
                f"{anc}_expression_hcp{args.k_low}_vs_hcp{args.k_high}"
                f"_exact_pairs.tsv.gz"))

    print(f"\n[{anc}] 29b done.")


if __name__ == "__main__":
    main()
