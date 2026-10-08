#!/usr/bin/env python3
"""
55_sensitivity_snpxcov.py — SNP x covariate sensitivity refits (Objective 2.1)

Aims: "As sensitivity, we will refit interaction models adjusting for
SNP x covariate interactions (maternal age, parity, smoking, gestational
age at delivery, ancestry PCs); stability of b_int would imply any residual
unmeasured confounder must interact with G orthogonally to all measured
factors."

For every tier-1 significant hit (qval <= --fdr), refit the pooled model
    Y ~ G + E + G:E + sum_j G:C_j + Z_cov
where C_j are the enabled sensitivity covariates (sensitivity_config.tsv;
new columns also enter as main effects). Reports b_int before/after, the
absolute and relative change, both p-values, and a stability flag
(|delta b_int / b_int| <= 0.2 and same-significance call at alpha=0.05).

STATUS-GATED: with no enabled covariates (the default — maternal age,
parity and smoking are not yet in the cohort metadata) the script writes
{sensitivity_dir}/status.tsv with status=OFF and exits 0.

Usage:
  python3 55_sensitivity_snpxcov.py --qtl-dir ... --results-dir ... \
      --gxe-dir ... [--exposures GA] [--fdr 0.05]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))


def log(msg):
    print(f"[{pd.Timestamp.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def ols_b_se(y, X, idx):
    """OLS via lstsq; returns (b, se, p) for column idx."""
    beta, res, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = max(len(y) - rank, 1)
    sigma2 = float(resid @ resid) / dof
    XtX_inv = np.linalg.pinv(X.T @ X)
    se = np.sqrt(max(XtX_inv[idx, idx] * sigma2, 1e-300))
    t = beta[idx] / se
    return float(beta[idx]), float(se), float(2 * stats.t.sf(abs(t), dof))


def main():
    p = argparse.ArgumentParser(description="SNP x covariate sensitivity refits")
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--gxe-dir", required=True)
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    p.add_argument("--exposures", nargs="+", default=None,
                   help="default: enabled rows of gxe_config.tsv")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--config", default=os.path.join(HERE, "sensitivity_config.tsv"))
    args = p.parse_args()

    out_dir = os.path.join(args.gxe_dir, "sensitivity")
    os.makedirs(out_dir, exist_ok=True)
    status_path = os.path.join(out_dir, "status.tsv")

    cfg = pd.read_csv(args.config, sep="\t")
    enabled = cfg[cfg["enabled"] == 1]
    if len(enabled) == 0:
        pd.DataFrame([dict(component="sensitivity_snpxcov", status="OFF",
                           detail="no enabled covariates in "
                                  "sensitivity_config.tsv (maternal age, "
                                  "parity, smoking pending metadata)")]
                     ).to_csv(status_path, sep="\t", index=False)
        log("no enabled sensitivity covariates — status=OFF (see status.tsv)")
        return

    # enabled path: refit tier-1 significant hits with GxC terms
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "gxe_scan_52", os.path.join(HERE, "52_gxe_scan.py"))
    scanner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scanner)

    meta = pd.read_csv(os.path.join(args.gxe_dir, "inputs",
                                    "pooled_sample_manifest.tsv"), sep="\t")
    meta = meta.set_index("array_id")
    missing_cols = [c for c in enabled["column"] if c and c not in meta.columns]
    if missing_cols:
        pd.DataFrame([dict(component="sensitivity_snpxcov", status="BLOCKED",
                           detail=f"metadata columns missing: {missing_cols}")]
                     ).to_csv(status_path, sep="\t", index=False)
        log(f"BLOCKED: metadata columns missing: {missing_cols}")
        return

    hits_path = os.path.join(args.gxe_dir, "aggregate",
                             "gxe_tier1_significant.tsv")
    if not os.path.exists(hits_path):
        sys.exit(f"ERROR: {hits_path} not found — run 57_aggregate_gxe.py first")
    hits = pd.read_csv(hits_path, sep="\t")
    log(f"{len(hits)} tier-1 significant hits to refit")

    # Build the pooled covariate main-effect matrix once per modality
    results = []
    for (mod, exp_id), h in hits.groupby(["modality", "exposure"]):
        cov = pd.read_csv(os.path.join(args.gxe_dir, "inputs",
                                       f"pooled_covariates_{mod}.tsv"),
                          sep="\t", index_col=0)
        if exp_id in cov.index:
            cov = cov.drop(index=exp_id)
        exp = pd.read_csv(os.path.join(args.gxe_dir, "inputs", "exposures.tsv"),
                          sep="\t", index_col=0).loc[exp_id]
        bed_path = os.path.join(args.gxe_dir, "inputs", f"pooled_{mod}.bed.gz")
        import tensorqtl
        pheno, pheno_pos = tensorqtl.read_phenotype_bed(bed_path)
        for _, hit in h.iterrows():
            pid, vid = hit["phenotype_id"], hit["variant_id"]
            if pid not in pheno.index:
                continue
            y_full = pheno.loc[pid]
            samples = [s for s in y_full.index
                       if s in cov.columns and pd.notna(exp.get(s))]
            y = y_full[samples].values.astype(float)
            e = exp[samples].values.astype(float)
            Z = cov[samples].T.values
            # genotype at the hit variant (pooled, per-ancestry blocks)
            chrom = str(hit.get("chrom", "")) or vid.split(":")[0]
            man = pd.read_csv(os.path.join(args.gxe_dir, "inputs",
                                           "pooled_sample_manifest.tsv"),
                              sep="\t")
            anc_of = dict(zip(man["array_id"], man["ancestry"]))
            sba = {a: [s for s in samples if anc_of.get(s) == a]
                   for a in args.ancestries}
            pgens = {a: os.path.join(args.qtl_dir, f"{a}_qtl")
                     for a in args.ancestries}
            G, vids, _, _ = scanner.load_chromosome_genotypes(
                pgens, sba, chrom, maf_threshold=0.0)
            gi = np.where(vids == vid)[0]
            if len(gi) == 0:
                continue
            g = G[gi[0]]
            # extra covariate main effects + GxC terms
            extra_main, gxc = [], []
            for _, crow in enabled.iterrows():
                cname, ccol = crow["covariate"], crow["column"]
                if cname == "ancestry_pcs":
                    pc_rows = [i for i in cov.index if "_PC" in i]
                    C = cov.loc[pc_rows, samples].T.values
                else:
                    cv = pd.to_numeric(meta.loc[samples, ccol],
                                       errors="coerce")
                    cv = ((cv - cv.mean()) / cv.std(ddof=0)).fillna(0.0)
                    C = cv.values.reshape(-1, 1)
                    extra_main.append(C)
                gxc.append(g.reshape(-1, 1) * C)
            X0 = np.column_stack([np.ones(len(samples)), g, e, g * e, Z])
            X1 = np.column_stack([X0] + extra_main + gxc)
            b0, se0, p0 = ols_b_se(y, X0, 3)
            b1, se1, p1 = ols_b_se(y, X1, 3)
            results.append(dict(
                modality=mod, exposure=exp_id, phenotype_id=pid,
                variant_id=vid, b_int_base=b0, se_base=se0, p_base=p0,
                b_int_adj=b1, se_adj=se1, p_adj=p1,
                delta_b=b1 - b0,
                rel_change=abs(b1 - b0) / abs(b0) if b0 != 0 else np.nan,
                stable=bool((abs(b1 - b0) <= 0.2 * abs(b0))
                            and ((p0 < 0.05) == (p1 < 0.05)))))
    out = pd.DataFrame(results)
    out_path = os.path.join(out_dir, "sensitivity_snpxcov.tsv")
    out.to_csv(out_path, sep="\t", index=False)
    pd.DataFrame([dict(component="sensitivity_snpxcov", status="DONE",
                       detail=f"{len(out)} refits; "
                              f"{int(out['stable'].sum()) if len(out) else 0} stable")]
                 ).to_csv(status_path, sep="\t", index=False)
    log(f"wrote {out_path} ({len(out)} refits)")


if __name__ == "__main__":
    main()
