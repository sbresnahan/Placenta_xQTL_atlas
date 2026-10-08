#!/usr/bin/env python3
"""
56_transmitted_nontransmitted.py — transmitted/non-transmitted allele
decomposition at top GxE loci (Objective 2.1)

Genetic-nurture guard: a fetal SNP x exposure interaction could reflect the
MATERNAL genotype (which influences both the exposure and, via transmission,
the fetal genotype) rather than direct fetal-genetic regulation. In
mother-child duos (the aims name 5 cohorts with maternal genotypes: GUSTO,
Gen3G, ITU, ELGAN, RICHS), the maternal genotype g_m decomposes into the
transmitted allele T and the non-transmitted allele NT = g_m - T
(gxe_core.transmitted_dosage; T is deterministic from Mendelian rules
except mother-het/child-het, where T = 0.5 in expectation). The model

    Y ~ T + NT + E + T:E + NT:E + Z_cov

separates the direct fetal GxE (T:E) from the maternal-mediated component
(NT:E); b(T:E) stable with b(NT:E) ~ 0 supports direct fetal regulation.

STATUS-GATED: maternal genotypes are not yet processed in this project, so
unless --maternal-pgen-dir and --pairs are provided (and exist), the script
writes {tnt_dir}/status.tsv with status=BLOCKED and exits 0. The
decomposition logic itself is fixture-tested in test_gxe_scan.py.

Input contract when unblocked:
  --maternal-pgen-dir : directory with {COHORT}_maternal.pgen/.pvar/.psam
                        (ALT-count dosages matching the child variant IDs)
  --pairs             : TSV with columns child_array_id, mother_id, cohort

Usage:
  python3 56_transmitted_nontransmitted.py --qtl-dir ... --results-dir ... \
      --gxe-dir ... [--maternal-pgen-dir D --pairs P] [--fdr 0.05]
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gxe_core


def log(msg):
    print(f"[{pd.Timestamp.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def load_tier1_hits(gxe_dir, fdr, ancestries):
    """Load significant ancestry-specific merged tier-1 parquets directly."""
    rows = []
    for path in sorted(glob.glob(os.path.join(gxe_dir, "tier1", "*.gxe_cis.parquet"))):
        name = os.path.basename(path)
        m = re.match(r"^([^_]+)_(.+)_([^_]+)\.gxe_cis\.parquet$", name)
        if not m:
            continue
        anc, mod, exp = m.groups()
        if anc not in set(ancestries):
            continue
        df = pd.read_parquet(path).reset_index()
        if "qval" not in df.columns:
            continue
        df = df[df["qval"] <= fdr].copy()
        if len(df):
            df.insert(0, "exposure", exp)
            df.insert(0, "modality", mod)
            df.insert(0, "ancestry", anc)
            rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def ols(y, X):
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = max(len(y) - rank, 1)
    sigma2 = float(resid @ resid) / dof
    XtX_inv = np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.maximum(np.diag(XtX_inv) * sigma2, 1e-300))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    return beta, se, p


def main():
    p = argparse.ArgumentParser(description="T/NT decomposition at top GxE loci")
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--gxe-dir", required=True)
    p.add_argument("--maternal-pgen-dir", default=None)
    p.add_argument("--pairs", default=None,
                   help="TSV: child_array_id, mother_id, cohort")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    args = p.parse_args()

    out_dir = os.path.join(args.gxe_dir, "tnt")
    os.makedirs(out_dir, exist_ok=True)
    status_path = os.path.join(out_dir, "status.tsv")

    if (not args.maternal_pgen_dir or not args.pairs
            or not os.path.isdir(args.maternal_pgen_dir or "")
            or not os.path.exists(args.pairs or "")):
        pd.DataFrame([dict(
            component="transmitted_nontransmitted", status="BLOCKED",
            detail="maternal pgens (--maternal-pgen-dir) and/or mother-child "
                   "pairs file (--pairs) not provided; the decomposition is "
                   "implemented and fixture-tested — supply duo genotypes "
                   "for GUSTO/Gen3G/ITU/ELGAN/RICHS to enable")]
        ).to_csv(status_path, sep="\t", index=False)
        log("BLOCKED: maternal genotypes/pairs not available — see status.tsv")
        return

    hits = load_tier1_hits(args.gxe_dir, args.fdr, args.ancestries)
    if hits.empty:
        log(f"no ancestry-specific tier-1 hits at q<={args.fdr}; nothing to refit")
        hits = pd.DataFrame(columns=["ancestry", "modality", "exposure",
                                     "phenotype_id", "variant_id"])
    pairs = pd.read_csv(args.pairs, sep="\t")
    log(f"{len(hits)} hits x {len(pairs)} duos")

    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "gxe_scan_52", os.path.join(HERE, "52_gxe_scan.py"))
    scanner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scanner)
    from tensorqtl import pgen as tq_pgen
    import tensorqtl

    results = []
    if "ancestry" not in hits.columns:
        sys.exit("ERROR: gxe_tier1_significant.tsv lacks ancestry; rerun ancestry-first aggregation")
    for (anc, mod, exp_id), h in hits.groupby(["ancestry", "modality", "exposure"]):
        cov = pd.read_csv(os.path.join(args.gxe_dir, "inputs",
                                       f"{anc}_covariates_{mod}.tsv"),
                          sep="\t", index_col=0)
        if exp_id in cov.index:
            cov = cov.drop(index=exp_id)
        exp = pd.read_csv(os.path.join(args.gxe_dir, "inputs", "exposures.tsv"),
                          sep="\t", index_col=0).loc[exp_id]
        pheno, _ = tensorqtl.read_phenotype_bed(
            os.path.join(args.gxe_dir, "inputs", f"{anc}_{mod}.bed.gz"))
        for _, hit in h.iterrows():
            pid, vid = hit["phenotype_id"], hit["variant_id"]
            if pid not in pheno.index:
                continue
            chrom = vid.split(":")[0]
            duo_children = [s for s in pairs["child_array_id"] if s in pheno.columns]
            G, vids, _, _ = scanner.load_chromosome_genotypes(
                {anc: os.path.join(args.qtl_dir, f"{anc}_qtl")},
                {anc: duo_children}, chrom, maf_threshold=0.0)
            gi = np.where(vids == vid)[0]
            if len(gi) == 0:
                continue
            g_child = pd.Series(G[gi[0]], index=duo_children)
            # maternal genotypes at the same variant, per cohort pgen
            g_mother = {}
            for cohort, sub in pairs.groupby("cohort"):
                prefix = os.path.join(args.maternal_pgen_dir,
                                      f"{cohort}_maternal")
                if not os.path.exists(prefix + ".pgen"):
                    log(f"  WARNING: maternal pgen missing for {cohort}")
                    continue
                pvar = tq_pgen.read_pvar(prefix + ".pvar")
                mrow = pvar.index[pvar["id"] == vid]
                if len(mrow) == 0:
                    continue
                gm = tq_pgen.read_list(prefix + ".pgen",
                                       np.array(mrow, dtype=np.uint32))[0]
                iids = scanner.read_psam_iids(prefix + ".psam")
                gm = pd.Series(np.where(gm == -9, np.nan, gm).astype(float),
                               index=iids)
                for _, r in sub.iterrows():
                    if r["mother_id"] in gm.index:
                        g_mother[r["child_array_id"]] = gm[r["mother_id"]]
            duo = pairs[pairs["child_array_id"].isin(g_child.index)
                        & pairs["child_array_id"].isin(g_mother)]
            if len(duo) < 30:
                log(f"  {pid} x {vid}: only {len(duo)} duos — skipping")
                continue
            children = duo["child_array_id"].tolist()
            gc = g_child[children].values
            gm = np.array([g_mother[c] for c in children])
            T = gxe_core.transmitted_dosage(gm, gc)
            NT = gm - T
            ok = np.isfinite(T) & np.isfinite(NT)
            y = pheno.loc[pid, children].values.astype(float)
            e = exp[children].values.astype(float)
            Z = cov[children].T.values
            ok &= np.isfinite(y) & np.isfinite(e)
            if ok.sum() < 30:
                continue
            X = np.column_stack([np.ones(ok.sum()), T[ok], NT[ok], e[ok],
                                 T[ok] * e[ok], NT[ok] * e[ok], Z[ok]])
            beta, se, pv = ols(y[ok], X)
            results.append(dict(
                ancestry=anc, modality=mod, exposure=exp_id, phenotype_id=pid,
                variant_id=vid, n_duos=int(ok.sum()),
                b_T=beta[1], se_T=se[1], p_T=pv[1],
                b_NT=beta[2], se_NT=se[2], p_NT=pv[2],
                b_TxE=beta[4], se_TxE=se[4], p_TxE=pv[4],
                b_NTxE=beta[5], se_NTxE=se[5], p_NTxE=pv[5],
                direct_fetal=bool(pv[4] < 0.05 and pv[5] >= 0.05)))
    out = pd.DataFrame(results)
    out_path = os.path.join(out_dir, "tnt_results.tsv")
    out.to_csv(out_path, sep="\t", index=False)
    pd.DataFrame([dict(component="transmitted_nontransmitted", status="DONE",
                       detail=f"{len(out)} locus refits")]
                 ).to_csv(status_path, sep="\t", index=False)
    log(f"wrote {out_path} ({len(out)} refits)")


if __name__ == "__main__":
    main()
