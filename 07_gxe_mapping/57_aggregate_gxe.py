#!/usr/bin/env python3
"""57_aggregate_gxe.py — Aggregate ancestry-first Module-07 GxE results."""

import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

MODALITIES = ("expression isoforms isoform_expression splicing "
              "intron_retention alt_TSS alt_polyA RNA_editing stability").split()


def log(msg):
    print(f"[{pd.Timestamp.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def parse_ancestry_mod_exp(basename, allowed_ancestries):
    m = re.match(r"^([^_]+)_(.+)_([^_]+)\.gxe_cis\.parquet$", basename)
    if not m:
        return None, None, None
    anc, mod, exp = m.groups()
    if anc not in allowed_ancestries:
        return None, None, None
    if mod not in MODALITIES:
        return None, None, None
    return anc, mod, exp


def main():
    p = argparse.ArgumentParser(description="Aggregate ancestry-first GxE results")
    p.add_argument("--results-dir", required=True)
    p.add_argument("--gxe-dir", default=None)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    args = p.parse_args()

    gxe_dir = args.gxe_dir or os.path.join(args.results_dir, "gxe")
    out_dir = os.path.join(gxe_dir, "aggregate")
    os.makedirs(out_dir, exist_ok=True)
    status_rows = []

    # Primary discoveries: ancestry-specific full-feature scans.
    tier1 = []
    for path in sorted(glob.glob(os.path.join(gxe_dir, "tier1", "*.gxe_cis.parquet"))):
        anc, mod, exp = parse_ancestry_mod_exp(os.path.basename(path), set(args.ancestries))
        if mod is None:
            continue
        df = pd.read_parquet(path).reset_index()
        df.insert(0, "exposure", exp)
        df.insert(0, "modality", mod)
        df.insert(0, "ancestry", anc)
        tier1.append(df)
        n_sig = int((df["qval"] <= args.fdr).sum()) if "qval" in df else 0
        status_rows.append(dict(component="tier1", ancestry=anc, modality=mod,
                                exposure=exp, status="DONE", n_tested=len(df),
                                n_sig=n_sig))
        log(f"  tier1 {anc} {mod} x {exp}: {len(df)} phenotypes, {n_sig} at q<={args.fdr}")
    if tier1:
        allt = pd.concat(tier1, ignore_index=True)
        allt.to_csv(os.path.join(out_dir, "gxe_tier1_results.tsv.gz"),
                    sep="\t", index=False, compression="gzip")
        sig = allt[allt["qval"] <= args.fdr].sort_values(["ancestry", "pval_beta"])
        sig.to_csv(os.path.join(out_dir, "gxe_tier1_significant.tsv"),
                   sep="\t", index=False)
        log(f"  tier1 total: {len(allt)} ancestry-specific tests, {len(sig)} significant")
    else:
        allt = pd.DataFrame()
        log("  WARNING: no ancestry-specific tier-1 merged parquets found")

    # Cross-ancestry synthesis keeps ancestry-only features and exact-lead IVW.
    meta_frames = []
    for path in sorted(glob.glob(os.path.join(gxe_dir, "tier1_meta", "*.meta.tsv.gz"))):
        df = pd.read_csv(path, sep="\t")
        if len(df):
            meta_frames.append(df)
            status_rows.append(dict(
                component="tier1_meta", ancestry="cross", modality=df["modality"].iloc[0],
                exposure=df["exposure"].iloc[0], status="DONE", n_tested=len(df),
                n_sig=int((pd.to_numeric(df.get("q_feature_acat"), errors="coerce") <= args.fdr).sum())))
    if meta_frames:
        meta_all = pd.concat(meta_frames, ignore_index=True)
        meta_all.to_csv(os.path.join(out_dir, "gxe_tier1_meta.tsv.gz"),
                        sep="\t", index=False, compression="gzip")
        log(f"  tier1 meta/synthesis: {len(meta_all)} phenotype rows")

    # Aim-1 prioritized tier 2.
    tier2 = []
    for path in sorted(glob.glob(os.path.join(gxe_dir, "tier2", "*.tier2.tsv"))):
        df = pd.read_csv(path, sep="\t")
        if len(df):
            tier2.append(df)
            status_rows.append(dict(
                component="tier2", ancestry="cross", modality=df["modality"].iloc[0],
                exposure=df["exposure"].iloc[0], status="DONE", n_tested=len(df),
                n_sig=int((df["q_meta"] <= args.fdr).sum())))
    if tier2:
        all2 = pd.concat(tier2, ignore_index=True)
        all2.to_csv(os.path.join(out_dir, "gxe_tier2_results.tsv.gz"),
                    sep="\t", index=False, compression="gzip")
        log(f"  tier2 total: {len(all2)} prioritized phenotypes")

    for comp, fname in [("sensitivity_snpxcov", "sensitivity/status.tsv"),
                        ("transmitted_nontransmitted", "tnt/status.tsv")]:
        sp = os.path.join(gxe_dir, fname)
        if os.path.exists(sp):
            s = pd.read_csv(sp, sep="\t").iloc[0]
            status_rows.append(dict(component=comp, ancestry="", modality="", exposure="",
                                    status=s["status"], n_tested=np.nan, n_sig=np.nan))
        else:
            status_rows.append(dict(component=comp, ancestry="", modality="", exposure="",
                                    status="NOT_RUN", n_tested=np.nan, n_sig=np.nan))
    pd.DataFrame(status_rows).to_csv(
        os.path.join(out_dir, "gxe_run_status.tsv"), sep="\t", index=False)

    if len(allt):
        summ = (allt.groupby(["ancestry", "modality", "exposure"])
                .agg(n_tested=("pval_nominal", "size"),
                     n_sig=("qval", lambda q: int((q <= args.fdr).sum())),
                     min_qval=("qval", "min"),
                     min_pval_beta=("pval_beta", "min"))
                .reset_index())
        top = (allt.sort_values("pval_beta")
               .groupby(["ancestry", "modality", "exposure"]).head(1)
               [["ancestry", "modality", "exposure", "phenotype_id", "variant_id",
                 "b_gi", "pval_beta"]]
               .rename(columns={"phenotype_id": "top_phenotype",
                                "variant_id": "top_variant",
                                "b_gi": "top_b_gi",
                                "pval_beta": "top_pval_beta"}))
        summ = summ.merge(top, on=["ancestry", "modality", "exposure"], how="left")
        summ.to_csv(os.path.join(out_dir, "gxe_summary.tsv"), sep="\t", index=False)
        log("  wrote gxe_summary.tsv")

    log("done.")


if __name__ == "__main__":
    main()
