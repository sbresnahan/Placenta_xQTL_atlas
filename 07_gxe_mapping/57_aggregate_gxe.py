#!/usr/bin/env python3
"""
57_aggregate_gxe.py — Aggregate module-07 GxE results (Objective 2.1)

Collects, across modalities x exposures:
  * tier-1 pooled scans ({GXE_DIR}/tier1/pooled_{MOD}_{EXP}.gxe_cis.parquet,
    q-values already added at merge) -> gxe_tier1_results.tsv.gz (all
    phenotypes) and gxe_tier1_significant.tsv (qval <= --fdr);
  * tier-2 ancestry-stratified + IVW meta ({GXE_DIR}/tier2/*.tier2.tsv)
    -> gxe_tier2_results.tsv.gz;
  * component run status (tier1/tier2 done counts, sensitivity OFF/DONE,
    T/NT BLOCKED/DONE) -> gxe_run_status.tsv;
  * per-modality x exposure summary -> gxe_summary.tsv + gxe_summary.png.

Usage:
  python3 57_aggregate_gxe.py --results-dir <qtl_results> [--gxe-dir ...]
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

MODALITIES = ("expression isoforms isoform_expression splicing "
              "intron_retention alt_TSS alt_polyA RNA_editing stability").split()
EXPOSURES = ["GA", "gdm", "ogtt", "ppBMI"]


def log(msg):
    print(f"[{pd.Timestamp.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def parse_mod_exp(basename):
    """pooled_{MOD}_{EXP}.gxe_cis.parquet -> (MOD, EXP) using known lists."""
    m = re.match(r"^pooled_(.+)_([^_]+)\.gxe_cis\.parquet$", basename)
    if not m:
        return None, None
    mod, exp = m.group(1), m.group(2)
    if mod not in MODALITIES:
        return None, None
    return mod, exp


def main():
    p = argparse.ArgumentParser(description="Aggregate GxE results")
    p.add_argument("--results-dir", required=True)
    p.add_argument("--gxe-dir", default=None)
    p.add_argument("--fdr", type=float, default=0.05)
    args = p.parse_args()

    gxe_dir = args.gxe_dir or os.path.join(args.results_dir, "gxe")
    out_dir = os.path.join(gxe_dir, "aggregate")
    os.makedirs(out_dir, exist_ok=True)

    # ---- tier 1 ----
    tier1 = []
    status_rows = []
    for path in sorted(glob.glob(os.path.join(gxe_dir, "tier1",
                                              "pooled_*.gxe_cis.parquet"))):
        mod, exp = parse_mod_exp(os.path.basename(path))
        if mod is None:
            continue
        df = pd.read_parquet(path)
        df = df.reset_index()
        df.insert(0, "exposure", exp)
        df.insert(0, "modality", mod)
        tier1.append(df)
        n_sig = int((df["qval"] <= args.fdr).sum()) if "qval" in df else 0
        status_rows.append(dict(component="tier1", modality=mod, exposure=exp,
                                status="DONE", n_tested=len(df), n_sig=n_sig))
        log(f"  tier1 {mod} x {exp}: {len(df)} phenotypes, {n_sig} at "
            f"q<={args.fdr}")
    if tier1:
        allt = pd.concat(tier1, ignore_index=True)
        allt.to_csv(os.path.join(out_dir, "gxe_tier1_results.tsv.gz"),
                    sep="\t", index=False, compression="gzip")
        sig = allt[allt["qval"] <= args.fdr].sort_values("pval_beta")
        sig.to_csv(os.path.join(out_dir, "gxe_tier1_significant.tsv"),
                   sep="\t", index=False)
        log(f"  tier1 total: {len(allt)} tests, {len(sig)} significant")
    else:
        log("  WARNING: no tier-1 merged parquets found")

    # ---- tier 2 ----
    tier2 = []
    for path in sorted(glob.glob(os.path.join(gxe_dir, "tier2", "*.tier2.tsv"))):
        df = pd.read_csv(path, sep="\t")
        if len(df):
            tier2.append(df)
            status_rows.append(dict(
                component="tier2",
                modality=df["modality"].iloc[0], exposure=df["exposure"].iloc[0],
                status="DONE", n_tested=len(df),
                n_sig=int((df["q_meta"] <= args.fdr).sum())))
    if tier2:
        all2 = pd.concat(tier2, ignore_index=True)
        all2.to_csv(os.path.join(out_dir, "gxe_tier2_results.tsv.gz"),
                    sep="\t", index=False, compression="gzip")
        log(f"  tier2 total: {len(all2)} prioritized phenotypes")

    # ---- component status ----
    for comp, fname in [("sensitivity_snpxcov", "sensitivity/status.tsv"),
                        ("transmitted_nontransmitted", "tnt/status.tsv")]:
        sp = os.path.join(gxe_dir, fname)
        if os.path.exists(sp):
            s = pd.read_csv(sp, sep="\t").iloc[0]
            status_rows.append(dict(component=comp, modality="", exposure="",
                                    status=s["status"], n_tested=np.nan,
                                    n_sig=np.nan))
        else:
            status_rows.append(dict(component=comp, modality="", exposure="",
                                    status="NOT_RUN", n_tested=np.nan,
                                    n_sig=np.nan))
    pd.DataFrame(status_rows).to_csv(
        os.path.join(out_dir, "gxe_run_status.tsv"), sep="\t", index=False)

    # ---- summary table + figure ----
    if tier1:
        allt = pd.concat(tier1, ignore_index=True)
        summ = (allt.groupby(["modality", "exposure"])
                .agg(n_tested=("pval_nominal", "size"),
                     n_sig=("qval", lambda q: int((q <= 0.05).sum())),
                     min_qval=("qval", "min"),
                     min_pval_nominal=("pval_nominal", "min"))
                .reset_index())
        top = (allt.sort_values("pval_beta")
               .groupby(["modality", "exposure"]).head(1)
               [["modality", "exposure", "phenotype_id", "variant_id",
                 "b_gi", "pval_beta"]]
               .rename(columns={"phenotype_id": "top_phenotype",
                                "variant_id": "top_variant",
                                "b_gi": "top_b_gi",
                                "pval_beta": "top_pval_beta"}))
        summ = summ.merge(top, on=["modality", "exposure"], how="left")
        summ.to_csv(os.path.join(out_dir, "gxe_summary.tsv"),
                    sep="\t", index=False)

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["font.family"] = ["Liberation Sans", "Arimo",
                                       "DejaVu Sans"]
        fig, ax = plt.subplots(figsize=(7, 4))
        labels = [f"{m}\n{e}" for m, e in zip(summ["modality"], summ["exposure"])]
        ax.bar(range(len(summ)), summ["n_sig"], color="#0279EE")
        ax.set_xticks(range(len(summ)))
        ax.set_xticklabels(labels, rotation=90, fontsize=7)
        ax.set_ylabel("GxE phenotypes (q <= 0.05)")
        ax.set_xlabel("Modality x exposure")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "gxe_summary.png"), dpi=200)
        plt.close(fig)
        log(f"  wrote gxe_summary.tsv / gxe_summary.png")

    log("done.")


if __name__ == "__main__":
    main()
