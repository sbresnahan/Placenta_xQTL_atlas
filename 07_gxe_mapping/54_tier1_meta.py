#!/usr/bin/env python3
"""
54_tier1_meta.py — Cross-ancestry synthesis of ancestry-specific tier-1 GxE scans.

Primary discovery remains ancestry-specific and therefore retains each ancestry's
full feature space.  This script never filters those discovery results down to an
EAS∩EUR phenotype set.  Instead it outer-joins the merged tier-1 results by
canonical phenotype_id and reports:

  * ancestry-only phenotypes unchanged;
  * shared phenotypes with per-ancestry feature-level permutation/Beta p-values;
  * ACAT combination of pval_beta across ancestries as a feature-level
    cross-ancestry evidence score (the ancestry scans already account for the
    within-cis variant search);
  * fixed-effect IVW beta/se + Cochran Q/I2 only when every contributing
    ancestry selected the exact same lead variant_id.  If lead variants differ,
    no effect meta-analysis is claimed.

This is intentionally conservative: it does not pretend that two different
ancestry-specific lead variants are the same feature-variant test.  Exhaustive
variant-level cross-ancestry meta-analysis would require a second full nominal
scan of every shared cis variant and is not used as the primary discovery path.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gxe_core


def log(msg):
    print(f"[{pd.Timestamp.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def acat(pvalues):
    p = np.asarray(pvalues, dtype=float)
    p = p[np.isfinite(p)]
    if len(p) == 0:
        return np.nan
    p = np.clip(p, 1e-15, 1 - 1e-15)
    t = np.mean(np.tan((0.5 - p) * np.pi))
    if t > 1e15:
        return float(1.0 / (np.pi * t))
    return float(0.5 - np.arctan(t) / np.pi)


def bh_qvalues(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full(len(p), np.nan, dtype=float)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    x = p[ok]
    order = np.argsort(x)
    ranks = np.arange(1, len(x) + 1)
    q = x[order] * len(x) / ranks
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(len(x), dtype=float)
    tmp[order] = q
    out[ok] = tmp
    return out


def read_one(path, ancestry):
    if not os.path.exists(path):
        log(f"  WARNING: missing {path}")
        return None
    df = pd.read_parquet(path).reset_index()
    if "phenotype_id" not in df.columns:
        sys.exit(f"ERROR: {path} lacks phenotype_id")
    if df["phenotype_id"].duplicated().any():
        sys.exit(f"ERROR: duplicate phenotype_id rows in {path}")
    keep = [c for c in [
        "phenotype_id", "variant_id", "b_gi", "b_gi_se", "pval_nominal",
        "pval_perm", "pval_beta", "qval", "af", "ma_count", "nperm", "dof"
    ] if c in df.columns]
    out = df[keep].copy()
    rename = {c: f"{c}_{ancestry}" for c in keep if c != "phenotype_id"}
    return out.rename(columns=rename)


def main():
    p = argparse.ArgumentParser(description="Cross-ancestry synthesis of tier-1 GxE")
    p.add_argument("--gxe-dir", required=True)
    p.add_argument("--modality", required=True)
    p.add_argument("--exposure", required=True)
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    p.add_argument("--fdr", type=float, default=0.05)
    args = p.parse_args()

    tier1 = os.path.join(args.gxe_dir, "tier1")
    out_dir = os.path.join(args.gxe_dir, "tier1_meta")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{args.modality}_{args.exposure}.meta.tsv.gz")
    if os.path.exists(out_path) and os.environ.get("FORCE", "0") != "1":
        log(f"{out_path} exists — skipping (FORCE=1 to redo)")
        return

    frames = []
    present = []
    for anc in args.ancestries:
        path = os.path.join(tier1, f"{anc}_{args.modality}_{args.exposure}.gxe_cis.parquet")
        df = read_one(path, anc)
        if df is not None:
            frames.append(df)
            present.append(anc)
    if not frames:
        sys.exit("ERROR: no ancestry-specific tier-1 merged results found")

    merged = frames[0]
    for df in frames[1:]:
        merged = merged.merge(df, on="phenotype_id", how="outer", validate="one_to_one")

    rows = []
    for _, r in merged.iterrows():
        available = [a for a in present if pd.notna(r.get(f"pval_beta_{a}"))]
        variants = [str(r[f"variant_id_{a}"]) for a in available
                    if pd.notna(r.get(f"variant_id_{a}"))]
        pvals = [float(r[f"pval_beta_{a}"]) for a in available
                 if pd.notna(r.get(f"pval_beta_{a}"))]
        rec = dict(
            modality=args.modality,
            exposure=args.exposure,
            phenotype_id=r["phenotype_id"],
            n_ancestries=len(available),
            ancestries=";".join(available),
            p_feature_acat=acat(pvals),
            significant_any_ancestry=any(
                pd.notna(r.get(f"qval_{a}")) and float(r[f"qval_{a}"]) <= args.fdr
                for a in available),
            significant_all_ancestries=(len(available) >= 2 and all(
                pd.notna(r.get(f"qval_{a}")) and float(r[f"qval_{a}"]) <= args.fdr
                for a in available)),
            lead_concordant=(len(available) >= 2 and len(set(variants)) == 1),
        )
        if len(available) == 1:
            rec["scope"] = "ancestry_only"
        elif rec["lead_concordant"]:
            rec["scope"] = "shared_concordant_lead"
        else:
            rec["scope"] = "shared_discordant_lead"

        if rec["lead_concordant"]:
            beta = np.array([float(r[f"b_gi_{a}"]) for a in available])
            se = np.array([float(r[f"b_gi_se_{a}"]) for a in available])
            m = gxe_core.ivw_meta(beta, se)
            rec.update(meta_variant_id=variants[0], beta_meta=m["beta"],
                       se_meta=m["se"], z_meta=m["z"], p_meta=m["p"],
                       q_het=m["q_het"], p_het=m["p_het"], i2=m["i2"])
        else:
            rec.update(meta_variant_id=np.nan, beta_meta=np.nan, se_meta=np.nan,
                       z_meta=np.nan, p_meta=np.nan, q_het=np.nan,
                       p_het=np.nan, i2=np.nan)

        for anc in present:
            for c in ["variant_id", "b_gi", "b_gi_se", "pval_nominal",
                      "pval_perm", "pval_beta", "qval", "af", "ma_count",
                      "nperm", "dof"]:
                rec[f"{c}_{anc}"] = r.get(f"{c}_{anc}", np.nan)
        rows.append(rec)

    out = pd.DataFrame(rows)
    shared = out["n_ancestries"] >= 2
    out["q_feature_acat"] = np.nan
    if shared.any():
        out.loc[shared, "q_feature_acat"] = bh_qvalues(
            out.loc[shared, "p_feature_acat"].values)
    out = out.sort_values(["scope", "p_feature_acat"], na_position="last")
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")
    log(f"wrote {out_path}: {len(out)} phenotypes; "
        f"{int(shared.sum())} shared; {int(out['lead_concordant'].sum())} shared with concordant lead")


if __name__ == "__main__":
    main()
