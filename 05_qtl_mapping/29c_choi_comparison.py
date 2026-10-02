#!/usr/bin/env python3
"""
29c_choi_comparison.py — External-reference comparison (Choi 2024 SNUH)
for the within-cohort-INT validation gate (module 05; runs after 29).

Revision of the ad-hoc prepare_choi_vs_eas.py, generalized over ancestry
labels, plus the retention/catalog generator that feeds the retention
plotter and the 29d within-cohort scan of reference-defined pairs.

Inputs:
  --choi              Choi 2024 full cis summary stats (gzip TSV; columns
                      gene_id, variant_id [chr01_000778597_C_T], slope,
                      slope_se, mlogp, PIP, credible_set)
  --choi-significant  Choi significant eGene list (TSV/TXT; the gene column
                      is auto-detected among gene_id/gene/phenotype_id/
                      gene_name, else the first column; version suffixes
                      stripped)
  --top-table         {ANC}_expression_cisqtl_top.tsv (module 29 output;
                      needs phenotype_id, variant_id, slope, slope_se,
                      pval_nominal, pval_beta, qval)
  --out-prefix        e.g. $DIAG_DIR/EAS/Choi_vs_pooled_EAS

Outputs (<out-prefix>.*):
  .gene_top_comparison.tsv.gz       per-gene best cis variant, both studies
  .exact_{ANC}_lead_comparison.tsv.gz  our lead pairs evaluated in Choi
  .Choi_significant_retention.tsv.gz   mlog10_pbeta_Choi, mlog10_pbeta_{ANC},
                                       retention label (per Choi eGene
                                       tested in both)
  .common_tested_catalog.tsv.gz     per-gene category: Shared significant /
                                    Choi-only / Pooled-{ANC}-only / Neither
  .choi_defined_pairs.tsv           Choi lead SNP per Choi-significant gene
                                    tested in both (phenotype_id in our
                                    namespace; variant_id normalized to our
                                    chr:pos:ref:alt convention) — pair set B
                                    for 29d

Notes:
  * mlog10_pbeta_Choi is the per-gene max of Choi's per-variant mlogp (the
    full summary stats carry no permutation-adjusted gene-level p); it is a
    best-variant nominal signal, used only for the retention scatter.
  * mlog10_pbeta_{ANC} is -log10(pval_beta) of the gene lead from the top
    table (falls back to -log10(pval_nominal) with a warning if pval_beta
    is absent).
  * Gene significance on our side: qval <= --fdr (Storey q from 27/29;
    falls back to BH on pval_beta with a warning).

Usage:
  python3 29c_choi_comparison.py --ancestry EAS \
      --choi /path/to/choi_full_summary_stats.tsv.gz \
      --choi-significant /path/to/choi_significant_egenes.tsv \
      --top-table $RESULTS_DIR/EAS_expression_cisqtl_top.tsv \
      --out-prefix $DIAG_DIR/EAS/Choi_vs_pooled_EAS
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


def gene_base(x):
    return x.astype(str).str.replace(r"\.\d+$", "", regex=True)


def normalize_our_variant(v):
    # Ours: 1:63613970:T:A
    out = []
    for s in v.astype(str):
        p = s.split(":")
        if len(p) < 4:
            out.append(s)
            continue
        chrom = p[0].replace("chr", "")
        if chrom.isdigit():
            chrom = str(int(chrom))
        pos = str(int(p[1]))
        out.append(":".join([chrom, pos, p[2], p[3]]))
    return pd.Series(out, index=v.index)


def normalize_choi_variant(v):
    # Choi: chr01_000778597_C_T
    out = []
    for s in v.astype(str):
        p = s.split("_")
        if len(p) < 4:
            out.append(s)
            continue
        chrom = p[0].replace("chr", "")
        if chrom.isdigit():
            chrom = str(int(chrom))
        pos = str(int(p[1]))
        out.append(":".join([chrom, pos, p[2], p[3]]))
    return pd.Series(out, index=v.index)


def read_gene_list(path):
    # Detect the delimiter from the header line explicitly: pandas' sep=None
    # sniffer misparses single-column files (e.g. splits "gene_id" on 'e').
    with open(path) as fh:
        header = ""
        for line in fh:
            if line.strip() and not line.startswith("#"):
                header = line.rstrip("\n")
                break
    sep = "\t" if "\t" in header else ("," if "," in header else None)
    df = pd.read_csv(path, sep=sep or "\t", comment="#")
    if df.shape[1] == 0:
        sys.exit(f"ERROR: no columns parsed in {path}")
    col = next((c for c in df.columns
                if c.lower() in ("gene_id", "gene", "phenotype_id",
                                 "gene_name", "ensembl_id", "egene")),
               df.columns[0])
    genes = set(gene_base(df[col].dropna()))
    print(f"  Read {len(genes)} genes from {path} (column '{col}')")
    return genes


def bh_qvalues(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Choi 2024 vs pooled-ancestry eQTL comparison")
    ap.add_argument("--ancestry", required=True)
    ap.add_argument("--choi", required=True)
    ap.add_argument("--choi-significant", required=True)
    ap.add_argument("--top-table", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--chunksize", type=int, default=500000)
    args = ap.parse_args()

    anc = args.ancestry
    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)),
                exist_ok=True)

    # ---- our top table ----
    ours = pd.read_csv(args.top_table, sep="\t")
    for col in ("phenotype_id", "variant_id", "slope", "slope_se",
                "pval_nominal"):
        if col not in ours.columns:
            sys.exit(f"ERROR: {args.top_table} missing column '{col}'")
    ours["gene_base"] = gene_base(ours["phenotype_id"])
    ours["variant_norm"] = normalize_our_variant(ours["variant_id"])
    ours[f"t_{anc}"] = ours["slope"] / ours["slope_se"]
    ours[f"abs_t_{anc}"] = ours[f"t_{anc}"].abs()
    ours[f"mlogp_{anc}"] = -np.log10(pd.to_numeric(ours["pval_nominal"],
                                                 errors="coerce"))

    # Gene-level significance on our side
    if "qval" in ours.columns and ours["qval"].notna().any():
        ours["qval_gene"] = pd.to_numeric(ours["qval"], errors="coerce")
        sig_source = "qval"
    elif "pval_beta" in ours.columns:
        print("  WARN: no qval column; BH on pval_beta for gene significance")
        ours["qval_gene"] = bh_qvalues(pd.to_numeric(
            ours["pval_beta"], errors="coerce").fillna(1.0))
        sig_source = "BH(pval_beta)"
    else:
        sys.exit("ERROR: top table has neither qval nor pval_beta")
    ours[f"sig_{anc}"] = ours["qval_gene"] <= args.fdr
    if "pval_beta" in ours.columns:
        ours[f"mlog10_pbeta_{anc}"] = -np.log10(pd.to_numeric(
            ours["pval_beta"], errors="coerce").clip(lower=1e-300))
    else:
        print("  WARN: no pval_beta column; using pval_nominal for the "
              "retention scatter y-axis")
        ours[f"mlog10_pbeta_{anc}"] = ours[f"mlogp_{anc}"]
    print(f"  {anc}: {len(ours)} genes tested, "
          f"{int(ours[f'sig_{anc}'].sum())} significant ({sig_source} "
          f"q <= {args.fdr})")

    choi_sig = read_gene_list(args.choi_significant)

    # One expression row per gene in our top table.
    keys = ours[["gene_base", "variant_norm"]].drop_duplicates()
    genes = set(keys["gene_base"])

    usecols = ["gene_id", "variant_id", "slope", "slope_se", "mlogp",
               "PIP", "credible_set"]

    top_candidates = []
    exact_matches = []
    n = 0
    reader = pd.read_csv(args.choi, sep="\t", compression="gzip",
                         usecols=usecols, chunksize=args.chunksize,
                         low_memory=False, na_values=["."])
    for chunk in reader:
        n += len(chunk)
        chunk["gene_base"] = gene_base(chunk["gene_id"])
        chunk = chunk[chunk["gene_base"].isin(genes)].copy()
        if chunk.empty:
            continue
        chunk["variant_norm"] = normalize_choi_variant(chunk["variant_id"])
        chunk["mlogp"] = pd.to_numeric(chunk["mlogp"], errors="coerce")
        chunk["slope"] = pd.to_numeric(chunk["slope"], errors="coerce")
        chunk["slope_se"] = pd.to_numeric(chunk["slope_se"], errors="coerce")
        chunk["t_Choi"] = chunk["slope"] / chunk["slope_se"]
        chunk["abs_t_Choi"] = chunk["t_Choi"].abs()

        valid = chunk["mlogp"].notna()
        if valid.any():
            idx = chunk.loc[valid].groupby("gene_base")["mlogp"].idxmax()
            top_candidates.append(chunk.loc[idx].copy())

        m = chunk.merge(keys, on=["gene_base", "variant_norm"], how="inner")
        if len(m):
            exact_matches.append(m)

        if n % 5000000 < args.chunksize:
            print(f"  scanned {n:,} Choi rows")

    if not top_candidates:
        sys.exit("ERROR: no overlapping genes between the Choi summary "
                 "stats and our top table — check gene ID conventions")

    choi_top = pd.concat(top_candidates, ignore_index=True)
    idx = choi_top.groupby("gene_base")["mlogp"].idxmax()
    choi_top = choi_top.loc[idx].copy()

    choi_cols = ["gene_base", "gene_id", "variant_id", "variant_norm",
                 "slope", "slope_se", "t_Choi", "abs_t_Choi",
                 "mlogp", "PIP", "credible_set"]
    top_cmp = ours.merge(choi_top[choi_cols], on="gene_base", how="inner",
                         suffixes=(f"_{anc}", "_Choi"))
    top_cmp["same_lead_variant"] = (
        top_cmp[f"variant_norm_{anc}"] == top_cmp["variant_norm_Choi"])

    if exact_matches:
        exact = pd.concat(exact_matches, ignore_index=True).drop_duplicates(
            ["gene_base", "variant_norm"])
    else:
        # No exact lead-variant overlaps (possible for some ancestries);
        # keep an empty, correctly-shaped frame so downstream code works.
        exact = pd.DataFrame(columns=choi_cols)
    exact_cmp = ours.merge(exact[choi_cols], on=["gene_base", "variant_norm"],
                           how="inner", suffixes=(f"_{anc}", "_Choi"))
    if f"slope_{anc}" not in exact_cmp.columns:
        exact_cmp = exact_cmp.rename(columns={"slope": f"slope_{anc}",
                                              "slope_se": f"slope_se_{anc}"})

    # ---- retention table (Choi-significant genes tested in both) ----
    choi_top["sig_Choi"] = choi_top["gene_base"].isin(choi_sig)
    gene_level = ours[["gene_base", "phenotype_id", f"sig_{anc}",
                       f"mlog10_pbeta_{anc}", "qval_gene"]].copy()
    ret = choi_top.merge(gene_level, on="gene_base", how="inner")
    ret_sig = ret[ret["sig_Choi"]].copy()
    ret_sig["mlog10_pbeta_Choi"] = ret_sig["mlogp"].clip(lower=0)
    ret_sig[f"mlog10_pbeta_{anc}"] = ret_sig[f"mlog10_pbeta_{anc}"]
    ret_sig["retention"] = np.where(
        ret_sig[f"sig_{anc}"],
        f"Retained in pooled {anc}", f"Lost in pooled {anc}")
    retention = ret_sig[["gene_base", "phenotype_id", "mlog10_pbeta_Choi",
                         f"mlog10_pbeta_{anc}", "retention"]]

    # ---- common-tested catalog ----
    catalog = ret[["gene_base", "phenotype_id", "sig_Choi", f"sig_{anc}"]].copy()
    catalog["category"] = np.select(
        [catalog["sig_Choi"] & catalog[f"sig_{anc}"],
         catalog["sig_Choi"] & ~catalog[f"sig_{anc}"],
         ~catalog["sig_Choi"] & catalog[f"sig_{anc}"]],
        ["Shared significant", "Choi-only significant",
         f"Pooled-{anc}-only significant"],
        default="Neither significant")
    catalog = catalog[["gene_base", "phenotype_id", "category"]]

    # ---- Choi-defined pair set (pair set B for 29d) ----
    pairs_b = ret_sig[["phenotype_id", "variant_norm"]].rename(
        columns={"variant_norm": "variant_id"}).drop_duplicates()

    # ---- write ----
    top_out = args.out_prefix + ".gene_top_comparison.tsv.gz"
    exact_out = args.out_prefix + f".exact_{anc}_lead_comparison.tsv.gz"
    ret_out = args.out_prefix + ".Choi_significant_retention.tsv.gz"
    cat_out = args.out_prefix + ".common_tested_catalog.tsv.gz"
    pairs_out = args.out_prefix + ".choi_defined_pairs.tsv"
    top_cmp.to_csv(top_out, sep="\t", index=False, compression="gzip")
    exact_cmp.to_csv(exact_out, sep="\t", index=False, compression="gzip")
    retention.to_csv(ret_out, sep="\t", index=False, compression="gzip")
    catalog.to_csv(cat_out, sep="\t", index=False, compression="gzip")
    pairs_b.to_csv(pairs_out, sep="\t", index=False)

    n_ret = int((ret_sig["retention"] == f"Retained in pooled {anc}").sum())
    print("\n========================================")
    print(f"CHOI vs POOLED {anc}")
    print("========================================")
    print(f"{anc} genes:", ours["gene_base"].nunique())
    print("Genes found in Choi full summary:", top_cmp["gene_base"].nunique())
    print("Same lead variant:", int(top_cmp["same_lead_variant"].sum()),
          f"/ {len(top_cmp)} ({top_cmp['same_lead_variant'].mean():.3f})")
    print(f"Exact {anc} lead SNP-gene pairs found in Choi:", len(exact_cmp))
    print(f"Choi-significant eGenes tested in both: {len(ret_sig)}; "
          f"retained: {n_ret} ({n_ret / max(len(ret_sig), 1):.3f})")
    print("Catalog:", catalog["category"].value_counts().to_dict())

    z = exact_cmp.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[f"t_{anc}", "t_Choi"])
    if len(z):
        print("Exact-pair direction agreement:",
              np.mean(np.sign(z[f"t_{anc}"]) == np.sign(z["t_Choi"])))
        print("Exact-pair Spearman t:",
              z[[f"t_{anc}", "t_Choi"]].corr(method="spearman").iloc[0, 1])

    print("\nWritten:")
    for f in (top_out, exact_out, ret_out, cat_out, pairs_out):
        print(" ", f)


if __name__ == "__main__":
    main()
