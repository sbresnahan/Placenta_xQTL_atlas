#!/usr/bin/env python3
"""
54_tier2_stratified.py — Objective 2.1 tier 2: ancestry-stratified GxE at
Aim-1 prioritized loci + inverse-variance-weighted meta-analysis.

Prioritized phenotypes per modality (union):
  * colocalized xQTL phenotypes — module-06 coloc results with
    PP.H4.abf >= --pp-h4 (default 0.7), and
  * phenotypes whose gene has a significant isoTWAS association —
    module-06 gene-level TWAS results with acat_q <= --twas-q (default
    0.05). Gene -> phenotype mapping uses the module-05
    {ANC}_{MOD}.phenotype_groups.txt (expression: phenotype_id IS the gene).

Per ancestry, the nominal interaction model (same tensorQTL-mirroring math
as the tier-1 scanner, gxe_core) is fit over the cis window of each
prioritized phenotype: per-ancestry BED, per-ancestry covariates
({ANC}_covariates_{MOD}.tsv), per-ancestry pgen, and the pooled exposure
vector restricted to that ancestry's samples.

Per variant x phenotype, per-ancestry b_gi/se are meta-analyzed by fixed-
effect IVW (Cochran's Q / I^2 for heterogeneity). Per phenotype, the top
variant by meta p-value is reported; BH q-values across the prioritized set
within the modality x exposure (tier 2 is a small hypothesis-driven set —
Storey pi0 is unstable there, so BH; documented in the README).

Execution modes:
  full   — scan all autosomes and finalize in one process (backward compatible)
  scan   — scan one chromosome and write a shard for LSF array execution
  merge  — combine requested chromosome shards and run the IVW/BH finalization

Outputs ({GXE_DIR}/tier2/):
  {MOD}_{EXP}.tier2.tsv             — one row per prioritized phenotype
  {MOD}_{EXP}.tier2_variants.tsv.gz — every tested variant x phenotype row
  shards/{MOD}_{EXP}.chrN.tier2_variants.tsv.gz — per-chromosome scan shards
"""

import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gxe_core


VARIANT_COLUMNS = [
    "phenotype_id", "variant_id", "tss_distance", "ancestry",
    "b_gi", "b_gi_se", "tstat_gi", "pval_nominal", "dof",
]


def log(msg):
    print(f"[{pd.Timestamp.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def canonicalize_phenotype_ids(ids, modality):
    """Canonicalize ancestry-local splicing/IR phenotype IDs.

    Module 03 assigns LeafCutter meta-cluster numbers and MAJIQ IR unified
    indices separately within each ancestry.  Remove only those local numeric
    labels so equivalent genomic events share one ID for Tier-2 selection and
    cross-ancestry IVW meta-analysis.
    """
    s = pd.Index(ids).astype(str).to_series(index=range(len(ids)))
    if modality == "splicing":
        out = s.str.replace(r"_clu_\d+_([+-])$", r"_\1", regex=True)
    elif modality == "intron_retention":
        out = s.str.replace(r"_([+-])_\d+$", r"_\1", regex=True)
    else:
        return pd.Index(ids)
    return pd.Index(out.values)


def load_scanner_module():
    """Import 52_gxe_scan.py by path (digit-leading filename)."""
    spec = importlib.util.spec_from_file_location(
        "gxe_scan_52", os.path.join(HERE, "52_gxe_scan.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def bh_qvalues(p):
    p = np.asarray(p, float)
    p = np.where(np.isfinite(p), p, 1.0)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p, kind="stable")
    q = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def resolve_existing(candidates, label):
    for c in candidates:
        if c and os.path.exists(c):
            return c
    log(f"  WARNING: {label} not found (tried: {candidates})")
    return None


def prioritized_phenotypes(modality, ancestries, qtl_dir, results_dir,
                           pp_h4, twas_q):
    """Union of coloc-called phenotype_ids and TWAS-significant genes'
    phenotypes for this modality."""
    phenos = set()
    # Module 06 (49_aggregate_coloc_twas.py) writes to coloc/aggregated/;
    # keep the older coloc/aggregate/ and coloc/ layouts as fallbacks.
    coloc_path = resolve_existing([
        os.path.join(results_dir, "coloc", "aggregated", "coloc_results.tsv.gz"),
        os.path.join(results_dir, "coloc", "aggregate", "coloc_results.tsv.gz"),
        os.path.join(results_dir, "coloc", "coloc_results.tsv.gz"),
    ], "coloc results")
    if coloc_path:
        coloc = pd.read_csv(coloc_path, sep="\t")
        m = (coloc["modality"] == modality)
        if "coloc_call" in coloc.columns:
            m &= coloc["coloc_call"].astype(bool)
        elif "PP.H4.abf" in coloc.columns:
            m &= coloc["PP.H4.abf"] >= pp_h4
        coloc_ids = canonicalize_phenotype_ids(
            coloc.loc[m, "phenotype_id"].astype(str).tolist(), modality)
        phenos |= set(coloc_ids)
        log(f"  coloc PP.H4>={pp_h4}: {len(phenos)} phenotypes")

    twas_path = resolve_existing([
        os.path.join(results_dir, "coloc", "aggregated", "twas_gene_results.tsv.gz"),
        os.path.join(results_dir, "coloc", "aggregate", "twas_gene_results.tsv.gz"),
        os.path.join(results_dir, "coloc", "twas_gene_results.tsv.gz"),
    ], "TWAS gene results")
    if twas_path:
        twas = pd.read_csv(twas_path, sep="\t")
        qcol = "acat_q" if "acat_q" in twas.columns else None
        # module 06 writes the gene column as "GENE" (FUSION .pos convention)
        gcol = next((c for c in ["gene", "gene_id", "id", "GENE"]
                     if c in twas.columns), None)
        if qcol and gcol:
            genes = set(twas.loc[twas[qcol] <= twas_q, gcol])
            log(f"  TWAS acat_q<={twas_q}: {len(genes)} genes")
            if modality == "expression":
                phenos |= genes
            else:
                # Splicing/IR IDs contain ancestry-local numbering upstream,
                # and phenotype availability can differ by ancestry.  Read all
                # available ancestry group files and canonicalize before union.
                for anc in ancestries:
                    groups = resolve_existing([
                        os.path.join(qtl_dir, f"{anc}_{modality}.phenotype_groups.txt")
                    ], f"phenotype groups ({anc})")
                    if not groups:
                        continue
                    gdf = pd.read_csv(groups, sep="\t", header=None,
                                      names=["phenotype_id", "gene_id"])
                    ids = gdf.loc[gdf["gene_id"].isin(genes),
                                  "phenotype_id"].astype(str).tolist()
                    phenos |= set(canonicalize_phenotype_ids(ids, modality))
    log(f"  prioritized phenotypes ({modality}): {len(phenos)}")
    return phenos


def scan_ancestry(scanner, anc, modality, exposure_id, pheno_keep, args,
                  chrom_filter=None):
    """Nominal interaction scan at prioritized phenotypes for one ancestry."""
    import torch
    import tensorqtl

    inputs_dir = os.path.join(args.gxe_dir, "inputs")
    bed = os.path.join(inputs_dir, f"{anc}_{modality}.bed.gz")
    cov_path = os.path.join(inputs_dir, f"{anc}_covariates_{modality}.tsv")
    pgen = os.path.join(args.qtl_dir, f"{anc}_qtl")
    for f in [bed, cov_path, pgen + ".pgen"]:
        if not os.path.exists(f):
            log(f"  WARNING: missing {f} — skipping {anc}")
            return pd.DataFrame(columns=VARIANT_COLUMNS)

    pheno, pheno_pos = tensorqtl.read_phenotype_bed(bed)

    # Canonicalize ancestry-local splicing/IR IDs before prioritized-phenotype
    # filtering and before emitting rows that will be meta-analyzed.
    canonical = canonicalize_phenotype_ids(pheno.index.tolist(), modality)
    if canonical.duplicated().any():
        dup = canonical[canonical.duplicated(keep=False)][:5].tolist()
        sys.exit(f"ERROR: {anc} {modality} canonicalization creates duplicate "
                 f"phenotype IDs: {dup}")
    pheno.index = canonical
    pheno_pos = pheno_pos.copy()
    pheno_pos.index = canonical

    pheno = pheno[pheno.index.isin(pheno_keep)]
    pheno_pos = pheno_pos.loc[pheno.index]

    # Module 07 is autosomal only. Tier 2 reads the original ancestry BEDs
    # directly (rather than the pooled Stage-0 BED), so enforce chr1-22 here.
    pos_chrom = pheno_pos["chr"].astype(str).str.replace(
        r"^chr", "", case=False, regex=True)
    chrom_num = pd.to_numeric(pos_chrom, errors="coerce")
    autosomal = chrom_num.between(1, 22) & (chrom_num % 1 == 0)
    n_nonauto = int((~autosomal).sum())
    if n_nonauto:
        log(f"  {anc}: dropping {n_nonauto} prioritized non-autosomal phenotypes")
    pheno = pheno.loc[autosomal.values]
    pheno_pos = pheno_pos.loc[autosomal.values]

    if chrom_filter is not None:
        pos_chrom = pheno_pos["chr"].astype(str).str.replace(
            r"^chr", "", case=False, regex=True)
        keep_chrom = (pos_chrom == str(chrom_filter)).values
        pheno = pheno.loc[keep_chrom]
        pheno_pos = pheno_pos.loc[keep_chrom]

    if len(pheno) == 0:
        suffix = f" on chr{chrom_filter}" if chrom_filter is not None else ""
        log(f"  {anc}: no prioritized phenotypes{suffix} in the BED")
        return pd.DataFrame(columns=VARIANT_COLUMNS)

    exp = pd.read_csv(os.path.join(args.gxe_dir, "inputs", "exposures.tsv"),
                      sep="\t", index_col=0).loc[exposure_id].dropna()
    cov = pd.read_csv(cov_path, sep="\t", index_col=0)
    samples = [s for s in pheno.columns if s in exp.index and s in cov.columns]
    if len(samples) < 20:
        log(f"  WARNING: {anc} has {len(samples)} samples — skipping")
        return pd.DataFrame(columns=VARIANT_COLUMNS)
    pheno = pheno[samples]
    cov = cov[samples]
    # The exposure main effect is fitted explicitly in the GxE model.
    # If the exposure is also present in the covariate matrix, remove it
    # to avoid including the same variable twice.
    if exposure_id in cov.index:
        cov = cov.drop(index=exposure_id)
    e = exp[samples]

    device = gxe_core.get_device()
    Q_t, dof_base = gxe_core.make_residualizer(cov.T.values, device)
    dof = dof_base - 2
    e_t = torch.tensor(e.values.astype(np.float32), device=device)

    rows = []
    pos_chrom = pheno_pos["chr"].astype(str).str.replace(
        r"^chr", "", case=False, regex=True)
    for chrom in sorted(pos_chrom.unique(), key=int):
        m = (pos_chrom == chrom).values
        G, var_ids, var_pos, gsamples = scanner.load_chromosome_genotypes(
            {anc: pgen}, {anc: samples}, chrom,
            maf_threshold=args.maf_threshold,
            maf_threshold_interaction=args.maf_threshold_interaction,
            interaction=e.values, dosages=args.dosages)
        G_t = torch.tensor(G, dtype=torch.float32, device=device)
        g0_t, e0_t, ge0_t = gxe_core.prepare_interaction_terms(G_t, e_t, Q_t)
        sub_p = pheno[m]
        sub_pos = pheno_pos[m]
        for pid, prow in sub_p.iterrows():
            p0_t = gxe_core.prepare_phenotype(
                torch.tensor(prow.values.astype(np.float32), device=device), Q_t)
            start = int(sub_pos.loc[pid, "start"]) if "start" in sub_pos.columns \
                else int(sub_pos.loc[pid, "pos"])
            end = int(sub_pos.loc[pid, "end"]) if "end" in sub_pos.columns else start
            lo, hi = scanner.cis_slice(var_pos, start, end, args.cis_window)
            if hi - lo < 1:
                continue
            sl = slice(lo, hi)
            Xt_w, Xinv_w, _ = gxe_core.build_window_design(
                g0_t[sl], e0_t, ge0_t[sl])
            b, b_se, tstat = gxe_core.window_fit(Xt_w, Xinv_w, p0_t, dof)
            t_gi = tstat[:, 2].cpu().numpy()
            pval = gxe_core.r2_to_pval(gxe_core.t_to_r2(t_gi, dof), dof)
            df = pd.DataFrame(dict(
                phenotype_id=pid, variant_id=var_ids[sl],
                tss_distance=(var_pos[sl] - start).astype(int),
                ancestry=anc,
                b_gi=b[:, 2].cpu().numpy(),
                b_gi_se=b_se[:, 2].cpu().numpy(),
                tstat_gi=t_gi, pval_nominal=pval, dof=int(dof)))
            rows.append(df[np.isfinite(df["pval_nominal"])])
        log(f"  {anc} chr{chrom}: done ({sum(len(r) for r in rows)} rows so far)")
    return (pd.concat(rows, ignore_index=True) if rows
            else pd.DataFrame(columns=VARIANT_COLUMNS))


def shard_path(args, chrom):
    shard_dir = args.shard_dir or os.path.join(args.gxe_dir, "tier2", "shards")
    return os.path.join(
        shard_dir,
        f"{args.modality}_{args.exposure}.chr{chrom}.tier2_variants.tsv.gz")


def empty_top_columns(ancestries):
    cols = [
        "modality", "exposure", "phenotype_id", "variant_id", "tss_distance",
        "n_ancestries", "beta_meta", "se_meta", "z_meta", "p_meta",
        "q_het", "p_het", "i2",
    ]
    for anc in ancestries:
        cols.extend([f"b_gi_{anc}", f"b_gi_se_{anc}", f"pval_{anc}"])
    cols.append("q_meta")
    return cols


def write_empty_final(args, out_pheno, out_var):
    pd.DataFrame(columns=VARIANT_COLUMNS).to_csv(
        out_var, sep="\t", index=False, compression="gzip")
    pd.DataFrame(columns=empty_top_columns(args.ancestries)).to_csv(
        out_pheno, sep="\t", index=False)
    log(f"  wrote empty {out_var}")
    log(f"  wrote empty {out_pheno}")


def finalize_results(allv, args, out_pheno, out_var):
    """Write combined variant rows and the per-phenotype IVW/BH summary."""
    if allv.empty:
        write_empty_final(args, out_pheno, out_var)
        return

    # IVW meta per phenotype x variant across ancestries.
    metas = []
    for (pid, vid), g in allv.groupby(["phenotype_id", "variant_id"], sort=False):
        r = gxe_core.ivw_meta(g["b_gi"].values, g["b_gi_se"].values)
        metas.append(dict(phenotype_id=pid, variant_id=vid,
                          tss_distance=int(g["tss_distance"].iloc[0]),
                          n_ancestries=r["k"], beta_meta=r["beta"],
                          se_meta=r["se"], z_meta=r["z"], p_meta=r["p"],
                          q_het=r["q_het"], p_het=r["p_het"], i2=r["i2"]))
    meta_df = pd.DataFrame(metas)
    allv.to_csv(out_var, sep="\t", index=False, compression="gzip")
    log(f"  wrote {out_var} ({len(allv)} rows)")

    # Top variant per phenotype by meta p + per-ancestry stats at that variant.
    top = (meta_df.sort_values("p_meta")
           .groupby("phenotype_id", sort=False).head(1).copy())
    for anc in args.ancestries:
        sub = allv[allv["ancestry"] == anc][
            ["phenotype_id", "variant_id", "b_gi", "b_gi_se", "pval_nominal"]]
        sub = sub.rename(columns={"b_gi": f"b_gi_{anc}",
                                  "b_gi_se": f"b_gi_se_{anc}",
                                  "pval_nominal": f"pval_{anc}"})
        top = top.merge(sub, on=["phenotype_id", "variant_id"], how="left")
    top["q_meta"] = bh_qvalues(top["p_meta"].values)
    top.insert(0, "exposure", args.exposure)
    top.insert(0, "modality", args.modality)
    top = top.sort_values("p_meta")
    top.to_csv(out_pheno, sep="\t", index=False)
    log(f"  wrote {out_pheno} ({len(top)} phenotypes; "
        f"{(top['q_meta'] <= 0.05).sum()} at BH q<=0.05)")


def scan_requested(args, chrom_filter=None):
    scanner = load_scanner_module()
    pheno_keep = prioritized_phenotypes(args.modality, args.ancestries,
                                        args.qtl_dir, args.results_dir,
                                        args.pp_h4, args.twas_q)
    if not pheno_keep:
        log("  no prioritized phenotypes")
        return pd.DataFrame(columns=VARIANT_COLUMNS)

    per_anc = []
    for anc in args.ancestries:
        df = scan_ancestry(scanner, anc, args.modality, args.exposure,
                           pheno_keep, args, chrom_filter=chrom_filter)
        if len(df):
            per_anc.append(df)
    return (pd.concat(per_anc, ignore_index=True) if per_anc
            else pd.DataFrame(columns=VARIANT_COLUMNS))


def main():
    p = argparse.ArgumentParser(description="Tier-2 ancestry-stratified GxE + IVW meta")
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--gxe-dir", required=True)
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    p.add_argument("--modality", required=True)
    p.add_argument("--exposure", required=True)
    p.add_argument("--mode", choices=["full", "scan", "merge"], default="full",
                   help="full=legacy all-chromosome run; scan=one chromosome shard; "
                        "merge=combine shards and finalize")
    p.add_argument("--chrom", type=int,
                   help="chromosome for --mode scan")
    p.add_argument("--chroms", nargs="+", type=int,
                   help="expected chromosome shards for --mode merge")
    p.add_argument("--shard-dir",
                   help="override tier-2 shard directory")
    p.add_argument("--pp-h4", type=float, default=0.7)
    p.add_argument("--twas-q", type=float, default=0.05)
    p.add_argument("--cis-window", type=int, default=1000000)
    p.add_argument("--maf-threshold", type=float, default=0.01,
                   help="overall in-sample MAF prefilter (default: 0.01)")
    p.add_argument("--maf-threshold-interaction", type=float, default=0.05,
                   help="tensorQTL-style MAF floor required in both lower and "
                        "upper exposure halves (default: 0.05)")
    p.add_argument("--dosages", action="store_true")
    args = p.parse_args()

    for name, value in [("maf_threshold", args.maf_threshold),
                        ("maf_threshold_interaction", args.maf_threshold_interaction)]:
        if not (0.0 <= value <= 0.5):
            p.error(f"--{name.replace(chr(95), chr(45))} must be between 0 and 0.5")

    if args.mode == "scan":
        if args.chrom is None:
            p.error("--mode scan requires --chrom")
        if not 1 <= args.chrom <= 22:
            p.error("--chrom must be between 1 and 22")
    if args.mode == "merge":
        if not args.chroms:
            p.error("--mode merge requires --chroms")
        bad = [c for c in args.chroms if not 1 <= c <= 22]
        if bad:
            p.error(f"--chroms contains non-autosomal values: {bad}")

    out_dir = os.path.join(args.gxe_dir, "tier2")
    os.makedirs(out_dir, exist_ok=True)
    out_pheno = os.path.join(out_dir, f"{args.modality}_{args.exposure}.tier2.tsv")
    out_var = os.path.join(out_dir, f"{args.modality}_{args.exposure}.tier2_variants.tsv.gz")

    if args.mode == "scan":
        out_shard = shard_path(args, args.chrom)
        os.makedirs(os.path.dirname(out_shard), exist_ok=True)
        if os.path.exists(out_shard) and os.environ.get("FORCE", "0") != "1":
            log(f"{out_shard} exists — skipping (FORCE=1 to redo)")
            return
        allv = scan_requested(args, chrom_filter=args.chrom)
        allv.to_csv(out_shard, sep="\t", index=False, compression="gzip")
        log(f"  wrote {out_shard} ({len(allv)} rows)")
        return

    if os.path.exists(out_pheno) and os.environ.get("FORCE", "0") != "1":
        log(f"{out_pheno} exists — skipping (FORCE=1 to redo)")
        return

    if args.mode == "merge":
        paths = [shard_path(args, chrom) for chrom in args.chroms]
        missing = [path for path in paths if not os.path.exists(path)]
        if missing:
            preview = "\n    ".join(missing[:10])
            extra = "" if len(missing) <= 10 else f"\n    ... +{len(missing) - 10} more"
            sys.exit(f"ERROR: missing {len(missing)} tier-2 chromosome shard(s):\n"
                     f"    {preview}{extra}")
        frames = []
        for path in paths:
            df = pd.read_csv(path, sep="\t")
            if len(df):
                frames.append(df)
        allv = (pd.concat(frames, ignore_index=True) if frames
                else pd.DataFrame(columns=VARIANT_COLUMNS))
        log(f"  merging {len(paths)} chromosome shards ({len(allv)} rows)")
        finalize_results(allv, args, out_pheno, out_var)
        return

    # Backward-compatible single-process execution over all autosomes.
    allv = scan_requested(args)
    finalize_results(allv, args, out_pheno, out_var)


if __name__ == "__main__":
    main()
