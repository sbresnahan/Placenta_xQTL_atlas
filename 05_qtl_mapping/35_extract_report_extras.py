#!/usr/bin/env python3
"""35_extract_report_extras.py — targeted extractions for report figure panels.

The mapping parquets are cis.map_cis lead-only outputs (one row per
phenotype/group), so report panels that need per-variant statistics in
*other* contexts (the other ancestry, other modalities, or a whole locus)
require targeted re-computation. This script performs those extractions
with a lightweight OLS engine that mirrors the mapping model exactly:
phenotype ~ genotype + the same per-(ancestry, modality) covariates used
by 27_run_tensorqtl.py, on the same intersected {ANC}_qtl pgens.

Outputs land in {results_dir}/report_extras/ and are staged into the
report archive by make_report_archive.sh:

  cross_ancestry_lookup.tsv.gz   — each ancestry's FDR<=q lead (variant,
                                   phenotype) pairs re-tested in the other
                                   ancestry (all layers)
  cross_modality_lookup.tsv.gz   — each grouped-layer lead (variant, gene)
                                   re-tested against all phenotypes of the
                                   same gene in every other grouped modality
  locus_{GENE}_{ANC}.tsv.gz      — per-variant regional stats + r2 to the
                                   lead variant for showcase loci
  lead_variant_annotations.tsv.gz— fastVEP consequence class + cCRE/OCR
                                   overlap flags for every unique variant in
                                   any top table (foreground = q<=q lead)
  gene_constraint.tsv            — gnomAD v2.1.1 pLI/LOEUF for tested genes
  pcair_pcs.tsv                  — per-cohort PC-AiR PCs (staged, if found)

Each section is independent and skippable (--only / --skip). Sections that
need external resources (fastVEP binary, cCRE BEDs, gnomAD file, PC-AiR
outputs) degrade to instructions + warnings when those are absent.

Requires: pandas, numpy, scipy, tensorqtl (genotypeio) — the mapping env.
"""

import argparse
import glob
import gzip
import os
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sstats


def _strip_ver(ids):
    """Strip Ensembl version suffixes ('.13') from an id Series."""
    return ids.astype(str).str.replace(r"\.\d+$", "", regex=True)


def _unver(gid):
    """Strip the Ensembl version suffix from a single id."""
    return re.sub(r"\.\d+$", "", str(gid))


def _bed_pheno_id(pheno_id, mod, group_id=None):
    """Map a combined-layer compound phenotype id to the id used in the
    per-modality BED. Combined tables encode phenotype_id as
    '{modality}__{phenotype_id}' (ungrouped) or
    '{modality}__{group_id}__{phenotype_id}' (grouped); per-modality top
    tables already carry the bare BED id and pass through unchanged."""
    pid = str(pheno_id)
    if pid.startswith(f"{mod}__"):
        pid = pid[len(mod) + 2:]
        if group_id is not None and pid.startswith(f"{group_id}__"):
            pid = pid[len(str(group_id)) + 2:]
    return pid

MODALITIES = ["expression", "isoforms", "isoform_expression", "splicing",
              "intron_retention", "alt_TSS", "alt_polyA", "RNA_editing",
              "stability"]
# layers: grouped (all 9), ungrouped (8 non-expression), combined (1)
UNGROUPED_MODALITIES = [m for m in MODALITIES if m != "expression"]

GNOMAD_CONSTRAINT_URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
                         "release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")


# ---------------------------------------------------------------------------
# Core statistics (fixture-testable without tensorqtl / pgens)
# ---------------------------------------------------------------------------

def residualize(M, C):
    """Residualize features x samples matrix M on samples x k covariates C.

    An intercept column is added internally. Returns M - M_hat.
    """
    C1 = np.column_stack([np.ones(C.shape[0]), C])
    Q, _ = np.linalg.qr(C1)
    return M - M @ Q @ Q.T


def ols_stats(G, Y, n_cov):
    """Vectorized OLS for residualized matrices.

    G: variants x samples (residualized), Y: phenotypes x samples
    (residualized), n_cov: number of covariates INCLUDING intercept.
    Returns per-(variant, phenotype) slope, se, pval as broadcast against
    the diagonal pairing used by callers — here we compute the full
    cross-product only when requested; callers normally pass one phenotype
    at a time (Y shape 1 x n) or use pair_ols below.
    """
    GY = G @ Y.T                      # variants x phenotypes
    G2 = (G * G).sum(axis=1)          # variants
    Y2 = (Y * Y).sum(axis=1)          # phenotypes
    slope = GY / G2[:, None]
    n = G.shape[1]
    df = n - n_cov - 1
    # SSR = ||y||^2 - slope^2 ||g||^2
    ssr = Y2[None, :] - slope ** 2 * G2[:, None]
    ssr = np.clip(ssr, 1e-300, None)
    sigma2 = ssr / df
    se = np.sqrt(sigma2 / G2[:, None])
    with np.errstate(divide="ignore", invalid="ignore"):
        t = slope / se
    pval = 2 * sstats.t.sf(np.abs(t), df)
    return slope, se, pval


def pair_lookup(geno_res, pheno_res, pairs, var_index, phen_index, n_cov):
    """OLS for specific (variant_id, phenotype_id) pairs.

    geno_res / pheno_res: residualized matrices (features x samples).
    pairs: DataFrame with variant_id, phenotype_id columns.
    var_index / phen_index: {id: row_position}.
    Returns arrays slope, se, pval aligned to pairs (NaN when an id is
    absent from the target matrices).
    """
    n = len(pairs)
    slope = np.full(n, np.nan)
    se = np.full(n, np.nan)
    pval = np.full(n, np.nan)
    # group pairs by phenotype row for vectorization
    p_rows = pairs["phenotype_id"].map(phen_index)
    v_rows = pairs["variant_id"].map(var_index)
    ok = p_rows.notna() & v_rows.notna()
    if not ok.any():
        return slope, se, pval
    pairs_ok = np.where(ok.to_numpy())[0]
    pr = p_rows[ok].astype(int).to_numpy()
    vr = v_rows[ok].astype(int).to_numpy()
    df = geno_res.shape[1] - n_cov - 1
    # iterate over unique phenotypes (vectorized across their variants)
    for p in pd.unique(pr):
        mask = pr == p
        v_idx = vr[mask]
        G = geno_res[v_idx, :]
        y = pheno_res[p, :]
        G2 = (G * G).sum(axis=1)
        G2[G2 == 0] = np.nan
        b = (G @ y) / G2
        ssr = np.clip((y * y).sum() - b ** 2 * G2, 1e-300, None)
        s = np.sqrt(ssr / df / G2)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = b / s
        pv = 2 * sstats.t.sf(np.abs(t), df)
        tgt = pairs_ok[mask]
        slope[tgt] = b
        se[tgt] = s
        pval[tgt] = pv
    return slope, se, pval


# ---------------------------------------------------------------------------
# IO adapters (seadragon: tensorqtl env)
# ---------------------------------------------------------------------------

def read_bed(bed_path):
    """Read a tensorQTL phenotype BED (.bed.gz) -> (DataFrame phenotypes
    (phenotype_id x samples), DataFrame pos (phenotype_id, chr, start, end)).
    """
    bed = pd.read_csv(bed_path, sep="\t", dtype={"#chr": str})
    bed = bed.rename(columns={"#chr": "chr"})
    pos = bed[["chr", "start", "end", "phenotype_id"]].copy()
    pos["chr"] = pos["chr"].astype(str).str.replace("^chr", "", regex=True)
    pheno = bed.drop(columns=["chr", "start", "end"]).set_index("phenotype_id")
    pheno = pheno[~pheno.index.duplicated(keep="first")]
    return pheno, pos


def load_ancestry_genotypes(qtl_dir, anc):
    """Load the intersected {ANC}_qtl pgen -> (genotype_df, variant_df).

    Uses tensorqtl's genotypeio (same as 27_run_tensorqtl.py). Chrom names
    are normalized to bare ('1'-style) to match top-table variant_ids.
    """
    from tensorqtl import genotypeio
    prefix = os.path.join(qtl_dir, f"{anc}_qtl")
    genotype_df, variant_df = genotypeio.load_genotypes(prefix)
    variant_df = variant_df.copy()
    variant_df["chrom"] = variant_df["chrom"].astype(str).str.replace(
        "^chr", "", regex=True)
    return genotype_df, variant_df


def load_covariates(qtl_dir, anc, mod):
    """Per-modality covariates -> samples x covariates (with fallback)."""
    per_mod = os.path.join(qtl_dir, f"{anc}_covariates_{mod}.tsv")
    shared = os.path.join(qtl_dir, f"{anc}_covariates.tsv")
    path = per_mod if os.path.exists(per_mod) else shared
    if not os.path.exists(path):
        sys.exit(f"ERROR: no covariates for {anc}/{mod} ({per_mod} or {shared})")
    cov = pd.read_csv(path, sep="\t", index_col=0).T
    return cov


def load_groups(qtl_dir, anc, mod):
    """phenotype_id -> group_id Series for a grouped modality."""
    path = os.path.join(qtl_dir, f"{anc}_{mod}.phenotype_groups.txt")
    g = pd.read_csv(path, sep="\t", header=None,
                    names=["phenotype_id", "group_id"])
    g = g.drop_duplicates("phenotype_id")
    return pd.Series(g["group_id"].values, index=g["phenotype_id"])


def load_leads(results_dir, anc, q_max):
    """All lead tables for one ancestry -> long DataFrame with a 'layer'
    column in {grouped, ungrouped, combined} and modality label."""
    frames = []
    for mod in MODALITIES:
        p = os.path.join(results_dir, f"{anc}_{mod}_cisqtl_top.tsv")
        if os.path.exists(p):
            df = pd.read_csv(p, sep="\t")
            df["layer"] = "grouped"
            df["modality"] = mod
            frames.append(df)
    for mod in UNGROUPED_MODALITIES:
        p = os.path.join(results_dir, f"{anc}_{mod}_ungrouped_cisqtl_top.tsv")
        if os.path.exists(p):
            df = pd.read_csv(p, sep="\t")
            df["layer"] = "ungrouped"
            df["modality"] = mod
            frames.append(df)
    p = os.path.join(results_dir, f"{anc}_combined_cisqtl_top.tsv")
    if os.path.exists(p):
        df = pd.read_csv(p, sep="\t")
        df["layer"] = "combined"
        df["modality"] = "combined"
        frames.append(df)
    if not frames:
        sys.exit(f"ERROR: no top tables found for {anc} in {results_dir}")
    out = pd.concat(frames, ignore_index=True)
    out["ancestry"] = anc
    out["is_lead"] = out["qval"] <= q_max
    return out


class TargetScanner:
    """Per-(ancestry, modality) lookup engine: residualizes the needed
    phenotype rows and variant rows once, then serves pair lookups."""

    def __init__(self, qtl_dir, anc, mod, genotype_df):
        self.anc, self.mod = anc, mod
        self.pheno, self.pos = read_bed(
            os.path.join(qtl_dir, f"{anc}_{mod}.bed.gz"))
        cov = load_covariates(qtl_dir, anc, mod)
        common = [s for s in self.pheno.columns
                  if s in genotype_df.columns and s in cov.index]
        if len(common) < 10:
            sys.exit(f"ERROR: {anc}/{mod}: only {len(common)} common samples")
        self.samples = common
        self.C = cov.loc[common].to_numpy(dtype=np.float64)
        self.n_cov = self.C.shape[1] + 1  # + intercept
        self.genotype_df = genotype_df[common]
        self._pheno_res = None
        self._var_cache = {}

    def phenotypes_residualized(self, pheno_ids):
        M = self.pheno.loc[pheno_ids, self.samples].to_numpy(dtype=np.float64)
        return residualize(M, self.C)

    def variants_residualized(self, variant_ids):
        missing = [v for v in variant_ids if v not in self._var_cache]
        if missing:
            present = [v for v in missing if v in self.genotype_df.index]
            if present:
                G = self.genotype_df.loc[present].to_numpy(dtype=np.float64)
                G = residualize(G, self.C)
                for i, v in enumerate(present):
                    self._var_cache[v] = G[i, :]
            for v in missing:
                if v not in self._var_cache:
                    self._var_cache[v] = None
        return self._var_cache

    def lookup(self, pairs):
        """pairs: DataFrame(variant_id, phenotype_id) -> DataFrame with
        slope/se/pval columns (NaN where variant or phenotype absent)."""
        pheno_ids = sorted(pairs["phenotype_id"].unique())
        pheno_ids = [p for p in pheno_ids if p in self.pheno.index]
        P = self.phenotypes_residualized(pheno_ids)
        phen_index = {p: i for i, p in enumerate(pheno_ids)}
        var_ids = sorted(pairs["variant_id"].unique())
        cache = self.variants_residualized(var_ids)
        present = [v for v in var_ids if cache.get(v) is not None]
        G = np.array([cache[v] for v in present]) if present else \
            np.zeros((0, len(self.samples)))
        var_index = {v: i for i, v in enumerate(present)}
        s, e, pv = pair_lookup(G, P, pairs, var_index, phen_index, self.n_cov)
        out = pairs.copy()
        out["repl_slope"], out["repl_se"], out["repl_pval"] = s, e, pv
        return out


# ---------------------------------------------------------------------------
# Section 1: cross-ancestry lookup
# ---------------------------------------------------------------------------

def cross_ancestry_lookup(results_dir, qtl_dir, ancestries, q_max, out_dir):
    leads = {a: load_leads(results_dir, a, q_max) for a in ancestries}
    genos = {}
    frames = []
    for disc_anc in ancestries:
        repl_ancs = [a for a in ancestries if a != disc_anc]
        if not repl_ancs:
            continue
        repl_anc = repl_ancs[0]
        hits = leads[disc_anc][leads[disc_anc]["is_lead"]]
        for mod, grp in hits.groupby("modality"):
            key = (repl_anc, mod)
            if key not in genos:
                bed = os.path.join(qtl_dir, f"{repl_anc}_{mod}.bed.gz")
                if not os.path.exists(bed):
                    print(f"  WARN: no BED for {repl_anc}/{mod} — skipped")
                    genos[key] = None
                    continue
                if repl_anc not in genos:
                    genos[repl_anc] = load_ancestry_genotypes(qtl_dir, repl_anc)[0]
                genos[key] = TargetScanner(qtl_dir, repl_anc, mod, genos[repl_anc])
            sc = genos[key]
            if sc is None:
                continue
            pairs = grp[["variant_id", "phenotype_id"]].drop_duplicates()
            res = sc.lookup(pairs)
            res["discovery_ancestry"] = disc_anc
            res["replication_ancestry"] = repl_anc
            res["modality"] = mod
            res = res.merge(grp[["variant_id", "phenotype_id", "layer",
                                 "slope", "qval", "group_id"]].drop_duplicates(),
                            on=["variant_id", "phenotype_id"], how="left")
            res = res.rename(columns={"slope": "discovery_slope",
                                      "qval": "discovery_qval"})
            frames.append(res)
            print(f"  {disc_anc}->{repl_anc} {mod}: {len(res)} pairs "
                  f"({res['repl_pval'].notna().sum()} tested)")
    if not frames:
        print("  cross-ancestry lookup: nothing to write")
        return None
    out = pd.concat(frames, ignore_index=True)
    path = os.path.join(out_dir, "cross_ancestry_lookup.tsv.gz")
    out.to_csv(path, sep="\t", index=False, compression="gzip")
    print(f"  wrote {len(out)} rows -> {path}")
    return path


# ---------------------------------------------------------------------------
# Section 2: cross-modality lookup (grouped layers, within ancestry)
# ---------------------------------------------------------------------------

def cross_modality_lookup(results_dir, qtl_dir, ancestries, q_max, out_dir):
    frames = []
    for anc in ancestries:
        leads = load_leads(results_dir, anc, q_max)
        hits = leads[(leads["is_lead"]) & (leads["layer"] == "grouped")]
        if hits.empty:
            continue
        geno = load_ancestry_genotypes(qtl_dir, anc)[0]
        scanners = {}
        for disc_mod, grp in hits.groupby("modality"):
            pairs_by_gene = grp[["variant_id", "group_id"]].drop_duplicates()
            for repl_mod in MODALITIES:
                if repl_mod == disc_mod:
                    continue
                try:
                    groups = load_groups(qtl_dir, anc, repl_mod)
                except FileNotFoundError:
                    continue
                key = repl_mod
                if key not in scanners:
                    bed = os.path.join(qtl_dir, f"{anc}_{repl_mod}.bed.gz")
                    if not os.path.exists(bed):
                        scanners[key] = None
                        continue
                    scanners[key] = TargetScanner(qtl_dir, anc, repl_mod, geno)
                sc = scanners[key]
                if sc is None:
                    continue
                # expand: each (variant, gene) -> all phenotypes of gene in repl_mod
                gene2pheno = groups.groupby(groups.values)
                rows = []
                for _, r in pairs_by_gene.iterrows():
                    phenos = gene2pheno.groups.get(r["group_id"])
                    if phenos is None:
                        continue
                    for ph in phenos:
                        rows.append((r["variant_id"], r["group_id"], ph))
                if not rows:
                    continue
                pairs = pd.DataFrame(rows, columns=["variant_id", "group_id",
                                                    "phenotype_id"])
                res = sc.lookup(pairs[["variant_id", "phenotype_id"]])
                res["group_id"] = pairs["group_id"].to_numpy()
                res["ancestry"] = anc
                res["discovery_modality"] = disc_mod
                res["replication_modality"] = repl_mod
                frames.append(res)
                print(f"  {anc} {disc_mod}->{repl_mod}: {len(res)} lookups")
    if not frames:
        print("  cross-modality lookup: nothing to write")
        return None
    out = pd.concat(frames, ignore_index=True)
    path = os.path.join(out_dir, "cross_modality_lookup.tsv.gz")
    out.to_csv(path, sep="\t", index=False, compression="gzip")
    print(f"  wrote {len(out)} rows -> {path}")
    return path


# ---------------------------------------------------------------------------
# Section 3: showcase-locus regional scans + LD to lead
# ---------------------------------------------------------------------------

def pick_showcase_loci(results_dir, ancestries, q_max, extra_gene=None):
    """ERAP2 + top shared non-expression-driver combined gene."""
    loci = {"ERAP2": "ENSG00000164308"}
    comb = []
    for anc in ancestries:
        p = os.path.join(results_dir, f"{anc}_combined_cisqtl_top.tsv")
        if os.path.exists(p):
            df = pd.read_csv(p, sep="\t")
            df["ancestry"] = anc
            comb.append(df)
    if len(comb) == 2:
        c0, c1 = comb
        hits0 = c0[c0["qval"] <= q_max]
        hits1 = c1[c1["qval"] <= q_max]
        # group_ids may carry Ensembl version suffixes; match unversioned
        shared = set(_strip_ver(hits0["group_id"])) & \
            set(_strip_ver(hits1["group_id"]))
        cand = pd.concat([hits0, hits1])
        cand = cand[_strip_ver(cand["group_id"]).isin(shared)]
        cand = cand[~_strip_ver(cand["group_id"]).isin(
            ["ENSG00000164308", "ENSG00000164307"])]  # not ERAP2/ERAP1
        if "modality" in cand.columns:
            cand = cand[cand["modality"] != "expression"]
        if not cand.empty:
            best = cand.nsmallest(1, "qval").iloc[0]
            gid = _unver(best["group_id"])
            loci[str(best.get("symbol", gid))] = gid
    if extra_gene:
        loci["extra"] = _unver(extra_gene)
    return loci


def locus_scans(results_dir, qtl_dir, ancestries, q_max, out_dir,
                window=1_000_000, extra_gene=None):
    loci = pick_showcase_loci(results_dir, ancestries, q_max, extra_gene)
    print(f"  showcase loci: {loci}")
    paths = []
    for label, gene_id in loci.items():
        for anc in ancestries:
            # driver modality + lead phenotype from the combined top table
            comb_p = os.path.join(results_dir, f"{anc}_combined_cisqtl_top.tsv")
            if not os.path.exists(comb_p):
                continue
            comb = pd.read_csv(comb_p, sep="\t")
            if "group_id" in comb.columns:
                row = comb[_strip_ver(comb["group_id"]) == gene_id]
            else:
                row = comb.iloc[0:0]
            if row.empty:
                # fall back: any modality hit for this gene (grouped tables
                # key on group_id; ungrouped tables key on phenotype_id)
                for mod in MODALITIES:
                    p = os.path.join(results_dir, f"{anc}_{mod}_cisqtl_top.tsv")
                    if os.path.exists(p):
                        t = pd.read_csv(p, sep="\t")
                        idcol = ("group_id" if "group_id" in t.columns
                                 else "phenotype_id")
                        r = t[_strip_ver(t[idcol]) == gene_id]
                        if not r.empty:
                            row = r.assign(modality=mod)
                            break
            if row.empty:
                print(f"  WARN: {gene_id} not found in {anc} tables — skipped")
                continue
            row = row.sort_values("qval").iloc[0]
            mod = row["modality"] if "modality" in row.index else "expression"
            pheno_id = _bed_pheno_id(row["phenotype_id"], mod,
                                     row.get("group_id"))
            lead_var = row["variant_id"]
            geno, variant_df = load_ancestry_genotypes(qtl_dir, anc)
            sc = TargetScanner(qtl_dir, anc, mod, geno)
            if pheno_id not in sc.pheno.index:
                print(f"  WARN: {pheno_id} not in {anc}/{mod} BED — skipped")
                continue
            chrom = sc.pos.loc[sc.pos["phenotype_id"] == pheno_id, "chr"].iloc[0]
            tss = sc.pos.loc[sc.pos["phenotype_id"] == pheno_id, "start"].iloc[0]
            lo, hi = max(0, tss - window), tss + window
            vsub = variant_df[(variant_df["chrom"] == str(chrom)) &
                              (variant_df["pos"] >= lo) &
                              (variant_df["pos"] <= hi)]
            if vsub.empty:
                print(f"  WARN: no variants in {label} region for {anc}")
                continue
            pairs = pd.DataFrame({"variant_id": vsub.index,
                                  "phenotype_id": pheno_id})
            res = sc.lookup(pairs)
            res["chrom"] = chrom
            res["pos"] = vsub["pos"].to_numpy()
            # LD to lead variant (squared Pearson on residualized dosages)
            cache = sc.variants_residualized([lead_var])
            g_lead = cache.get(lead_var)
            if g_lead is not None:
                cache = sc.variants_residualized(list(vsub.index))
                present = [v for v in vsub.index if cache.get(v) is not None]
                G = np.array([cache[v] for v in present])
                num = G @ g_lead
                den = np.sqrt((G * G).sum(axis=1) * (g_lead * g_lead).sum())
                r2 = pd.Series((num / den) ** 2, index=present)
                res["r2_lead"] = res["variant_id"].map(r2)
            else:
                res["r2_lead"] = np.nan
            res["ancestry"] = anc
            res["modality"] = mod
            res["gene_id"] = gene_id
            res["locus"] = label
            res["lead_variant"] = lead_var
            path = os.path.join(out_dir, f"locus_{label}_{anc}.tsv.gz")
            res.to_csv(path, sep="\t", index=False, compression="gzip")
            print(f"  wrote {len(res)} variants -> {path}")
            paths.append(path)
    return paths


# ---------------------------------------------------------------------------
# Section 4: lead-variant annotations (fastVEP + cCRE/OCR)
# ---------------------------------------------------------------------------

def annotate_lead_variants(results_dir, ancestries, q_max, out_dir,
                           fastvep_output=None, ccre_beds=None,
                           placenta_ocr=None):
    """Build the per-variant annotation table for every unique variant in
    any top table. fastVEP is run externally (or pass --fastvep-output);
    cCRE/OCR intersections are computed here with module-34 helpers."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "mod34", str(Path(__file__).with_name("34_pip_annotation_enrichment.py")))
    mod34 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod34)

    leads = pd.concat([load_leads(results_dir, a, q_max) for a in ancestries],
                      ignore_index=True)
    variants = pd.DataFrame({"variant_id": sorted(leads["variant_id"].unique())})
    parts = variants["variant_id"].str.split(":", expand=True)
    variants["chrom"] = parts[0].str.replace("^chr", "", regex=True)
    variants["pos"] = parts[1].astype(int)
    # foreground membership per ancestry (any layer)
    fg = (leads[leads["is_lead"]]
          .groupby("variant_id")["ancestry"]
          .agg(lambda s: ",".join(sorted(set(s)))))
    variants["foreground_in"] = variants["variant_id"].map(fg).fillna("")

    # fastVEP: write input VCF; parse output if provided
    vcf_path = os.path.join(out_dir, "lead_variants_fastvep_input.vcf")
    snps = variants["variant_id"]
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in variants["chrom"].unique():
            f.write(f"##contig=<ID={c}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v, c, p, r, a in zip(snps, variants["chrom"], variants["pos"],
                                 parts[2], parts[3]):
            f.write(f"{c}\t{p}\t{v}\t{r}\t{a}\t.\t.\t.\n")
    print(f"  wrote {len(variants)} variants -> {vcf_path}")
    if fastvep_output and os.path.exists(fastvep_output):
        cons = mod34.parse_fastvep(fastvep_output)
        variants["consequence"] = variants["variant_id"].map(cons).fillna("other")
    else:
        variants["consequence"] = np.nan
        print("  NOTE: run fastVEP on the VCF and re-run with "
              "--fastvep-output to fill consequence classes:")
        print(f"    fastvep annotate -i {vcf_path} -o lead_variants_fastvep.txt "
              "--output-format tab --gff3 <GRCh38.115.gff3> --fasta <GRCh38.fa>")

    # cCRE / OCR flags
    for bed_path, label in (ccre_beds or []):
        intervals = mod34.load_bed_intervals(bed_path)
        variants[f"ccre_{label}"] = mod34.overlap_flags(
            variants["chrom"].to_numpy(), variants["pos"].to_numpy(), intervals)
    if placenta_ocr and os.path.exists(placenta_ocr):
        intervals = mod34.load_bed_intervals(placenta_ocr)
        variants["placenta_ocr"] = mod34.overlap_flags(
            variants["chrom"].to_numpy(), variants["pos"].to_numpy(), intervals)

    path = os.path.join(out_dir, "lead_variant_annotations.tsv.gz")
    variants.to_csv(path, sep="\t", index=False, compression="gzip")
    print(f"  wrote {len(variants)} variants -> {path}")
    return path


# ---------------------------------------------------------------------------
# Section 5: gnomAD constraint subset
# ---------------------------------------------------------------------------

def constraint_table(results_dir, qtl_dir, ancestries, out_dir,
                     gnomad_path=None):
    """gnomAD v2.1.1 per-gene pLI/LOEUF for all tested genes (+ background)."""
    if gnomad_path is None:
        gnomad_path = os.path.join(out_dir, os.path.basename(GNOMAD_CONSTRAINT_URL))
        if not os.path.exists(gnomad_path):
            print(f"  downloading {GNOMAD_CONSTRAINT_URL}")
            try:
                urllib.request.urlretrieve(GNOMAD_CONSTRAINT_URL, gnomad_path)
            except Exception as e:
                print(f"  WARN: constraint download failed ({e}); section skipped")
                return None
    keep = ["gene_id", "gene", "pLI", "oe_lof", "oe_lof_upper", "lof_z"]
    con = pd.read_csv(gnomad_path, sep="\t", compression="gzip",
                      usecols=lambda c: c in keep, low_memory=False)
    con["gene_id_base"] = con["gene_id"].str.split(".").str[0]
    # tested genes: union of gene ids across top tables (grouped tables key
    # on group_id; ungrouped gene-level tables on phenotype_id). Ungrouped
    # tables lack group_id, so select columns by predicate; keep only
    # Ensembl gene ids (skips junction/peak ids) and strip versions to
    # match gnomAD's gene_id_base.
    tested = set()
    for anc in ancestries:
        for mod in MODALITIES:
            p = os.path.join(results_dir, f"{anc}_{mod}_cisqtl_top.tsv")
            if os.path.exists(p):
                t = pd.read_csv(p, sep="\t",
                                usecols=lambda c: c in ("group_id",
                                                        "phenotype_id"))
                col = "group_id" if "group_id" in t.columns else "phenotype_id"
                ids = t[col].dropna().astype(str)
                ids = ids[ids.str.startswith("ENSG")]
                tested |= set(ids.str.split(".").str[0])
    con["tested"] = con["gene_id_base"].isin(tested)
    out = con.drop(columns=["gene_id"]).rename(columns={"gene_id_base": "gene_id"})
    path = os.path.join(out_dir, "gene_constraint.tsv")
    out.to_csv(path, sep="\t", index=False)
    print(f"  wrote {len(out)} genes ({out['tested'].sum()} tested) -> {path}")
    return path


# ---------------------------------------------------------------------------
# Section 6: PC-AiR PC staging
# ---------------------------------------------------------------------------

def stage_pcair(pcair_dir, out_dir):
    """Concatenate per-cohort {COHORT}_pcair_pcs.tsv into pcair_pcs.tsv."""
    files = sorted(Path(pcair_dir).glob("*_pcair_pcs.tsv"))
    if not files:
        print(f"  WARN: no *_pcair_pcs.tsv under {pcair_dir} — PCA panel "
              "falls back to within-ancestry covariate PCs")
        return None
    frames = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        df["cohort"] = f.name.replace("_pcair_pcs.tsv", "")
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    path = os.path.join(out_dir, "pcair_pcs.tsv")
    out.to_csv(path, sep="\t", index=False)
    print(f"  wrote {len(out)} samples x {len(frames)} cohorts -> {path}")
    return path


# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Section 7: QC summaries — phenotype before/after processing, genotype
# stage counts, Picard sequencing metrics
# ---------------------------------------------------------------------------

def _cohort_labels(samples, ancestry_map, sidecar_dir, anc, mod):
    """Series sample_id -> cohort. Priority: per-modality sidecar
    ({anc}_{mod}_cohort_labels.tsv under sidecar_dir), then the ancestry map
    (sample_id, assigned_ancestry, cohort), then the namespaced-id prefix
    ({cohort}_{sample}; pool_modalities_within_ancestry.py convention)."""
    m = None
    if sidecar_dir:
        p = os.path.join(sidecar_dir, f"{anc}_{mod}_cohort_labels.tsv")
        if os.path.exists(p):
            d = pd.read_csv(p, sep="\t")
            m = dict(zip(d["sample_id"].astype(str), d["cohort"].astype(str)))
    if m is None and ancestry_map is not None:
        d = ancestry_map
        if "assigned_ancestry" in d.columns:
            d = d[d["assigned_ancestry"] == anc]
        m = dict(zip(d["sample_id"].astype(str), d["cohort"].astype(str)))
    out = {}
    for s in samples:
        s = str(s)
        if m and s in m:
            out[s] = m[s]
        elif "_" in s:
            out[s] = s.partition("_")[0]
        else:
            out[s] = "unknown"
    return pd.Series(out)


def _qc_pca(pheno, n_features, rng):
    """PCA scores (PC1-3) + per-PC variance explained for a phenotype matrix
    (features x samples). Features are subsampled to n_features,
    zero-variance dropped, per-feature mean-imputed and z-scored so the
    geometry is comparable across processing stages. Returns
    (DataFrame samples x PC1..PCk, array pct_var)."""
    X = pheno
    if len(X) > n_features:
        idx = np.sort(rng.choice(len(X), n_features, replace=False))
        X = X.iloc[idx]
    X = X.astype(float)
    sd = X.std(axis=1)
    X = X[sd > 0]
    X = X.sub(X.mean(axis=1), axis=0)
    X = X.fillna(0.0)  # NA -> feature mean (0 after centering)
    X = X.div(X.std(axis=1).replace(0, 1), axis=0)
    M = X.to_numpy().T  # samples x features
    k = min(3, min(M.shape))
    U, S, _ = np.linalg.svd(M, full_matrices=False)
    scores = pd.DataFrame(U[:, :k] * S[:k], index=X.columns,
                          columns=[f"PC{i + 1}" for i in range(k)])
    ev = S ** 2
    pct = ev / ev.sum()
    return scores, pct[:k]


def qc_phenotype_summaries(qtl_dir, ancestries, out_dir,
                           pooled_bed_dir=None, ancestry_map_path=None,
                           cohort_labels_dir=None, n_values=5000,
                           n_features=5000, seed=42):
    """Before/after processing summaries per ancestry x modality.

    Stages: 'before' = pooled unnormalized BED
    ({pooled_bed_dir}/{ANC}_{MOD}_pooled.bed, input to
    20_combat_modalities.sh), 'after' = final mapping BED
    ({qtl_dir}/{ANC}_{MOD}.bed.gz, post QN+INT+ComBat). Writes:
      qc_pheno_values.tsv.gz    — up to n_values sampled values per
                                  anc x mod x stage x cohort cell
      qc_pheno_quantiles.tsv    — exact 1..99% quantiles per cell
      qc_pheno_pca.tsv.gz       — sample PC1-3 scores + pct variance explained
      qc_pheno_feature_counts.tsv — n_features, n_samples per cell
    Missing stage files are skipped (cell simply absent from outputs)."""
    os.makedirs(out_dir, exist_ok=True)
    rng = np.random.default_rng(seed)
    amap = None
    if ancestry_map_path and os.path.exists(ancestry_map_path):
        amap = pd.read_csv(ancestry_map_path, sep="\t")
    val_rows, qtl_rows, pca_rows, cnt_rows = [], [], [], []
    for anc in ancestries:
        for mod in MODALITIES:
            stages = {
                "before": (os.path.join(pooled_bed_dir,
                                        f"{anc}_{mod}_pooled.bed")
                           if pooled_bed_dir else None),
                "after": os.path.join(qtl_dir, f"{anc}_{mod}.bed.gz"),
            }
            for stage, path in stages.items():
                if not path or not os.path.exists(path):
                    continue
                print(f"  {anc}/{mod}/{stage}: {os.path.basename(path)}")
                pheno, _ = read_bed(path)
                pheno = pheno.loc[:, pheno.notna().any(axis=0)]
                cohorts = _cohort_labels(pheno.columns, amap,
                                         cohort_labels_dir, anc, mod)
                cnt_rows.append({"ancestry": anc, "modality": mod,
                                 "stage": stage, "n_features": len(pheno),
                                 "n_samples": pheno.shape[1]})
                for cohort, grp in cohorts.groupby(cohorts):
                    V = pheno[grp.index].to_numpy(dtype=float).ravel()
                    V = V[~np.isnan(V)]
                    if V.size == 0:
                        continue
                    vs = (V[rng.choice(V.size, n_values, replace=False)]
                          if V.size > n_values else V)
                    val_rows.append(pd.DataFrame({
                        "ancestry": anc, "modality": mod, "stage": stage,
                        "cohort": cohort, "value": vs}))
                    qtl_rows.append(pd.DataFrame({
                        "ancestry": anc, "modality": mod, "stage": stage,
                        "cohort": cohort, "quantile": np.arange(1, 100),
                        "value": np.percentile(V, np.arange(1, 100))}))
                try:
                    scores, pct = _qc_pca(pheno, n_features, rng)
                    sdf = scores.rename_axis("sample_id").reset_index()
                    sdf["ancestry"], sdf["modality"], sdf["stage"] = \
                        anc, mod, stage
                    sdf["cohort"] = \
                        cohorts.reindex(sdf["sample_id"]).to_numpy()
                    for i, c in enumerate(scores.columns):
                        sdf[f"pct_var_{c}"] = pct[i]
                    pca_rows.append(sdf)
                except Exception as e:
                    print(f"  WARN: PCA failed for {anc}/{mod}/{stage}: {e}")
                del pheno
    paths = {}
    if val_rows:
        p = os.path.join(out_dir, "qc_pheno_values.tsv.gz")
        pd.concat(val_rows).to_csv(p, sep="\t", index=False,
                                   compression="gzip")
        paths["values"] = p
        p = os.path.join(out_dir, "qc_pheno_quantiles.tsv")
        pd.concat(qtl_rows).to_csv(p, sep="\t", index=False)
        paths["quantiles"] = p
    if pca_rows:
        p = os.path.join(out_dir, "qc_pheno_pca.tsv.gz")
        pd.concat(pca_rows).to_csv(p, sep="\t", index=False,
                                   compression="gzip")
        paths["pca"] = p
    if cnt_rows:
        p = os.path.join(out_dir, "qc_pheno_feature_counts.tsv")
        pd.DataFrame(cnt_rows).to_csv(p, sep="\t", index=False)
        paths["counts"] = p
    for k, p in paths.items():
        print(f"  wrote {k} -> {p}")
    if not paths:
        print("  WARN: no phenotype matrices found — QC outputs empty")
    return paths


def _count_data_lines(path):
    """Non-comment (non-'#') line count of a plink2 pvar/psam."""
    try:
        out = subprocess.check_output(["grep", "-c", "-v", "^#", path])
        return int(out.strip())
    except Exception:
        with open(path) as f:
            return sum(1 for line in f if not line.startswith("#"))


def qc_genotype_stage_counts(out_dir, cohort_qc_glob=None,
                             geno_pooled_dir=None, ancestries=("EAS", "EUR")):
    """Variant/sample counts per genotype processing stage. Per-cohort
    post-imputation counts from {cohort_qc_glob}/*.rsq_pass.pvar/.psam;
    pooled per-ancestry counts from {geno_pooled_dir}/{ANC}_pooled.pvar/.psam.
    Writes qc_genotype_stage_counts.tsv (entity, level, stage, n_variants,
    n_samples)."""
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    if cohort_qc_glob:
        pat = os.path.join(cohort_qc_glob, "*.rsq_pass.pvar")
        for pvar in sorted(glob.glob(pat)):
            cohort = os.path.basename(pvar).replace(".rsq_pass.pvar", "")
            psam = pvar[:-len(".pvar")] + ".psam"
            rows.append({"entity": cohort, "level": "cohort",
                         "stage": "imputed_rsq_pass",
                         "n_variants": _count_data_lines(pvar),
                         "n_samples": (_count_data_lines(psam)
                                       if os.path.exists(psam) else None)})
    if geno_pooled_dir:
        for anc in ancestries:
            pvar = os.path.join(geno_pooled_dir, f"{anc}_pooled.pvar")
            psam = os.path.join(geno_pooled_dir, f"{anc}_pooled.psam")
            if os.path.exists(pvar):
                rows.append({"entity": anc, "level": "ancestry",
                             "stage": "pooled_qc_pass",
                             "n_variants": _count_data_lines(pvar),
                             "n_samples": (_count_data_lines(psam)
                                           if os.path.exists(psam)
                                           else None)})
    if not rows:
        print("  WARN: no genotype stage inputs found — skipped")
        return None
    out = pd.DataFrame(rows)
    p = os.path.join(out_dir, "qc_genotype_stage_counts.tsv")
    out.to_csv(p, sep="\t", index=False)
    print(f"  wrote {len(out)} stage rows -> {p}")
    return p


def qc_picard_metrics(picard_glob, out_dir):
    """Concatenate per-cohort *_qc_metrics.tsv (module 19a outputs), cohort
    parsed from the filename prefix. Writes qc_picard_metrics.tsv."""
    os.makedirs(out_dir, exist_ok=True)
    files = sorted(glob.glob(picard_glob)) if picard_glob else []
    frames = []
    for f in files:
        d = pd.read_csv(f, sep="\t")
        d["cohort"] = os.path.basename(f).replace("_qc_metrics.tsv", "")
        frames.append(d)
    if not frames:
        print("  WARN: no Picard qc_metrics files matched — skipped")
        return None
    out = pd.concat(frames, ignore_index=True)
    p = os.path.join(out_dir, "qc_picard_metrics.tsv")
    out.to_csv(p, sep="\t", index=False)
    print(f"  wrote {len(out)} samples x {len(frames)} cohorts -> {p}")
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", required=True, help="qtl_results dir")
    ap.add_argument("--qtl-dir", required=True, help="qtl_inputs dir")
    ap.add_argument("--out-dir", default=None,
                    help="default: {results-dir}/report_extras")
    ap.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    ap.add_argument("--qvalue", type=float, default=0.05)
    ap.add_argument("--only", nargs="+", default=None,
                    choices=["ancestry", "modality", "locus", "annotate",
                             "constraint", "pcair", "qc"])
    ap.add_argument("--skip", nargs="+", default=[],
                    choices=["ancestry", "modality", "locus", "annotate",
                             "constraint", "pcair", "qc"])
    ap.add_argument("--locus-gene", default=None,
                    help="extra showcase locus gene_id")
    ap.add_argument("--fastvep-output", default=None)
    ap.add_argument("--ccre", nargs="*", default=[],
                    help="BED:LABEL pairs (as in 34_pip_annotation_enrichment.py)")
    ap.add_argument("--placenta-ocr", default=None)
    ap.add_argument("--gnomad-constraint", default=None,
                    help="local gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz "
                         "(downloaded if absent)")
    ap.add_argument("--pcair-dir", default=None,
                    help="dir holding per-cohort *_pcair_pcs.tsv")
    ap.add_argument("--pooled-bed-dir", default=None,
                    help="dir with {ANC}_{MOD}_pooled.bed (unnormalized; "
                         "combat_modalities/pooled) for before-stage QC")
    ap.add_argument("--ancestry-map", default=None,
                    help="pooled_sample_ancestry_RNAseq.tsv (sample_id, "
                         "assigned_ancestry, cohort) for cohort labels")
    ap.add_argument("--cohort-labels-dir", default=None,
                    help="dir with {ANC}_{MOD}_cohort_labels.tsv sidecars")
    ap.add_argument("--cohort-qc-glob", default=None,
                    help="glob matching per-cohort post-imputation QC dirs "
                         "(containing {COHORT}.rsq_pass.pvar/.psam)")
    ap.add_argument("--geno-pooled-dir", default=None,
                    help="dir with {ANC}_pooled.pvar/.psam")
    ap.add_argument("--picard-glob", default=None,
                    help="glob matching per-cohort *_qc_metrics.tsv "
                         "(module 19a outputs)")
    ap.add_argument("--n-violin-values", type=int, default=5000,
                    help="sampled values per anc x mod x stage x cohort")
    ap.add_argument("--n-pca-features", type=int, default=5000,
                    help="feature subsample size for QC PCA")
    args = ap.parse_args()

    out_dir = args.out_dir or os.path.join(args.results_dir, "report_extras")
    os.makedirs(out_dir, exist_ok=True)

    sections = args.only or ["ancestry", "modality", "locus", "annotate",
                             "constraint", "pcair", "qc"]
    sections = [s for s in sections if s not in args.skip]
    print(f"== 35_extract_report_extras == sections: {sections}")

    if "ancestry" in sections:
        print("\n[1] cross-ancestry lookup")
        cross_ancestry_lookup(args.results_dir, args.qtl_dir,
                              args.ancestries, args.qvalue, out_dir)
    if "modality" in sections:
        print("\n[2] cross-modality lookup")
        cross_modality_lookup(args.results_dir, args.qtl_dir,
                              args.ancestries, args.qvalue, out_dir)
    if "locus" in sections:
        print("\n[3] showcase-locus scans")
        locus_scans(args.results_dir, args.qtl_dir, args.ancestries,
                    args.qvalue, out_dir, extra_gene=args.locus_gene)
    if "annotate" in sections:
        print("\n[4] lead-variant annotations")
        ccre = [tuple(x.rsplit(":", 1)) for x in args.ccre]
        annotate_lead_variants(args.results_dir, args.ancestries, args.qvalue,
                               out_dir, fastvep_output=args.fastvep_output,
                               ccre_beds=ccre, placenta_ocr=args.placenta_ocr)
    if "constraint" in sections:
        print("\n[5] gnomAD constraint")
        constraint_table(args.results_dir, args.qtl_dir, args.ancestries,
                         out_dir, gnomad_path=args.gnomad_constraint)
    if "pcair" in sections:
        print("\n[6] PC-AiR staging")
        if args.pcair_dir:
            stage_pcair(args.pcair_dir, out_dir)
        else:
            print("  --pcair-dir not given — skipped")
    if "qc" in sections:
        print("\n[7] QC summaries")
        qc_phenotype_summaries(args.qtl_dir, args.ancestries, out_dir,
                               pooled_bed_dir=args.pooled_bed_dir,
                               ancestry_map_path=args.ancestry_map,
                               cohort_labels_dir=args.cohort_labels_dir,
                               n_values=args.n_violin_values,
                               n_features=args.n_pca_features)
        qc_genotype_stage_counts(out_dir,
                                 cohort_qc_glob=args.cohort_qc_glob,
                                 geno_pooled_dir=args.geno_pooled_dir,
                                 ancestries=args.ancestries)
        if args.picard_glob:
            qc_picard_metrics(args.picard_glob, out_dir)
    print("\nDone.")


if __name__ == "__main__":
    main()
