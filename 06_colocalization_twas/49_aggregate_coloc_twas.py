#!/usr/bin/env python3
"""49_aggregate_coloc_twas.py — aggregate module 06 outputs.

Objective 1.6 roll-up:

  1. SuSiE-coloc: all credible-set-pair rows -> coloc_results.tsv.gz; the
     best pair per task -> coloc_best.tsv.gz; shard diagnostics rollup ->
     coloc_diagnostics_summary.tsv.
  2. colocBoost: all region cluster rows -> colocboost_clusters.tsv.gz.
  3. isoTWAS/TWAS: all FUSION per-model rows -> twas_results.tsv.gz with
     per-(weight set, trait) BH q-values; gene-level ACAT combination across
     each gene's models -> twas_gene_results.tsv.gz.
  4. Cell-type annotation: for every gene with a coloc call (PP.H4 >=
     --pp-h4) or a gene-level TWAS q < 0.05, the primary cell type is the
     MuSiC deconvolution fraction with the largest absolute Spearman
     correlation with the gene's covariate-residualized expression
     (|rho| >= 0.2, else "unassigned"), computed per ancestry ->
     gene_celltype_annotation.tsv.

Usage:
  python3 49_aggregate_coloc_twas.py --results-dir $RESULTS_DIR \
      --qtl-dir $QTL_DIR
"""

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

MOD_COLS = ["modality", "phenotype_id", "trait_id", "ancestry"]


def bh_qvalues(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = p[order[i]] * n / rank
        prev = min(prev, val)
        q[order[i]] = prev
    return np.minimum(q, 1.0)


def acat(pvals):
    """ACAT (Cauchy) combination of p-values (equal weights)."""
    p = np.asarray(pvals, dtype=float)
    p = p[np.isfinite(p)]
    if len(p) == 0:
        return np.nan
    p = np.clip(p, 1e-300, 1 - 1e-16)
    t = np.mean(np.tan((0.5 - p) * np.pi))
    return 0.5 - np.arctan(t) / np.pi


def read_many(pattern, cols_dtype=None):
    files = sorted(Path("/").glob(str(pattern).lstrip("/")))
    frames = []
    for f in files:
        try:
            frames.append(pd.read_csv(f, sep="\t"))
        except Exception as e:
            print(f"  WARNING: could not read {f}: {e}")
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def aggregate_coloc(coloc_dir, out_dir, pp_h4):
    res = read_many(coloc_dir / "results" / "*" / "*.coloc.tsv")
    if len(res) == 0:
        print("  coloc: no result files")
        return pd.DataFrame()
    res.to_csv(out_dir / "coloc_results.tsv.gz", sep="\t", index=False)
    pp_col = next((c for c in res.columns if c.startswith("PP.H4")), None)
    best = (res.sort_values(pp_col, ascending=False)
               .groupby(MOD_COLS, as_index=False).first())
    best.to_csv(out_dir / "coloc_best.tsv.gz", sep="\t", index=False)
    n_calls = int(best[pp_col].ge(pp_h4).sum())
    print(f"  coloc: {len(res)} CS-pair rows, {len(best)} tasks, "
          f"{n_calls} colocalization calls (PP.H4 >= {pp_h4})")

    diag = read_many(coloc_dir / "diagnostics" / "*.diagnostics.tsv")
    if len(diag):
        summ = (diag.groupby(["modality", "status"]).size()
                    .rename("n").reset_index())
        summ.to_csv(out_dir / "coloc_diagnostics_summary.tsv",
                    sep="\t", index=False)
    return best


def aggregate_colocboost(coloc_dir, out_dir):
    clu = read_many(coloc_dir / "colocboost" / "results" / "*" /
                    "*.clusters.tsv")
    clu = clu[clu["cos_id"].notna()] if len(clu) else clu
    if len(clu) == 0:
        print("  colocboost: no clusters")
        return
    clu.to_csv(out_dir / "colocboost_clusters.tsv.gz", sep="\t", index=False)
    print(f"  colocboost: {len(clu)} cluster rows -> colocboost_clusters.tsv.gz")


def aggregate_twas(results_dir, out_dir):
    fusion_dir = results_dir / "isotwas" / "fusion"
    frames = []
    for f in sorted(fusion_dir.glob("*/*.twas.tsv")):
        ws = f.parent.name
        trait = f.name.replace(f"{ws}_", "").replace(".twas.tsv", "")
        df = pd.read_csv(f, sep=r"\s+")
        df["weight_set"] = ws
        df["trait_id"] = trait
        frames.append(df)
    if not frames:
        print("  twas: no FUSION outputs")
        return pd.DataFrame()
    twas = pd.concat(frames, ignore_index=True)
    twas["TWAS.P"] = pd.to_numeric(twas["TWAS.P"], errors="coerce")
    # model-level BH within weight set x trait
    twas["TWAS.Q"] = np.nan
    for (_, _), idx in twas.groupby(["weight_set", "trait_id"]).groups.items():
        p = twas.loc[idx, "TWAS.P"]
        ok = p.notna()
        twas.loc[idx[ok.to_numpy()], "TWAS.Q"] = bh_qvalues(p[ok].to_numpy())
    twas.to_csv(out_dir / "twas_results.tsv.gz", sep="\t", index=False)

    # gene level: ACAT across the gene's models
    if "GENE" not in twas.columns:
        # FUSION output lacks GENE; recover from the .pos files
        pos_frames = []
        for pos in sorted((results_dir / "isotwas" / "weights").glob("*.pos")):
            pf = pd.read_csv(pos, sep="\t")
            pf["weight_set"] = pos.stem
            pos_frames.append(pf[["ID", "GENE", "weight_set"]])
        pos_df = pd.concat(pos_frames, ignore_index=True)
        twas = twas.merge(pos_df, on=["ID", "weight_set"], how="left")
    rows = []
    for (ws, trait, gene), sub in twas.groupby(["weight_set", "trait_id",
                                                "GENE"]):
        p = sub["TWAS.P"].dropna().to_numpy()
        rows.append({"weight_set": ws, "trait_id": trait, "GENE": gene,
                     "n_models": len(p),
                     "min_p": p.min() if len(p) else np.nan,
                     "acat_p": acat(p)})
    genes = pd.DataFrame(rows)
    genes["acat_q"] = np.nan
    for (_, _), idx in genes.groupby(["weight_set", "trait_id"]).groups.items():
        p = genes.loc[idx, "acat_p"]
        ok = p.notna()
        genes.loc[idx[ok.to_numpy()], "acat_q"] = bh_qvalues(p[ok].to_numpy())
    genes.to_csv(out_dir / "twas_gene_results.tsv.gz", sep="\t", index=False)
    print(f"  twas: {len(twas)} model rows, {len(genes)} gene rows "
          f"({int((genes['acat_q'] < 0.05).sum())} genes at ACAT q < 0.05)")
    return genes


# --- cell-type annotation -----------------------------------------------------

def load_residual_expression(bed_path, cov_path, genes):
    """Return genes x samples matrix of covariate-residualized expression."""
    bed = pd.read_csv(bed_path, sep="\t", compression="gzip")
    bed = bed.rename(columns={bed.columns[0]: "chr"})
    bed = bed[bed["phenotype_id"].isin(genes)]
    if len(bed) == 0:
        return None
    sample_cols = [c for c in bed.columns
                   if c not in ("chr", "start", "end", "phenotype_id")]
    expr = bed.set_index("phenotype_id")[sample_cols].T  # samples x genes

    cov = pd.read_csv(cov_path, sep="\t", index_col=0).T  # samples x covs
    common = expr.index.intersection(cov.index)
    if len(common) < 20:
        return None
    expr = expr.loc[common]
    C = np.column_stack([np.ones(len(common)), cov.loc[common].to_numpy()])
    # OLS residuals for all genes at once
    beta, *_ = np.linalg.lstsq(C, expr.to_numpy(), rcond=None)
    resid = expr.to_numpy() - C @ beta
    return pd.DataFrame(resid, index=common, columns=expr.columns)


def annotate_cell_types(genes, qtl_dir, ancestries, out_dir):
    from scipy.stats import spearmanr
    rows = {}
    for anc in ancestries:
        deconv_path = qtl_dir / f"{anc}_deconvolution_harmonized.tsv"
        bed_path = qtl_dir / f"{anc}_expression.bed.gz"
        cov_path = qtl_dir / f"{anc}_covariates_expression.tsv"
        if not all(p.exists() for p in (deconv_path, bed_path, cov_path)):
            print(f"  celltype: skipping {anc} (missing inputs)")
            continue
        deconv = pd.read_csv(deconv_path, sep="\t")
        id_col = deconv.columns[0]
        deconv = deconv.set_index(id_col)
        resid = load_residual_expression(bed_path, cov_path, set(genes))
        if resid is None:
            print(f"  celltype: no residual expression for {anc}")
            continue
        common = resid.index.intersection(deconv.index)
        resid = resid.loc[common]
        deconv = deconv.loc[common]
        for gene in genes:
            if gene not in resid.columns:
                continue
            x = resid[gene].to_numpy()
            best_ct, best_rho = "unassigned", 0.0
            for ct in deconv.columns:
                y = deconv[ct].to_numpy()
                ok = np.isfinite(x) & np.isfinite(y)
                if ok.sum() < 20 or np.std(y[ok]) == 0:
                    continue
                rho = spearmanr(x[ok], y[ok])[0]
                if np.isfinite(rho) and abs(rho) > abs(best_rho):
                    best_rho, best_ct = rho, ct
            if abs(best_rho) < 0.2:
                best_ct = "unassigned"
            rows.setdefault(gene, {})[f"primary_cell_type_{anc}"] = best_ct
            rows[gene][f"max_abs_spearman_{anc}"] = round(best_rho, 4)
    if not rows:
        print("  celltype: nothing annotated")
        return
    ann = pd.DataFrame(rows).T.rename_axis("gene").reset_index()
    ann.to_csv(out_dir / "gene_celltype_annotation.tsv", sep="\t", index=False)
    print(f"  celltype: {len(ann)} genes -> gene_celltype_annotation.tsv")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    p.add_argument("--pp-h4", type=float, default=0.7)
    p.add_argument("--twas-q", type=float, default=0.05)
    p.add_argument("--skip-celltype", action="store_true")
    args = p.parse_args()

    results_dir = Path(args.results_dir)
    qtl_dir = Path(args.qtl_dir)
    coloc_dir = results_dir / "coloc"
    out_dir = coloc_dir / "aggregated"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[1/4] SuSiE-coloc")
    best = aggregate_coloc(coloc_dir, out_dir, args.pp_h4)
    print("[2/4] colocBoost")
    aggregate_colocboost(coloc_dir, out_dir)
    print("[3/4] isoTWAS/TWAS")
    genes_twas = aggregate_twas(results_dir, out_dir)

    print("[4/4] cell-type annotation")
    if args.skip_celltype:
        print("  skipped (--skip-celltype)")
        return
    genes = set()
    pp_col = next((c for c in best.columns if c.startswith("PP.H4")), None)
    if len(best) and pp_col:
        genes |= set(best.loc[best[pp_col] >= args.pp_h4, "phenotype_id"])
    if isinstance(genes_twas, pd.DataFrame) and len(genes_twas):
        genes |= set(genes_twas.loc[genes_twas["acat_q"] < args.twas_q, "GENE"])
    # coloc phenotype_ids may be non-gene modalities; keep gene-like ids
    genes = sorted(g for g in genes if isinstance(g, str) and "__" not in g)
    if genes:
        annotate_cell_types(genes, qtl_dir, args.ancestries, out_dir)
    else:
        print("  no genes to annotate")


if __name__ == "__main__":
    main()
