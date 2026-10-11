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
  4. Cell-type attribution (R21 Objective 1.4): for every ancestry-specific
     FDR-significant xQTL phenotype in module 05, take its exact Module-05
     discovery lead cis variant and test genotype x cell-type-proportion interactions using the
     same normalized phenotype and optimized covariates as the xQTL scan.
     Feature-level interaction statistics and FDR-gated primary cell types are
     written first; gene-level annotations for coloc/TWAS hits are then derived
     from those genetic interaction results.

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


# --- cell-type attribution ----------------------------------------------------


def _norm_name(value):
    return "".join(ch.lower() for ch in str(value) if ch.isalnum())


def _parse_sig_in(value):
    if pd.isna(value):
        return set()
    return {x.strip().upper() for x in str(value).split(",") if x.strip()}


def _open_text(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def _load_selected_bed(path, wanted):
    """Stream a BED and return samples x requested phenotypes only."""
    import csv

    wanted = set(wanted)
    if not wanted:
        return pd.DataFrame()
    rows = {}
    with _open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 5:
            raise ValueError(f"invalid/empty phenotype BED: {path}")
        samples = header[4:]
        for row in reader:
            if len(row) < 4 or row[3] not in wanted:
                continue
            vals = pd.to_numeric(pd.Series(row[4:]), errors="coerce").to_numpy()
            rows[row[3]] = vals
    if not rows:
        return pd.DataFrame(index=samples)
    return pd.DataFrame(rows, index=samples)


def _load_covariates(qtl_dir, anc, mod):
    path = qtl_dir / f"{anc}_covariates_{mod}.tsv"
    if not path.exists():
        fallback = qtl_dir / f"{anc}_covariates.tsv"
        if not fallback.exists():
            return None, path
        print(f"  celltype: WARNING {path.name} missing; using {fallback.name}")
        path = fallback
    cov = pd.read_csv(path, sep="\t", index_col=0).T
    cov.index = cov.index.astype(str)
    return cov.apply(pd.to_numeric, errors="coerce"), path


def _load_deconvolution(path):
    df = pd.read_csv(path, sep="\t")
    if len(df.columns) < 2:
        raise ValueError(f"deconvolution table has <2 columns: {path}")
    id_col = "sample_id" if "sample_id" in df.columns else df.columns[0]
    df[id_col] = df[id_col].astype(str)
    df = df.set_index(id_col)
    numeric = {}
    for col in df.columns:
        if str(col).lower() == "cohort":
            continue
        vals = pd.to_numeric(df[col], errors="coerce")
        if vals.notna().any():
            numeric[col] = vals
    out = pd.DataFrame(numeric, index=df.index)
    return out


def _resolve_cell_types(deconv, requested):
    available = list(deconv.columns)
    if requested:
        by_norm = {_norm_name(c): c for c in available}
        resolved, missing = [], []
        for name in requested:
            key = _norm_name(name)
            if key in by_norm:
                resolved.append(by_norm[key])
            else:
                missing.append(name)
        if missing:
            raise SystemExit(
                "ERROR: requested --cell-types not found in deconvolution "
                f"table: {missing}; available={available}")
        return resolved
    # The harmonized Module-05 table is already the collapsed placental cell
    # composition. Maternal is a contamination/QC fraction, not a fetal target.
    return [c for c in available if "maternal" not in _norm_name(c)]


def _maternal_column(deconv):
    return next((c for c in deconv.columns if "maternal" in _norm_name(c)), None)


def _candidate_xqtl_leads(results_dir, ancestries):
    """Return the exact Module-05 discovery lead per significant phenotype.

    ``finemap/loci/*.loci.tsv`` is the union of phenotype-level FDR discoveries
    and records the significant ancestry strata in ``sig_in``. The lead variant
    is read from the same Module-05 discovery table used to build those loci:
    expression uses ``{ANC}_expression_cisqtl_top.tsv`` and the multi-phenotype
    modalities use ``{ANC}_{MOD}_ungrouped_cisqtl_top.tsv``. A stage-39 nominal
    lead is accepted only as a backward-compatible fallback when the discovery
    table is absent.
    """
    loci_dir = results_dir / "finemap" / "loci"
    locus_files = sorted(loci_dir.glob("*.loci.tsv"))
    if not locus_files:
        print(f"  celltype: no module-05 locus files under {loci_dir}")
        return pd.DataFrame()

    sig_rows = []
    anc_set = {a.upper() for a in ancestries}
    for path in locus_files:
        try:
            loci = pd.read_csv(path, sep="\t")
        except Exception as exc:
            print(f"  celltype: WARNING could not read {path}: {exc}")
            continue
        if "phenotype_id" not in loci.columns:
            print(f"  celltype: WARNING {path.name} lacks phenotype_id")
            continue
        file_mod = path.name.replace(".loci.tsv", "")
        for _, row in loci.iterrows():
            mod = str(row.get("modality", file_mod))
            pid = str(row["phenotype_id"])
            sig_in = _parse_sig_in(row.get("sig_in", ""))
            for anc in sorted(sig_in & anc_set):
                sig_rows.append({"ancestry": anc, "modality": mod,
                                 "phenotype_id": pid})
    sig = pd.DataFrame(sig_rows).drop_duplicates() if sig_rows else pd.DataFrame()
    if sig.empty:
        print("  celltype: no ancestry-specific significant xQTL phenotypes")
        return sig

    merged = []
    for (anc, mod), sub in sig.groupby(["ancestry", "modality"]):
        result_label = mod if mod == "expression" else f"{mod}_ungrouped"
        discovery_path = results_dir / f"{anc}_{result_label}_cisqtl_top.tsv"
        fallback_path = results_dir / "nominal" / anc / f"{anc}_{mod}.nominal.top.tsv"
        top_path = discovery_path
        source = "module05_discovery"
        if not top_path.exists():
            if fallback_path.exists():
                top_path = fallback_path
                source = "stage39_fallback"
                print(f"  celltype: WARNING {discovery_path.name} missing; "
                      f"using stage-39 lead table {fallback_path.name}")
            else:
                print(f"  celltype: WARNING missing Module-05 lead table "
                      f"{discovery_path} (and no stage-39 fallback)")
                continue
        try:
            top = pd.read_csv(top_path, sep="\t")
        except Exception as exc:
            print(f"  celltype: WARNING could not read {top_path}: {exc}")
            continue
        required = {"phenotype_id", "variant_id"}
        if not required <= set(top.columns):
            print(f"  celltype: WARNING malformed {top_path.name}; "
                  f"missing {sorted(required - set(top.columns))}")
            continue
        top = top[top["phenotype_id"].astype(str).isin(set(sub["phenotype_id"]))].copy()
        if "pval_nominal" in top.columns:
            top["pval_nominal"] = pd.to_numeric(top["pval_nominal"], errors="coerce")
            top = (top.sort_values("pval_nominal", na_position="last")
                      .drop_duplicates("phenotype_id", keep="first"))
        else:
            top = top.drop_duplicates("phenotype_id", keep="first")
            top["pval_nominal"] = np.nan
        top["phenotype_id"] = top["phenotype_id"].astype(str)
        keep = top[["phenotype_id", "variant_id", "pval_nominal"]].copy()
        keep["lead_source"] = source
        keep.insert(0, "modality", mod)
        keep.insert(0, "ancestry", anc)
        merged.append(keep)
        n_missing = len(sub) - len(keep)
        if n_missing:
            print(f"  celltype: WARNING {anc}/{mod}: {n_missing} significant "
                  "phenotypes lack an ancestry-specific discovery lead")
    if not merged:
        return pd.DataFrame()
    return pd.concat(merged, ignore_index=True).drop_duplicates(
        ["ancestry", "modality", "phenotype_id"])


def _export_lead_dosages(qtl_dir, anc, variants, plink2, work_dir):
    """Export selected in-sample pgen dosages with plink2 --export A."""
    import subprocess

    variants = list(dict.fromkeys(str(v) for v in variants if pd.notna(v)))
    if not variants:
        return pd.DataFrame()
    prefix = qtl_dir / f"{anc}_qtl"
    for ext in (".pgen", ".pvar", ".psam"):
        if not Path(str(prefix) + ext).exists():
            raise FileNotFoundError(f"missing genotype input: {prefix}{ext}")
    extract = Path(work_dir) / f"{anc}.lead_variants.txt"
    out_prefix = Path(work_dir) / f"{anc}.lead_dosage"
    extract.write_text("\n".join(variants) + "\n")
    cmd = [plink2, "--pfile", str(prefix), "--extract", str(extract),
           "--export", "A", "--out", str(out_prefix), "--silent"]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    raw_path = Path(str(out_prefix) + ".raw")
    if proc.returncode != 0 or not raw_path.exists():
        msg = (proc.stderr or proc.stdout).strip()
        raise RuntimeError(f"plink2 dosage export failed for {anc}: {msg[-1200:]}")

    raw = pd.read_csv(raw_path, sep=r"\s+")
    iid_col = next((c for c in raw.columns if str(c).lstrip("#").upper() == "IID"), None)
    if iid_col is None:
        raise RuntimeError(f"cannot identify IID column in {raw_path}")
    raw[iid_col] = raw[iid_col].astype(str)
    meta_norm = {"fid", "iid", "pat", "mat", "sex", "phenotype"}
    geno_cols = [c for c in raw.columns if str(c).lstrip("#").lower() not in meta_norm]
    out = pd.DataFrame(index=raw[iid_col])
    missing = []
    for vid in variants:
        matches = [c for c in geno_cols if str(c) == vid or str(c).startswith(vid + "_")]
        if not matches:
            missing.append(vid)
            continue
        # Variant IDs are chr:pos:ref:alt and therefore do not contain '_';
        # PLINK appends '_<counted allele>' to the dosage column name.
        col = matches[0]
        out[vid] = pd.to_numeric(raw[col], errors="coerce").to_numpy()
    if missing:
        print(f"  celltype: WARNING {anc}: {len(missing)}/{len(variants)} lead "
              "variants absent from dosage export")
    return out


def _drop_target_covariate(cov, cell_type):
    if cov is None:
        return pd.DataFrame()
    target = _norm_name("ct_" + str(cell_type))
    drop = [c for c in cov.columns if _norm_name(c) == target]
    return cov.drop(columns=drop, errors="ignore")


def _fit_gxc(y, g, cell, cov, min_n=20):
    """OLS for Y ~ G + C + G:C + covariates; return interaction statistics."""
    from scipy import stats

    common = y.index.intersection(g.index).intersection(cell.index)
    if cov is not None and len(cov.columns):
        common = common.intersection(cov.index)
    if len(common) == 0:
        return {"status": "no_common_samples"}

    yv = pd.to_numeric(y.loc[common], errors="coerce").to_numpy(dtype=float)
    gv = pd.to_numeric(g.loc[common], errors="coerce").to_numpy(dtype=float)
    cv = pd.to_numeric(cell.loc[common], errors="coerce").to_numpy(dtype=float)
    if cov is not None and len(cov.columns):
        zdf = cov.loc[common].apply(pd.to_numeric, errors="coerce")
        # Subset-specific constant covariates can arise after sample alignment.
        keep_cov = [c for c in zdf.columns
                    if np.nanstd(zdf[c].to_numpy(dtype=float)) > 1e-12]
        z = zdf[keep_cov].to_numpy(dtype=float) if keep_cov else np.empty((len(common), 0))
    else:
        z = np.empty((len(common), 0))

    finite = np.isfinite(yv) & np.isfinite(gv) & np.isfinite(cv)
    if z.shape[1]:
        finite &= np.isfinite(z).all(axis=1)
    yv, gv, cv, z = yv[finite], gv[finite], cv[finite], z[finite]
    n = len(yv)
    if n < min_n:
        return {"status": "insufficient_n", "n": n}
    if np.std(gv) <= 1e-12:
        return {"status": "monomorphic", "n": n}
    if np.std(cv) <= 1e-12:
        return {"status": "no_celltype_variation", "n": n}

    af = float(np.mean(gv) / 2.0)
    maf = min(af, 1.0 - af)
    g0 = gv - np.mean(gv)
    c0 = cv - np.mean(cv)
    X = np.column_stack([np.ones(n), g0, c0, g0 * c0, z])
    rank = int(np.linalg.matrix_rank(X))
    p = X.shape[1]
    if rank < p:
        return {"status": "rank_deficient", "n": n, "af": af, "maf": maf,
                "design_rank": rank, "design_cols": p}
    df = n - p
    if df <= 2:
        return {"status": "insufficient_df", "n": n, "af": af, "maf": maf,
                "design_rank": rank, "design_cols": p}

    beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
    resid = yv - X @ beta
    sigma2 = float(np.dot(resid, resid) / df)
    xtx_inv = np.linalg.inv(X.T @ X)
    var_b = sigma2 * xtx_inv[3, 3]
    if not np.isfinite(var_b) or var_b <= 0:
        return {"status": "invalid_se", "n": n, "af": af, "maf": maf}
    se = float(np.sqrt(var_b))
    tstat = float(beta[3] / se)
    pval = float(2 * stats.t.sf(abs(tstat), df))
    q25, q50, q75 = np.quantile(cv, [0.25, 0.50, 0.75])
    iqr = float(q75 - q25)
    # cv is already arcsinh-transformed; center at the analyzed-sample mean for
    # interpretation of the genotype main effect and genetic slopes by cell mix.
    cmean = float(np.mean(cv))
    cq25, cq50, cq75 = q25 - cmean, q50 - cmean, q75 - cmean
    return {
        "status": "ok", "n": n, "af": af, "maf": maf,
        "design_rank": rank, "design_cols": p,
        "beta_G": float(beta[1]), "beta_C": float(beta[2]),
        "beta_GxC": float(beta[3]), "se_GxC": se, "p_GxC": pval,
        "cell_transformed_q25": float(q25),
        "cell_transformed_q50": float(q50),
        "cell_transformed_q75": float(q75),
        "cell_transformed_iqr": iqr,
        "beta_GxC_iqr": float(beta[3] * iqr),
        "beta_G_at_cell_q25": float(beta[1] + beta[3] * cq25),
        "beta_G_at_cell_q50": float(beta[1] + beta[3] * cq50),
        "beta_G_at_cell_q75": float(beta[1] + beta[3] * cq75),
    }


def _feature_gene(modality, phenotype_id):
    """Map TWAS-relevant xQTL phenotypes to gene IDs without guessing.

    The module-06 TWAS models are gene expression and isoform expression.
    Other modalities retain feature-level cell-type annotations but are not
    collapsed into the compatibility gene summary.
    """
    pid = str(phenotype_id)
    mod = str(modality)
    if mod == "expression":
        return pid
    if mod == "isoform_expression" and "__" in pid:
        return pid.split("__", 1)[0]
    return None


def _genes_needing_annotation(best, genes_twas, pp_h4, twas_q):
    genes = set()
    pp_col = next((c for c in best.columns if c.startswith("PP.H4")), None)
    if len(best) and pp_col:
        called = best.loc[best[pp_col] >= pp_h4]
        for _, row in called.iterrows():
            gene = _feature_gene(row.get("modality", ""), row.get("phenotype_id", ""))
            if gene:
                genes.add(gene)
    if isinstance(genes_twas, pd.DataFrame) and len(genes_twas):
        genes |= set(genes_twas.loc[genes_twas["acat_q"] < twas_q, "GENE"].astype(str))
    return sorted(genes)


def annotate_cell_types(best, genes_twas, results_dir, qtl_dir, ancestries,
                        out_dir, pp_h4=0.7, twas_q=0.05, celltype_q=0.05,
                        cell_types=None, min_n=20, maternal_threshold=0.10,
                        plink2="plink2"):
    """R21 cell-type attribution via genotype x cell-proportion interactions."""
    import tempfile

    leads = _candidate_xqtl_leads(results_dir, ancestries)
    if leads.empty:
        print("  celltype: no significant xQTL leads to test")
        return

    interaction_rows = []
    with tempfile.TemporaryDirectory(prefix="celltype_xqtl_", dir=str(out_dir)) as tmp:
        for anc in ancestries:
            anc_leads = leads[leads["ancestry"] == anc]
            if anc_leads.empty:
                continue
            deconv_path = qtl_dir / f"{anc}_deconvolution_harmonized.tsv"
            if not deconv_path.exists():
                print(f"  celltype: skipping {anc} (missing {deconv_path.name})")
                continue
            deconv_raw = _load_deconvolution(deconv_path)
            targets = _resolve_cell_types(deconv_raw, cell_types)
            if not targets:
                print(f"  celltype: skipping {anc} (no fetal cell types)")
                continue
            # Match Module 05: arcsinh-transform then mean-center each retained
            # cell-type proportion. Maternal remains raw for the >10% sensitivity.
            deconv_t = np.arcsinh(deconv_raw[targets])
            deconv_t = deconv_t - deconv_t.mean(axis=0)
            maternal_col = _maternal_column(deconv_raw)
            maternal = deconv_raw[maternal_col] if maternal_col else None

            print(f"  celltype: {anc}: {len(anc_leads)} significant xQTL phenotypes, "
                  f"{len(targets)} fetal cell types")
            try:
                dosages = _export_lead_dosages(
                    qtl_dir, anc, anc_leads["variant_id"].unique(), plink2,
                    Path(tmp))
            except Exception as exc:
                raise SystemExit(f"ERROR: cell-type genotype export failed: {exc}")

            for mod, mod_leads in anc_leads.groupby("modality"):
                bed_path = qtl_dir / f"{anc}_{mod}.bed.gz"
                cov, cov_path = _load_covariates(qtl_dir, anc, mod)
                if not bed_path.exists() or cov is None:
                    missing = []
                    if not bed_path.exists():
                        missing.append(str(bed_path))
                    if cov is None:
                        missing.append(str(cov_path))
                    print(f"  celltype: WARNING skipping {anc}/{mod}; missing "
                          + ", ".join(missing))
                    continue
                phen = _load_selected_bed(bed_path, set(mod_leads["phenotype_id"]))
                if phen.empty or phen.shape[1] == 0:
                    print(f"  celltype: WARNING no selected phenotypes in {bed_path.name}")
                    continue

                for _, lead in mod_leads.iterrows():
                    pid = str(lead["phenotype_id"])
                    vid = str(lead["variant_id"])
                    if pid not in phen.columns or vid not in dosages.columns:
                        continue
                    y = phen[pid]
                    g = dosages[vid]
                    for ct in targets:
                        cov_test = _drop_target_covariate(cov, ct)
                        fit = _fit_gxc(y, g, deconv_t[ct], cov_test, min_n=min_n)
                        rec = {
                            "ancestry": anc, "modality": mod,
                            "phenotype_id": pid, "variant_id": vid,
                            "lead_pval_nominal": lead.get("pval_nominal", np.nan),
                            "lead_source": lead.get("lead_source", ""),
                            "cell_type": ct,
                            "cell_prop_raw_q25": float(deconv_raw[ct].quantile(0.25)),
                            "cell_prop_raw_q50": float(deconv_raw[ct].quantile(0.50)),
                            "cell_prop_raw_q75": float(deconv_raw[ct].quantile(0.75)),
                        }
                        rec.update(fit)

                        # R21 sensitivity: exclude samples with >10% inferred
                        # maternal contribution, when a maternal fraction exists.
                        rec.update({
                            "maternal_sensitivity_n": np.nan,
                            "maternal_sensitivity_beta_GxC": np.nan,
                            "maternal_sensitivity_se_GxC": np.nan,
                            "maternal_sensitivity_p_GxC": np.nan,
                            "maternal_sensitivity_status": "no_maternal_fraction",
                        })
                        if maternal is not None:
                            keep = maternal <= maternal_threshold
                            y_s = y.where(keep.reindex(y.index).fillna(False))
                            sens = _fit_gxc(y_s, g, deconv_t[ct], cov_test,
                                            min_n=min_n)
                            rec["maternal_sensitivity_status"] = sens.get("status")
                            rec["maternal_sensitivity_n"] = sens.get("n", np.nan)
                            rec["maternal_sensitivity_beta_GxC"] = sens.get(
                                "beta_GxC", np.nan)
                            rec["maternal_sensitivity_se_GxC"] = sens.get(
                                "se_GxC", np.nan)
                            rec["maternal_sensitivity_p_GxC"] = sens.get(
                                "p_GxC", np.nan)
                        interaction_rows.append(rec)

    if not interaction_rows:
        print("  celltype: no GxC interaction models could be fit")
        return

    inter = pd.DataFrame(interaction_rows)
    inter["q_GxC"] = np.nan
    for (_, _), idx in inter.groupby(["ancestry", "modality"]).groups.items():
        pvals = pd.to_numeric(inter.loc[idx, "p_GxC"], errors="coerce")
        ok = pvals.notna() & (inter.loc[idx, "status"] == "ok")
        if ok.any():
            good_idx = pvals.index[ok]
            inter.loc[good_idx, "q_GxC"] = bh_qvalues(pvals.loc[good_idx].to_numpy())
    inter_path = out_dir / "xqtl_celltype_interactions.tsv.gz"
    inter.to_csv(inter_path, sep="\t", index=False)

    # One primary cell type per significant xQTL feature. Assignment is
    # significance-gated; best_cell_type is retained even when unassigned.
    feat_rows = []
    for key, sub in inter.groupby(["ancestry", "modality", "phenotype_id"], sort=False):
        valid = sub[(sub["status"] == "ok") & sub["p_GxC"].notna()].copy()
        base = {"ancestry": key[0], "modality": key[1], "phenotype_id": key[2]}
        if valid.empty:
            feat_rows.append({**base, "primary_cell_type": "unassigned",
                              "best_cell_type": np.nan, "primary_beta_GxC": np.nan,
                              "primary_se_GxC": np.nan, "primary_p_GxC": np.nan,
                              "primary_q_GxC": np.nan, "variant_id": np.nan,
                              "n_celltypes_tested": 0})
            continue
        valid = valid.sort_values(["q_GxC", "p_GxC"], na_position="last")
        b = valid.iloc[0]
        q = b["q_GxC"]
        assigned = str(b["cell_type"]) if pd.notna(q) and q <= celltype_q else "unassigned"
        feat_rows.append({
            **base, "variant_id": b["variant_id"],
            "best_cell_type": b["cell_type"], "primary_cell_type": assigned,
            "primary_beta_GxC": b["beta_GxC"], "primary_se_GxC": b["se_GxC"],
            "primary_p_GxC": b["p_GxC"], "primary_q_GxC": q,
            "n_celltypes_tested": len(valid),
            "maternal_sensitivity_beta_GxC": b["maternal_sensitivity_beta_GxC"],
            "maternal_sensitivity_p_GxC": b["maternal_sensitivity_p_GxC"],
        })
    feat = pd.DataFrame(feat_rows)
    feat_path = out_dir / "xqtl_celltype_annotation.tsv.gz"
    feat.to_csv(feat_path, sep="\t", index=False)

    # Compatibility/output layer for Objective 1.6: annotate genes implicated by
    # coloc/TWAS using their strongest feature-level genetic cell-type signal.
    genes = _genes_needing_annotation(best, genes_twas, pp_h4, twas_q)
    gene_rows = {g: {"gene": g} for g in genes}
    if genes:
        feat = feat.copy()
        feat["gene"] = [
            _feature_gene(m, p) for m, p in zip(feat["modality"], feat["phenotype_id"])
        ]
        feat = feat[feat["gene"].isin(genes)]
        for anc in ancestries:
            for gene in genes:
                sub = feat[(feat["ancestry"] == anc) & (feat["gene"] == gene)].copy()
                if sub.empty:
                    gene_rows[gene][f"primary_cell_type_{anc}"] = "unassigned"
                    continue
                sub = sub.sort_values(["primary_q_GxC", "primary_p_GxC"],
                                      na_position="last")
                b = sub.iloc[0]
                gene_rows[gene].update({
                    f"primary_cell_type_{anc}": b["primary_cell_type"],
                    f"best_cell_type_{anc}": b["best_cell_type"],
                    f"primary_GxC_beta_{anc}": b["primary_beta_GxC"],
                    f"primary_GxC_se_{anc}": b["primary_se_GxC"],
                    f"primary_GxC_p_{anc}": b["primary_p_GxC"],
                    f"primary_GxC_q_{anc}": b["primary_q_GxC"],
                    f"source_modality_{anc}": b["modality"],
                    f"source_phenotype_id_{anc}": b["phenotype_id"],
                    f"lead_variant_{anc}": b["variant_id"],
                })
        ann = pd.DataFrame(gene_rows.values())
        ann.to_csv(out_dir / "gene_celltype_annotation.tsv", sep="\t", index=False)
        print(f"  celltype: {len(inter)} GxC tests, {len(feat_rows)} feature "
              f"annotations, {len(ann)} coloc/TWAS genes")
    else:
        print(f"  celltype: {len(inter)} GxC tests, {len(feat_rows)} feature "
              "annotations; no coloc/TWAS genes met annotation thresholds")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    p.add_argument("--pp-h4", type=float, default=0.7)
    p.add_argument("--twas-q", type=float, default=0.05)
    p.add_argument("--celltype-q", type=float, default=0.05,
                   help="BH FDR threshold for assigning a primary cell type")
    p.add_argument("--cell-types", nargs="+", default=None,
                   help="fetal cell-type columns to test; default: all numeric "
                        "deconvolution columns except Maternal")
    p.add_argument("--celltype-min-n", type=int, default=20,
                   help="minimum complete samples per GxC model")
    p.add_argument("--maternal-threshold", type=float, default=0.10,
                   help="maternal-fraction exclusion threshold for sensitivity")
    p.add_argument("--plink2", default="plink2",
                   help="plink2 executable for extracting lead-variant dosages")
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

    print("[4/4] cell-type attribution (genotype x cell proportion)")
    if args.skip_celltype:
        print("  skipped (--skip-celltype)")
        return
    annotate_cell_types(
        best, genes_twas, results_dir, qtl_dir, args.ancestries, out_dir,
        pp_h4=args.pp_h4, twas_q=args.twas_q, celltype_q=args.celltype_q,
        cell_types=args.cell_types, min_n=args.celltype_min_n,
        maternal_threshold=args.maternal_threshold, plink2=args.plink2)


if __name__ == "__main__":
    main()
