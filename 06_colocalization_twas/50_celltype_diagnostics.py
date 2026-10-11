#!/usr/bin/env python3
"""50_celltype_diagnostics.py — QC/power diagnostics for Module-06 GxC attribution.

This is a diagnostic follow-up to 49_aggregate_coloc_twas.py.  It does not
change cell-type assignments.  It is designed to explain ancestry-specific
power/calibration differences when ancestries have very different cohort sizes.
EUR N=136 is treated as the known full EUR cohort size, not as unexplained
sample attrition.  The primary diagnostic repeatedly downsamples EAS to the EUR
cohort size while preserving the production model and multiple-testing scheme.

Core outputs
------------
  sample_overlap_summary.tsv             (QC only; EUR N=136 is expected)
  sample_overlap_by_sample.tsv.gz
  celltype_distribution.tsv
  gxc_support.tsv.gz
  gxc_hc3.tsv.gz                         (unless --skip-hc3)
  gxc_permutation.tsv.gz                 (when --permutations > 0)
  gxc_downsample_replicates.tsv          (EAS -> EUR-sized repeated downsampling)
  gxc_downsample_tests.tsv.gz
  gxc_downsample_features.tsv.gz
  gxc_matched_ancestry.tsv.gz
  gxc_matched_ancestry_summary.tsv
  diagnostic_summary.tsv
  run_metadata.json
  celltype_diagnostics_bundle.zip        (upload this bundle for review)

The script reuses the input readers and cell-type transformations from the
co-located 49_aggregate_coloc_twas.py so diagnostics stay synchronized with the
production attribution model.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import math
import multiprocessing as mp
import os
import platform
import socket
import sys
import tempfile
import time
import warnings
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

KEYS = ["ancestry", "modality", "phenotype_id", "variant_id", "cell_type"]
_DIAG_STATE = None


def _norm_name(value):
    return "".join(ch.lower() for ch in str(value) if ch.isalnum())


def _nonempty(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _auto_workers(requested=None):
    if requested is not None:
        return max(1, int(requested))
    for key in ("LSB_DJOB_NUMPROC", "SLURM_CPUS_PER_TASK", "NSLOTS"):
        value = os.environ.get(key)
        if value:
            try:
                return max(1, int(value))
            except ValueError:
                pass
    return 1


def _load_stage49(path: Path):
    path = Path(path)
    if not path.exists():
        raise SystemExit(f"ERROR: stage-49 script not found: {path}")
    spec = importlib.util.spec_from_file_location("module06_stage49", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"ERROR: cannot import stage-49 helpers from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    required = [
        "_load_selected_bed", "_load_covariates", "_load_deconvolution",
        "_export_lead_dosages", "_drop_target_covariate", "bh_qvalues",
    ]
    missing = [x for x in required if not hasattr(mod, x)]
    if missing:
        raise SystemExit(
            "ERROR: stage-49 script is missing helpers required by diagnostics: "
            + ", ".join(missing)
        )
    return mod


def _open_text(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def _bed_samples(path: Path):
    with _open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if len(header) < 5:
        raise ValueError(f"invalid BED header: {path}")
    return [str(x) for x in header[4:]]


def _psam_samples(prefix: Path):
    path = Path(str(prefix) + ".psam")
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep=r"\s+", dtype=str)
    iid = next((c for c in df.columns if str(c).lstrip("#").upper() == "IID"), None)
    if iid is None:
        raise ValueError(f"cannot identify IID in {path}")
    return df[iid].astype(str).tolist()


def _load_deconv_full(path: Path):
    """Return numeric deconvolution matrix plus optional cohort labels."""
    df = pd.read_csv(path, sep="\t")
    if len(df.columns) < 2:
        raise ValueError(f"deconvolution table has <2 columns: {path}")
    id_col = "sample_id" if "sample_id" in df.columns else df.columns[0]
    df[id_col] = df[id_col].astype(str)
    df = df.set_index(id_col)

    cohort_col = next((c for c in df.columns if _norm_name(c) == "cohort"), None)
    cohort = None
    if cohort_col is not None:
        cohort = df[cohort_col].astype(str)
        cohort.name = "cohort"

    numeric = {}
    for col in df.columns:
        if col == cohort_col:
            continue
        vals = pd.to_numeric(df[col], errors="coerce")
        if vals.notna().any():
            numeric[col] = vals
    return pd.DataFrame(numeric, index=df.index), cohort


def _finite_numeric_rows(df: pd.DataFrame):
    if df is None or df.shape[1] == 0:
        return pd.Series(True, index=(df.index if df is not None else pd.Index([])))
    z = df.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    return pd.Series(np.isfinite(z).all(axis=1), index=df.index)


def _describe_values(values):
    x = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return {
            "n": 0, "mean": np.nan, "sd": np.nan, "min": np.nan,
            "q25": np.nan, "q50": np.nan, "q75": np.nan, "iqr": np.nan,
            "max": np.nan, "zero_fraction": np.nan,
        }
    q25, q50, q75 = np.quantile(x, [0.25, 0.50, 0.75])
    return {
        "n": int(len(x)), "mean": float(np.mean(x)),
        "sd": float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
        "min": float(np.min(x)), "q25": float(q25), "q50": float(q50),
        "q75": float(q75), "iqr": float(q75 - q25), "max": float(np.max(x)),
        "zero_fraction": float(np.mean(np.isclose(x, 0.0, atol=1e-15))),
    }


def _sample_overlap_tables(qtl_dir, interactions, ancestries, stage49):
    summaries = []
    sample_rows = []
    distribution_rows = []

    for anc in ancestries:
        anc_i = interactions[interactions["ancestry"] == anc]
        if anc_i.empty:
            continue
        prefix = qtl_dir / f"{anc}_qtl"
        geno_ids = set(_psam_samples(prefix))
        deconv_path = qtl_dir / f"{anc}_deconvolution_harmonized.tsv"
        deconv_raw, cohort = _load_deconv_full(deconv_path)
        deconv_ids = set(deconv_raw.index.astype(str))
        targets = sorted(set(anc_i["cell_type"].dropna().astype(str)))
        targets = [c for c in targets if c in deconv_raw.columns]

        # All deconvolution samples: useful to distinguish cohort-wide cell-range
        # differences from range restriction introduced by the QTL intersection.
        for ct in targets:
            raw_desc = _describe_values(deconv_raw[ct])
            trans = np.arcsinh(pd.to_numeric(deconv_raw[ct], errors="coerce"))
            trans = trans - np.nanmean(trans)
            tr_desc = _describe_values(trans)
            distribution_rows.append({
                "ancestry": anc, "modality": "ALL", "scope": "deconvolution_all",
                "cell_type": ct,
                **{f"raw_{k}": v for k, v in raw_desc.items()},
                **{f"transformed_{k}": v for k, v in tr_desc.items()},
            })

        for mod in sorted(anc_i["modality"].dropna().astype(str).unique()):
            bed_path = qtl_dir / f"{anc}_{mod}.bed.gz"
            if not bed_path.exists():
                warnings.warn(f"missing BED for overlap diagnostics: {bed_path}")
                continue
            bed_ids = set(_bed_samples(bed_path))
            cov, cov_path = stage49._load_covariates(qtl_dir, anc, mod)
            if cov is None:
                warnings.warn(f"missing covariates for overlap diagnostics: {cov_path}")
                continue
            cov_ids = set(cov.index.astype(str))
            cov_complete = _finite_numeric_rows(cov)
            cov_complete_ids = set(cov_complete.index[cov_complete].astype(str))
            universe = sorted(geno_ids | bed_ids | cov_ids | deconv_ids)
            common = geno_ids & bed_ids & cov_ids & deconv_ids
            common_cov_complete = common & cov_complete_ids

            observed = pd.to_numeric(
                anc_i.loc[anc_i["modality"] == mod, "n"], errors="coerce"
            ).dropna()
            summaries.append({
                "ancestry": anc, "modality": mod,
                "n_genotype_psam": len(geno_ids),
                "n_phenotype_bed": len(bed_ids),
                "n_covariate_ids": len(cov_ids),
                "n_covariate_complete_all": len(cov_complete_ids),
                "n_deconvolution": len(deconv_ids),
                "n_id_intersection": len(common),
                "n_id_intersection_cov_complete_all": len(common_cov_complete),
                "model_n_min": float(observed.min()) if len(observed) else np.nan,
                "model_n_median": float(observed.median()) if len(observed) else np.nan,
                "model_n_max": float(observed.max()) if len(observed) else np.nan,
            })

            for sid in universe:
                sample_rows.append({
                    "ancestry": anc, "modality": mod, "sample_id": sid,
                    "in_genotype_psam": sid in geno_ids,
                    "in_phenotype_bed": sid in bed_ids,
                    "in_covariates": sid in cov_ids,
                    "covariates_complete_all": sid in cov_complete_ids,
                    "in_deconvolution": sid in deconv_ids,
                    "in_all_id_sets": sid in common,
                    "in_all_and_cov_complete_all": sid in common_cov_complete,
                    "cohort": cohort.get(sid, np.nan) if cohort is not None else np.nan,
                })

            overlap_ids = sorted(common)
            for ct in targets:
                raw = deconv_raw.loc[deconv_raw.index.intersection(overlap_ids), ct]
                raw_desc = _describe_values(raw)
                trans = np.arcsinh(pd.to_numeric(raw, errors="coerce"))
                # Production stage 49 centers using the ancestry-wide
                # deconvolution mean, not the modality-specific overlap mean.
                all_trans = np.arcsinh(pd.to_numeric(deconv_raw[ct], errors="coerce"))
                trans = trans - np.nanmean(all_trans)
                tr_desc = _describe_values(trans)
                distribution_rows.append({
                    "ancestry": anc, "modality": mod, "scope": "analysis_id_overlap",
                    "cell_type": ct,
                    **{f"raw_{k}": v for k, v in raw_desc.items()},
                    **{f"transformed_{k}": v for k, v in tr_desc.items()},
                })

    return (pd.DataFrame(summaries), pd.DataFrame(sample_rows),
            pd.DataFrame(distribution_rows))


def _prepare_design(y, g, cell_t, cell_raw, cov, cohort=None):
    """Build the exact GxC design and return arrays plus attrition metrics."""
    common = y.index.intersection(g.index).intersection(cell_t.index)
    common = common.intersection(cell_raw.index)
    if cov is not None and cov.shape[1]:
        common = common.intersection(cov.index)
    n_id = len(common)
    if n_id == 0:
        return None, {"n_id_intersection": 0, "n_complete": 0}

    yv = pd.to_numeric(y.loc[common], errors="coerce").to_numpy(dtype=float)
    gv = pd.to_numeric(g.loc[common], errors="coerce").to_numpy(dtype=float)
    cv = pd.to_numeric(cell_t.loc[common], errors="coerce").to_numpy(dtype=float)
    cr = pd.to_numeric(cell_raw.loc[common], errors="coerce").to_numpy(dtype=float)

    zdf = None
    z = np.empty((n_id, 0))
    keep_cov = []
    if cov is not None and cov.shape[1]:
        zdf = cov.loc[common].apply(pd.to_numeric, errors="coerce")
        keep_cov = [
            c for c in zdf.columns
            if np.nanstd(zdf[c].to_numpy(dtype=float)) > 1e-12
        ]
        if keep_cov:
            z = zdf[keep_cov].to_numpy(dtype=float)

    miss_y = ~np.isfinite(yv)
    miss_g = ~np.isfinite(gv)
    miss_c = ~np.isfinite(cv) | ~np.isfinite(cr)
    miss_z = ~np.isfinite(z).all(axis=1) if z.shape[1] else np.zeros(n_id, dtype=bool)
    finite = ~(miss_y | miss_g | miss_c | miss_z)

    attr = {
        "n_id_intersection": int(n_id),
        "n_missing_y": int(miss_y.sum()),
        "n_missing_g": int(miss_g.sum()),
        "n_missing_cell": int(miss_c.sum()),
        "n_missing_covariate_any": int(miss_z.sum()),
        "n_nonfinite_any": int((~finite).sum()),
        "n_complete": int(finite.sum()),
        "n_covariates_retained": int(len(keep_cov)),
    }
    if not finite.any():
        return None, attr

    yv, gv, cv, cr, z = yv[finite], gv[finite], cv[finite], cr[finite], z[finite]
    ids = np.asarray(common.astype(str))[finite]
    if cohort is None:
        cohorts = np.array(["ALL"] * len(yv), dtype=object)
    else:
        cohorts = cohort.reindex(ids).fillna("MISSING").astype(str).to_numpy()

    g0 = gv - np.mean(gv)
    c0 = cv - np.mean(cv)
    w = g0 * c0
    X0 = np.column_stack([np.ones(len(yv)), g0, c0, z])
    X = np.column_stack([np.ones(len(yv)), g0, c0, w, z])
    return {
        "y": yv, "g": gv, "cell_t": cv, "cell_raw": cr, "z": z,
        "g0": g0, "c0": c0, "w": w, "X0": X0, "X": X,
        "ids": ids, "cohort": cohorts,
    }, attr


def _support_metrics(d):
    g = d["g"]
    c = d["cell_t"]
    cr = d["cell_raw"]
    X0 = d["X0"]
    w = d["w"]
    n = len(g)

    af = float(np.mean(g) / 2.0)
    maf = min(af, 1.0 - af)
    minor = g if af <= 0.5 else 2.0 - g
    hard = np.clip(np.rint(minor), 0, 2).astype(int)
    n0, n1, n2 = (int(np.sum(hard == k)) for k in (0, 1, 2))
    carrier = hard >= 1

    # Cell-composition tertiles. arcsinh is monotone, so raw/transformed
    # tertiles induce the same ordering; use raw fractions for readability.
    q33, q67 = np.quantile(cr, [1 / 3, 2 / 3])
    bins = np.where(cr <= q33, 0, np.where(cr <= q67, 1, 2))
    carrier_counts = [int(np.sum(carrier & (bins == k))) for k in range(3)]
    bin_counts = [int(np.sum(bins == k)) for k in range(3)]

    coef_w, *_ = np.linalg.lstsq(X0, w, rcond=None)
    wr = w - X0 @ coef_w
    interaction_info = float(np.dot(wr, wr))

    corr = np.corrcoef(g, c)[0, 1] if np.std(g) > 0 and np.std(c) > 0 else np.nan
    return {
        "recalc_af": af, "recalc_maf": maf,
        "n_minor_dosage0": n0, "n_minor_dosage1": n1, "n_minor_dosage2": n2,
        "n_minor_carriers": int(carrier.sum()),
        "minor_allele_count_dosage": float(np.sum(minor)),
        "minor_carriers_low_cell_tertile": carrier_counts[0],
        "minor_carriers_mid_cell_tertile": carrier_counts[1],
        "minor_carriers_high_cell_tertile": carrier_counts[2],
        "min_minor_carriers_across_cell_tertiles": int(min(carrier_counts)),
        "n_low_cell_tertile": bin_counts[0],
        "n_mid_cell_tertile": bin_counts[1],
        "n_high_cell_tertile": bin_counts[2],
        "corr_G_cell_transformed": float(corr) if np.isfinite(corr) else np.nan,
        "cell_raw_sd_complete": float(np.std(cr, ddof=1)) if n > 1 else 0.0,
        "cell_raw_iqr_complete": float(np.quantile(cr, .75) - np.quantile(cr, .25)),
        "cell_transformed_sd_complete": float(np.std(c, ddof=1)) if n > 1 else 0.0,
        "cell_transformed_iqr_complete": float(np.quantile(c, .75) - np.quantile(c, .25)),
        "interaction_raw_sd": float(np.std(w, ddof=1)) if n > 1 else 0.0,
        "interaction_resid_sd": float(np.std(wr, ddof=1)) if n > 1 else 0.0,
        "interaction_information": interaction_info,
    }


def _ols_hc3_metrics(d):
    X, y = d["X"], d["y"]
    n, p = X.shape
    beta, _, rank, singular = np.linalg.lstsq(X, y, rcond=None)
    rank = int(rank)
    out = {"diag_design_rank": rank, "diag_design_cols": int(p)}
    if rank < p:
        out["diag_status"] = "rank_deficient"
        return out
    df = n - p
    if df <= 2:
        out["diag_status"] = "insufficient_df"
        return out

    resid = y - X @ beta
    xtx_inv = np.linalg.inv(X.T @ X)
    sigma2 = float(np.dot(resid, resid) / df)
    classical_var = sigma2 * xtx_inv[3, 3]
    if classical_var <= 0 or not np.isfinite(classical_var):
        out["diag_status"] = "invalid_classical_se"
        return out
    classical_se = float(np.sqrt(classical_var))
    classical_t = float(beta[3] / classical_se)
    classical_p = float(2 * stats.t.sf(abs(classical_t), df))

    h = np.einsum("ij,jk,ik->i", X, xtx_inv, X)
    one_minus_h = np.clip(1.0 - h, 1e-8, None)
    omega = (resid / one_minus_h) ** 2
    meat = X.T @ (X * omega[:, None])
    hc3_cov = xtx_inv @ meat @ xtx_inv
    hc3_var = float(hc3_cov[3, 3])
    if hc3_var <= 0 or not np.isfinite(hc3_var):
        out["diag_status"] = "invalid_hc3_se"
        return out
    hc3_se = float(np.sqrt(hc3_var))
    hc3_t = float(beta[3] / hc3_se)
    hc3_p = float(2 * stats.t.sf(abs(hc3_t), df))
    cond = float(singular[0] / singular[-1]) if len(singular) and singular[-1] > 0 else np.inf
    out.update({
        "diag_status": "ok", "diag_df": int(df),
        "recalc_beta_GxC": float(beta[3]),
        "recalc_classical_se_GxC": classical_se,
        "recalc_classical_p_GxC": classical_p,
        "hc3_se_GxC": hc3_se, "hc3_t_GxC": hc3_t, "hc3_p_GxC": hc3_p,
        "hc3_se_ratio": hc3_se / classical_se,
        "max_leverage": float(np.max(h)), "median_leverage": float(np.median(h)),
        "design_condition_number": cond,
    })
    return out


def _stable_seed(base_seed, key):
    text = "|".join(map(str, key)).encode("utf-8")
    digest = hashlib.blake2b(text, digest_size=8).digest()
    return (int.from_bytes(digest, "little") + int(base_seed)) % (2**32 - 1)

def _make_downsample_sets(ids, n_target, n_reps, seed, key):
    """Return reproducible common-subject subsets for one ancestry x modality.

    The sampling frame is the analysis-ID intersection for that modality, not
    the genotype PSAM alone.  This isolates the consequence of reducing sample
    size while avoiding avoidable losses from samples that could never enter a
    stage-49 model.
    """
    ids = np.asarray(sorted(set(map(str, ids))), dtype=object)
    if n_reps <= 0 or n_target is None or n_target <= 0:
        return []
    if len(ids) < n_target:
        return []
    rng = np.random.default_rng(_stable_seed(seed, key))
    return [set(rng.choice(ids, size=n_target, replace=False).tolist())
            for _ in range(int(n_reps))]


def _fit_classical_subset(d, keep_ids, min_n=20):
    """Fast classical GxC fit after restricting a reconstructed design to IDs."""
    ids = d["ids"]
    mask = np.fromiter((sid in keep_ids for sid in ids), dtype=bool, count=len(ids))
    X = d["X"][mask]
    y = d["y"][mask]
    n, p = X.shape
    if n < min_n or n <= p + 2:
        return np.nan, np.nan, int(n), "insufficient_n"
    xtx = X.T @ X
    xty = X.T @ y
    try:
        beta = np.linalg.solve(xtx, xty)
        e3 = np.zeros(p)
        e3[3] = 1.0
        inv_col = np.linalg.solve(xtx, e3)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, int(n), "rank_deficient"
    resid = y - X @ beta
    df = n - p
    sigma2 = float(np.dot(resid, resid) / df)
    var_b = sigma2 * float(inv_col[3])
    if not np.isfinite(var_b) or var_b <= 0:
        return np.nan, np.nan, int(n), "invalid_se"
    se = math.sqrt(var_b)
    tval = float(beta[3] / se)
    pval = float(2 * stats.t.sf(abs(tval), df))
    return float(beta[3]), pval, int(n), "ok"


def _summarize_downsample_modality(records, stage49, source_anc, mod,
                                   n_target, reference_anc):
    """Apply BH per replicate and reduce a large downsampling experiment.

    Returns compact replicate-, test-, and feature-level summaries rather than
    a multi-million-row replicate-by-test table.
    """
    if not records:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    n_reps = len(records[0]["pvals"])
    n_tests = len(records)
    pmat = np.vstack([r["pvals"] for r in records]).astype(float)
    bmat = np.vstack([r["betas"] for r in records]).astype(float)
    nmat = np.vstack([r["ns"] for r in records]).astype(float)
    qmat = np.full_like(pmat, np.nan, dtype=float)
    for r in range(n_reps):
        ok = np.isfinite(pmat[:, r])
        if ok.any():
            qmat[ok, r] = stage49.bh_qvalues(pmat[ok, r])

    feature_ids = np.asarray([str(x["phenotype_id"]) for x in records], dtype=object)
    cell_types = np.asarray([str(x["cell_type"]) for x in records], dtype=object)
    variants = np.asarray([str(x["variant_id"]) for x in records], dtype=object)
    unique_features = pd.unique(feature_ids)
    feat_to_idx = {f: np.where(feature_ids == f)[0] for f in unique_features}

    replicate_rows = []
    primary_by_rep = []
    for r in range(n_reps):
        p = pmat[:, r]
        q = qmat[:, r]
        ok = np.isfinite(p)
        sig = np.isfinite(q) & (q < 0.05)
        assigned = 0
        valid_features = 0
        prim = {}
        for f, idx in feat_to_idx.items():
            valid_idx = idx[np.isfinite(p[idx])]
            if len(valid_idx):
                valid_features += 1
            sig_idx = idx[sig[idx]]
            if len(sig_idx):
                assigned += 1
                # Stage 49 chooses the smallest q, then p.
                best = min(sig_idx, key=lambda j: (q[j], p[j]))
                prim[f] = cell_types[best]
        primary_by_rep.append(prim)
        replicate_rows.append({
            "source_ancestry": source_anc,
            "reference_ancestry": reference_anc,
            "modality": mod,
            "replicate": r + 1,
            "target_n": int(n_target),
            "n_tests_total": int(n_tests),
            "n_tests_valid": int(ok.sum()),
            "median_model_n": float(np.nanmedian(nmat[:, r])) if np.isfinite(nmat[:, r]).any() else np.nan,
            "nominal_p_lt_0.05_rate": float(np.mean(p[ok] < .05)) if ok.any() else np.nan,
            "fdr_q_lt_0.05_rate": float(np.mean(q[ok] < .05)) if ok.any() else np.nan,
            "n_features_total": int(len(unique_features)),
            "n_features_with_valid_test": int(valid_features),
            "n_features_assigned": int(assigned),
            "feature_assignment_rate": float(assigned / len(unique_features)) if len(unique_features) else np.nan,
        })

    test_rows = []
    for i, rec in enumerate(records):
        valid = np.isfinite(pmat[i])
        qs = qmat[i]
        bs = bmat[i]
        full_beta = rec.get("full_beta", np.nan)
        sign_ok = np.isfinite(bs) & np.isfinite(full_beta) & (full_beta != 0)
        sign_conc = np.mean(np.sign(bs[sign_ok]) == np.sign(full_beta)) if sign_ok.any() else np.nan
        test_rows.append({
            "source_ancestry": source_anc,
            "reference_ancestry": reference_anc,
            "modality": mod,
            "phenotype_id": rec["phenotype_id"],
            "variant_id": rec["variant_id"],
            "cell_type": rec["cell_type"],
            "target_n": int(n_target),
            "full_n": rec.get("full_n", np.nan),
            "full_beta_GxC": full_beta,
            "full_p_GxC": rec.get("full_p", np.nan),
            "full_q_GxC": rec.get("full_q", np.nan),
            "downsample_valid_reps": int(valid.sum()),
            "downsample_median_n": float(np.nanmedian(nmat[i])) if np.isfinite(nmat[i]).any() else np.nan,
            "downsample_median_beta_GxC": float(np.nanmedian(bs)) if np.isfinite(bs).any() else np.nan,
            "downsample_beta_sign_concordance": float(sign_conc) if np.isfinite(sign_conc) else np.nan,
            "downsample_median_p_GxC": float(np.nanmedian(pmat[i])) if valid.any() else np.nan,
            "downsample_nominal_p_lt_0.05_rate": float(np.mean(pmat[i, valid] < .05)) if valid.any() else np.nan,
            "downsample_median_q_GxC": float(np.nanmedian(qs)) if np.isfinite(qs).any() else np.nan,
            "downsample_fdr_q_lt_0.05_rate": float(np.mean(qs[np.isfinite(qs)] < .05)) if np.isfinite(qs).any() else np.nan,
        })

    feature_rows = []
    for f in unique_features:
        idx = feat_to_idx[f]
        assigned_flags = []
        primary = []
        for r in range(n_reps):
            sig_idx = idx[np.isfinite(qmat[idx, r]) & (qmat[idx, r] < .05)]
            assigned_flags.append(len(sig_idx) > 0)
            if len(sig_idx):
                best = min(sig_idx, key=lambda j: (qmat[j, r], pmat[j, r]))
                primary.append(cell_types[best])
        mode_ct = np.nan
        mode_frac = np.nan
        if primary:
            vc = pd.Series(primary).value_counts()
            mode_ct = vc.index[0]
            mode_frac = float(vc.iloc[0] / len(primary))
        feature_rows.append({
            "source_ancestry": source_anc,
            "reference_ancestry": reference_anc,
            "modality": mod,
            "phenotype_id": f,
            "variant_id": variants[idx[0]],
            "target_n": int(n_target),
            "n_replicates": int(n_reps),
            "downsample_assignment_rate": float(np.mean(assigned_flags)),
            "downsample_n_assigned_reps": int(np.sum(assigned_flags)),
            "downsample_modal_primary_cell_type": mode_ct,
            "downsample_modal_primary_fraction_among_assigned": mode_frac,
        })
    return pd.DataFrame(replicate_rows), pd.DataFrame(test_rows), pd.DataFrame(feature_rows)


def _matched_ancestry_tables(merged, annotations, anc_a="EAS", anc_b="EUR"):
    """Compare the same modality/phenotype/cell-type tests across ancestries.

    Lead variants can differ by ancestry.  The detailed table therefore carries
    both variant IDs and a same_variant flag; direct beta concordance is only
    summarized for same-variant pairs.
    """
    base_cols = [
        "modality", "phenotype_id", "cell_type", "variant_id", "n",
        "beta_GxC", "se_GxC", "p_GxC", "q_GxC", "recalc_maf",
        "cell_raw_iqr_complete", "interaction_information",
    ]
    have = [c for c in base_cols if c in merged.columns]
    a = merged[merged["ancestry"] == anc_a][have].copy()
    b = merged[merged["ancestry"] == anc_b][have].copy()
    keys = ["modality", "phenotype_id", "cell_type"]
    matched = a.merge(b, on=keys, how="inner", suffixes=(f"_{anc_a}", f"_{anc_b}"))
    if matched.empty:
        return matched, pd.DataFrame()
    va, vb = f"variant_id_{anc_a}", f"variant_id_{anc_b}"
    matched["same_variant"] = matched[va].astype(str) == matched[vb].astype(str)
    for metric in ("n", "recalc_maf", "cell_raw_iqr_complete", "interaction_information",
                   "se_GxC", "p_GxC", "q_GxC"):
        ca, cb = f"{metric}_{anc_a}", f"{metric}_{anc_b}"
        if ca in matched.columns and cb in matched.columns:
            matched[f"{metric}_ratio_{anc_a}_over_{anc_b}"] = (
                pd.to_numeric(matched[ca], errors="coerce") /
                pd.to_numeric(matched[cb], errors="coerce").replace(0, np.nan)
            )

    rows = []
    for mod, sub in matched.groupby("modality", sort=True):
        same = sub[sub["same_variant"]]
        rec = {
            "modality": mod,
            "n_matched_tests": int(len(sub)),
            "n_matched_features": int(sub["phenotype_id"].nunique()),
            "same_variant_rate": float(sub["same_variant"].mean()),
        }
        for anc in (anc_a, anc_b):
            rec[f"median_n_{anc}"] = float(pd.to_numeric(sub.get(f"n_{anc}"), errors="coerce").median())
            rec[f"nominal_p_lt_0.05_rate_{anc}"] = float((pd.to_numeric(sub.get(f"p_GxC_{anc}"), errors="coerce") < .05).mean())
            rec[f"fdr_q_lt_0.05_rate_{anc}"] = float((pd.to_numeric(sub.get(f"q_GxC_{anc}"), errors="coerce") < .05).mean())
            if f"se_GxC_{anc}" in sub:
                rec[f"median_se_GxC_{anc}"] = float(pd.to_numeric(sub[f"se_GxC_{anc}"], errors="coerce").median())
            if f"interaction_information_{anc}" in sub:
                rec[f"median_interaction_information_{anc}"] = float(pd.to_numeric(sub[f"interaction_information_{anc}"], errors="coerce").median())
        if len(same) and f"beta_GxC_{anc_a}" in same and f"beta_GxC_{anc_b}" in same:
            ba = pd.to_numeric(same[f"beta_GxC_{anc_a}"], errors="coerce")
            bb = pd.to_numeric(same[f"beta_GxC_{anc_b}"], errors="coerce")
            ok = ba.notna() & bb.notna()
            rec["same_variant_n_beta_comparable"] = int(ok.sum())
            rec["same_variant_effect_sign_concordance"] = float((np.sign(ba[ok]) == np.sign(bb[ok])).mean()) if ok.any() else np.nan
            rec["same_variant_beta_spearman"] = float(stats.spearmanr(ba[ok], bb[ok]).statistic) if ok.sum() >= 3 else np.nan
        rows.append(rec)
    return matched, pd.DataFrame(rows)


def _permute_residual_test(d, n_perm, seed):
    """Freedman-Lane-style reduced-model residual permutation within cohort.

    The reduced model contains G, C, and the covariates; the tested term is GxC.
    Residuals are shuffled within cohort when cohort labels are available.
    """
    if n_perm <= 0:
        return {}
    X0, y, w = d["X0"], d["y"], d["w"]
    cohort = np.asarray(d["cohort"], dtype=object)
    n = len(y)
    p_full = X0.shape[1] + 1
    df = n - p_full
    if df <= 2:
        return {"perm_status": "insufficient_df"}

    b0, _, rank0, _ = np.linalg.lstsq(X0, y, rcond=None)
    if int(rank0) < X0.shape[1]:
        return {"perm_status": "reduced_rank_deficient"}
    resid0 = y - X0 @ b0
    bw, *_ = np.linalg.lstsq(X0, w, rcond=None)
    wr = w - X0 @ bw
    info = float(np.dot(wr, wr))
    if not np.isfinite(info) or info <= 1e-12:
        return {"perm_status": "no_interaction_information"}

    # FWL observed statistic.
    yr = resid0
    beta_obs = float(np.dot(wr, yr) / info)
    full_resid = yr - beta_obs * wr
    se_obs = math.sqrt(float(np.dot(full_resid, full_resid) / df) / info)
    if not np.isfinite(se_obs) or se_obs <= 0:
        return {"perm_status": "invalid_observed_se"}
    t_obs = beta_obs / se_obs

    groups = [np.where(cohort == label)[0] for label in pd.unique(cohort)]
    rng = np.random.default_rng(seed)
    exceed = 0
    t_perm_abs = np.empty(n_perm, dtype=float)
    for b in range(n_perm):
        ep = resid0.copy()
        for idx in groups:
            if len(idx) > 1:
                ep[idx] = ep[rng.permutation(idx)]
        # Shuffling can reintroduce projection onto X0; residualize again.
        bp, *_ = np.linalg.lstsq(X0, ep, rcond=None)
        epr = ep - X0 @ bp
        beta = float(np.dot(wr, epr) / info)
        rr = epr - beta * wr
        se = math.sqrt(float(np.dot(rr, rr) / df) / info)
        tval = beta / se if se > 0 and np.isfinite(se) else 0.0
        t_perm_abs[b] = abs(tval)
        exceed += abs(tval) >= abs(t_obs)
    empirical = (exceed + 1.0) / (n_perm + 1.0)
    return {
        "perm_status": "ok", "permutations": int(n_perm),
        "permutation_strata": "within_cohort" if len(groups) > 1 else "global",
        "n_permutation_strata": int(len(groups)),
        "perm_observed_t": float(t_obs),
        "perm_empirical_p": float(empirical),
        "perm_abs_t_q95": float(np.quantile(t_perm_abs, 0.95)),
        "perm_abs_t_q99": float(np.quantile(t_perm_abs, 0.99)),
    }


def _diagnose_feature(task):
    state = _DIAG_STATE
    if state is None:
        raise RuntimeError("diagnostic worker state not initialized")
    pid, vid = task
    phen = state["phen"]
    dosages = state["dosages"]
    if pid not in phen.columns or vid not in dosages.columns:
        return []
    y, g = phen[pid], dosages[vid]
    rows = []
    target_rows = state["tests_by_feature"].get((pid, vid), [])
    for ct in target_rows:
        if ct not in state["deconv_t"].columns:
            continue
        cov_test = state["drop_target"](state["cov"], ct)
        d, attr = _prepare_design(
            y, g, state["deconv_t"][ct], state["deconv_raw"][ct],
            cov_test, state.get("cohort")
        )
        rec = {"phenotype_id": pid, "variant_id": vid, "cell_type": ct, **attr}
        if d is None or attr.get("n_complete", 0) < state["min_n"]:
            rec["diag_status"] = "insufficient_n"
            rows.append(rec)
            continue
        if np.std(d["g"]) <= 1e-12:
            rec["diag_status"] = "monomorphic"
            rows.append(rec)
            continue
        if np.std(d["cell_t"]) <= 1e-12:
            rec["diag_status"] = "no_celltype_variation"
            rows.append(rec)
            continue
        rec.update(_support_metrics(d))
        if state["do_hc3"]:
            rec.update(_ols_hc3_metrics(d))
        else:
            # We still need condition/rank diagnostics cheaply enough to assess
            # near-collinearity even when HC3 is skipped.
            _, _, rank, singular = np.linalg.lstsq(d["X"], d["y"], rcond=None)
            rec["diag_design_rank"] = int(rank)
            rec["diag_design_cols"] = int(d["X"].shape[1])
            rec["diag_status"] = "ok" if int(rank) == d["X"].shape[1] else "rank_deficient"
            rec["design_condition_number"] = (
                float(singular[0] / singular[-1])
                if len(singular) and singular[-1] > 0 else np.inf
            )
        key = (state["ancestry"], state["modality"], pid, vid, ct)
        label = state["permute_keys"].get(key)
        if label and state["n_perm"] > 0 and rec.get("diag_status") == "ok":
            rec["permutation_sample"] = label
            rec.update(_permute_residual_test(
                d, state["n_perm"], _stable_seed(state["seed"], key)
            ))

        # Exact repeated EAS -> EUR-sized fits.  The same participant subset is
        # used for every feature/cell type within a modality and replicate.
        if state.get("downsample_sets") and rec.get("diag_status") == "ok":
            betas, pvals, ns, statuses = [], [], [], []
            for keep_ids in state["downsample_sets"]:
                b, pv, nn, st = _fit_classical_subset(
                    d, keep_ids, min_n=state["min_n"]
                )
                betas.append(b); pvals.append(pv); ns.append(nn); statuses.append(st)
            rec["_downsample"] = {
                "phenotype_id": pid, "variant_id": vid, "cell_type": ct,
                "full_n": rec.get("n_complete", np.nan),
                "full_beta": rec.get("recalc_beta_GxC", np.nan),
                "full_p": rec.get("recalc_classical_p_GxC", np.nan),
                "betas": np.asarray(betas, dtype=float),
                "pvals": np.asarray(pvals, dtype=float),
                "ns": np.asarray(ns, dtype=float),
                "statuses": statuses,
            }
        rows.append(rec)
    return rows


def _select_permutation_tests(interactions, n_top, n_random, seed):
    selected = {}
    rng = np.random.default_rng(seed)
    if n_top <= 0 and n_random <= 0:
        return selected
    valid = interactions[(interactions["status"] == "ok") & interactions["p_GxC"].notna()].copy()
    for (anc, mod), sub in valid.groupby(["ancestry", "modality"], sort=False):
        sub = sub.sort_values("p_GxC")
        top = sub.head(n_top)
        for row in top.itertuples(index=False):
            key = (str(anc), str(mod), str(row.phenotype_id), str(row.variant_id), str(row.cell_type))
            selected[key] = "top_nominal"
        remaining = sub.iloc[len(top):]
        if len(remaining) and n_random > 0:
            k = min(n_random, len(remaining))
            idx = rng.choice(len(remaining), size=k, replace=False)
            for row in remaining.iloc[np.sort(idx)].itertuples(index=False):
                key = (str(anc), str(mod), str(row.phenotype_id), str(row.variant_id), str(row.cell_type))
                selected.setdefault(key, "random")
    return selected


def _run_model_diagnostics(qtl_dir, interactions, ancestries, stage49, args, tmpdir):
    global _DIAG_STATE
    all_rows = []
    down_rep_frames, down_test_frames, down_feat_frames = [], [], []
    permute_keys = _select_permutation_tests(
        interactions, args.permute_top, args.permute_random, args.seed
    ) if args.permutations > 0 else {}

    for anc in ancestries:
        anc_i = interactions[interactions["ancestry"] == anc].copy()
        if anc_i.empty:
            continue
        deconv_path = qtl_dir / f"{anc}_deconvolution_harmonized.tsv"
        deconv_raw, cohort = _load_deconv_full(deconv_path)
        targets = sorted(set(anc_i["cell_type"].dropna().astype(str)))
        targets = [c for c in targets if c in deconv_raw.columns]
        deconv_t = np.arcsinh(deconv_raw[targets].apply(pd.to_numeric, errors="coerce"))
        deconv_t = deconv_t - deconv_t.mean(axis=0)

        variants = anc_i["variant_id"].dropna().astype(str).unique()
        print(f"  {anc}: exporting {len(variants)} lead variants")
        dosages = stage49._export_lead_dosages(
            qtl_dir, anc, variants, args.plink2, Path(tmpdir)
        )

        for mod in sorted(anc_i["modality"].dropna().astype(str).unique()):
            sub = anc_i[anc_i["modality"] == mod].copy()
            bed_path = qtl_dir / f"{anc}_{mod}.bed.gz"
            cov, cov_path = stage49._load_covariates(qtl_dir, anc, mod)
            if not bed_path.exists() or cov is None:
                print(f"  WARNING: skipping {anc}/{mod}; missing BED/covariates")
                continue
            leads = sub[["phenotype_id", "variant_id"]].drop_duplicates().copy()
            phen = stage49._load_selected_bed(bed_path, set(leads["phenotype_id"].astype(str)))
            if phen.empty or phen.shape[1] == 0:
                print(f"  WARNING: skipping {anc}/{mod}; no selected phenotypes")
                continue

            tests_by_feature = {}
            for (pid, vid), ss in sub.groupby(["phenotype_id", "variant_id"], sort=False):
                tests_by_feature[(str(pid), str(vid))] = list(ss["cell_type"].astype(str))
            tasks = list(tests_by_feature.keys())

            downsample_sets = []
            if (anc == args.downsample_ancestry and args.downsample_reps > 0
                    and args.downsample_n is not None):
                eligible = (set(phen.index.astype(str)) & set(dosages.index.astype(str)) &
                            set(deconv_raw.index.astype(str)) & set(cov.index.astype(str)))
                downsample_sets = _make_downsample_sets(
                    eligible, args.downsample_n, args.downsample_reps, args.seed,
                    ("downsample", anc, mod, args.downsample_n)
                )
                if downsample_sets:
                    print(f"    downsampling {anc}/{mod}: {len(eligible)} eligible -> "
                          f"N={args.downsample_n}, {len(downsample_sets)} replicate(s)")
                else:
                    print(f"    WARNING: cannot downsample {anc}/{mod}: "
                          f"eligible N={len(eligible)} < target N={args.downsample_n}")

            _DIAG_STATE = {
                "ancestry": anc, "modality": mod,
                "phen": phen, "dosages": dosages,
                "deconv_raw": deconv_raw, "deconv_t": deconv_t,
                "cohort": cohort, "cov": cov,
                "tests_by_feature": tests_by_feature,
                "drop_target": stage49._drop_target_covariate,
                "min_n": args.min_n, "do_hc3": not args.skip_hc3,
                "permute_keys": permute_keys, "n_perm": args.permutations,
                "seed": args.seed, "downsample_sets": downsample_sets,
            }
            workers = args.workers
            print(
                f"  {anc}/{mod}: {len(tasks)} features, {len(sub)} GxC tests, "
                f"{workers} worker(s)"
            )
            mod_rows = []
            if workers > 1 and sys.platform.startswith("linux"):
                ctx = mp.get_context("fork")
                chunksize = max(1, len(tasks) // (workers * 8))
                with ctx.Pool(workers) as pool:
                    for i, rows in enumerate(
                        pool.imap_unordered(_diagnose_feature, tasks, chunksize=chunksize), 1
                    ):
                        mod_rows.extend(rows)
                        if i % 250 == 0 or i == len(tasks):
                            print(f"    completed {i}/{len(tasks)} features")
            else:
                if workers > 1:
                    print("    multiprocessing fork unavailable; using serial fits")
                for i, task in enumerate(tasks, 1):
                    mod_rows.extend(_diagnose_feature(task))
                    if i % 250 == 0 or i == len(tasks):
                        print(f"    completed {i}/{len(tasks)} features")
            down_records = []
            for rec in mod_rows:
                rec["ancestry"] = anc
                rec["modality"] = mod
                ds = rec.pop("_downsample", None)
                if ds is not None:
                    full_key = ((sub["phenotype_id"].astype(str) == str(ds["phenotype_id"])) &
                                (sub["variant_id"].astype(str) == str(ds["variant_id"])) &
                                (sub["cell_type"].astype(str) == str(ds["cell_type"])))
                    ss = sub.loc[full_key]
                    if len(ss):
                        ds["full_p"] = pd.to_numeric(ss.iloc[0].get("p_GxC"), errors="coerce")
                        ds["full_q"] = pd.to_numeric(ss.iloc[0].get("q_GxC"), errors="coerce")
                        ds["full_beta"] = pd.to_numeric(ss.iloc[0].get("beta_GxC"), errors="coerce")
                        ds["full_n"] = pd.to_numeric(ss.iloc[0].get("n"), errors="coerce")
                    down_records.append(ds)
            all_rows.extend(mod_rows)
            if down_records:
                rr, tt, ff = _summarize_downsample_modality(
                    down_records, stage49, anc, mod, args.downsample_n,
                    args.downsample_reference_ancestry
                )
                if len(rr): down_rep_frames.append(rr)
                if len(tt): down_test_frames.append(tt)
                if len(ff): down_feat_frames.append(ff)
            _DIAG_STATE = None

    def _cat(frames):
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return (pd.DataFrame(all_rows), permute_keys,
            _cat(down_rep_frames), _cat(down_test_frames), _cat(down_feat_frames))


def _add_hc3_qvalues(hc3, stage49):
    hc3 = hc3.copy()
    hc3["hc3_q_GxC"] = np.nan
    for (_, _), idx in hc3.groupby(["ancestry", "modality"]).groups.items():
        p = pd.to_numeric(hc3.loc[idx, "hc3_p_GxC"], errors="coerce")
        ok = p.notna() & (hc3.loc[idx, "diag_status"] == "ok")
        if ok.any():
            good = p.index[ok]
            hc3.loc[good, "hc3_q_GxC"] = stage49.bh_qvalues(p.loc[good].to_numpy())
    return hc3


def _summary_table(merged, annotations, hc3=None, perm=None):
    rows = []
    groups = list(merged.groupby(["ancestry", "modality"], sort=True))
    # Add ancestry-wide totals to make the EAS/EUR contrast immediately visible.
    groups += [((anc, "ALL"), sub) for anc, sub in merged.groupby("ancestry", sort=True)]
    for (anc, mod), sub in groups:
        valid = sub[sub["diag_status"] == "ok"]
        original_ok = sub[sub["status"] == "ok"] if "status" in sub.columns else sub
        ann = annotations[annotations["ancestry"] == anc]
        if mod != "ALL":
            ann = ann[ann["modality"] == mod]
        rec = {
            "ancestry": anc, "modality": mod,
            "n_features": int(sub[["phenotype_id", "variant_id"]].drop_duplicates().shape[0]),
            "n_tests": int(len(sub)), "n_diag_ok": int(len(valid)),
            "model_n_median": float(pd.to_numeric(original_ok["n"], errors="coerce").median()),
            "complete_n_median_recalc": float(pd.to_numeric(valid["n_complete"], errors="coerce").median()) if len(valid) else np.nan,
            "maf_median": float(pd.to_numeric(valid["recalc_maf"], errors="coerce").median()) if len(valid) else np.nan,
            "cell_raw_iqr_median": float(pd.to_numeric(valid["cell_raw_iqr_complete"], errors="coerce").median()) if len(valid) else np.nan,
            "interaction_resid_sd_median": float(pd.to_numeric(valid["interaction_resid_sd"], errors="coerce").median()) if len(valid) else np.nan,
            "interaction_information_median": float(pd.to_numeric(valid["interaction_information"], errors="coerce").median()) if len(valid) else np.nan,
            "min_minor_carriers_tertile_median": float(pd.to_numeric(valid["min_minor_carriers_across_cell_tertiles"], errors="coerce").median()) if len(valid) else np.nan,
            "design_condition_median": float(pd.to_numeric(valid["design_condition_number"], errors="coerce").replace([np.inf, -np.inf], np.nan).median()) if len(valid) else np.nan,
            "original_p_lt_0.05_rate": float((pd.to_numeric(original_ok["p_GxC"], errors="coerce") < .05).mean()) if len(original_ok) else np.nan,
            "original_q_lt_0.05_rate": float((pd.to_numeric(original_ok["q_GxC"], errors="coerce") < .05).mean()) if len(original_ok) else np.nan,
            "assigned_feature_rate": float((ann["primary_cell_type"].astype(str) != "unassigned").mean()) if len(ann) else np.nan,
        }
        if hc3 is not None and len(hc3):
            hh = hc3[(hc3["ancestry"] == anc)]
            if mod != "ALL":
                hh = hh[hh["modality"] == mod]
            hh = hh[hh["diag_status"] == "ok"]
            rec.update({
                "hc3_se_ratio_median": float(pd.to_numeric(hh["hc3_se_ratio"], errors="coerce").median()) if len(hh) else np.nan,
                "hc3_p_lt_0.05_rate": float((pd.to_numeric(hh["hc3_p_GxC"], errors="coerce") < .05).mean()) if len(hh) else np.nan,
                "hc3_q_lt_0.05_rate": float((pd.to_numeric(hh["hc3_q_GxC"], errors="coerce") < .05).mean()) if len(hh) else np.nan,
                "original_sig_lost_with_hc3": int(((pd.to_numeric(hh["q_GxC"], errors="coerce") < .05) & (pd.to_numeric(hh["hc3_q_GxC"], errors="coerce") >= .05)).sum()) if len(hh) else 0,
                "hc3_sig_gained": int(((pd.to_numeric(hh["q_GxC"], errors="coerce") >= .05) & (pd.to_numeric(hh["hc3_q_GxC"], errors="coerce") < .05)).sum()) if len(hh) else 0,
            })
        if perm is not None and len(perm):
            pp = perm[perm["ancestry"] == anc]
            if mod != "ALL":
                pp = pp[pp["modality"] == mod]
            pp = pp[pp["perm_status"] == "ok"]
            random = pp[pp["permutation_sample"] == "random"]
            top = pp[pp["permutation_sample"] == "top_nominal"]
            rec.update({
                "perm_random_n": int(len(random)),
                "perm_random_empirical_p_lt_0.05_rate": float((pd.to_numeric(random["perm_empirical_p"], errors="coerce") < .05).mean()) if len(random) else np.nan,
                "perm_top_n": int(len(top)),
                "perm_top_empirical_p_lt_0.05_rate": float((pd.to_numeric(top["perm_empirical_p"], errors="coerce") < .05).mean()) if len(top) else np.nan,
            })
        rows.append(rec)
    return pd.DataFrame(rows)


def _write_readme(out_dir, args, paths):
    text = f"""Module 06 cell-type attribution diagnostics
Generated: {datetime.now(timezone.utc).isoformat()}

Purpose
-------
This directory diagnoses ancestry-specific power/calibration of the Module-06
GxC attribution model Y ~ G + C + GxC + Z.  It does not modify assignments.

Files
-----
sample_overlap_summary.tsv
    Input/sample overlap QC only.  EUR N=136 is the known full cohort size, so
    this table is used only to detect additional losses below the cohort N.

sample_overlap_by_sample.tsv.gz
    Sample-level membership flags for each ancestry x modality.  This identifies
    exactly which samples disappear at each input intersection.

celltype_distribution.tsv
    Raw and arcsinh-transformed cell-fraction distributions for all deconvolution
    samples and for each modality's analysis ID intersection.

gxc_support.tsv.gz
    One row per GxC test.  Includes exact complete-case attrition, MAF/genotype
    counts, minor-allele carrier counts across cell-fraction tertiles, G-C
    correlation, residualized GxC variance/information, leverage, and design
    condition number.

gxc_hc3.tsv.gz
    Classical vs HC3 heteroskedasticity-robust interaction SE/P/FDR comparison.
    Present only when HC3 was requested.

gxc_permutation.tsv.gz
    Reduced-model residual permutation diagnostic, shuffled within cohort when
    a cohort label is available.  Contains top-nominal and random tests from each
    ancestry x modality stratum.  Present only when --permutations > 0.

gxc_downsample_replicates.tsv
    Repeated exact {args.downsample_ancestry} -> {args.downsample_reference_ancestry}-sized
    refits, summarized per replicate x modality.  This is the primary test of
    whether N alone explains the ancestry difference.

gxc_downsample_tests.tsv.gz
    Per-test retention of nominal/FDR significance and effect direction across
    the repeated downsampled fits.

gxc_downsample_features.tsv.gz
    Per-feature probability of receiving any FDR-significant cell-type assignment
    after downsampling, plus the modal assigned cell type.

gxc_matched_ancestry.tsv.gz / gxc_matched_ancestry_summary.tsv
    Same modality/phenotype/cell-type tests present in EAS and EUR, with both
    ancestry-specific lead variants and a same_variant flag.

diagnostic_summary.tsv
    Compact ancestry x modality and ancestry-wide summary of the above metrics.

run_metadata.json
    Command/provenance and diagnostic settings.

celltype_diagnostics_bundle.zip
    Convenience bundle containing all reportable diagnostics.  Upload this ZIP
    for downstream summarization.

Settings
--------
workers: {args.workers}
minimum N: {args.min_n}
HC3: {not args.skip_hc3}
permutations per selected test: {args.permutations}
permutation top tests / stratum: {args.permute_top}
permutation random tests / stratum: {args.permute_random}
seed: {args.seed}
downsample ancestry: {args.downsample_ancestry}
downsample reference ancestry: {args.downsample_reference_ancestry}
downsample target N: {args.downsample_n}
downsample replicates: {args.downsample_reps}
"""
    (out_dir / "README.txt").write_text(text)


def _bundle_outputs(out_dir):
    bundle = out_dir / "celltype_diagnostics_bundle.zip"
    members = [
        "README.txt", "sample_overlap_summary.tsv", "sample_overlap_by_sample.tsv.gz",
        "celltype_distribution.tsv", "gxc_support.tsv.gz", "gxc_hc3.tsv.gz",
        "gxc_permutation.tsv.gz", "gxc_downsample_replicates.tsv",
        "gxc_downsample_tests.tsv.gz", "gxc_downsample_features.tsv.gz",
        "gxc_matched_ancestry.tsv.gz", "gxc_matched_ancestry_summary.tsv",
        "diagnostic_summary.tsv", "run_metadata.json",
    ]
    with zipfile.ZipFile(bundle, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name in members:
            path = out_dir / name
            if _nonempty(path):
                zf.write(path, arcname=name)
    return bundle


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    p.add_argument("--out-dir", default=None,
                   help="default: RESULTS_DIR/coloc/aggregated/celltype_diagnostics")
    p.add_argument("--aggregate-script", default=None,
                   help="stage-49 script to reuse for exact input readers/transforms")
    p.add_argument("--plink2", default="plink2")
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--min-n", type=int, default=20)
    p.add_argument("--skip-hc3", action="store_true")
    p.add_argument("--permutations", type=int, default=0,
                   help="residual permutations per selected GxC test; 0 disables")
    p.add_argument("--permute-top", type=int, default=20,
                   help="smallest nominal P tests per ancestry x modality")
    p.add_argument("--permute-random", type=int, default=40,
                   help="random additional tests per ancestry x modality")
    p.add_argument("--downsample-ancestry", default="EAS",
                   help="ancestry to repeatedly downsample (default: EAS)")
    p.add_argument("--downsample-reference-ancestry", default="EUR",
                   help="ancestry whose cohort N defines the default target (default: EUR)")
    p.add_argument("--downsample-n", type=int, default=None,
                   help="target N; default is reference ancestry QTL PSAM sample count")
    p.add_argument("--downsample-reps", type=int, default=100,
                   help="exact repeated downsample refits per modality; 0 disables")
    p.add_argument("--seed", type=int, default=20261010)
    p.add_argument("--force", action="store_true")
    args = p.parse_args()
    args.workers = _auto_workers(args.workers)
    if args.min_n < 3:
        p.error("--min-n must be >= 3")
    if args.permutations < 0 or args.permute_top < 0 or args.permute_random < 0:
        p.error("permutation counts must be >= 0")
    if args.downsample_reps < 0:
        p.error("--downsample-reps must be >= 0")

    results_dir = Path(args.results_dir).resolve()
    qtl_dir = Path(args.qtl_dir).resolve()
    if args.downsample_reps > 0 and args.downsample_n is None:
        ref_prefix = qtl_dir / f"{args.downsample_reference_ancestry}_qtl"
        args.downsample_n = len(_psam_samples(ref_prefix))
    if args.downsample_n is not None and args.downsample_n < args.min_n:
        p.error("--downsample-n must be >= --min-n")
    if args.downsample_reps > 0:
        print(f"Downsample experiment: {args.downsample_ancestry} -> "
              f"{args.downsample_reference_ancestry}-sized N={args.downsample_n}, "
              f"{args.downsample_reps} replicate(s)")
    agg_dir = results_dir / "coloc" / "aggregated"
    out_dir = Path(args.out_dir).resolve() if args.out_dir else agg_dir / "celltype_diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)
    done = out_dir / ".done"
    bundle = out_dir / "celltype_diagnostics_bundle.zip"
    if done.exists() and _nonempty(bundle) and not args.force:
        print(f"SKIP: completed diagnostics exist: {bundle}")
        return

    interaction_path = agg_dir / "xqtl_celltype_interactions.tsv.gz"
    annotation_path = agg_dir / "xqtl_celltype_annotation.tsv.gz"
    if not _nonempty(interaction_path) or not _nonempty(annotation_path):
        raise SystemExit(
            "ERROR: stage-49 cell-type outputs are required under " + str(agg_dir)
        )
    interactions = pd.read_csv(interaction_path, sep="\t")
    annotations = pd.read_csv(annotation_path, sep="\t")
    for col in ("ancestry", "modality", "phenotype_id", "variant_id", "cell_type", "status"):
        if col not in interactions.columns:
            raise SystemExit(f"ERROR: {interaction_path.name} lacks required column {col}")
    for col in ("ancestry", "modality", "phenotype_id", "primary_cell_type"):
        if col not in annotations.columns:
            raise SystemExit(f"ERROR: {annotation_path.name} lacks required column {col}")
    interactions["p_GxC"] = pd.to_numeric(interactions.get("p_GxC"), errors="coerce")
    interactions["q_GxC"] = pd.to_numeric(interactions.get("q_GxC"), errors="coerce")

    aggregate_script = Path(args.aggregate_script) if args.aggregate_script else Path(__file__).with_name("49_aggregate_coloc_twas.py")
    stage49 = _load_stage49(aggregate_script)

    print("[1/5] sample overlap QC + cell-composition range")
    sample_summary, sample_detail, cell_dist = _sample_overlap_tables(
        qtl_dir, interactions, args.ancestries, stage49
    )
    sample_summary.to_csv(out_dir / "sample_overlap_summary.tsv", sep="\t", index=False)
    sample_detail.to_csv(out_dir / "sample_overlap_by_sample.tsv.gz", sep="\t", index=False)
    cell_dist.to_csv(out_dir / "celltype_distribution.tsv", sep="\t", index=False)

    print("[2/5] GxC support / robust-SE + EAS->EUR-sized downsampling")
    with tempfile.TemporaryDirectory(prefix="celltype_diag_", dir=str(out_dir)) as tmp:
        diag, permute_keys, down_rep, down_test, down_feat = _run_model_diagnostics(
            qtl_dir, interactions, args.ancestries, stage49, args, tmp
        )
    if diag.empty:
        raise SystemExit("ERROR: no diagnostic GxC models could be reconstructed")

    merged = interactions.merge(diag, on=KEYS, how="left", validate="one_to_one")
    if "recalc_classical_p_GxC" in merged.columns:
        merged["classical_p_abs_diff"] = (
            pd.to_numeric(merged["p_GxC"], errors="coerce") -
            pd.to_numeric(merged["recalc_classical_p_GxC"], errors="coerce")
        ).abs()
    support_cols = KEYS + [
        c for c in [
            "status", "n", "maf", "p_GxC", "q_GxC", "n_id_intersection",
            "n_missing_y", "n_missing_g", "n_missing_cell", "n_missing_covariate_any",
            "n_nonfinite_any", "n_complete", "n_covariates_retained", "diag_status",
            "diag_design_rank", "diag_design_cols", "recalc_af", "recalc_maf",
            "n_minor_dosage0", "n_minor_dosage1", "n_minor_dosage2", "n_minor_carriers",
            "minor_allele_count_dosage", "minor_carriers_low_cell_tertile",
            "minor_carriers_mid_cell_tertile", "minor_carriers_high_cell_tertile",
            "min_minor_carriers_across_cell_tertiles", "n_low_cell_tertile",
            "n_mid_cell_tertile", "n_high_cell_tertile", "corr_G_cell_transformed",
            "cell_raw_sd_complete", "cell_raw_iqr_complete", "cell_transformed_sd_complete",
            "cell_transformed_iqr_complete", "interaction_raw_sd", "interaction_resid_sd",
            "interaction_information", "max_leverage", "median_leverage",
            "design_condition_number",
        ] if c in merged.columns
    ]
    merged[support_cols].to_csv(out_dir / "gxc_support.tsv.gz", sep="\t", index=False)

    hc3 = None
    if not args.skip_hc3 and "hc3_p_GxC" in merged.columns:
        hc3_cols = KEYS + [
            c for c in [
                "status", "n", "beta_GxC", "se_GxC", "p_GxC", "q_GxC", "diag_status",
                "recalc_beta_GxC", "recalc_classical_se_GxC", "recalc_classical_p_GxC",
                "classical_p_abs_diff", "hc3_se_GxC", "hc3_t_GxC", "hc3_p_GxC",
                "hc3_se_ratio", "max_leverage", "median_leverage", "design_condition_number",
            ] if c in merged.columns
        ]
        hc3 = _add_hc3_qvalues(merged[hc3_cols], stage49)
        hc3.to_csv(out_dir / "gxc_hc3.tsv.gz", sep="\t", index=False)

    if len(down_rep):
        down_rep.to_csv(out_dir / "gxc_downsample_replicates.tsv", sep="\t", index=False)
        down_test.to_csv(out_dir / "gxc_downsample_tests.tsv.gz", sep="\t", index=False)
        # Attach the production full-N assignment to feature summaries.
        src_ann = annotations[annotations["ancestry"] == args.downsample_ancestry][
            ["modality", "phenotype_id", "primary_cell_type"]
        ].copy()
        src_ann = src_ann.rename(columns={"primary_cell_type": "full_primary_cell_type"})
        down_feat = down_feat.merge(src_ann, on=["modality", "phenotype_id"], how="left")
        down_feat["full_assigned"] = down_feat["full_primary_cell_type"].astype(str) != "unassigned"
        down_feat.to_csv(out_dir / "gxc_downsample_features.tsv.gz", sep="\t", index=False)
        print(f"  downsampling: {len(down_rep)} replicate summaries, "
              f"{len(down_test)} tests, {len(down_feat)} features")
    elif args.downsample_reps > 0:
        print("  WARNING: downsampling requested but no downsample outputs were generated")

    print("[3/5] permutation calibration")
    perm = None
    if args.permutations > 0 and "perm_status" in merged.columns:
        perm_cols = KEYS + [
            c for c in [
                "status", "n", "p_GxC", "q_GxC", "permutation_sample", "perm_status",
                "permutations", "permutation_strata", "n_permutation_strata",
                "perm_observed_t", "perm_empirical_p", "perm_abs_t_q95", "perm_abs_t_q99",
            ] if c in merged.columns
        ]
        perm = merged.loc[merged["permutation_sample"].notna(), perm_cols].copy()
        perm.to_csv(out_dir / "gxc_permutation.tsv.gz", sep="\t", index=False)
        print(f"  permutation diagnostics: {len(perm)} selected tests")
    else:
        print("  skipped (--permutations 0)")

    print("[4/5] matched-feature ancestry comparison")
    matched, matched_summary = _matched_ancestry_tables(merged, annotations, "EAS", "EUR")
    if len(matched):
        matched.to_csv(out_dir / "gxc_matched_ancestry.tsv.gz", sep="\t", index=False)
        matched_summary.to_csv(out_dir / "gxc_matched_ancestry_summary.tsv", sep="\t", index=False)
        print(f"  matched ancestry: {len(matched)} cell-type tests across "
              f"{matched['phenotype_id'].nunique()} features")
    else:
        print("  no EAS/EUR matched feature-cell tests")

    print("[5/5] summaries + bundle")
    summary = _summary_table(merged, annotations, hc3=hc3, perm=perm)
    summary.to_csv(out_dir / "diagnostic_summary.tsv", sep="\t", index=False)

    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version,
        "command": " ".join(sys.argv),
        "results_dir": str(results_dir), "qtl_dir": str(qtl_dir),
        "aggregate_script": str(aggregate_script.resolve()),
        "ancestries": args.ancestries, "workers": args.workers,
        "min_n": args.min_n, "hc3": not args.skip_hc3,
        "permutations": args.permutations, "permute_top": args.permute_top,
        "permute_random": args.permute_random, "seed": args.seed,
        "downsample_ancestry": args.downsample_ancestry,
        "downsample_reference_ancestry": args.downsample_reference_ancestry,
        "downsample_n": args.downsample_n, "downsample_reps": args.downsample_reps,
        "n_input_interaction_rows": int(len(interactions)),
        "n_diagnostic_rows": int(len(diag)),
        "n_permutation_selected": int(len(permute_keys)),
        "n_downsample_replicate_rows": int(len(down_rep)),
        "n_downsample_test_rows": int(len(down_test)),
        "n_downsample_feature_rows": int(len(down_feat)),
        "n_matched_ancestry_rows": int(len(matched)),
    }
    (out_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    _write_readme(out_dir, args, None)
    bundle = _bundle_outputs(out_dir)
    done.write_text(datetime.now(timezone.utc).isoformat() + "\n")
    print(f"\nDiagnostics complete: {out_dir}")
    print(f"Upload for summarization: {bundle}")


if __name__ == "__main__":
    main()
