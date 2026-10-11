#!/usr/bin/env python3
"""Prepare sharded Module-06 cell-type diagnostics.

Creates deterministic per-modality shard plans and, in build-cache mode,
materializes reusable phenotype/dosage/covariate/deconvolution caches so every
LSF shard does not rescan BEDs or rerun PLINK exports.
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import json
import math
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def load_core(script_dir: Path):
    path = script_dir / "50_celltype_diagnostics.py"
    spec = importlib.util.spec_from_file_location("celltype_diag_core", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"ERROR: cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    # multiprocessing pickles worker callables by fully qualified module name.
    # Dynamic modules are not registered automatically by module_from_spec(),
    # so Pool.imap_unordered() cannot pickle core._diagnose_feature /
    # core._downsample_replicate unless the module is present in sys.modules.
    # Register before exec so forked workers inherit an import-resolvable module.
    sys.modules[spec.name] = mod
    try:
        spec.loader.exec_module(mod)
    except Exception:
        sys.modules.pop(spec.name, None)
        raise
    return mod


def nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def sanitize(value: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in str(value))


def build_plan(interactions: pd.DataFrame, ancestries, work_dir: Path,
               shard_features: int, force: bool):
    plan_dir = work_dir / "plans"
    feature_dir = work_dir / "features"
    plan_dir.mkdir(parents=True, exist_ok=True)
    feature_dir.mkdir(parents=True, exist_ok=True)
    rows_all = []
    modalities = sorted(interactions.loc[
        interactions["ancestry"].astype(str).isin(ancestries), "modality"
    ].dropna().astype(str).unique())
    for mod in modalities:
        mod_rows = []
        array_index = 0
        for anc in ancestries:
            sub = interactions[(interactions["ancestry"].astype(str) == anc) &
                               (interactions["modality"].astype(str) == mod)]
            feats = (sub[["phenotype_id", "variant_id"]]
                     .drop_duplicates()
                     .astype(str)
                     .sort_values(["phenotype_id", "variant_id"])
                     .reset_index(drop=True))
            if feats.empty:
                continue
            n_shards = int(math.ceil(len(feats) / shard_features))
            for shard0 in range(n_shards):
                array_index += 1
                lo = shard0 * shard_features
                hi = min(len(feats), (shard0 + 1) * shard_features)
                chunk = feats.iloc[lo:hi].copy()
                feat_path = feature_dir / f"{sanitize(mod)}.{anc}.shard{shard0+1:04d}.tsv"
                if force or not nonempty(feat_path):
                    chunk.to_csv(feat_path, sep="\t", index=False)
                rec = {
                    "array_index": array_index,
                    "ancestry": anc,
                    "modality": mod,
                    "shard_index": shard0 + 1,
                    "n_shards_ancestry": n_shards,
                    "n_features": len(chunk),
                    "feature_file": str(feat_path.resolve()),
                }
                mod_rows.append(rec)
                rows_all.append(rec)
        if mod_rows:
            pd.DataFrame(mod_rows).to_csv(plan_dir / f"{sanitize(mod)}.tsv",
                                          sep="\t", index=False)
    plan = pd.DataFrame(rows_all)
    plan.to_csv(work_dir / "plan_all.tsv", sep="\t", index=False)
    return plan


def build_cache(core, interactions, ancestries, results_dir: Path, qtl_dir: Path,
                out_dir: Path, force: bool, plink2: str):
    cache_dir = out_dir / "work" / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    stage49 = core._load_stage49(Path(__file__).with_name("49_aggregate_coloc_twas.py"))

    # QC tables are cheap and global; write them once here.
    ann_path = results_dir / "coloc" / "aggregated" / "xqtl_celltype_annotation.tsv.gz"
    annotations = pd.read_csv(ann_path, sep="\t") if ann_path.exists() else pd.DataFrame()
    sample_summary, sample_detail, cell_dist = core._sample_overlap_tables(
        qtl_dir, interactions, ancestries, stage49
    )
    sample_summary.to_csv(out_dir / "sample_overlap_summary.tsv", sep="\t", index=False)
    sample_detail.to_csv(out_dir / "sample_overlap_by_sample.tsv.gz", sep="\t", index=False)
    cell_dist.to_csv(out_dir / "celltype_distribution.tsv", sep="\t", index=False)

    with tempfile.TemporaryDirectory(prefix="celltype_diag_prepare_", dir=str(out_dir)) as tmp:
        for anc in ancestries:
            anc_i = interactions[interactions["ancestry"].astype(str) == anc].copy()
            if anc_i.empty:
                continue
            deconv_path = qtl_dir / f"{anc}_deconvolution_harmonized.tsv"
            deconv_raw, cohort = core._load_deconv_full(deconv_path)
            targets = sorted(set(anc_i["cell_type"].dropna().astype(str)))
            targets = [c for c in targets if c in deconv_raw.columns]
            deconv_t = np.arcsinh(deconv_raw[targets].apply(pd.to_numeric, errors="coerce"))
            deconv_t = deconv_t - deconv_t.mean(axis=0)
            deconv_raw.to_pickle(cache_dir / f"{anc}.deconv_raw.pkl.gz", compression="gzip")
            deconv_t.to_pickle(cache_dir / f"{anc}.deconv_t.pkl.gz", compression="gzip")
            cohort_cache = (cohort if cohort is not None else pd.Series(dtype=object, name="cohort"))
            cohort_cache.to_pickle(cache_dir / f"{anc}.cohort.pkl.gz", compression="gzip")

            dose_path = cache_dir / f"{anc}.dosages.pkl.gz"
            if force or not nonempty(dose_path):
                variants = anc_i["variant_id"].dropna().astype(str).unique()
                print(f"prepare: {anc}: exporting {len(variants)} lead variants once")
                dosages = stage49._export_lead_dosages(qtl_dir, anc, variants, plink2, Path(tmp))
                dosages.to_pickle(dose_path, compression="gzip")
            else:
                print(f"prepare: {anc}: reusing {dose_path.name}")

            for mod in sorted(anc_i["modality"].dropna().astype(str).unique()):
                sub = anc_i[anc_i["modality"].astype(str) == mod]
                phen_path = cache_dir / f"{anc}.{sanitize(mod)}.phen.pkl.gz"
                cov_path = cache_dir / f"{anc}.{sanitize(mod)}.cov.pkl.gz"
                if force or not nonempty(phen_path):
                    bed = qtl_dir / f"{anc}_{mod}.bed.gz"
                    wanted = set(sub["phenotype_id"].dropna().astype(str))
                    print(f"prepare: {anc}/{mod}: extracting {len(wanted)} phenotypes once")
                    phen = stage49._load_selected_bed(bed, wanted)
                    phen.to_pickle(phen_path, compression="gzip")
                if force or not nonempty(cov_path):
                    cov, source = stage49._load_covariates(qtl_dir, anc, mod)
                    if cov is None:
                        raise SystemExit(f"ERROR: missing covariates for {anc}/{mod}: {source}")
                    cov.to_pickle(cov_path, compression="gzip")

    meta = {
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "ancestries": list(ancestries),
        "n_interaction_rows": int(len(interactions)),
        "n_annotation_rows": int(len(annotations)),
    }
    (out_dir / "work" / "prepare_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    (out_dir / "work" / ".prepare_done").write_text(meta["prepared_utc"] + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    p.add_argument("--shard-features", type=int, default=100)
    p.add_argument("--plink2", default="plink2")
    p.add_argument("--plan-only", action="store_true")
    p.add_argument("--force", action="store_true")
    args = p.parse_args()
    if args.shard_features < 1:
        p.error("--shard-features must be >=1")

    results_dir = Path(args.results_dir).resolve()
    qtl_dir = Path(args.qtl_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    work_dir = out_dir / "work"
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    inter_path = results_dir / "coloc" / "aggregated" / "xqtl_celltype_interactions.tsv.gz"
    if not inter_path.exists():
        raise SystemExit(f"ERROR: missing {inter_path}")
    interactions = pd.read_csv(inter_path, sep="\t")
    need = {"ancestry", "modality", "phenotype_id", "variant_id", "cell_type"}
    if not need <= set(interactions.columns):
        raise SystemExit(f"ERROR: interactions missing columns {sorted(need-set(interactions.columns))}")
    interactions["ancestry"] = interactions["ancestry"].astype(str)
    interactions["modality"] = interactions["modality"].astype(str)
    interactions["phenotype_id"] = interactions["phenotype_id"].astype(str)
    interactions["variant_id"] = interactions["variant_id"].astype(str)

    plan = build_plan(interactions, args.ancestries, work_dir,
                      args.shard_features, args.force)
    print(f"plan: {len(plan)} shard task(s) across {plan['modality'].nunique() if len(plan) else 0} modalities")
    if len(plan):
        for mod, sub in plan.groupby("modality", sort=True):
            print(f"PLAN\t{mod}\t{len(sub)}")
    if args.plan_only:
        return
    core = load_core(Path(__file__).resolve().parent)
    build_cache(core, interactions, args.ancestries, results_dir, qtl_dir,
                out_dir, args.force, args.plink2)
    print(f"prepare complete: {out_dir / 'work'}")


if __name__ == "__main__":
    main()
