#!/usr/bin/env python3
"""Finalize sharded Module-06 cell-type diagnostics into reportable outputs."""
from __future__ import annotations

import argparse
import importlib.util
import json
import platform
import socket
import sys
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
    spec.loader.exec_module(mod)
    return mod


def nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def sanitize(value: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in str(value))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--plan", required=True)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--min-n", type=int, default=20)
    p.add_argument("--skip-hc3", action="store_true")
    p.add_argument("--permutations", type=int, default=100)
    p.add_argument("--permute-top", type=int, default=20)
    p.add_argument("--permute-random", type=int, default=40)
    p.add_argument("--downsample-reps", type=int, default=100)
    p.add_argument("--downsample-n", type=int, default=None)
    p.add_argument("--downsample-ancestry", default="EAS")
    p.add_argument("--downsample-reference-ancestry", default="EUR")
    p.add_argument("--seed", type=int, default=20261010)
    p.add_argument("--force", action="store_true")
    args = p.parse_args()

    script_dir = Path(__file__).resolve().parent
    core = load_core(script_dir)
    results_dir = Path(args.results_dir).resolve()
    qtl_dir = Path(args.qtl_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    bundle = out_dir / "celltype_diagnostics_bundle.zip"
    done = out_dir / ".done"
    if done.exists() and nonempty(bundle) and not args.force:
        print(f"SKIP finalization: {bundle}")
        return
    if args.downsample_reps > 0 and args.downsample_n is None:
        args.downsample_n = len(core._psam_samples(qtl_dir / f"{args.downsample_reference_ancestry}_qtl"))

    plan = pd.read_csv(args.plan, sep="\t")
    modalities = sorted(plan["modality"].astype(str).unique())
    diag_frames = []
    down_rep_frames, down_test_frames, down_feat_frames = [], [], []
    for mod in modalities:
        mdir = out_dir / "work" / "modality" / sanitize(mod)
        dpath = mdir / "diag.tsv.gz"
        if not (mdir / ".done").exists() or not nonempty(dpath):
            raise SystemExit(f"ERROR: modality aggregate incomplete: {mod}")
        diag_frames.append(pd.read_csv(dpath, sep="\t"))
        for name, frames in [
            ("downsample_replicates.tsv", down_rep_frames),
            ("downsample_tests.tsv.gz", down_test_frames),
            ("downsample_features.tsv.gz", down_feat_frames),
        ]:
            path = mdir / name
            if nonempty(path):
                frames.append(pd.read_csv(path, sep="\t"))
    diag = pd.concat(diag_frames, ignore_index=True)

    agg_dir = results_dir / "coloc" / "aggregated"
    interactions = pd.read_csv(agg_dir / "xqtl_celltype_interactions.tsv.gz", sep="\t")
    annotations = pd.read_csv(agg_dir / "xqtl_celltype_annotation.tsv.gz", sep="\t")
    for col in ("ancestry", "modality", "phenotype_id", "variant_id", "cell_type"):
        interactions[col] = interactions[col].astype(str)
        if col in diag.columns:
            diag[col] = diag[col].astype(str)
    interactions["p_GxC"] = pd.to_numeric(interactions.get("p_GxC"), errors="coerce")
    interactions["q_GxC"] = pd.to_numeric(interactions.get("q_GxC"), errors="coerce")
    merged = interactions.merge(diag, on=core.KEYS, how="left", validate="one_to_one")
    if "recalc_classical_p_GxC" in merged.columns:
        merged["classical_p_abs_diff"] = (
            pd.to_numeric(merged["p_GxC"], errors="coerce") -
            pd.to_numeric(merged["recalc_classical_p_GxC"], errors="coerce")
        ).abs()

    support_cols = core.KEYS + [c for c in [
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
    ] if c in merged.columns]
    merged[support_cols].to_csv(out_dir / "gxc_support.tsv.gz", sep="\t", index=False)

    stage49 = core._load_stage49(script_dir / "49_aggregate_coloc_twas.py")
    hc3 = None
    if not args.skip_hc3 and "hc3_p_GxC" in merged.columns:
        hc3_cols = core.KEYS + [c for c in [
            "status", "n", "beta_GxC", "se_GxC", "p_GxC", "q_GxC", "diag_status",
            "recalc_beta_GxC", "recalc_classical_se_GxC", "recalc_classical_p_GxC",
            "classical_p_abs_diff", "hc3_se_GxC", "hc3_t_GxC", "hc3_p_GxC",
            "hc3_se_ratio", "max_leverage", "median_leverage", "design_condition_number",
        ] if c in merged.columns]
        hc3 = core._add_hc3_qvalues(merged[hc3_cols], stage49)
        hc3.to_csv(out_dir / "gxc_hc3.tsv.gz", sep="\t", index=False)

    down_rep = pd.concat(down_rep_frames, ignore_index=True) if down_rep_frames else pd.DataFrame()
    down_test = pd.concat(down_test_frames, ignore_index=True) if down_test_frames else pd.DataFrame()
    down_feat = pd.concat(down_feat_frames, ignore_index=True) if down_feat_frames else pd.DataFrame()
    if len(down_rep):
        down_rep.to_csv(out_dir / "gxc_downsample_replicates.tsv", sep="\t", index=False)
        down_test.to_csv(out_dir / "gxc_downsample_tests.tsv.gz", sep="\t", index=False)
        src_ann = annotations[annotations["ancestry"].astype(str) == args.downsample_ancestry][
            ["modality", "phenotype_id", "primary_cell_type"]
        ].copy()
        src_ann["modality"] = src_ann["modality"].astype(str)
        src_ann["phenotype_id"] = src_ann["phenotype_id"].astype(str)
        src_ann = src_ann.rename(columns={"primary_cell_type": "full_primary_cell_type"})
        down_feat["modality"] = down_feat["modality"].astype(str)
        down_feat["phenotype_id"] = down_feat["phenotype_id"].astype(str)
        down_feat = down_feat.merge(src_ann, on=["modality", "phenotype_id"], how="left")
        down_feat["full_assigned"] = down_feat["full_primary_cell_type"].astype(str) != "unassigned"
        down_feat.to_csv(out_dir / "gxc_downsample_features.tsv.gz", sep="\t", index=False)

    perm = None
    if args.permutations > 0 and "perm_status" in merged.columns:
        perm_cols = core.KEYS + [c for c in [
            "status", "n", "p_GxC", "q_GxC", "permutation_sample", "perm_status",
            "permutations", "permutation_strata", "n_permutation_strata",
            "perm_observed_t", "perm_empirical_p", "perm_abs_t_q95", "perm_abs_t_q99",
        ] if c in merged.columns]
        mask = merged.get("permutation_sample", pd.Series(index=merged.index, dtype=object)).notna()
        perm = merged.loc[mask, perm_cols].copy()
        perm.to_csv(out_dir / "gxc_permutation.tsv.gz", sep="\t", index=False)

    matched, matched_summary = core._matched_ancestry_tables(merged, annotations, "EAS", "EUR")
    if len(matched):
        matched.to_csv(out_dir / "gxc_matched_ancestry.tsv.gz", sep="\t", index=False)
        matched_summary.to_csv(out_dir / "gxc_matched_ancestry_summary.tsv", sep="\t", index=False)
    summary = core._summary_table(merged, annotations, hc3=hc3, perm=perm)
    summary.to_csv(out_dir / "diagnostic_summary.tsv", sep="\t", index=False)

    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": socket.gethostname(), "platform": platform.platform(),
        "python": sys.version, "command": " ".join(sys.argv),
        "execution_engine": "lsf_modality_arrays_phenotype_shards_v4",
        "results_dir": str(results_dir), "qtl_dir": str(qtl_dir),
        "workers_per_shard": args.workers, "min_n": args.min_n,
        "hc3": not args.skip_hc3, "permutations": args.permutations,
        "permute_top": args.permute_top, "permute_random": args.permute_random,
        "seed": args.seed, "downsample_ancestry": args.downsample_ancestry,
        "downsample_reference_ancestry": args.downsample_reference_ancestry,
        "downsample_n": args.downsample_n, "downsample_reps": args.downsample_reps,
        "n_plan_tasks": int(len(plan)), "n_modalities": int(len(modalities)),
        "n_input_interaction_rows": int(len(interactions)),
        "n_diagnostic_rows": int(len(diag)),
        "n_downsample_replicate_rows": int(len(down_rep)),
        "n_downsample_test_rows": int(len(down_test)),
        "n_downsample_feature_rows": int(len(down_feat)),
        "n_matched_ancestry_rows": int(len(matched)),
    }
    (out_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    core._write_readme(out_dir, args, None)
    core._bundle_outputs(out_dir)
    done.write_text(datetime.now(timezone.utc).isoformat() + "\n")
    print(f"Diagnostics complete: {out_dir}")
    print(f"Upload: {bundle}")


if __name__ == "__main__":
    main()
