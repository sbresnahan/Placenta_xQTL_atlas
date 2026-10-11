#!/usr/bin/env python3
"""Aggregate all phenotype shards for one diagnostic modality.

This is the first point where downsampled BH/FDR is calculated: raw P values
from every phenotype shard in the modality are concatenated before correction,
so sharding does not change the production multiple-testing unit.
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
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


def sanitize(value: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in str(value))


def nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--plan", required=True)
    p.add_argument("--modality", required=True)
    p.add_argument("--downsample-ancestry", default="EAS")
    p.add_argument("--downsample-reference-ancestry", default="EUR")
    p.add_argument("--force", action="store_true")
    args = p.parse_args()

    out_dir = Path(args.out_dir).resolve()
    work = out_dir / "work"
    mod = args.modality
    mod_safe = sanitize(mod)
    plan = pd.read_csv(args.plan, sep="\t")
    plan = plan[plan["modality"].astype(str) == mod]
    if plan.empty:
        raise SystemExit(f"ERROR: no plan rows for modality {mod}")
    shard_dir = work / "shards" / mod_safe
    mod_dir = work / "modality" / mod_safe
    mod_dir.mkdir(parents=True, exist_ok=True)
    done = mod_dir / ".done"
    diag_out = mod_dir / "diag.tsv.gz"
    if done.exists() and nonempty(diag_out) and not args.force:
        print(f"SKIP aggregate {mod}: existing {diag_out}")
        return

    diag_frames = []
    down_payloads = []
    missing = []
    for row in plan.itertuples(index=False):
        anc = str(row.ancestry)
        shard = int(row.shard_index)
        stem = f"{anc}.shard{shard:04d}"
        diag_path = shard_dir / f"{stem}.diag.tsv.gz"
        shard_done = shard_dir / f"{stem}.done"
        if not shard_done.exists() or not nonempty(diag_path):
            missing.append(str(diag_path))
            continue
        diag_frames.append(pd.read_csv(diag_path, sep="\t"))
        down_path = shard_dir / f"{stem}.downsample.pkl.gz"
        if anc == args.downsample_ancestry and nonempty(down_path):
            down_payloads.append(pd.read_pickle(down_path, compression="gzip"))
    if missing:
        raise SystemExit("ERROR: missing shard outputs:\n  " + "\n  ".join(missing[:20]))
    diag = pd.concat(diag_frames, ignore_index=True) if diag_frames else pd.DataFrame()
    if diag.empty:
        raise SystemExit(f"ERROR: no diagnostic rows for {mod}")
    diag.to_csv(diag_out, sep="\t", index=False)
    print(f"aggregate {mod}: {len(diag)} diagnostic rows from {len(diag_frames)} shards")

    if down_payloads:
        core = load_core(Path(__file__).resolve().parent)
        stage49 = core._load_stage49(Path(__file__).with_name("49_aggregate_coloc_twas.py"))
        n_reps = int(down_payloads[0]["n_reps"])
        n_target = int(down_payloads[0]["target_n"])
        records = []
        for payload in down_payloads:
            if int(payload["n_reps"]) != n_reps or int(payload["target_n"]) != n_target:
                raise SystemExit(f"ERROR: inconsistent downsample settings across {mod} shards")
            meta = payload["meta"].reset_index(drop=True)
            pmat = np.asarray(payload["pvals"], dtype=float)
            bmat = np.asarray(payload["betas"], dtype=float)
            nmat = np.asarray(payload["ns"], dtype=float)
            if len(meta) != pmat.shape[0] or pmat.shape != bmat.shape or pmat.shape != nmat.shape:
                raise SystemExit(f"ERROR: malformed downsample payload for {mod}")
            for i, row in meta.iterrows():
                records.append({
                    "phenotype_id": str(row["phenotype_id"]),
                    "variant_id": str(row["variant_id"]),
                    "cell_type": str(row["cell_type"]),
                    "full_n": pd.to_numeric(row.get("n", np.nan), errors="coerce"),
                    "full_beta": pd.to_numeric(row.get("beta_GxC", np.nan), errors="coerce"),
                    "full_p": pd.to_numeric(row.get("p_GxC", np.nan), errors="coerce"),
                    "full_q": pd.to_numeric(row.get("q_GxC", np.nan), errors="coerce"),
                    "betas": bmat[i], "pvals": pmat[i], "ns": nmat[i],
                })
        rr, tt, ff = core._summarize_downsample_modality(
            records, stage49, args.downsample_ancestry, mod, n_target,
            args.downsample_reference_ancestry
        )
        rr.to_csv(mod_dir / "downsample_replicates.tsv", sep="\t", index=False)
        tt.to_csv(mod_dir / "downsample_tests.tsv.gz", sep="\t", index=False)
        ff.to_csv(mod_dir / "downsample_features.tsv.gz", sep="\t", index=False)
        print(f"aggregate {mod}: downsample {len(records)} tests x {n_reps} reps -> "
              f"{len(rr)} replicate rows")
    done.write_text("ok\n")


if __name__ == "__main__":
    main()
