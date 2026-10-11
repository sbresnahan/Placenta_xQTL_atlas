#!/usr/bin/env python3
"""Run one phenotype shard of Module-06 GxC diagnostics.

The shard consumes caches created by 50a_prepare_celltype_diagnostics.py.  It
writes full-N support/HC3/permutation diagnostics and, for the downsampled
ancestry, raw replicate matrices.  BH/FDR is intentionally deferred until all
shards of a modality are aggregated.
"""
from __future__ import annotations

import argparse
import importlib.util
import multiprocessing as mp
import os
import sys
import time
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


def sanitize(value: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in str(value))


def nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def load_plan_row(plan_path: Path, array_index: int):
    plan = pd.read_csv(plan_path, sep="\t")
    hit = plan[plan["array_index"] == int(array_index)]
    if len(hit) != 1:
        raise SystemExit(f"ERROR: array index {array_index} maps to {len(hit)} rows in {plan_path}")
    return hit.iloc[0]


def run_full_diagnostics(core, sub, phen, dosages, deconv_raw, deconv_t,
                         cohort, cov, stage49, anc, mod, min_n, do_hc3,
                         permute_keys, n_perm, seed, workers):
    tests_by_feature = {}
    for (pid, vid), ss in sub.groupby(["phenotype_id", "variant_id"], sort=False):
        tests_by_feature[(str(pid), str(vid))] = list(ss["cell_type"].astype(str))
    tasks = list(tests_by_feature.keys())
    core._DIAG_STATE = {
        "ancestry": anc, "modality": mod,
        "phen": phen, "dosages": dosages,
        "deconv_raw": deconv_raw, "deconv_t": deconv_t,
        "cohort": cohort, "cov": cov,
        "tests_by_feature": tests_by_feature,
        "drop_target": stage49._drop_target_covariate,
        "min_n": min_n, "do_hc3": do_hc3,
        "permute_keys": permute_keys, "n_perm": n_perm,
        "seed": seed,
    }
    rows = []
    workers = max(1, int(workers))
    if workers > 1 and sys.platform.startswith("linux") and len(tasks) > 1:
        ctx = mp.get_context("fork")
        nproc = min(workers, len(tasks))
        chunksize = max(1, len(tasks) // (nproc * 4))
        with ctx.Pool(nproc) as pool:
            for done, rr in enumerate(pool.imap_unordered(core._diagnose_feature, tasks,
                                                          chunksize=chunksize), 1):
                rows.extend(rr)
                if done % 100 == 0 or done == len(tasks):
                    print(f"    full-N {done}/{len(tasks)} features")
    else:
        for done, task in enumerate(tasks, 1):
            rows.extend(core._diagnose_feature(task))
            if done % 100 == 0 or done == len(tasks):
                print(f"    full-N {done}/{len(tasks)} features")
    core._DIAG_STATE = None
    out = pd.DataFrame(rows)
    if len(out):
        out["ancestry"] = anc
        out["modality"] = mod
    return out


def run_downsample_raw(core, sub, phen, dosages, deconv_t, cov, stage49,
                       anc, mod, n_target, n_reps, seed, min_n, workers):
    """Exact batched downsampling for one phenotype shard.

    The sample subsets are deterministic from ancestry/modality/replicate only,
    so every phenotype shard uses the same N=target subjects for a replicate.
    Raw P/beta/N matrices are returned; modality-wide BH is done later.
    """
    if n_reps <= 0 or n_target is None or sub.empty:
        return None
    base_ids = sorted(set(phen.index.astype(str)) & set(dosages.index.astype(str)) &
                      set(deconv_t.index.astype(str)) & set(cov.index.astype(str)))
    if len(base_ids) < n_target:
        print(f"    WARNING: eligible N={len(base_ids)} < target N={n_target}; no downsampling")
        return None

    feature_table = (sub[["phenotype_id", "variant_id"]].drop_duplicates()
                     .astype(str).sort_values(["phenotype_id", "variant_id"])
                     .reset_index(drop=True))
    feat_map = {(r.phenotype_id, r.variant_id): i
                for i, r in feature_table.iterrows()}
    pids = feature_table["phenotype_id"].tolist()
    vids = feature_table["variant_id"].tolist()
    Y = phen.reindex(base_ids)[pids].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    G = dosages.reindex(base_ids)[vids].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    sub2 = (sub.copy().assign(
        phenotype_id=lambda d: d["phenotype_id"].astype(str),
        variant_id=lambda d: d["variant_id"].astype(str),
        cell_type=lambda d: d["cell_type"].astype(str),
    ).sort_values(["phenotype_id", "variant_id", "cell_type"]).reset_index(drop=True))
    ct_groups, cell_by_ct, cov_by_ct = {}, {}, {}
    for ct, ss in sub2.groupby("cell_type", sort=False):
        if ct not in deconv_t.columns:
            continue
        test_idx = ss.index.to_numpy(dtype=int)
        feat_idx = np.asarray([feat_map[(p, v)] for p, v in
                               zip(ss["phenotype_id"], ss["variant_id"])], dtype=int)
        ct_groups[ct] = {"test_idx": test_idx, "feat_idx": feat_idx}
        cell_by_ct[ct] = pd.to_numeric(deconv_t[ct].reindex(base_ids), errors="coerce").to_numpy(dtype=float)
        ctest = stage49._drop_target_covariate(cov, ct).reindex(base_ids)
        cov_by_ct[ct] = ctest.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    sample_indices = core._make_downsample_indices(
        len(base_ids), n_target, n_reps, seed, ("downsample", anc, mod, n_target)
    )
    if len(sample_indices) == 0:
        return None
    n_tests = len(sub2)
    pmat = np.full((n_tests, n_reps), np.nan, dtype=float)
    bmat = np.full((n_tests, n_reps), np.nan, dtype=float)
    nmat = np.full((n_tests, n_reps), np.nan, dtype=np.float32)
    core._DOWNSAMPLE_STATE = {
        "Y": Y, "G": G,
        "cell_by_ct": cell_by_ct, "cov_by_ct": cov_by_ct,
        "ct_groups": ct_groups, "sample_indices": sample_indices,
        "n_tests": n_tests, "min_n": min_n,
    }
    workers = max(1, int(workers))
    print(f"    downsample raw: {n_tests} tests x {n_reps} reps; "
          f"{len(ct_groups)} cell-type batches/rep; {min(workers,n_reps)} worker(s)")
    if workers > 1 and sys.platform.startswith("linux") and n_reps > 1:
        ctx = mp.get_context("fork")
        nproc = min(workers, n_reps)
        chunksize = max(1, n_reps // (nproc * 4))
        with ctx.Pool(nproc) as pool:
            for done, (r, b, pv, nn) in enumerate(
                pool.imap_unordered(core._downsample_replicate, range(n_reps),
                                    chunksize=chunksize), 1
            ):
                bmat[:, r] = b; pmat[:, r] = pv; nmat[:, r] = nn
                if done % max(1, n_reps // 4) == 0 or done == n_reps:
                    print(f"      downsample reps {done}/{n_reps}")
    else:
        for done, r in enumerate(range(n_reps), 1):
            rr, b, pv, nn = core._downsample_replicate(r)
            bmat[:, rr] = b; pmat[:, rr] = pv; nmat[:, rr] = nn
            if done % max(1, n_reps // 4) == 0 or done == n_reps:
                print(f"      downsample reps {done}/{n_reps}")
    core._DOWNSAMPLE_STATE = None
    meta_cols = [c for c in ["phenotype_id", "variant_id", "cell_type",
                             "n", "beta_GxC", "p_GxC", "q_GxC"] if c in sub2.columns]
    meta = sub2[meta_cols].copy()
    return {
        "ancestry": anc, "modality": mod, "target_n": int(n_target),
        "n_reps": int(n_reps), "base_n": int(len(base_ids)),
        "meta": meta, "betas": bmat, "pvals": pmat, "ns": nmat,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--plan", required=True)
    p.add_argument("--array-index", type=int, required=True)
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--min-n", type=int, default=20)
    p.add_argument("--permutations", type=int, default=100)
    p.add_argument("--permute-top", type=int, default=20)
    p.add_argument("--permute-random", type=int, default=40)
    p.add_argument("--downsample-reps", type=int, default=100)
    p.add_argument("--downsample-n", type=int, default=None)
    p.add_argument("--downsample-ancestry", default="EAS")
    p.add_argument("--downsample-reference-ancestry", default="EUR")
    p.add_argument("--seed", type=int, default=20261010)
    p.add_argument("--skip-hc3", action="store_true")
    p.add_argument("--force", action="store_true")
    args = p.parse_args()

    script_dir = Path(__file__).resolve().parent
    core = load_core(script_dir)
    args.workers = core._auto_workers(args.workers)
    results_dir = Path(args.results_dir).resolve()
    qtl_dir = Path(args.qtl_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    work_dir = out_dir / "work"
    row = load_plan_row(Path(args.plan), args.array_index)
    anc, mod = str(row["ancestry"]), str(row["modality"])
    shard = int(row["shard_index"])
    feature_file = Path(str(row["feature_file"]))
    shard_dir = work_dir / "shards" / sanitize(mod)
    shard_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{anc}.shard{shard:04d}"
    diag_path = shard_dir / f"{stem}.diag.tsv.gz"
    down_path = shard_dir / f"{stem}.downsample.pkl.gz"
    done_path = shard_dir / f"{stem}.done"
    if done_path.exists() and nonempty(diag_path) and not args.force:
        print(f"SKIP completed shard: {anc}/{mod} shard {shard}")
        return

    if args.downsample_reps > 0 and args.downsample_n is None:
        args.downsample_n = len(core._psam_samples(qtl_dir / f"{args.downsample_reference_ancestry}_qtl"))

    inter_path = results_dir / "coloc" / "aggregated" / "xqtl_celltype_interactions.tsv.gz"
    interactions = pd.read_csv(inter_path, sep="\t")
    for c in ["ancestry", "modality", "phenotype_id", "variant_id", "cell_type"]:
        interactions[c] = interactions[c].astype(str)
    pairs = pd.read_csv(feature_file, sep="\t", dtype=str)
    base = interactions[(interactions["ancestry"] == anc) & (interactions["modality"] == mod)]
    sub = base.merge(pairs, on=["phenotype_id", "variant_id"], how="inner")
    sub = sub.sort_values(["phenotype_id", "variant_id", "cell_type"]).reset_index(drop=True)
    if sub.empty:
        raise SystemExit(f"ERROR: empty shard {anc}/{mod} {shard}")

    cache = work_dir / "cache"
    phen_all = pd.read_pickle(cache / f"{anc}.{sanitize(mod)}.phen.pkl.gz", compression="gzip")
    dose_all = pd.read_pickle(cache / f"{anc}.dosages.pkl.gz", compression="gzip")
    cov = pd.read_pickle(cache / f"{anc}.{sanitize(mod)}.cov.pkl.gz", compression="gzip")
    deconv_raw = pd.read_pickle(cache / f"{anc}.deconv_raw.pkl.gz", compression="gzip")
    deconv_t = pd.read_pickle(cache / f"{anc}.deconv_t.pkl.gz", compression="gzip")
    cohort = pd.read_pickle(cache / f"{anc}.cohort.pkl.gz", compression="gzip")
    if len(cohort) == 0:
        cohort = None
    pids = pairs["phenotype_id"].astype(str).tolist()
    vids = pairs["variant_id"].astype(str).tolist()
    phen = phen_all[[x for x in pids if x in phen_all.columns]].copy()
    dose_cols = list(dict.fromkeys(x for x in vids if x in dose_all.columns))
    dosages = dose_all[dose_cols].copy()
    stage49 = core._load_stage49(script_dir / "49_aggregate_coloc_twas.py")

    permute_keys = core._select_permutation_tests(
        interactions, args.permute_top, args.permute_random, args.seed
    ) if args.permutations > 0 else {}
    t0 = time.time()
    print(f"shard {args.array_index}: {anc}/{mod} shard {shard}/{int(row['n_shards_ancestry'])}; "
          f"{len(pairs)} features, {len(sub)} tests, {args.workers} CPU(s)")
    diag = run_full_diagnostics(
        core, sub, phen, dosages, deconv_raw, deconv_t, cohort, cov, stage49,
        anc, mod, args.min_n, not args.skip_hc3, permute_keys,
        args.permutations, args.seed, args.workers
    )
    diag.to_csv(diag_path, sep="\t", index=False)

    if anc == args.downsample_ancestry and args.downsample_reps > 0:
        payload = run_downsample_raw(
            core, sub, phen, dosages, deconv_t, cov, stage49,
            anc, mod, args.downsample_n, args.downsample_reps,
            args.seed, args.min_n, args.workers
        )
        if payload is not None:
            pd.to_pickle(payload, down_path, compression="gzip")
    done_path.write_text(f"completed {time.time()-t0:.1f}s\n")
    print(f"completed shard in {time.time()-t0:.1f}s: {diag_path}")


if __name__ == "__main__":
    main()
