#!/usr/bin/env python3
"""Finalize per-k HCP optimization jobs.

This helper combines two sources of completed grid points:
  1) legacy serial 25a/25b staging directories; and
  2) isolated per-k sandboxes created by 25a_submit_hcp_k_jobs.sh / 25b_submit_hcp_k_jobs.sh.

It selects k* by maximum eGene count (ties -> smaller k), writes the canonical
optimization TSV/PNG, and installs the winning HCP file into qtl_inputs.
"""

import argparse
import os
import shutil
import sys
from pathlib import Path

import pandas as pd


def valid_parquet(path: Path) -> bool:
    """Cheap corruption guard: Parquet files begin and end with PAR1."""
    try:
        if not path.is_file() or path.stat().st_size < 8:
            return False
        with path.open("rb") as fh:
            head = fh.read(4)
            fh.seek(-4, os.SEEK_END)
            tail = fh.read(4)
        return head == b"PAR1" and tail == b"PAR1"
    except OSError:
        return False


def count_egenes(parquet_path: Path, fdr: float) -> tuple[int, int]:
    df = pd.read_parquet(parquet_path)
    if "qval" not in df.columns:
        raise RuntimeError(f"no qval column in {parquet_path}")
    pheno = df["phenotype_id"] if "phenotype_id" in df.columns else df.index.to_series()
    hits = pheno[df["qval"] <= fdr]
    return int(hits.nunique()), int(df.shape[0])


def hcp_pruning_counts(pruning_path: Path, k: int) -> tuple[int, int]:
    if k == 0 or not pruning_path.is_file():
        return k, 0
    pruned = pd.read_csv(pruning_path, sep="\t")
    if "covariate" not in pruned.columns:
        return k, 0
    dropped = int(pruned["covariate"].astype(str).str.startswith("HCP_").sum())
    return k - dropped, dropped


def atomic_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + f".tmp.{os.getpid()}")
    shutil.copy2(src, tmp)
    os.replace(tmp, dst)


def atomic_tsv(df: pd.DataFrame, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + f".tmp.{os.getpid()}")
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, dst)


def sandbox_expression(qtl_dir: Path, anc: str, k: int, sandbox_suffix: str = ""):
    root = qtl_dir / f"hcp_optimization_per_k{sandbox_suffix}" / anc / f"k{k}"
    summary = root / "work" / f"{anc}_optimal_hcp.tsv"
    hcp = root / "qtl_inputs" / f"{anc}_hcp_factors_harmonized.tsv"
    if not (summary.is_file() and summary.stat().st_size > 0 and hcp.is_file() and hcp.stat().st_size > 0):
        return None
    df = pd.read_csv(summary, sep="\t")
    if "k" not in df.columns:
        return None
    row = df[df["k"].astype(int) == k]
    if row.empty:
        return None
    rec = row.iloc[0].to_dict()
    rec.update({"ancestry": anc, "k": k})
    return rec, hcp, "sandbox"


def legacy_expression(qtl_dir: Path, anc: str, k: int, fdr: float):
    staging = qtl_dir / "hcp_optimization" / anc
    hcp = staging / f"{anc}_hcp_k{k}_harmonized.tsv"
    parquet = staging / f"results_k{k}" / f"{anc}_expression_cisqtl.parquet"
    if not (hcp.is_file() and hcp.stat().st_size > 0 and valid_parquet(parquet)):
        return None
    n_eg, n_tested = count_egenes(parquet, fdr)
    n_used, n_dropped = hcp_pruning_counts(staging / f"{anc}_covariate_pruning_k{k}.tsv", k)
    rec = {
        "ancestry": anc,
        "k": k,
        "n_egenes": n_eg,
        "n_tested": n_tested,
        "n_hcp_used": n_used,
        "n_hcp_dropped": n_dropped,
    }
    return rec, hcp, "legacy"


def sandbox_modality(qtl_dir: Path, anc: str, mod: str, k: int, sandbox_suffix: str = ""):
    root = qtl_dir / f"hcp_optimization_modalities_per_k{sandbox_suffix}" / anc / mod / f"k{k}"
    summary = root / "work" / f"{anc}_{mod}_optimal_hcp.tsv"
    hcp = root / "qtl_inputs" / f"{anc}_{mod}_hcp_factors_optimized.tsv"
    if not (summary.is_file() and summary.stat().st_size > 0 and hcp.is_file() and hcp.stat().st_size > 0):
        return None
    df = pd.read_csv(summary, sep="\t")
    if "k" not in df.columns:
        return None
    row = df[df["k"].astype(int) == k]
    if row.empty:
        return None
    rec = row.iloc[0].to_dict()
    rec.update({"ancestry": anc, "modality": mod, "k": k})
    return rec, hcp, "sandbox"


def infer_scope(qtl_dir: Path, anc: str, mod: str, chr1_min: int) -> str:
    bed = qtl_dir / f"{anc}_{mod}_harmonized.bed"
    if not bed.is_file():
        raise RuntimeError(f"missing harmonized BED needed to infer mapping scope: {bed}")
    chrom = pd.read_csv(bed, sep="\t", usecols=["#chr"])["#chr"]
    n_chr1 = int((chrom == "chr1").sum())
    return "chr1" if n_chr1 >= chr1_min else "genome-wide"


def legacy_modality(qtl_dir: Path, anc: str, mod: str, k: int, fdr: float, chr1_min: int):
    staging = qtl_dir / "hcp_optimization_modalities" / anc / mod
    hcp = staging / f"{anc}_{mod}_hcp_k{k}.tsv"
    parquet = staging / f"results_k{k}" / f"{anc}_{mod}_cisqtl.parquet"
    if not (hcp.is_file() and hcp.stat().st_size > 0 and valid_parquet(parquet)):
        return None
    n_eg, n_tested = count_egenes(parquet, fdr)
    n_used, n_dropped = hcp_pruning_counts(staging / f"{anc}_covariate_pruning_k{k}.tsv", k)
    rec = {
        "ancestry": anc,
        "modality": mod,
        "k": k,
        "n_egenes": n_eg,
        "n_tested": n_tested,
        "mapping_scope": infer_scope(qtl_dir, anc, mod, chr1_min),
        "n_hcp_used": n_used,
        "n_hcp_dropped": n_dropped,
    }
    return rec, hcp, "legacy"


def load_plotters(scripts_dir: Path):
    sys.path.insert(0, str(scripts_dir))
    import optimize_hcp_chr1 as expression_driver  # noqa: E402
    import optimize_hcp_modalities as modality_driver  # noqa: E402
    return expression_driver, modality_driver


def finalize_expression(args, k_grid, expr_driver):
    qtl = Path(args.qtl_dir)
    rows = []
    hcp_by_k = {}
    source_by_k = {}
    for k in k_grid:
        item = sandbox_expression(qtl, args.ancestry, k, args.sandbox_suffix)
        if item is None and not args.sandbox_suffix:
            # legacy trees are lambda1=0.5 artifacts; skip for suffixed runs
            item = legacy_expression(qtl, args.ancestry, k, args.fdr)
        if item is None:
            raise RuntimeError(f"missing completed grid point: {args.ancestry} expression k={k}")
        rec, hcp, source = item
        rec.pop("chosen", None)
        rows.append(rec)
        hcp_by_k[k] = hcp
        source_by_k[k] = source

    cols = ["ancestry", "k", "n_egenes", "n_tested", "n_hcp_used", "n_hcp_dropped"]
    df = pd.DataFrame(rows)
    for c in cols:
        if c not in df.columns:
            raise RuntimeError(f"missing required column {c} while finalizing expression")
    df = df[cols].sort_values("k").reset_index(drop=True)
    best = int(df["n_egenes"].max())
    k_star = int(df.loc[df["n_egenes"] == best, "k"].min())
    df["chosen"] = df["k"] == k_star

    out_root = qtl / "hcp_optimization"
    out_tsv = out_root / f"{args.ancestry}_optimal_hcp.tsv"
    out_png = out_root / f"{args.ancestry}_optimal_hcp.png"
    atomic_tsv(df, out_tsv)
    expr_driver.plot_optimization(df, args.ancestry, k_star, args.fdr, str(out_png))

    canonical = qtl / f"{args.ancestry}_hcp_factors_harmonized.tsv"
    backup = qtl / f"{args.ancestry}_hcp_factors_harmonized.pre25a_backup.tsv"
    if canonical.is_file() and not backup.exists():
        atomic_copy(canonical, backup)
    atomic_copy(hcp_by_k[k_star], canonical)

    print(f"Finalized {args.ancestry} expression: k*={k_star}, eGenes={best}")
    print(f"  winner source: {source_by_k[k_star]}")
    print(f"  summary: {out_tsv}")
    print(f"  installed: {canonical}")


def finalize_modality(args, k_grid, mod_driver):
    qtl = Path(args.qtl_dir)
    rows = []
    hcp_by_k = {}
    source_by_k = {}
    for k in k_grid:
        item = sandbox_modality(qtl, args.ancestry, args.modality, k, args.sandbox_suffix)
        if item is None and not args.sandbox_suffix:
            # legacy trees are lambda1=0.5 artifacts; skip for suffixed runs
            item = legacy_modality(qtl, args.ancestry, args.modality, k, args.fdr, args.chr1_min)
        if item is None:
            raise RuntimeError(f"missing completed grid point: {args.ancestry} {args.modality} k={k}")
        rec, hcp, source = item
        rec.pop("chosen", None)
        rows.append(rec)
        hcp_by_k[k] = hcp
        source_by_k[k] = source

    cols = ["ancestry", "modality", "k", "n_egenes", "n_tested", "mapping_scope",
            "n_hcp_used", "n_hcp_dropped"]
    df = pd.DataFrame(rows)
    for c in cols:
        if c not in df.columns:
            raise RuntimeError(f"missing required column {c} while finalizing modality")
    df = df[cols].sort_values("k").reset_index(drop=True)
    scopes = sorted(df["mapping_scope"].astype(str).unique())
    if len(scopes) != 1:
        raise RuntimeError(f"inconsistent mapping scopes across k for {args.ancestry} {args.modality}: {scopes}")
    scope = scopes[0]
    best = int(df["n_egenes"].max())
    k_star = int(df.loc[df["n_egenes"] == best, "k"].min())
    df["chosen"] = df["k"] == k_star

    out_root = qtl / "hcp_optimization_modalities"
    out_tsv = out_root / f"{args.ancestry}_{args.modality}_optimal_hcp.tsv"
    out_png = out_root / f"{args.ancestry}_{args.modality}_optimal_hcp.png"
    atomic_tsv(df, out_tsv)
    mod_driver.plot_optimization(df, args.ancestry, args.modality, k_star,
                                 args.fdr, scope, str(out_png))

    canonical = qtl / f"{args.ancestry}_{args.modality}_hcp_factors_optimized.tsv"
    atomic_copy(hcp_by_k[k_star], canonical)

    print(f"Finalized {args.ancestry} {args.modality}: k*={k_star}, eGenes={best}, scope={scope}")
    print(f"  winner source: {source_by_k[k_star]}")
    print(f"  summary: {out_tsv}")
    print(f"  installed: {canonical}")


def main():
    p = argparse.ArgumentParser(description="Finalize isolated per-k HCP optimization jobs")
    p.add_argument("--mode", choices=["expression", "modality"], required=True)
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--scripts-dir", required=True)
    p.add_argument("--ancestry", required=True)
    p.add_argument("--modality")
    p.add_argument("--k-grid", default="0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--chr1-min", type=int, default=300)
    p.add_argument("--sandbox-suffix", default="",
                   help="Suffix for the per-k sandbox tree dir names "
                        "(e.g. '_lam5' reads hcp_optimization_per_k_lam5/). "
                        "Used by non-default-lambda1 sensitivity/production runs.")
    args = p.parse_args()

    if args.mode == "modality" and not args.modality:
        p.error("--modality is required with --mode modality")
    k_grid = [int(x) for x in args.k_grid.split()]
    if not k_grid:
        p.error("--k-grid is empty")
    if len(k_grid) != len(set(k_grid)):
        p.error("--k-grid contains duplicates")

    expr_driver, mod_driver = load_plotters(Path(args.scripts_dir))
    if args.mode == "expression":
        finalize_expression(args, k_grid, expr_driver)
    else:
        finalize_modality(args, k_grid, mod_driver)


if __name__ == "__main__":
    main()
