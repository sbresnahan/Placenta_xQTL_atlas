#!/usr/bin/env python3
"""
check_coloc_completeness.py

Audit completeness of Module 06 SuSiE-coloc (stages 40/41/42) and colocBoost
(stages 43/44/45) from the manifests that define the expected work.

The checker is intentionally stdlib-only.

Completeness semantics
----------------------
SuSiE-coloc:
  Expected units are rows of:
    {COLOC_DIR}/loci/{modality}.tasks.tsv
  A task is normally terminal when its per-task .done sentinel exists.
  Successful tasks also have both .coloc.tsv and .variants.tsv. Terminal
  scientific/numerical outcomes (e.g. no credible set, too few variants,
  SuSiE/coloc failure) may legitimately have a .done sentinel and no result
  TSVs. A partial result pair is always treated as needing rerun.

colocBoost:
  Expected units are rows of:
    {COLOC_DIR}/colocboost/{ancestry}.regions.tsv
  A region is normally terminal when its .done sentinel exists. Successful
  regions also have .clusters.tsv; .vcp.tsv is optional. Terminal outcomes
  such as too few variants/outcomes or colocBoost failure may legitimately
  have .done but no clusters file.

Diagnostics are used for status/message annotation, but .done sentinels and
output consistency determine whether an expected unit needs rerun. This makes
the checker robust to reruns where diagnostics say "ok / already done".

Outputs
-------
  <out-dir>/summary.tsv
  <out-dir>/susie_tasks.tsv
  <out-dir>/colocboost_regions.tsv
  <out-dir>/recovery_shards.tsv

The recovery_shards.tsv file reports array indices containing at least one unit
that needs rerun, using the same round-robin assignment as 42/45. The shard
sizes must match the submission run (defaults: SuSiE 25, colocBoost 10).

Exit status is 0 by default even when incomplete. Use --fail-on-incomplete to
exit 2 when any expected unit needs rerun or a required manifest is missing.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


DEFAULT_MODALITIES = [
    "expression",
    "isoforms",
    "isoform_expression",
    "splicing",
    "intron_retention",
    "alt_TSS",
    "alt_polyA",
    "RNA_editing",
    "stability",
]

DEFAULT_ANCESTRIES = ["EAS", "EUR"]


def sanitize_id(x: str) -> str:
    return re.sub(r"[^A-Za-z0-9._+-]", "_", str(x))


def read_tsv(path: Path) -> List[dict]:
    with path.open("r", newline="", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: Sequence[dict], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=list(fieldnames),
            extrasaction="ignore",
            lineterminator="\n",
        )
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})


def nonempty(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def yn(value: bool) -> str:
    return "1" if value else "0"


def compress_indices(indices: Iterable[int]) -> str:
    vals = sorted(set(int(x) for x in indices))
    if not vals:
        return ""
    out: List[str] = []
    start = prev = vals[0]
    for x in vals[1:]:
        if x == prev + 1:
            prev = x
            continue
        out.append(str(start) if start == prev else f"{start}-{prev}")
        start = prev = x
    out.append(str(start) if start == prev else f"{start}-{prev}")
    return ",".join(out)


def latest_diagnostics(
    files: Iterable[Path],
    key_fields: Sequence[str],
) -> Tuple[Dict[Tuple[str, ...], dict], List[str]]:
    """Index diagnostics by key, preferring the most recently modified file."""
    index: Dict[Tuple[str, ...], dict] = {}
    source_mtime: Dict[Tuple[str, ...], float] = {}
    unreadable: List[str] = []

    paths = sorted(set(Path(p) for p in files), key=lambda p: (p.stat().st_mtime if p.exists() else -1, str(p)))
    for path in paths:
        try:
            rows = read_tsv(path)
            mtime = path.stat().st_mtime
        except Exception as exc:
            unreadable.append(f"{path}: {exc}")
            continue
        for row in rows:
            if not all(k in row for k in key_fields):
                continue
            key = tuple(str(row.get(k, "")) for k in key_fields)
            if key not in source_mtime or mtime >= source_mtime[key]:
                rec = dict(row)
                rec["_diagnostic_file"] = str(path)
                rec["_diagnostic_mtime"] = str(mtime)
                index[key] = rec
                source_mtime[key] = mtime
    return index, unreadable


def expected_shard(row_index_zero_based: int, n_items: int, shard_size: int) -> Tuple[int, int]:
    if shard_size < 1:
        raise ValueError("shard size must be >= 1")
    n_shards = max(1, math.ceil(n_items / shard_size))
    shard = (row_index_zero_based % n_shards) + 1
    return shard, n_shards


def diagnostic_is_fresh_ok(diag_status: str, diag_message: str) -> bool:
    """True when diag explicitly says this run succeeded, not just 'already done'."""
    return (
        diag_status == "ok"
        and "already done" not in (diag_message or "").lower()
    )


def audit_susie(
    coloc_dir: Path,
    modalities: Sequence[str],
    shard_size: int,
) -> Tuple[List[dict], List[dict], List[dict], List[str]]:
    out_rows: List[dict] = []
    summary_rows: List[dict] = []
    recovery_rows: List[dict] = []
    problems: List[str] = []

    diag_dir = coloc_dir / "diagnostics"
    results_root = coloc_dir / "results"

    for mod in modalities:
        tasks_path = coloc_dir / "loci" / f"{mod}.tasks.tsv"
        if not nonempty(tasks_path):
            summary_rows.append({
                "stage": "susie_coloc",
                "group": mod,
                "manifest": str(tasks_path),
                "manifest_present": "0",
                "expected": "0",
                "terminal": "0",
                "with_primary_result": "0",
                "terminal_no_primary_result": "0",
                "partial_output": "0",
                "needs_rerun": "0",
                "missing_done": "0",
                "diagnostic_rows_seen": "0",
                "missing_diagnostic": "0",
                "expected_n_shards": "0",
                "rerun_shards": "",
                "missing_shard_diagnostics": "",
                "status_counts": "",
            })
            problems.append(f"SuSiE-coloc manifest missing/empty: {tasks_path}")
            continue

        tasks = read_tsv(tasks_path)
        if not tasks:
            summary_rows.append({
                "stage": "susie_coloc",
                "group": mod,
                "manifest": str(tasks_path),
                "manifest_present": "1",
                "expected": "0",
                "terminal": "0",
                "with_primary_result": "0",
                "terminal_no_primary_result": "0",
                "partial_output": "0",
                "needs_rerun": "0",
                "missing_done": "0",
                "diagnostic_rows_seen": "0",
                "missing_diagnostic": "0",
                "expected_n_shards": "0",
                "rerun_shards": "",
                "missing_shard_diagnostics": "",
                "status_counts": "",
            })
            continue

        required = {"modality", "phenotype_id", "trait_id", "ancestry"}
        missing_cols = required - set(tasks[0].keys())
        if missing_cols:
            problems.append(
                f"SuSiE-coloc manifest {tasks_path} missing columns: {sorted(missing_cols)}"
            )
            continue

        diag_files = list(diag_dir.glob(f"{mod}.shard-*.diagnostics.tsv"))
        diag_index, unreadable = latest_diagnostics(
            diag_files, ["modality", "phenotype_id", "trait_id", "ancestry"]
        )
        problems.extend(f"Unreadable SuSiE diagnostic: {x}" for x in unreadable)

        n_items = len(tasks)
        n_shards = max(1, math.ceil(n_items / shard_size))
        seen_bases: Dict[str, int] = defaultdict(int)
        status_counts: Counter = Counter()
        rerun_indices: set[int] = set()
        n_terminal = 0
        n_primary = 0
        n_terminal_no_primary = 0
        n_partial = 0
        n_rerun = 0
        n_missing_done = 0
        n_missing_diag = 0

        for i, task in enumerate(tasks):
            modality = str(task.get("modality", mod))
            phenotype = str(task.get("phenotype_id", ""))
            trait = str(task.get("trait_id", ""))
            anc = str(task.get("ancestry", ""))
            base = f"{anc}_{sanitize_id(phenotype)}_{trait}"
            seen_bases[base] += 1

            res_dir = results_root / modality
            coloc_path = res_dir / f"{base}.coloc.tsv"
            variants_path = res_dir / f"{base}.variants.tsv"
            done_path = res_dir / f"{base}.done"

            has_coloc = nonempty(coloc_path)
            has_variants = nonempty(variants_path)
            has_done = done_path.is_file()
            partial_output = has_coloc != has_variants

            dkey = (modality, phenotype, trait, anc)
            diag = diag_index.get(dkey, {})
            diag_status = str(diag.get("status", ""))
            diag_message = str(diag.get("message", ""))
            diag_file = str(diag.get("_diagnostic_file", ""))
            if diag_status:
                status_counts[diag_status] += 1
            else:
                n_missing_diag += 1

            fresh_ok = diagnostic_is_fresh_ok(diag_status, diag_message)

            # A missing sentinel always needs rerun. A partial output pair is
            # inconsistent. If the latest diagnostic says this run really was
            # successful, both primary outputs must exist.
            needs_rerun = (
                (not has_done)
                or partial_output
                or (has_done and fresh_ok and not (has_coloc and has_variants))
            )

            if has_done and not needs_rerun:
                n_terminal += 1
            if has_coloc and has_variants:
                n_primary += 1
            if has_done and not has_coloc and not has_variants and not needs_rerun:
                n_terminal_no_primary += 1
            if partial_output:
                n_partial += 1
            if not has_done:
                n_missing_done += 1
            if needs_rerun:
                n_rerun += 1

            shard, _ = expected_shard(i, n_items, shard_size)
            if needs_rerun:
                rerun_indices.add(shard)

            if needs_rerun:
                state = "needs_rerun"
            elif has_coloc and has_variants:
                state = "complete_with_result"
            elif has_done:
                state = "complete_terminal_no_result"
            else:
                state = "unknown"

            out_rows.append({
                "stage": "susie_coloc",
                "group": mod,
                "row_index": str(i + 1),
                "expected_shard": str(shard),
                "modality": modality,
                "phenotype_id": phenotype,
                "trait_id": trait,
                "ancestry": anc,
                "base": base,
                "state": state,
                "needs_rerun": yn(needs_rerun),
                "done": yn(has_done),
                "coloc_tsv": yn(has_coloc),
                "variants_tsv": yn(has_variants),
                "partial_output": yn(partial_output),
                "diag_status": diag_status,
                "diag_message": diag_message,
                "diagnostic_file": diag_file,
                "done_path": str(done_path),
                "primary_output": str(coloc_path),
                "secondary_output": str(variants_path),
            })

        collisions = sorted(k for k, n in seen_bases.items() if n > 1)
        if collisions:
            problems.append(
                f"SuSiE-coloc {mod}: {len(collisions)} sanitized output-name collisions; "
                f"examples: {', '.join(collisions[:5])}"
            )

        expected_diag = [
            diag_dir / f"{mod}.shard-{j:04d}.diagnostics.tsv"
            for j in range(1, n_shards + 1)
        ]
        missing_diag_shards = [
            j for j, p in enumerate(expected_diag, start=1) if not nonempty(p)
        ]

        summary_rows.append({
            "stage": "susie_coloc",
            "group": mod,
            "manifest": str(tasks_path),
            "manifest_present": "1",
            "expected": str(n_items),
            "terminal": str(n_terminal),
            "with_primary_result": str(n_primary),
            "terminal_no_primary_result": str(n_terminal_no_primary),
            "partial_output": str(n_partial),
            "needs_rerun": str(n_rerun),
            "missing_done": str(n_missing_done),
            "diagnostic_rows_seen": str(n_items - n_missing_diag),
            "missing_diagnostic": str(n_missing_diag),
            "expected_n_shards": str(n_shards),
            "rerun_shards": compress_indices(rerun_indices),
            "missing_shard_diagnostics": compress_indices(missing_diag_shards),
            "status_counts": ";".join(
                f"{k}:{v}" for k, v in sorted(status_counts.items())
            ),
        })

        for shard in sorted(rerun_indices):
            count = sum(
                1 for r in out_rows
                if r["stage"] == "susie_coloc"
                and r["group"] == mod
                and r["needs_rerun"] == "1"
                and int(r["expected_shard"]) == shard
            )
            recovery_rows.append({
                "stage": "susie_coloc",
                "group": mod,
                "array_job_name": f"coloc_{mod}",
                "shard_index": str(shard),
                "n_items_needing_rerun": str(count),
                "shard_size_assumed": str(shard_size),
                "expected_n_shards": str(n_shards),
            })

    return out_rows, summary_rows, recovery_rows, problems


def audit_colocboost(
    coloc_dir: Path,
    ancestries: Sequence[str],
    shard_size: int,
) -> Tuple[List[dict], List[dict], List[dict], List[str]]:
    cb_dir = coloc_dir / "colocboost"
    out_rows: List[dict] = []
    summary_rows: List[dict] = []
    recovery_rows: List[dict] = []
    problems: List[str] = []

    for anc in ancestries:
        regions_path = cb_dir / f"{anc}.regions.tsv"
        outcomes_path = cb_dir / f"{anc}.outcomes.tsv"

        if not nonempty(regions_path):
            summary_rows.append({
                "stage": "colocboost",
                "group": anc,
                "manifest": str(regions_path),
                "manifest_present": "0",
                "expected": "0",
                "terminal": "0",
                "with_primary_result": "0",
                "terminal_no_primary_result": "0",
                "partial_output": "0",
                "needs_rerun": "0",
                "missing_done": "0",
                "diagnostic_rows_seen": "0",
                "missing_diagnostic": "0",
                "expected_n_shards": "0",
                "rerun_shards": "",
                "missing_shard_diagnostics": "",
                "status_counts": "",
            })
            problems.append(f"colocBoost regions manifest missing/empty: {regions_path}")
            continue
        if not nonempty(outcomes_path):
            problems.append(f"colocBoost outcomes manifest missing/empty: {outcomes_path}")

        regions = read_tsv(regions_path)
        if not regions:
            continue

        required = {"region_id", "ancestry"}
        missing_cols = required - set(regions[0].keys())
        if missing_cols:
            problems.append(
                f"colocBoost manifest {regions_path} missing columns: {sorted(missing_cols)}"
            )
            continue

        diag_files = list((cb_dir / "diagnostics").glob(f"{anc}.shard-*.diagnostics.tsv"))
        diag_index, unreadable = latest_diagnostics(
            diag_files, ["region_id", "ancestry"]
        )
        problems.extend(f"Unreadable colocBoost diagnostic: {x}" for x in unreadable)

        n_items = len(regions)
        n_shards = max(1, math.ceil(n_items / shard_size))
        status_counts: Counter = Counter()
        rerun_indices: set[int] = set()
        n_terminal = 0
        n_primary = 0
        n_terminal_no_primary = 0
        n_partial = 0
        n_rerun = 0
        n_missing_done = 0
        n_missing_diag = 0

        for i, reg in enumerate(regions):
            region_id = str(reg.get("region_id", ""))
            ancestry = str(reg.get("ancestry", anc))
            res_dir = cb_dir / "results" / ancestry
            clusters_path = res_dir / f"{region_id}.clusters.tsv"
            vcp_path = res_dir / f"{region_id}.vcp.tsv"
            done_path = res_dir / f"{region_id}.done"

            has_clusters = nonempty(clusters_path)
            has_vcp = nonempty(vcp_path)
            has_done = done_path.is_file()

            dkey = (region_id, ancestry)
            diag = diag_index.get(dkey, {})
            diag_status = str(diag.get("status", ""))
            diag_message = str(diag.get("message", ""))
            diag_file = str(diag.get("_diagnostic_file", ""))
            if diag_status:
                status_counts[diag_status] += 1
            else:
                n_missing_diag += 1

            fresh_ok = diagnostic_is_fresh_ok(diag_status, diag_message)

            # vcp without clusters is inconsistent because clusters are written
            # first on every successful colocBoost run.
            partial_output = has_vcp and not has_clusters
            needs_rerun = (
                (not has_done)
                or partial_output
                or (has_done and fresh_ok and not has_clusters)
            )

            if has_done and not needs_rerun:
                n_terminal += 1
            if has_clusters:
                n_primary += 1
            if has_done and not has_clusters and not needs_rerun:
                n_terminal_no_primary += 1
            if partial_output:
                n_partial += 1
            if not has_done:
                n_missing_done += 1
            if needs_rerun:
                n_rerun += 1

            shard, _ = expected_shard(i, n_items, shard_size)
            if needs_rerun:
                rerun_indices.add(shard)

            if needs_rerun:
                state = "needs_rerun"
            elif has_clusters:
                state = "complete_with_result"
            elif has_done:
                state = "complete_terminal_no_result"
            else:
                state = "unknown"

            out_rows.append({
                "stage": "colocboost",
                "group": anc,
                "row_index": str(i + 1),
                "expected_shard": str(shard),
                "region_id": region_id,
                "ancestry": ancestry,
                "state": state,
                "needs_rerun": yn(needs_rerun),
                "done": yn(has_done),
                "clusters_tsv": yn(has_clusters),
                "vcp_tsv": yn(has_vcp),
                "partial_output": yn(partial_output),
                "diag_status": diag_status,
                "diag_message": diag_message,
                "diagnostic_file": diag_file,
                "done_path": str(done_path),
                "primary_output": str(clusters_path),
                "secondary_output": str(vcp_path),
            })

        expected_diag = [
            cb_dir / "diagnostics" / f"{anc}.shard-{j:04d}.diagnostics.tsv"
            for j in range(1, n_shards + 1)
        ]
        missing_diag_shards = [
            j for j, p in enumerate(expected_diag, start=1) if not nonempty(p)
        ]

        summary_rows.append({
            "stage": "colocboost",
            "group": anc,
            "manifest": str(regions_path),
            "manifest_present": "1",
            "expected": str(n_items),
            "terminal": str(n_terminal),
            "with_primary_result": str(n_primary),
            "terminal_no_primary_result": str(n_terminal_no_primary),
            "partial_output": str(n_partial),
            "needs_rerun": str(n_rerun),
            "missing_done": str(n_missing_done),
            "diagnostic_rows_seen": str(n_items - n_missing_diag),
            "missing_diagnostic": str(n_missing_diag),
            "expected_n_shards": str(n_shards),
            "rerun_shards": compress_indices(rerun_indices),
            "missing_shard_diagnostics": compress_indices(missing_diag_shards),
            "status_counts": ";".join(
                f"{k}:{v}" for k, v in sorted(status_counts.items())
            ),
        })

        for shard in sorted(rerun_indices):
            count = sum(
                1 for r in out_rows
                if r["stage"] == "colocboost"
                and r["group"] == anc
                and r["needs_rerun"] == "1"
                and int(r["expected_shard"]) == shard
            )
            recovery_rows.append({
                "stage": "colocboost",
                "group": anc,
                "array_job_name": f"colocboost_{anc}",
                "shard_index": str(shard),
                "n_items_needing_rerun": str(count),
                "shard_size_assumed": str(shard_size),
                "expected_n_shards": str(n_shards),
            })

    return out_rows, summary_rows, recovery_rows, problems


def print_summary(summary_rows: Sequence[dict], problems: Sequence[str]) -> None:
    print("=== Module 06 coloc completeness ===")
    print(
        f"{'stage':<13} {'group':<20} {'expected':>9} {'terminal':>9} "
        f"{'result':>9} {'rerun':>7} {'missing .done':>13} {'rerun shards'}"
    )
    print("-" * 104)
    for r in summary_rows:
        print(
            f"{r.get('stage',''):<13} {r.get('group',''):<20} "
            f"{r.get('expected','0'):>9} {r.get('terminal','0'):>9} "
            f"{r.get('with_primary_result','0'):>9} {r.get('needs_rerun','0'):>7} "
            f"{r.get('missing_done','0'):>13} {r.get('rerun_shards','')}"
        )

    total_expected = sum(int(r.get("expected", 0) or 0) for r in summary_rows)
    total_rerun = sum(int(r.get("needs_rerun", 0) or 0) for r in summary_rows)
    print()
    print(f"Expected units: {total_expected}")
    print(f"Need rerun:     {total_rerun}")

    if problems:
        print()
        print("Warnings / manifest problems:")
        for p in problems:
            print(f"  - {p}")


def main() -> int:
    p = argparse.ArgumentParser(
        description="Check completeness of Module 06 SuSiE-coloc and colocBoost."
    )
    p.add_argument(
        "--results-dir",
        required=True,
        help="Module output qtl_results directory",
    )
    p.add_argument(
        "--coloc-dir",
        default=None,
        help="override coloc directory [default: <results-dir>/coloc]",
    )
    p.add_argument(
        "--modalities",
        nargs="+",
        default=DEFAULT_MODALITIES,
        help="SuSiE-coloc modalities to audit",
    )
    p.add_argument(
        "--ancestries",
        nargs="+",
        default=DEFAULT_ANCESTRIES,
        help="colocBoost ancestries to audit",
    )
    p.add_argument(
        "--susie-shard-size",
        type=int,
        default=25,
        help="SHARD_SIZE used by 42_submit_coloc.sh [default: 25]",
    )
    p.add_argument(
        "--colocboost-shard-size",
        type=int,
        default=10,
        help="SHARD_SIZE used by 45_submit_colocboost.sh [default: 10]",
    )
    p.add_argument(
        "--out-dir",
        default=None,
        help="audit output directory [default: <coloc-dir>/completeness]",
    )
    p.add_argument(
        "--fail-on-incomplete",
        action="store_true",
        help="exit 2 if any expected unit needs rerun or a manifest problem is detected",
    )
    args = p.parse_args()

    if args.susie_shard_size < 1 or args.colocboost_shard_size < 1:
        p.error("shard sizes must be >= 1")

    results_dir = Path(args.results_dir)
    coloc_dir = Path(args.coloc_dir) if args.coloc_dir else results_dir / "coloc"
    out_dir = Path(args.out_dir) if args.out_dir else coloc_dir / "completeness"
    out_dir.mkdir(parents=True, exist_ok=True)

    susie_rows, susie_summary, susie_recovery, susie_problems = audit_susie(
        coloc_dir, args.modalities, args.susie_shard_size
    )
    cb_rows, cb_summary, cb_recovery, cb_problems = audit_colocboost(
        coloc_dir, args.ancestries, args.colocboost_shard_size
    )

    summary_rows = susie_summary + cb_summary
    recovery_rows = susie_recovery + cb_recovery
    problems = susie_problems + cb_problems

    summary_fields = [
        "stage", "group", "manifest", "manifest_present", "expected", "terminal",
        "with_primary_result", "terminal_no_primary_result", "partial_output",
        "needs_rerun", "missing_done", "diagnostic_rows_seen",
        "missing_diagnostic", "expected_n_shards", "rerun_shards",
        "missing_shard_diagnostics", "status_counts",
    ]
    susie_fields = [
        "stage", "group", "row_index", "expected_shard", "modality",
        "phenotype_id", "trait_id", "ancestry", "base", "state",
        "needs_rerun", "done", "coloc_tsv", "variants_tsv", "partial_output",
        "diag_status", "diag_message", "diagnostic_file", "done_path",
        "primary_output", "secondary_output",
    ]
    cb_fields = [
        "stage", "group", "row_index", "expected_shard", "region_id", "ancestry",
        "state", "needs_rerun", "done", "clusters_tsv", "vcp_tsv",
        "partial_output", "diag_status", "diag_message", "diagnostic_file",
        "done_path", "primary_output", "secondary_output",
    ]
    recovery_fields = [
        "stage", "group", "array_job_name", "shard_index",
        "n_items_needing_rerun", "shard_size_assumed", "expected_n_shards",
    ]

    write_tsv(out_dir / "summary.tsv", summary_rows, summary_fields)
    write_tsv(out_dir / "susie_tasks.tsv", susie_rows, susie_fields)
    write_tsv(out_dir / "colocboost_regions.tsv", cb_rows, cb_fields)
    write_tsv(out_dir / "recovery_shards.tsv", recovery_rows, recovery_fields)

    print_summary(summary_rows, problems)
    print()
    print(f"Audit files -> {out_dir}")
    print(f"  {out_dir / 'summary.tsv'}")
    print(f"  {out_dir / 'susie_tasks.tsv'}")
    print(f"  {out_dir / 'colocboost_regions.tsv'}")
    print(f"  {out_dir / 'recovery_shards.tsv'}")
    print()
    print(
        "NOTE: recovery shard indices assume "
        f"SuSiE SHARD_SIZE={args.susie_shard_size} and "
        f"colocBoost SHARD_SIZE={args.colocboost_shard_size}."
    )

    total_rerun = sum(int(r.get("needs_rerun", 0) or 0) for r in summary_rows)
    if args.fail_on_incomplete and (total_rerun > 0 or problems):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
