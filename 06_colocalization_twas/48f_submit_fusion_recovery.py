#!/usr/bin/env python3
"""
48f_submit_fusion_recovery.py

Read the TSV produced by 48e_check_fusion_incomplete.py and submit only the
still-incomplete adaptive FUSION Step-48 array indices, followed by one
dependent merge job per recovered pair.

Safety / restart behavior
-------------------------
* Re-checks the filesystem at submission time. If shards completed after 48e
  wrote the TSV, those indices are NOT resubmitted.
* Skips pairs whose final merged TWAS output has appeared.
* Refuses to duplicate a pair while an unfinished LSF job with the normal
  Step-48 array or merge job name is still present.
* Uses the original array indices from SHARD_MANIFEST, so 48a_run_fusion.sh
  reads exactly the intended manifest rows.
* If all shards are already complete but the merged file is absent, submits
  only 48c_merge_fusion.sh.
* Writes a submission log TSV containing the exact indices and LSF job IDs.

This script assumes the adaptive feature-shard Step-48 implementation:
  48a_run_fusion.sh
  48c_merge_fusion.sh
  48d_prepare_fusion_shards.py

Python standard library only.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable


JOB_ID_RE = re.compile(r"Job <(\d+)>")


def nonempty(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def compress_indices(indices: Iterable[int]) -> str:
    vals = sorted(set(int(x) for x in indices))
    if not vals:
        return ""
    out: list[str] = []
    start = prev = vals[0]
    for x in vals[1:]:
        if x == prev + 1:
            prev = x
            continue
        out.append(str(start) if start == prev else f"{start}-{prev}")
        start = prev = x
    out.append(str(start) if start == prev else f"{start}-{prev}")
    return ",".join(out)


def positive_int(text: str) -> int:
    x = int(text)
    if x < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return x


def env_int(name: str, default: int) -> int:
    text = os.environ.get(name)
    if text is None or text == "":
        return default
    try:
        x = int(text)
    except ValueError as exc:
        raise SystemExit(f"ERROR: ${name} must be an integer, got {text!r}") from exc
    if x < 1:
        raise SystemExit(f"ERROR: ${name} must be >= 1, got {x}")
    return x


def safe_env_string(values: dict[str, str]) -> str:
    for key, value in values.items():
        if "," in value or "\n" in value:
            raise ValueError(
                f"cannot pass {key} through LSF -env because its value contains "
                f"a comma/newline: {value!r}"
            )
    return ",".join(f"{k}={v}" for k, v in values.items())


def read_recovery(path: Path) -> list[dict[str, str]]:
    if not nonempty(path):
        raise FileNotFoundError(f"recovery TSV missing/empty: {path}")
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {
            "action",
            "weight_set",
            "trait_id",
            "ld_ancestry",
            "manifest",
            "job_name",
            "merge_job_name",
            "final_output",
            "per_shard_dir",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        rows = list(reader)
    return rows


def read_manifest(path: Path) -> dict[int, dict[str, str]]:
    if not nonempty(path):
        raise FileNotFoundError(f"shard manifest missing/empty: {path}")
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {
            "array_index",
            "chromosome",
            "shard_index",
            "n_features",
            "weights_file",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")

        rows: dict[int, dict[str, str]] = {}
        for row in reader:
            idx = int(row["array_index"])
            if idx in rows:
                raise ValueError(f"{path}: duplicate array_index {idx}")
            rows[idx] = row

    if not rows:
        raise ValueError(f"{path}: no shard rows")
    expected = list(range(1, len(rows) + 1))
    if sorted(rows) != expected:
        raise ValueError(
            f"{path}: array_index must be contiguous 1..N; "
            f"observed {compress_indices(rows)}"
        )
    return rows


def shard_stem(row: dict[str, str]) -> str:
    return f"chr{int(row['chromosome']):02d}.shard{int(row['shard_index']):04d}"


def still_missing_indices(
    manifest_rows: dict[int, dict[str, str]], per_shard_dir: Path
) -> list[int]:
    missing: list[int] = []
    for idx in sorted(manifest_rows):
        stem = shard_stem(manifest_rows[idx])
        out = per_shard_dir / f"{stem}.twas.tsv"
        done = Path(str(out) + ".done")
        if not (nonempty(out) and done.is_file()):
            missing.append(idx)
    return missing


def active_named_job(name: str) -> bool:
    """Use the same unfinished-job semantics as the Step-48 driver."""
    proc = subprocess.run(
        ["bjobs", "-J", name],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    text = (proc.stdout + "\n" + proc.stderr).strip()
    if proc.returncode != 0:
        return False
    if not text:
        return False
    if "No unfinished job found" in text or "not found" in text.lower():
        return False
    # A successful bjobs -J query with nonempty output means at least one
    # unfinished job matched that name.
    return True


def run_bsub(cmd: list[str], script_path: Path, dry_run: bool) -> str | None:
    print("  " + " ".join(cmd))
    print(f"    < {script_path}")
    if dry_run:
        return None

    script_text = script_path.read_text(encoding="utf-8")
    proc = subprocess.run(
        cmd,
        input=script_text,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.stdout:
        print(proc.stdout.rstrip())
    if proc.stderr:
        print(proc.stderr.rstrip(), file=sys.stderr)
    if proc.returncode != 0:
        raise RuntimeError(
            f"bsub failed ({proc.returncode}) for {script_path.name}"
        )
    match = JOB_ID_RE.search(proc.stdout + "\n" + proc.stderr)
    if not match:
        raise RuntimeError(
            f"could not parse LSF job ID from bsub output for {script_path.name}"
        )
    return match.group(1)


def derive_isotwas_dir(final_output: Path) -> Path:
    # .../<isotwas>/fusion/<weight_set>/<weight_set>_<trait>.twas.tsv
    try:
        return final_output.resolve().parents[2]
    except IndexError as exc:
        raise ValueError(f"cannot infer ISOTWAS_DIR from {final_output}") from exc


def preflight_missing(
    missing: list[int],
    manifest_rows: dict[int, dict[str, str]],
    sumstats_prefix: Path,
    ldref_prefix: Path,
) -> None:
    chroms: set[int] = set()
    for idx in missing:
        row = manifest_rows[idx]
        wgt = Path(row["weights_file"])
        if not nonempty(wgt):
            raise FileNotFoundError(
                f"array index {idx}: shard weights missing/empty: {wgt}"
            )
        chroms.add(int(row["chromosome"]))

    for chrom in sorted(chroms):
        sumstats = Path(f"{sumstats_prefix}{chrom}.fusion.tsv")
        if not nonempty(sumstats):
            raise FileNotFoundError(
                f"chr{chrom}: chromosome FUSION GWAS missing/empty: {sumstats}"
            )

        prefix = Path(f"{ldref_prefix}{chrom}")
        for suffix in (".bed", ".bim", ".fam"):
            p = Path(str(prefix) + suffix)
            if not p.is_file():
                raise FileNotFoundError(
                    f"chr{chrom}: FUSION LD reference component missing: {p}"
                )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Submit only missing adaptive FUSION Step-48 shard indices from "
            "48e recovery output, plus dependent merges."
        )
    )
    parser.add_argument(
        "--recovery-tsv",
        type=Path,
        default=Path("fusion_recovery.tsv"),
        help="Output from 48e_check_fusion_incomplete.py (default: fusion_recovery.tsv)",
    )
    parser.add_argument(
        "--scripts-dir",
        type=Path,
        default=Path(os.environ["SCRIPTS_DIR"]) if os.environ.get("SCRIPTS_DIR") else Path(__file__).resolve().parent,
        help="06_colocalization_twas directory (default: $SCRIPTS_DIR or script directory)",
    )
    parser.add_argument(
        "--gwas-dir",
        type=Path,
        default=Path(os.environ["GWAS_DIR"]) if os.environ.get("GWAS_DIR") else None,
        help="GWAS root (default: $GWAS_DIR; otherwise derived from $CONFIG)",
    )
    parser.add_argument(
        "--ldref-dir",
        type=Path,
        default=Path(os.environ["LDREF_DIR"]) if os.environ.get("LDREF_DIR") else None,
        help="FUSION LD reference directory (default: $LDREF_DIR or GWAS_DIR/fusion_ldref)",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=Path(os.environ["LOG_DIR"]) if os.environ.get("LOG_DIR") else None,
        help="LSF log directory (default: $LOG_DIR; otherwise derived from $CONFIG)",
    )
    parser.add_argument(
        "--fusion-dir",
        type=Path,
        default=Path(
            os.environ.get(
                "FUSION_DIR",
                "/rsrch5/home/epi/bhattacharya_lab/software/fusion_twas",
            )
        ),
        help="FUSION clone directory",
    )
    parser.add_argument(
        "--rscript",
        type=Path,
        default=Path(os.environ["RSCRIPT"]) if os.environ.get("RSCRIPT") else None,
        help="Rscript wrapper (default: $RSCRIPT or <repo>/bin/Rscript_sif)",
    )
    parser.add_argument(
        "--queue",
        default=os.environ.get("QUEUE", "medium"),
    )
    parser.add_argument(
        "--walltime",
        default=os.environ.get("WALLTIME", "04:00"),
        help="per-shard walltime (default: $WALLTIME or 04:00)",
    )
    parser.add_argument(
        "--mem-gb",
        type=positive_int,
        default=env_int("MEM_GB", 128),
        help="memory per shard array element in GB (default: $MEM_GB or 128)",
    )
    parser.add_argument(
        "--max-parallel",
        type=positive_int,
        default=env_int("FUSION_MAX_PARALLEL", 4),
        help="maximum concurrent recovery indices per pair",
    )
    parser.add_argument(
        "--merge-walltime",
        default=os.environ.get("MERGE_WALLTIME", "00:30"),
    )
    parser.add_argument(
        "--merge-mem-gb",
        type=positive_int,
        default=env_int("MERGE_MEM_GB", 4),
    )
    parser.add_argument(
        "--min-r2pred",
        default=os.environ.get("MIN_R2PRED", "0.7"),
    )
    parser.add_argument(
        "--submission-log",
        type=Path,
        default=None,
        help="TSV recording submitted indices/job IDs (default: <recovery>.submitted.tsv)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="submit at most N recovery rows; 0 means all",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="show exact bsub commands without submitting",
    )
    args = parser.parse_args()

    if args.limit < 0:
        parser.error("--limit must be >= 0")

    scripts_dir = args.scripts_dir.resolve()
    worker = scripts_dir / "48a_run_fusion.sh"
    merger = scripts_dir / "48c_merge_fusion.sh"
    if not worker.is_file():
        parser.error(f"missing worker: {worker}")
    if not merger.is_file():
        parser.error(f"missing merger: {merger}")

    config = Path(os.environ["CONFIG"]).resolve() if os.environ.get("CONFIG") else None
    output_base = (
        Path(os.environ["OUTPUT_BASE"]).resolve()
        if os.environ.get("OUTPUT_BASE")
        else (config.parent if config is not None else None)
    )

    gwas_dir = (
        args.gwas_dir.resolve()
        if args.gwas_dir is not None
        else ((output_base / "gwas") if output_base is not None else None)
    )
    if gwas_dir is None:
        parser.error("--gwas-dir is required unless GWAS_DIR or CONFIG is set")

    ldref_dir = (
        args.ldref_dir.resolve()
        if args.ldref_dir is not None
        else (gwas_dir / "fusion_ldref")
    )
    log_dir = (
        args.log_dir.resolve()
        if args.log_dir is not None
        else ((output_base / "logs") if output_base is not None else None)
    )
    if log_dir is None:
        parser.error("--log-dir is required unless LOG_DIR or CONFIG is set")
    log_dir.mkdir(parents=True, exist_ok=True)

    repo_root = scripts_dir.parent.resolve()
    rscript = (
        args.rscript.resolve()
        if args.rscript is not None
        else (repo_root / "bin" / "Rscript_sif")
    )
    fusion_dir = args.fusion_dir.resolve()

    if not rscript.is_file():
        parser.error(f"Rscript wrapper not found: {rscript}")
    if not (fusion_dir / "FUSION.assoc_test.R").is_file():
        parser.error(f"FUSION.assoc_test.R not found under: {fusion_dir}")

    recovery = args.recovery_tsv.resolve()
    try:
        rows = read_recovery(recovery)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    submission_log = (
        args.submission_log.resolve()
        if args.submission_log is not None
        else recovery.with_name(recovery.stem + ".submitted.tsv")
    )

    log_fields = [
        "timestamp",
        "weight_set",
        "trait_id",
        "submit_action",
        "array_indices",
        "n_array_indices",
        "array_job_id",
        "merge_job_id",
        "manifest",
        "final_output",
        "status",
        "message",
    ]
    submission_records: list[dict[str, str]] = []

    n_seen = 0
    n_submitted_pairs = 0
    n_skipped_complete = 0
    n_skipped_active = 0
    n_failed = 0

    print("=== 48f_submit_fusion_recovery.py ===")
    print(f"  recovery TSV:       {recovery}")
    print(f"  GWAS_DIR:           {gwas_dir}")
    print(f"  LDREF_DIR:          {ldref_dir}")
    print(f"  queue:              {args.queue}")
    print(f"  shard resources:    1 core, {args.mem_gb}G, {args.walltime}")
    print(f"  max parallel/pair:  {args.max_parallel}")
    print(f"  merge resources:    1 core, {args.merge_mem_gb}G, {args.merge_walltime}")
    print(f"  dry run:            {args.dry_run}")
    print()

    for rec in rows:
        if args.limit and n_seen >= args.limit:
            break
        n_seen += 1

        ws = rec["weight_set"].strip()
        trait = rec["trait_id"].strip()
        ld_anc = rec["ld_ancestry"].strip()
        manifest = Path(rec["manifest"]).resolve()
        final = Path(rec["final_output"]).resolve()
        per_shard_dir = Path(rec["per_shard_dir"]).resolve()
        job_name = rec["job_name"].strip() or f"fusion_{ws}_{trait}"
        merge_job_name = rec["merge_job_name"].strip() or f"fusionmerge_{ws}_{trait}"

        timestamp = dt.datetime.now().astimezone().isoformat(timespec="seconds")
        base_record = {
            "timestamp": timestamp,
            "weight_set": ws,
            "trait_id": trait,
            "manifest": str(manifest),
            "final_output": str(final),
        }

        print(f"{ws} x {trait}")

        try:
            if nonempty(final):
                print(f"  SKIP: merged output now exists: {final}")
                n_skipped_complete += 1
                submission_records.append(
                    {
                        **base_record,
                        "submit_action": "none",
                        "array_indices": "",
                        "n_array_indices": "0",
                        "array_job_id": "",
                        "merge_job_id": "",
                        "status": "skipped_complete",
                        "message": "merged output appeared before recovery submission",
                    }
                )
                continue

            manifest_rows = read_manifest(manifest)
            missing = still_missing_indices(manifest_rows, per_shard_dir)

            # Do not create duplicate workers/merges while a normal Step-48
            # job with the same pair-level name is unfinished.
            active_array = active_named_job(job_name)
            active_merge = active_named_job(merge_job_name)
            if active_array or active_merge:
                which = []
                if active_array:
                    which.append(job_name)
                if active_merge:
                    which.append(merge_job_name)
                msg = "unfinished LSF job(s) still active: " + ", ".join(which)
                print(f"  SKIP: {msg}")
                n_skipped_active += 1
                submission_records.append(
                    {
                        **base_record,
                        "submit_action": "none",
                        "array_indices": compress_indices(missing),
                        "n_array_indices": str(len(missing)),
                        "array_job_id": "",
                        "merge_job_id": "",
                        "status": "skipped_active",
                        "message": msg,
                    }
                )
                continue

            isotwas_dir = derive_isotwas_dir(final)
            expected_per_shard = (
                isotwas_dir
                / "fusion"
                / ws
                / "per_shard"
                / trait
                / manifest.parent.name
            ).resolve()
            if expected_per_shard != per_shard_dir:
                raise ValueError(
                    "recovery row is internally inconsistent: "
                    f"per_shard_dir={per_shard_dir}, expected={expected_per_shard}"
                )

            if missing:
                indices = compress_indices(missing)
                parallel = min(args.max_parallel, len(missing))
                sumstats_prefix = gwas_dir / "fusion" / trait / "chr"
                ldref_prefix = ldref_dir / f"{ld_anc}.1kg.chr"

                preflight_missing(
                    missing,
                    manifest_rows,
                    sumstats_prefix=sumstats_prefix,
                    ldref_prefix=ldref_prefix,
                )

                array_env = safe_env_string(
                    {
                        "SCRIPTS_DIR": str(scripts_dir),
                        "ISOTWAS_DIR": str(isotwas_dir),
                        "WEIGHT_SET": ws,
                        "TRAIT_ID": trait,
                        "SUMSTATS_PREFIX": str(sumstats_prefix),
                        "LDREF_PREFIX": str(ldref_prefix),
                        "SHARD_MANIFEST": str(manifest),
                        "MIN_R2PRED": str(args.min_r2pred),
                        "FUSION_DIR": str(fusion_dir),
                        "RSCRIPT": str(rscript),
                        # Missing indices are already selected. FORCE_RUN=0
                        # retains 48a's sentinel-based restart protection.
                        "FORCE_RUN": "0",
                    }
                )

                array_cmd = [
                    "bsub",
                    "-J",
                    f"{job_name}[{indices}]%{parallel}",
                    "-q",
                    args.queue,
                    "-n",
                    "1",
                    "-W",
                    args.walltime,
                    "-M",
                    str(args.mem_gb),
                    "-R",
                    f"rusage[mem={args.mem_gb}]",
                    "-o",
                    str(log_dir / f"{job_name}.%J_%I.out"),
                    "-e",
                    str(log_dir / f"{job_name}.%J_%I.err"),
                    "-env",
                    array_env,
                ]

                print(
                    f"  recover {len(missing)}/{len(manifest_rows)} shards: "
                    f"[{indices}] (max {parallel} concurrent)"
                )
                array_job_id = run_bsub(array_cmd, worker, args.dry_run)

                merge_env = safe_env_string(
                    {
                        "ISOTWAS_DIR": str(isotwas_dir),
                        "WEIGHT_SET": ws,
                        "TRAIT_ID": trait,
                        "SHARD_MANIFEST": str(manifest),
                    }
                )
                merge_cmd = [
                    "bsub",
                    "-J",
                    merge_job_name,
                    "-q",
                    args.queue,
                    "-n",
                    "1",
                    "-W",
                    args.merge_walltime,
                    "-M",
                    str(args.merge_mem_gb),
                    "-R",
                    f"rusage[mem={args.merge_mem_gb}]",
                ]
                if args.dry_run:
                    merge_cmd += ["-w", "ended(<RECOVERY_ARRAY_JOB_ID>)"]
                else:
                    assert array_job_id is not None
                    merge_cmd += ["-w", f"ended({array_job_id})"]
                merge_cmd += [
                    "-o",
                    str(log_dir / f"{merge_job_name}.%J.out"),
                    "-e",
                    str(log_dir / f"{merge_job_name}.%J.err"),
                    "-env",
                    merge_env,
                ]
                merge_job_id = run_bsub(merge_cmd, merger, args.dry_run)

                submission_records.append(
                    {
                        **base_record,
                        "submit_action": "resubmit_shards",
                        "array_indices": indices,
                        "n_array_indices": str(len(missing)),
                        "array_job_id": array_job_id or "",
                        "merge_job_id": merge_job_id or "",
                        "status": "dry_run" if args.dry_run else "submitted",
                        "message": "",
                    }
                )
                n_submitted_pairs += 1

            else:
                # Filesystem says every shard succeeded; only the final merge is absent.
                merge_env = safe_env_string(
                    {
                        "ISOTWAS_DIR": str(isotwas_dir),
                        "WEIGHT_SET": ws,
                        "TRAIT_ID": trait,
                        "SHARD_MANIFEST": str(manifest),
                    }
                )
                merge_cmd = [
                    "bsub",
                    "-J",
                    merge_job_name,
                    "-q",
                    args.queue,
                    "-n",
                    "1",
                    "-W",
                    args.merge_walltime,
                    "-M",
                    str(args.merge_mem_gb),
                    "-R",
                    f"rusage[mem={args.merge_mem_gb}]",
                    "-o",
                    str(log_dir / f"{merge_job_name}.%J.out"),
                    "-e",
                    str(log_dir / f"{merge_job_name}.%J.err"),
                    "-env",
                    merge_env,
                ]
                print("  all shards complete; submit merge only")
                merge_job_id = run_bsub(merge_cmd, merger, args.dry_run)
                submission_records.append(
                    {
                        **base_record,
                        "submit_action": "merge_only",
                        "array_indices": "",
                        "n_array_indices": "0",
                        "array_job_id": "",
                        "merge_job_id": merge_job_id or "",
                        "status": "dry_run" if args.dry_run else "submitted",
                        "message": "",
                    }
                )
                n_submitted_pairs += 1

        except (OSError, ValueError, RuntimeError) as exc:
            n_failed += 1
            print(f"  ERROR: {exc}", file=sys.stderr)
            submission_records.append(
                {
                    **base_record,
                    "submit_action": "error",
                    "array_indices": "",
                    "n_array_indices": "0",
                    "array_job_id": "",
                    "merge_job_id": "",
                    "status": "error",
                    "message": str(exc),
                }
            )

    # Write/overwrite the log for this invocation. It is intentionally separate
    # from the 48e recovery TSV so the latter remains an immutable status snapshot.
    submission_log.parent.mkdir(parents=True, exist_ok=True)
    with submission_log.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=log_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(submission_records)

    print()
    print(f"Submission log:        {submission_log}")
    print(f"Recovery rows checked: {n_seen}")
    print(f"Pairs submitted:       {n_submitted_pairs}")
    print(f"Skipped complete:      {n_skipped_complete}")
    print(f"Skipped active:        {n_skipped_active}")
    print(f"Errors:                {n_failed}")

    return 1 if n_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
