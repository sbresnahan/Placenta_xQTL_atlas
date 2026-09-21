#!/usr/bin/env python3
"""
run_pipeline.py — Driver for PANTRY rewrite on MDACC seadragon HPC (LSF/BSUB)

Reads the master config.yml, then for a single cohort:
  0. Stages cohort inputs: symlinks fastq_map, samples_file, and fastq_dir
     from their config-specified source paths into ${COHORT_DIR} so the stage
     scripts (which hardcode ${COHORT_DIR}/fastq_map.txt, /samples.txt, /fastq)
     can find them.
  1. Submits per-sample STAR alignment jobs (fan-out)
  2. Submits per-sample quantification jobs (Salmon, RegTools, featureCounts,
     RNA editing) with LSF dependencies on each sample's BAM job
  3. Submits per-modality aggregation jobs with LSF dependencies on all
     relevant per-sample jobs
  4. Submits final indexing job
  5. (Default) Monitors all submitted jobs by polling `bjobs` every 60s; on
     detecting any EXIT job, prints a failure report (which job, its
     stage/sample/modality, and which downstream jobs are stranded by the
     unsatisfiable done() dependency) and exits non-zero. Exits 0 once the
     final index job reaches DONE.

Each stage script receives --config <path> and --scripts-dir <path> as its
first two arguments, so the script can load config.yml via config_get.py and
locate its own directory. The remaining arguments are cohort/sample/etc as
before.

LSF dependency model: bsub -w "done(jobA)" or "done(jobA) && done(jobB) && ..."

Resource model: the driver passes -n, -M, -R, -W explicitly on the bsub
command line (per-stage, from STAGE_RESOURCES below). This is REQUIRED because
the driver invokes bsub with the script as an argument (bsub ... script.sh),
and LSF only parses embedded #BSUB directives when the script is fed on stdin
(bsub < script.sh). The argument form treats #BSUB lines as plain comments,
so the per-stage resources defined in each script's #BSUB block are silently
ignored. The driver therefore re-declares them on the command line.

Usage:
    python3 run_pipeline.py --config config.yml --cohort cohort1
    python3 run_pipeline.py --config config.yml --cohort cohort1 --dry-run
    python3 run_pipeline.py --config config.yml --cohort cohort1 --no-monitor
    python3 run_pipeline.py --config config.yml --cohort cohort1 --monitor-all
    python3 run_pipeline.py --config config.yml --cohort cohort1 --scripts-dir /path/to/scripts

The driver does NOT run Tier 1 (shared reference prep). Run 001_refprep_pre.sh,
002_txrevise_array.sh, 003_refprep_post.sh separately first. The driver assumes
Tier 1 outputs exist in reference_dir.

Python compatibility: written for Python 3.6+ (seadragon's system python3 is
3.6). Avoids subprocess.run(capture_output=..., text=...) which require 3.7+;
uses the _run_capture() helper instead.
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML not installed. Install with: pip install pyyaml", file=sys.stderr)
    sys.exit(1)


# Strandedness -> featureCounts -s value
STRANDEDNESS_MAP = {
    "unstranded": "0",
    "stranded": "1",
    "reverse_stranded": "2",
}


# =============================================================================
# Per-stage LSF resources — mirrors the #BSUB directives in each script
# =============================================================================
# The driver invokes bsub with the script as an ARGUMENT (bsub ... script.sh),
# which means LSF does NOT parse the embedded #BSUB directives in the script
# (those are only read when the script is piped via stdin: bsub < script.sh).
# So the per-stage resources must be re-declared on the bsub command line.
#
# Values below are transcribed EXACTLY from each script's #BSUB block. If you
# change a script's #BSUB -n/-M/-R/-W/-q, update the corresponding entry here
# too (and vice versa) so the two stay in sync.
#
# Keys are script basenames. Each value is a dict with:
#   queue: LSF queue
#   n:     cores (-n)
#   M:     memory in GB (-M)
#   R:     rusage string (-R "rusage[mem=N]")  [N matches M]
#   W:     walltime HH:MM (-W)
STAGE_RESOURCES = {
    # Per-sample stages
    "01_star_align.sh":            {"queue": "medium", "n": 16, "M": 60, "R": "rusage[mem=60]", "W": "24:00"},
    "02_salmon_expression.sh":     {"queue": "medium", "n": 16, "M": 32, "R": "rusage[mem=32]", "W": "16:00"},
    "03_salmon_alt_tss_polya.sh":  {"queue": "medium", "n": 16, "M": 32, "R": "rusage[mem=32]", "W": "16:00"},
    "04_regtools_junctions.sh":    {"queue": "medium", "n": 4,  "M": 16, "R": "rusage[mem=16]", "W": "8:00"},
    "05_featureCounts.sh":         {"queue": "medium", "n": 8,  "M": 16, "R": "rusage[mem=16]", "W": "8:00"},
    "06_rna_editing_pileup.sh":    {"queue": "medium", "n": 4,  "M": 16, "R": "rusage[mem=16]", "W": "8:00"},
    # Aggregation stages
    "10_aggregate_expression.sh":  {"queue": "medium", "n": 8,  "M": 32, "R": "rusage[mem=32]", "W": "8:00"},
    "11_aggregate_alt_tss_polya.sh": {"queue": "medium", "n": 8, "M": 32, "R": "rusage[mem=32]", "W": "8:00"},
    "12_aggregate_splicing.sh":    {"queue": "medium", "n": 8,  "M": 32, "R": "rusage[mem=32]", "W": "8:00"},
    "13_aggregate_intron_retention.sh": {"queue": "long", "n": 16, "M": 64, "R": "rusage[mem=64]", "W": "25:00"},
    "14_aggregate_rna_editing.sh": {"queue": "medium", "n": 8,  "M": 32, "R": "rusage[mem=32]", "W": "8:00"},
    "15_aggregate_stability.sh":   {"queue": "medium", "n": 8,  "M": 32, "R": "rusage[mem=32]", "W": "8:00"},
    # Final indexing
    "16_index_outputs.sh":         {"queue": "short",  "n": 2,  "M": 8,  "R": "rusage[mem=8]",  "W": "3:00"},
}


# =============================================================================
# subprocess helper — Python 3.6 compatible (no capture_output= / text= kwargs)
# =============================================================================
def _run_capture(cmd, timeout=None):
    """Run cmd and return (returncode, stdout_str, stderr_str).

    Equivalent to subprocess.run(cmd, capture_output=True, text=True,
    timeout=...) but works on Python 3.6 where capture_output and text are
    not available. Decodes output as UTF-8 with errors='replace'.
    """
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        # Re-raised as FileNotFoundError so callers' except clauses match.
        raise
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        out, err = proc.communicate()
    return (
        proc.returncode,
        out.decode("utf-8", errors="replace") if out else "",
        err.decode("utf-8", errors="replace") if err else "",
    )


# =============================================================================
# Job registry — tracks every submitted job with its semantic role + deps
# =============================================================================
class JobRegistry:
    """Tracks submitted jobs: job_id -> metadata (name, stage, sample,
    modality, depends_on). Used by the monitor to report which job failed and
    which downstream jobs are stranded."""

    def __init__(self):
        self.jobs = []  # list of dicts: {job_id, name, stage, sample, modality, depends_on}

    def add(self, job_id, name, stage, sample=None, modality=None, depends_on=None):
        if not job_id:  # dry-run or skipped
            return
        self.jobs.append({
            "job_id": job_id,
            "name": name,
            "stage": stage,
            "sample": sample,
            "modality": modality,
            "depends_on": list(depends_on) if depends_on else [],
        })

    @property
    def all_ids(self):
        return [j["job_id"] for j in self.jobs]

    def find(self, job_id):
        for j in self.jobs:
            if j["job_id"] == job_id:
                return j
        return None

    def transitive_dependents(self, failed_id):
        """All jobs that (transitively) depend on failed_id — i.e., their
        done() chain includes failed_id, so they can never start if failed_id
        EXITed. Returns list of job dicts."""
        # Build child map: parent_id -> [children]
        children = {}
        for j in self.jobs:
            for dep in j["depends_on"]:
                children.setdefault(dep, []).append(j)
        # BFS from failed_id
        stranded = []
        seen = set()
        queue = [failed_id]
        while queue:
            cur = queue.pop(0)
            for child in children.get(cur, []):
                if child["job_id"] not in seen:
                    seen.add(child["job_id"])
                    stranded.append(child)
                    queue.append(child["job_id"])
        return stranded


# =============================================================================
# Config + samples
# =============================================================================
def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_samples(samples_file: str) -> list:
    with open(samples_file) as f:
        return [line.strip() for line in f if line.strip()]


# =============================================================================
# Cohort input staging — symlink source data into ${COHORT_DIR}
# =============================================================================
def stage_cohort_inputs(cohort_cfg: dict, cohort_dir: str, dry_run: bool = False):
    """Symlink fastq_map, samples_file, and fastq_dir from their config-specified
    source paths into ${COHORT_DIR} so the stage scripts (which hardcode
    ${COHORT_DIR}/fastq_map.txt, /samples.txt, /fastq) can find them.

    The stage scripts do NOT read these paths from config or env — they hardcode
        FASTQ_MAP="${COHORT_DIR}/fastq_map.txt"
        FASTQ_DIR="${COHORT_DIR}/fastq"
        SAMPLES_FILE="${COHORT_DIR}/samples.txt"
    so the driver must stage the source data into those locations before
    submitting any jobs.

    For each target:
      - If the target already exists and is a symlink pointing to the source,
        skip (idempotent re-runs).
      - If the target exists but is NOT the expected symlink, error out (don't
        clobber real files/dirs the user may have placed deliberately).
      - If the target does not exist, create the symlink.
      - If the source path does not exist, error out.

    In dry-run mode, prints what it would do without creating anything.

    Args:
        cohort_cfg: The per-cohort config dict (must contain fastq_map,
                    fastq_dir, samples_file).
        cohort_dir: Absolute path to the cohort output directory
                    (${OUTPUT_BASE}/${cohort}).
        dry_run: If True, only print the planned symlinks.
    """
    # Map: config key -> (target name under cohort_dir, human label)
    # fastq_dir is a directory symlink; the other two are file symlinks.
    staging = [
        ("fastq_map",    "fastq_map.txt", "fastq map"),
        ("samples_file", "samples.txt",   "samples file"),
        ("fastq_dir",    "fastq",         "FASTQ directory"),
    ]

    print("=== Stage 0: Staging cohort inputs ===")
    for cfg_key, target_name, label in staging:
        source = cohort_cfg.get(cfg_key)
        if not source:
            print(f"  ERROR: config key '{cfg_key}' not set for this cohort.",
                  file=sys.stderr)
            sys.exit(1)

        target = os.path.join(cohort_dir, target_name)

        # Resolve to absolute paths for comparison (symlinks store what you
        # give them, but we compare resolved real paths to detect existing
        # correct links).
        source_abs = os.path.abspath(source)

        if dry_run:
            print(f"  [DRY-RUN] would symlink {target} -> {source_abs}")
            continue

        # Check source exists.
        if not os.path.exists(source_abs):
            print(f"  ERROR: {label} source does not exist: {source_abs}",
                  file=sys.stderr)
            print(f"         Check the '{cfg_key}' key in config.yml for this cohort.",
                  file=sys.stderr)
            sys.exit(1)

        # Target already exists?
        if os.path.islink(target) or os.path.exists(target):
            # Is it already the right symlink?
            if os.path.islink(target):
                existing_dest = os.readlink(target)
                existing_resolved = os.path.abspath(
                    os.path.join(os.path.dirname(target), existing_dest)
                )
                if existing_resolved == source_abs or existing_dest == source_abs:
                    print(f"  {label}: already linked -> {source_abs}")
                    continue
                # Wrong symlink — error, don't silently re-point.
                print(f"  ERROR: {target} exists as a symlink to "
                      f"{existing_dest}, which does not match config source "
                      f"{source_abs}.", file=sys.stderr)
                print(f"         Remove the existing link first if you want to "
                      f"re-stage.", file=sys.stderr)
                sys.exit(1)
            else:
                # Real file/dir, not a symlink — don't clobber.
                print(f"  ERROR: {target} already exists (not a symlink). "
                      f"Refusing to overwrite.", file=sys.stderr)
                print(f"         If this is stale, remove it first: "
                      f"rm -r {target}", file=sys.stderr)
                sys.exit(1)
        else:
            # Create the symlink.
            try:
                os.symlink(source_abs, target)
                print(f"  {label}: linked {target} -> {source_abs}")
            except OSError as e:
                print(f"  ERROR: could not create symlink {target} -> "
                      f"{source_abs}: {e}", file=sys.stderr)
                sys.exit(1)
    print()


# =============================================================================
# Job submission
# =============================================================================
def submit_job(script_path: str, args: list, queue: str, logs_dir: str,
               job_name: str, config_path: str, scripts_dir: str,
               dependency_jobs: list = None, dry_run: bool = False,
               registry: JobRegistry = None, stage: str = None,
               sample: str = None, modality: str = None) -> str:
    """Submit a BSUB job and return the job ID.

    Resources (-n, -M, -R, -W) and queue are taken from STAGE_RESOURCES keyed
    by the script basename. This is REQUIRED because the driver invokes bsub
    with the script as an argument, so LSF does not parse the embedded #BSUB
    directives in the script. The `queue` argument is used as a fallback only
    if the script is not in STAGE_RESOURCES.

    The script receives --config <config_path> and --scripts-dir <scripts_dir>
    as its first two arguments (before `args`), so it can load config.yml via
    config_get.py and locate its own directory.

    Args:
        script_path: Path to the .sh script to run.
        args: List of string arguments to pass to the script (cohort, sample,
              etc.). --config and --scripts-dir are prepended automatically.
        queue: Fallback LSF queue name (used only if script not in STAGE_RESOURCES).
        logs_dir: Directory for LSF logs.
        job_name: Job name for LSF.
        config_path: Path to config.yml (passed to script as --config).
        scripts_dir: Path to scripts directory (passed to script as --scripts-dir).
        dependency_jobs: List of job IDs this job depends on (must all be 'done').
        dry_run: If True, print the command instead of submitting.
        registry: JobRegistry to record this job into (with metadata).
        stage/sample/modality: Metadata for the registry / failure reports.

    Returns:
        Job ID string (e.g. "12345"), or empty string if dry_run.
    """
    script_name = os.path.basename(script_path)
    res = STAGE_RESOURCES.get(script_name)
    if res:
        queue = res["queue"]

    cmd = ["bsub"]

    if dependency_jobs:
        dep_expr = " && ".join(f"done({j})" for j in dependency_jobs)
        cmd.extend(["-w", dep_expr])

    cmd.extend([
        "-q", queue,
        "-J", job_name,
    ])

    # Per-stage resources from STAGE_RESOURCES (mirrors each script's #BSUB
    # block, which LSF ignores under argument-style invocation).
    if res:
        cmd.extend([
            "-n", str(res["n"]),
            "-M", str(res["M"]),
            "-R", res["R"],
            "-W", res["W"],
        ])

    # Prepend --config and --scripts-dir to the script args so each script can
    # load config.yml via config_get.py.
    full_args = [config_path, scripts_dir] + args

    cmd.extend([
        "-o", f"{logs_dir}/{job_name}.%J.out",
        "-e", f"{logs_dir}/{job_name}.%J.err",
        script_path,
    ] + full_args)

    if dry_run:
        dep_str = f" (depends on: {dependency_jobs})" if dependency_jobs else ""
        print(f"  [DRY-RUN] {' '.join(cmd)}{dep_str}")
        return ""

    rc, stdout, stderr = _run_capture(cmd)
    if rc != 0:
        print(f"ERROR submitting job: {' '.join(cmd)}", file=sys.stderr)
        print(f"  stdout: {stdout}", file=sys.stderr)
        print(f"  stderr: {stderr}", file=sys.stderr)
        sys.exit(1)

    output = stdout.strip()
    try:
        job_id = output.split("<")[1].split(">")[0]
    except IndexError:
        print(f"WARNING: Could not parse job ID from bsub output: {output}", file=sys.stderr)
        return ""

    print(f"  Submitted {job_name} -> job {job_id}")
    if registry is not None:
        registry.add(job_id, job_name, stage, sample, modality, dependency_jobs)
    return job_id


# =============================================================================
# Monitor — polls bjobs, detects failures, reports stranded dependents
# =============================================================================
def poll_bjobs(job_ids: list) -> dict:
    """Query LSF for the status of the given job IDs.

    Returns dict {job_id: STAT} **only for job_ids that were requested** —
    jobs from other runs or unrelated work are filtered out so the monitor
    does not report them as "not in registry". STAT is RUN/PEND/DONE/EXIT/
    PSUSP/USUSP/SSUSP. A registered job absent from bjobs output is treated
    as DONE (LSF clears finished jobs after a retention period).

    Uses `bjobs -a -o "jobid stat job_name"` (wide, parseable, one line per
    job). Falls back to parsing `bjobs -a` default columns by position if the
    -o form is unsupported.
    """
    if not job_ids:
        return {}
    wanted = set(job_ids)
    statuses = {}

    # Try the parseable -o form first.
    try:
        rc, stdout, _ = _run_capture(["bjobs", "-a", "-o", "jobid stat job_name"], timeout=30)
        if rc == 0:
            lines = stdout.strip().split("\n")
            # Header: JOBID STAT JOB_NAME
            for line in lines[1:]:
                parts = line.split(None, 2)
                if len(parts) >= 2:
                    jid, stat = parts[0], parts[1]
                    if jid in wanted:
                        statuses[jid] = stat
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        pass

    # Fallback: default bjobs -a columns (JOBID USER STAT QUEUE ...).
    if not statuses:
        try:
            rc, stdout, _ = _run_capture(["bjobs", "-a"], timeout=30)
            if rc == 0:
                lines = stdout.strip().split("\n")
                for line in lines[1:]:
                    parts = line.split()
                    if len(parts) >= 3:
                        jid, stat = parts[0], parts[2]
                        if jid in wanted:
                            statuses[jid] = stat
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            pass

    # Absent registered jobs -> DONE (cleared from system).
    for jid in job_ids:
        if jid not in statuses:
            statuses[jid] = "DONE"
    return statuses


def print_failure_report(failed_job: dict, stranded: list, statuses: dict, logs_dir: str):
    """Print a detailed failure report for an EXIT job + its stranded dependents."""
    bar = "!" * 64
    print("\n" + bar, file=sys.stderr)
    print(f"MONITOR: job {failed_job['job_id']} FAILED (EXIT)", file=sys.stderr)
    print(f"  name:    {failed_job['name']}", file=sys.stderr)
    print(f"  stage:   {failed_job['stage']}", file=sys.stderr)
    if failed_job.get("sample"):
        print(f"  sample:  {failed_job['sample']}", file=sys.stderr)
    if failed_job.get("modality"):
        print(f"  modality:{failed_job['modality']}", file=sys.stderr)
    # Log path: the submit_job -o pattern was {logs_dir}/{job_name}.%J.out
    print(f"  log:     {logs_dir}/{failed_job['name']}.{failed_job['job_id']}.err",
          file=sys.stderr)
    if stranded:
        print("Stranded downstream jobs (their done() dependency can never be "
              "satisfied):", file=sys.stderr)
        for s in stranded:
            stat = statuses.get(s["job_id"], "?")
            print(f"  - {s['name']:<45} (job {s['job_id']}, {stat})", file=sys.stderr)
    else:
        print("No downstream jobs depend on this job.", file=sys.stderr)
    print(bar, file=sys.stderr)
    print("Monitor exiting (first failure detected). Fix the failed job and "
          "re-run,", file=sys.stderr)
    print("or cancel stranded jobs manually with: bkill <jobid> ...",
          file=sys.stderr)


def run_monitor(registry: JobRegistry, index_job_id: str, logs_dir: str,
                poll_interval: int, monitor_all: bool) -> int:
    """Poll bjobs until the index job terminates. On failure, print report and
    exit (immediately unless --monitor-all). Returns exit code (0 success, 1
    failure)."""
    if not index_job_id:
        print("\n[MONITOR] No index job to monitor (dry-run or no jobs).")
        return 0

    print(f"\n[MONITOR] Polling bjobs every {poll_interval}s until index job "
          f"{index_job_id} terminates.")
    print(f"[MONITOR] {'Reporting all failures' if monitor_all else 'Exiting on first failure'}.")
    print(f"[MONITOR] Ctrl-C to stop monitoring (jobs continue running via LSF deps).")

    failed = set()
    final_status = None

    try:
        while True:
            statuses = poll_bjobs(registry.all_ids)

            # Detect newly-failed jobs (only our own — poll_bjobs filters).
            for jid, stat in statuses.items():
                if stat == "EXIT" and jid not in failed:
                    failed.add(jid)
                    fj = registry.find(jid)
                    if fj is None:
                        # Should not happen now that poll_bjobs filters to our
                        # job IDs, but guard defensively.
                        print(f"\n[MONITOR] job {jid} FAILED (EXIT) — not in registry "
                              f"(unexpected)", file=sys.stderr)
                        continue
                    stranded = registry.transitive_dependents(jid)
                    print_failure_report(fj, stranded, statuses, logs_dir)
                    if not monitor_all:
                        return 1

            # Check terminal state of the index job.
            idx_stat = statuses.get(index_job_id)
            if idx_stat in ("DONE", "EXIT") or index_job_id not in statuses:
                final_status = idx_stat if idx_stat in ("DONE", "EXIT") else "DONE"
                break

            time.sleep(poll_interval)
    except KeyboardInterrupt:
        print("\n[MONITOR] Interrupted by user. Submitted jobs continue running "
              "via LSF dependencies.", file=sys.stderr)
        return 130

    # Index job reached terminal state.
    print(f"\n[MONITOR] Index job {index_job_id} reached {final_status}.")
    if failed or final_status == "EXIT":
        print(f"[MONITOR] {len(failed)} job(s) failed during the run.", file=sys.stderr)
        return 1
    print(f"[MONITOR] All jobs completed successfully.")
    return 0


# =============================================================================
# Main
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="PANTRY pipeline driver for MDACC seadragon HPC (LSF/BSUB)"
    )
    parser.add_argument("--config", required=True, help="Path to config.yml")
    parser.add_argument("--cohort", required=True, help="Cohort name (must be in config)")
    parser.add_argument("--scripts-dir", default=None,
                        help="Directory containing the BSUB .sh scripts "
                             "(defaults to the directory of this script)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print jobs that would be submitted without submitting")
    parser.add_argument("--skip-alignment", action="store_true",
                        help="Skip STAR alignment (assume BAMs already exist)")
    parser.add_argument("--skip-rna-editing", action="store_true",
                        help="Skip RNA editing modality")
    # Monitoring flags
    parser.add_argument("--no-monitor", action="store_true",
                        help="Submit all jobs and exit immediately (do not poll bjobs)")
    parser.add_argument("--monitor-all", action="store_true",
                        help="Keep polling until the index job terminates even after "
                             "a failure is detected (default: exit on first failure)")
    parser.add_argument("--poll-interval", type=int, default=60,
                        help="Seconds between bjobs polls (default 60)")
    args = parser.parse_args()

    config = load_config(args.config)

    if args.cohort not in config.get("cohorts", {}):
        print(f"ERROR: Cohort '{args.cohort}' not found in config. "
              f"Available: {list(config.get('cohorts', {}).keys())}", file=sys.stderr)
        sys.exit(1)

    cohort_cfg = config["cohorts"][args.cohort]
    cohort = args.cohort

    scripts_dir = args.scripts_dir or str(Path(__file__).parent.resolve())

    output_base = config["output_base"]
    cohort_dir = f"{output_base}/{cohort}"
    logs_dir = f"{cohort_dir}/logs"
    queue = config.get("lsf", {}).get("default_queue", "medium")

    Path(logs_dir).mkdir(parents=True, exist_ok=True)

    # ---- Stage 0: stage cohort inputs (symlink source data into cohort_dir) ----
    # The stage scripts hardcode ${COHORT_DIR}/fastq_map.txt, /samples.txt, and
    # /fastq. The config points to the real source paths. Symlink them into
    # place before submitting any jobs. In dry-run, just print the plan.
    stage_cohort_inputs(cohort_cfg, cohort_dir, dry_run=args.dry_run)

    samples = load_samples(cohort_cfg["samples_file"])
    print(f"\nCohort: {cohort}")
    print(f"  Samples: {len(samples)}")
    print(f"  Strandedness: {cohort_cfg['strandedness']}")
    print(f"  Scripts dir: {scripts_dir}")
    print(f"  Logs dir: {logs_dir}")
    print(f"  Queue: {queue}")
    print(f"  Dry run: {args.dry_run}")
    print(f"  Monitor: {'off (--no-monitor)' if args.no_monitor else 'on'}")
    print()

    strandedness = cohort_cfg["strandedness"]
    if strandedness not in STRANDEDNESS_MAP:
        print(f"ERROR: Invalid strandedness '{strandedness}'. "
              f"Must be one of: {list(STRANDEDNESS_MAP.keys())}", file=sys.stderr)
        sys.exit(1)
    fc_strand = STRANDEDNESS_MAP[strandedness]

    registry = JobRegistry()

    # Common kwargs for all submit_job calls: config + scripts dir forwarding.
    submit_kwargs = dict(
        config_path=args.config,
        scripts_dir=scripts_dir,
        dry_run=args.dry_run,
        registry=registry,
    )

    # =========================================================================
    # Stage 1: Per-sample STAR alignment
    # =========================================================================
    print("=== Stage 1: STAR alignment ===")
    bam_jobs = {}  # sample -> job_id
    align_script = f"{scripts_dir}/01_star_align.sh"

    if args.skip_alignment:
        print("  [SKIPPED] --skip-alignment flag set")
    else:
        for sample in samples:
            job_id = submit_job(
                align_script, [cohort, sample], queue, logs_dir,
                f"align_{cohort}_{sample}", stage="alignment", sample=sample,
                **submit_kwargs,
            )
            bam_jobs[sample] = job_id
    print()

    # =========================================================================
    # Stage 2: Per-sample quantification (depends on BAM)
    # =========================================================================
    print("=== Stage 2: Per-sample quantification ===")

    # Salmon expression (per sample)
    salmon_expr_jobs = {}
    expr_script = f"{scripts_dir}/02_salmon_expression.sh"
    for sample in samples:
        deps = [bam_jobs[sample]] if bam_jobs.get(sample) else None
        job_id = submit_job(
            expr_script, [cohort, sample], queue, logs_dir,
            f"salmon_expr_{cohort}_{sample}", dependency_jobs=deps,
            stage="salmon_expression", sample=sample, modality="expression",
            **submit_kwargs,
        )
        salmon_expr_jobs[sample] = job_id

    # Salmon alt_TSS_polyA (per sample x 6 indices)
    salmon_alt_jobs = []
    alt_script = f"{scripts_dir}/03_salmon_alt_tss_polya.sh"
    for sample in samples:
        deps = [bam_jobs[sample]] if bam_jobs.get(sample) else None
        for group in ["grp_1", "grp_2"]:
            for position in ["upstream", "contained", "downstream"]:
                job_id = submit_job(
                    alt_script, [cohort, sample, group, position], queue, logs_dir,
                    f"salmon_alt_{cohort}_{sample}_{group}_{position}",
                    dependency_jobs=deps, stage="salmon_alt_tss_polya",
                    sample=sample, modality=f"alt_{group}_{position}",
                    **submit_kwargs,
                )
                salmon_alt_jobs.append(job_id)

    # RegTools junctions (per sample, depends on BAM)
    regtools_jobs = {}
    regtools_script = f"{scripts_dir}/04_regtools_junctions.sh"
    for sample in samples:
        deps = [bam_jobs[sample]] if bam_jobs.get(sample) else None
        job_id = submit_job(
            regtools_script, [cohort, sample], queue, logs_dir,
            f"regtools_{cohort}_{sample}", dependency_jobs=deps,
            stage="regtools", sample=sample, modality="splicing",
            **submit_kwargs,
        )
        regtools_jobs[sample] = job_id

    # featureCounts (per sample x 2 feature types, depends on BAM)
    featcounts_jobs = []
    fc_script = f"{scripts_dir}/05_featureCounts.sh"
    for sample in samples:
        deps = [bam_jobs[sample]] if bam_jobs.get(sample) else None
        for feature_type in ["exonic", "intronic"]:
            job_id = submit_job(
                fc_script, [cohort, sample, feature_type, fc_strand], queue, logs_dir,
                f"featcounts_{cohort}_{sample}_{feature_type}",
                dependency_jobs=deps, stage="featurecounts",
                sample=sample, modality=f"stability_{feature_type}",
                **submit_kwargs,
            )
            featcounts_jobs.append(job_id)

    # RNA editing pileup (per sample, depends on BAM)
    rnaedit_jobs = {}
    rnaedit_script = f"{scripts_dir}/06_rna_editing_pileup.sh"
    if args.skip_rna_editing:
        print("  [SKIPPED] RNA editing (--skip-rna-editing)")
    else:
        edit_sites = config.get("edit_sites_bed", "")
        if not edit_sites or not Path(edit_sites).exists():
            print(f"  [SKIPPED] Edit sites BED not found: {edit_sites}")
            print("           Run Tier 1 reference prep with REDIPORTAL_INPUT set.")
        else:
            for sample in samples:
                deps = [bam_jobs[sample]] if bam_jobs.get(sample) else None
                job_id = submit_job(
                    rnaedit_script, [cohort, sample], queue, logs_dir,
                    f"rnaedit_{cohort}_{sample}", dependency_jobs=deps,
                    stage="rna_editing", sample=sample, modality="RNA_editing",
                    **submit_kwargs,
                )
                rnaedit_jobs[sample] = job_id
    print()

    # =========================================================================
    # Stage 3: Per-modality aggregation (depends on all samples for that modality)
    # =========================================================================
    print("=== Stage 3: Per-modality aggregation ===")

    expr_dep_ids = [j for j in salmon_expr_jobs.values() if j]
    agg_expr_job = submit_job(
        f"{scripts_dir}/10_aggregate_expression.sh", [cohort], queue, logs_dir,
        f"agg_expr_{cohort}", dependency_jobs=expr_dep_ids,
        stage="aggregation", modality="expression", **submit_kwargs,
    )

    alt_dep_ids = [j for j in salmon_alt_jobs if j]
    agg_alt_job = submit_job(
        f"{scripts_dir}/11_aggregate_alt_tss_polya.sh", [cohort], queue, logs_dir,
        f"agg_alt_{cohort}", dependency_jobs=alt_dep_ids,
        stage="aggregation", modality="alt_TSS_polyA", **submit_kwargs,
    )

    splice_dep_ids = [j for j in regtools_jobs.values() if j]
    agg_splice_job = submit_job(
        f"{scripts_dir}/12_aggregate_splicing.sh", [cohort], queue, logs_dir,
        f"agg_splice_{cohort}", dependency_jobs=splice_dep_ids,
        stage="aggregation", modality="splicing", **submit_kwargs,
    )

    ir_dep_ids = [j for j in bam_jobs.values() if j]
    agg_ir_job = submit_job(
        f"{scripts_dir}/13_aggregate_intron_retention.sh", [cohort], "long", logs_dir,
        f"agg_ir_{cohort}", dependency_jobs=ir_dep_ids,
        stage="aggregation", modality="intron_retention", **submit_kwargs,
    )

    if rnaedit_jobs:
        rnaedit_dep_ids = [j for j in rnaedit_jobs.values() if j]
        agg_rnaedit_job = submit_job(
            f"{scripts_dir}/14_aggregate_rna_editing.sh", [cohort], queue, logs_dir,
            f"agg_rnaedit_{cohort}", dependency_jobs=rnaedit_dep_ids,
            stage="aggregation", modality="RNA_editing", **submit_kwargs,
        )
    else:
        agg_rnaedit_job = ""

    stab_dep_ids = [j for j in featcounts_jobs if j]
    agg_stab_job = submit_job(
        f"{scripts_dir}/15_aggregate_stability.sh", [cohort], queue, logs_dir,
        f"agg_stab_{cohort}", dependency_jobs=stab_dep_ids,
        stage="aggregation", modality="stability", **submit_kwargs,
    )
    print()

    # =========================================================================
    # Stage 4: Final indexing (depends on all aggregation jobs)
    # =========================================================================
    print("=== Stage 4: Final indexing ===")
    all_agg_ids = [j for j in [agg_expr_job, agg_alt_job, agg_splice_job,
                               agg_ir_job, agg_rnaedit_job, agg_stab_job] if j]
    index_job = submit_job(
        f"{scripts_dir}/16_index_outputs.sh", [cohort], "short", logs_dir,
        f"index_{cohort}", dependency_jobs=all_agg_ids,
        stage="indexing", modality="all", **submit_kwargs,
    )
    print()

    # =========================================================================
    # Summary
    # =========================================================================
    print("=" * 60)
    print(f"Pipeline submitted for cohort: {cohort}")
    print(f"  Samples: {len(samples)}")
    print(f"  Alignment jobs: {len(bam_jobs)}")
    print(f"  Salmon expression jobs: {len(salmon_expr_jobs)}")
    print(f"  Salmon alt_TSS_polyA jobs: {len(salmon_alt_jobs)}")
    print(f"  RegTools jobs: {len(regtools_jobs)}")
    print(f"  featureCounts jobs: {len(featcounts_jobs)}")
    print(f"  RNA editing jobs: {len(rnaedit_jobs)}")
    print(f"  Aggregation jobs: {len(all_agg_ids)}")
    print(f"  Final indexing job: {index_job or '(dry-run)'}")
    print(f"  Total jobs in registry: {len(registry.jobs)}")
    if not args.dry_run and index_job:
        print(f"\n  Final job ID: {index_job}")
        print(f"  Monitor: bjobs {index_job}")
    print("=" * 60)

    # =========================================================================
    # Stage 5: Monitoring (default) — poll bjobs until index terminates
    # =========================================================================
    if args.dry_run or args.no_monitor:
        if args.no_monitor and not args.dry_run:
            print("\n[MONITOR] --no-monitor set; driver exiting. Jobs will run "
                  "via LSF dependencies. Check with: bjobs -J \"*_<cohort>\"")
        return 0

    exit_code = run_monitor(
        registry, index_job, logs_dir,
        poll_interval=args.poll_interval, monitor_all=args.monitor_all,
    )
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
