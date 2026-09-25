#!/bin/bash
# 00_run_pipeline.sh — bsub wrapper that runs the PANTRY driver as a single
# LSF job. The driver submits all per-sample and aggregation jobs (with LSF
# -w done() dependencies) and, by default, monitors them by polling bjobs
# until the final index job terminates or a failure is detected.
#
# Submit with (inline env — relies on LSF default env propagation):
#   COHORT=cohort1 bsub < 00_run_pipeline.sh
#
# Or with explicit -env (robust if your site restricts default propagation):
#   bsub -env "COHORT=\"cohort1\"" < 00_run_pipeline.sh
#
# Required environment:
#   COHORT          Cohort name (must exist in config.yml)
#
# Optional environment:
#   CONFIG          Path to config.yml (default: /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/config.yml)
#   DRIVER_ARGS     Extra args passed to run_pipeline.py (e.g. "--no-monitor --poll-interval 30")
#   SCRIPTS_DIR     Directory containing the BSUB .sh scripts (default: dir of run_pipeline.py)
#   PYTHON          Python interpreter (default: python3)
#
# Wrapper resources: long queue, 2 cores, 8 GB, 120h walltime. The driver
# itself is lightweight (bsub submissions + bjobs polling) but must stay alive
# for the full pipeline duration when monitoring is on. 120h is the seadragon
# long-queue maximum and covers the critical path (alignment 8h + intron
# retention 24h + overhead).
#
# Wrapper log paths are hardcoded below. seadragon's LSF does not auto-create
# log directories, so create this dir once before first use:
#   mkdir -p /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs

#BSUB -J pantry_run
#BSUB -q long
#BSUB -n 2
#BSUB -M 8
#BSUB -W 120:00
#BSUB -o /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs/run_pipeline.%J.out
#BSUB -e /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs/run_pipeline.%J.err

# --- Required: COHORT ---
if [ -z "${COHORT:-}" ]; then
    echo "ERROR: COHORT environment variable is not set." >&2
    echo "       Submit with: COHORT=cohort1 bsub < 00_run_pipeline.sh" >&2
    echo "       Or:          bsub -env \"COHORT=\\\"cohort1\\\"\" < 00_run_pipeline.sh" >&2
    exit 1
fi

# --- Defaults ---
CONFIG="${CONFIG:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/config.yml}"
# Default SCRIPTS_DIR to this script's own directory (repo-clone layout);
# explicit SCRIPTS_DIR env var overrides.
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PYTHON="${PYTHON:-python3}"
DRIVER_ARGS="${DRIVER_ARGS:-}"

# Resolve the driver script relative to this wrapper.
DRIVER="${SCRIPTS_DIR}/run_pipeline.py"
if [ ! -f "${DRIVER}" ]; then
    echo "ERROR: run_pipeline.py not found at ${DRIVER}" >&2
    echo "       Set SCRIPTS_DIR to the directory containing run_pipeline.py" >&2
    exit 1
fi

if [ ! -f "${CONFIG}" ]; then
    echo "ERROR: config.yml not found at ${CONFIG}" >&2
    echo "       Set CONFIG to the path of your config.yml" >&2
    exit 1
fi

echo "============================================================"
echo "PANTRY pipeline driver (bsub wrapper)"
echo "  Cohort:      ${COHORT}"
echo "  Config:      ${CONFIG}"
echo "  Scripts dir: ${SCRIPTS_DIR}"
echo "  Driver:      ${DRIVER}"
echo "  Driver args: ${DRIVER_ARGS:-<none>}"
echo "  Python:      ${PYTHON}"
echo "  Started:     $(date)"
echo "  Host:        $(hostname)"
echo "  LSF job:     ${LSB_JOBID:-<not under LSF>}"
echo "============================================================"

# Run the driver. LSF propagates the parent job's environment to this wrapper,
# so COHORT/CONFIG/DRIVER_ARGS are visible here. The driver submits child jobs
# with their own -w done() dependencies; those run independently of this
# wrapper's lifetime, but the monitor (default) keeps this wrapper alive until
# the index job terminates so failures are reported.
"${PYTHON}" "${DRIVER}" \
    --config "${CONFIG}" \
    --cohort "${COHORT}" \
    --scripts-dir "${SCRIPTS_DIR}" \
    ${DRIVER_ARGS}

DRIVER_EXIT=$?

echo "============================================================"
echo "Driver finished with exit code ${DRIVER_EXIT}"
echo "  Ended: $(date)"
echo "============================================================"
exit ${DRIVER_EXIT}
