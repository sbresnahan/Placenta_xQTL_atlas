#!/bin/bash
# =============================================================================
# 02_txrevise_array.sh — Reference prep STAGE 2 of 3 (txrevise event array)
# =============================================================================
# Runs txrevise constructEvents.R as an LSF job array: one element per batch.
# Each element processes batch $LSB_JOBINDEX of $TXREVISE_N_BATCHES, i.e. it
# calls `constructEvents.R --batch "$LSB_JOBINDEX $TXREVISE_N_BATCHES"`. Because
# genes are processed independently, this replaces the old single-core serial
# run (which processed only batch 1) with N batches running in parallel.
#
# Depends on STAGE 1 (refprep_pre) via LSF job-name dependency. Feeds STAGE 3
# (03_refprep_post.sh), which merges the per-batch GFF3s and builds Salmon
# indices.
#
# Config loading: CONFIG and SCRIPTS_DIR env vars must be set via LSF -env at
# submission time:
#   bsub -env "CONFIG=\"config.yml\",SCRIPTS_DIR=\"/path/to/scripts\"" < 02_txrevise_array.sh
#
# IMPORTANT — array range must match the batch count:
#   The #BSUB -J "txrevise[1-100]" range below and txrevise_n_batches in
#   config.yml must be the same number. The #BSUB range is parsed at LSF
#   submission time (before the script body runs) so it CANNOT read from
#   config — it must be hardcoded here. If you change the batch count, edit
#   BOTH the #BSUB range here AND txrevise_n_batches in config.yml.
#
# No concurrency throttle: all elements are eligible to run at once ([1-100],
# no % limit). To throttle later, change the range to e.g. "txrevise[1-100]%25".
# =============================================================================

# ---- LSF directives ----
#BSUB -J "txrevise[1-100]"
#BSUB -q long
#BSUB -n 1
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 25:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/txrevise.%J.%I.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/txrevise.%J.%I.err

set -eo pipefail

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# ---- Config loading (Tier 1: CONFIG and SCRIPTS_DIR via LSF -env) ----
CONFIG="${CONFIG:?Set CONFIG via bsub -env}"
# Default SCRIPTS_DIR to this script's own directory (repo-clone layout).
# Explicit SCRIPTS_DIR overrides — REQUIRED for `bsub < script` submission
# (LSF executes a spool copy; self-location would resolve to the spool dir).
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Load all config values as shell vars
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}")"

# ---- Paths from config ----
REFERENCE_DIR="${REFERENCE_DIR}"
PANTRY_SCRIPTS="${PANTRY_SCRIPTS}"

# Scratch log dir (infra — hardcoded, used by #BSUB directives above)
LOG_DIR="/rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs"

# Singularity R invocation (infra — hardcoded)
SING_R="singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript"
export R_LIBS_USER="/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"

# Batch count from config (must match #BSUB -J "txrevise[1-N]" range above)
TXREVISE_N_BATCHES="${TXREVISE_N_BATCHES}"
TXREVISE_DIR="$REFERENCE_DIR/txrevise"

mkdir -p "$LOG_DIR" "$TXREVISE_DIR/batch"
cd "$REFERENCE_DIR"

# ---- Array index (this element's batch number) ----
BATCH_INDEX="${LSB_JOBINDEX:?This script must be submitted as an LSF job array (LSB_JOBINDEX is unset)}"

echo "[$(date)] STAGE 2 txrevise event construction — batch $BATCH_INDEX of $TXREVISE_N_BATCHES"

# ---- Require the STAGE 1 annotations ----
if [ ! -f "$TXREVISE_DIR/txrevise_annotations.rds" ]; then
    echo "ERROR: $TXREVISE_DIR/txrevise_annotations.rds not found."
    echo "  STAGE 1 (01_refprep_pre.sh) must complete before this array runs."
    exit 1
fi

# ---- Per-element resume: skip if this batch's 6 output GFF3s already exist ----
batch_done=true
for group in grp_1 grp_2; do
    for position in upstream contained downstream; do
        [ -f "$TXREVISE_DIR/batch/txrevise.${group}.${position}.${BATCH_INDEX}_${TXREVISE_N_BATCHES}.gff3" ] || batch_done=false
    done
done

if [ "$batch_done" = true ]; then
    echo "[$(date)] Batch $BATCH_INDEX SKIPPED — outputs already exist in $TXREVISE_DIR/batch"
else
    $SING_R "$PANTRY_SCRIPTS/txrevise/constructEvents.R" \
        --annot "$TXREVISE_DIR/txrevise_annotations.rds" \
        --batch "$BATCH_INDEX $TXREVISE_N_BATCHES" \
        --out "$TXREVISE_DIR/batch"
    echo "[$(date)] Batch $BATCH_INDEX done"
fi
