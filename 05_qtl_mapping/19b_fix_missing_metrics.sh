#!/bin/bash
# =============================================================================
# 19b_fix_missing_metrics.sh — Re-run only the 3 broken Picard tools and merge
#                                results into existing sharded QC chunk files
# =============================================================================
# Fixes:
#   1. CollectAlignmentSummaryMetrics — empty READS_ALIGNED / PCT_READS_ALIGNED
#      (Picard 2.27 renamed to PF_READS_ALIGNED / PCT_PF_READS_ALIGNED)
#   2. CollectRnaSeqMetrics — was failing entirely (now with
#      VALIDATION_STRINGENCY=LENIENT + full stderr capture)
#   3. CollectGcBiasMetrics — was failing (needed R=<fasta>; now parses summary)
#
# Existing InsertSize and DupMetrics columns are preserved unchanged.
# Each array index processes one chunk file in-place.
#
# Usage:
#   bsub -J "fix_c1[1-16]" -env "CONFIG=...,SCRIPTS_DIR=...,COHORT=cohort1,NCHUNKS=16" < 19b_fix_missing_metrics.sh
#   bsub -J "fix_c2[1-16]" -env "...,COHORT=cohort2,NCHUNKS=16" < 19b_fix_missing_metrics.sh
#   bsub -J "fix_c3[1-16]" -env "...,COHORT=cohort3,NCHUNKS=16" < 19b_fix_missing_metrics.sh
#   bsub -J "fix_c4[1-16]" -env "...,COHORT=cohort4,NCHUNKS=16" < 19b_fix_missing_metrics.sh
#
# Required env vars: CONFIG, SCRIPTS_DIR, COHORT
# Optional env vars: NCHUNKS (default 16), QC_DIR, REFFLAT, FASTA
# =============================================================================

#BSUB -q medium
#BSUB -n 2
#BSUB -M 16
#BSUB -R "rusage[mem=16]"
#BSUB -W 12:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/fix.%J.%I.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/fix.%J.%I.err

set -eo pipefail

# ---- Required env vars ----
CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
# Default SCRIPTS_DIR to this script's own directory, so the pipeline runs
# directly from the git clone. An explicit SCRIPTS_DIR env var overrides —
# and is REQUIRED when submitting via `bsub < script` (LSF executes a spool
# copy of the script; self-location would resolve to the spool directory).
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# config_get.py lives in ../03_phenotyping in the repo layout; a flat copy in
# SCRIPTS_DIR (legacy deployment) takes precedence.
CONFIG_GET="${SCRIPTS_DIR}/config_get.py"
[ -f "$CONFIG_GET" ] || CONFIG_GET="${SCRIPTS_DIR}/../03_phenotyping/config_get.py"
if [ ! -f "$CONFIG_GET" ]; then
    echo "ERROR: config_get.py not found in $SCRIPTS_DIR or $SCRIPTS_DIR/../03_phenotyping" >&2
    echo "  Submitting via 'bsub <'? Export SCRIPTS_DIR=<repo>/05_qtl_mapping first." >&2
    exit 1
fi
COHORT="${COHORT:?ERROR: COHORT env var required}"
NCHUNKS="${NCHUNKS:-16}"
IDX="${LSB_JOBINDEX:-1}"

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load config values
eval "$(python3 "$CONFIG_GET" "${CONFIG}")"

OUTPUT_BASE="${OUTPUT_BASE}"
REFERENCE_DIR="${REFERENCE_DIR}"
NORMALIZED_GTF="${NORMALIZED_GTF}"
REF_GENOME="${REF_GENOME}"

QC_DIR="${QC_DIR:-${OUTPUT_BASE}/hcp/qc_metrics}"
REFFLAT="${REFFLAT:-${REFERENCE_DIR}/HPLRv2.refFlat}"
FASTA="${FASTA:-${REF_GENOME}}"

# ---- Environment: picard binary + java from picard env, python3 from MAJIQ venv ----
module load picard
conda activate picard-2.27.4
conda activate --stack samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

# ---- Locate this chunk's file ----
CHUNK_FILE="${QC_DIR}/${COHORT}.chunk${IDX}.qcmetrics.tsv"

if [ ! -f "$CHUNK_FILE" ]; then
    echo "[$(date)] Chunk file not found: $CHUNK_FILE — skipping (empty chunk from 19a)"
    exit 0
fi

if [ ! -s "$CHUNK_FILE" ]; then
    echo "[$(date)] Chunk file is empty: $CHUNK_FILE — skipping"
    exit 0
fi

BAM_DIR="${OUTPUT_BASE}/${COHORT}/intermediate/star_out"

echo "[$(date)] Fixing missing metrics: cohort=${COHORT} chunk=${IDX}/${NCHUNKS}"
echo "  CHUNK_FILE: ${CHUNK_FILE}"
echo "  BAM_DIR:    ${BAM_DIR}"
echo "  REFFLAT:    ${REFFLAT}"
echo "  FASTA:      ${FASTA}"

# ---- Run fix_missing_metrics.py ----
python3 "${SCRIPTS_DIR}/fix_missing_metrics.py" \
    --chunk-file "$CHUNK_FILE" \
    --bam-dir "$BAM_DIR" \
    --refflat "$REFFLAT" \
    --fasta "$FASTA" \
    --picard-cmd picard

echo "[$(date)] Done: cohort=${COHORT} chunk=${IDX}"
