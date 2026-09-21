#!/bin/bash
# =============================================================================
# regenerate_refflat.sh — Regenerate HPLRv2.refFlat with the corrected 11-field
#                          format (adds transcript_id as the 'name' field).
#
# Usage:
#   bsub -q short -W 0:30 -env "CONFIG=...,SCRIPTS_DIR=..." bash regenerate_refflat.sh
#   # or interactively:
#   CONFIG=... SCRIPTS_DIR=... bash regenerate_refflat.sh
# =============================================================================

#BSUB -q short
#BSUB -n 2
#BSUB -M 8
#BSUB -R "rusage[mem=8]"
#BSUB -W 3:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/refflat.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/refflat.%J.err

set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"

source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}")"

REFERENCE_DIR="${REFERENCE_DIR}"
NORMALIZED_GTF="${NORMALIZED_GTF}"
REFFLAT="${REFERENCE_DIR}/HPLRv2.refFlat"

# Back up the old refFlat
if [ -f "$REFFLAT" ]; then
    cp "$REFFLAT" "${REFFLAT}.bak10field"
    echo "Backed up old refFlat to ${REFFLAT}.bak10field"
fi

# Use samtools env + MAJIQ venv for python3 with pandas
conda activate samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

echo "Regenerating refFlat from $NORMALIZED_GTF"
python3 "${SCRIPTS_DIR}/picard_qc.py" \
    --generate-refflat \
    --gtf "$NORMALIZED_GTF" \
    --output "$REFFLAT"

# Verify field count
NFIELDS=$(head -1 "$REFFLAT" | awk -F'\t' '{print NF}')
echo "Verification: $NFIELDS fields per line (expected 11)"
if [ "$NFIELDS" != "11" ]; then
    echo "ERROR: expected 11 fields, got $NFIELDS" >&2
    exit 1
fi

echo "Done: $REFFLAT"
