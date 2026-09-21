#!/bin/bash
# =============================================================================
# 04_regtools_junctions.sh — Per-sample splice junction extraction
# =============================================================================
# Extracts junctions from the shrunk BAM using RegTools, for leafCutter.
# Uses -s XS (strand from STAR's XS tag — works for all library preps).
#
# Usage: 04_regtools_junctions.sh --config <config.yml> --scripts-dir <dir> <cohort> <sample_id>
# Outputs:
#   <cohort_dir>/intermediate/splicing/<sample>.junc
# =============================================================================

#BSUB -q medium
#BSUB -n 4
#BSUB -M 16
#BSUB -R "rusage[mem=16]"
#BSUB -W 4:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/regtools.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/regtools.%J.err

set -eo pipefail

# ---- Config loading (Tier 2: --config and --scripts-dir as args 1-2) ----
CONFIG="${1:?Usage: --config <config.yml> --scripts-dir <dir> <cohort> <sample>}"
SCRIPTS_DIR="${2:?}"
shift 2
COHORT="${1:?Usage: <cohort> <sample_id>}"
SAMPLE="${2:?Usage: <cohort> <sample_id>}"

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load all config values as shell vars
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}" --cohort "${COHORT}")"

# ---- Paths from config ----
OUTPUT_BASE="${OUTPUT_BASE}"
COHORT_DIR="${OUTPUT_BASE}/${COHORT}"
INTERM_DIR="${COHORT_DIR}/intermediate"

# RegTools binary (infra — hardcoded)
REGTOOLS="/rsrch5/home/epi/bhattacharya_lab/software/regtools/regtools/build/regtools"

BAM="${INTERM_DIR}/bam/${SAMPLE}.bam"
BAI="${INTERM_DIR}/bam/${SAMPLE}.bam.bai"
SPLICE_DIR="${INTERM_DIR}/splicing"
JUNC_OUT="${SPLICE_DIR}/${SAMPLE}.junc"
mkdir -p "$SPLICE_DIR"

if [ ! -f "$BAM" ] || [ ! -f "$BAI" ]; then
    echo "ERROR: BAM or BAI not found for $SAMPLE: $BAM" >&2
    exit 1
fi

echo "[$(date)] RegTools junctions for $SAMPLE (cohort $COHORT)"

# RegTools is a direct binary (no activation needed)
"$REGTOOLS" junctions extract \
    -a 8 \
    -m 50 \
    -M 500000 \
    -s XS \
    -o "$JUNC_OUT" \
    "$BAM"

echo "[$(date)] Done: $SAMPLE"
echo "  Junctions: $JUNC_OUT ($(wc -l < "$JUNC_OUT") lines)"
