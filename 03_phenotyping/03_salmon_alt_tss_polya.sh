#!/bin/bash
# =============================================================================
# 03_salmon_alt_tss_polya.sh — Per-sample Salmon quant for alt TSS/polyA
# =============================================================================
# Quantifies against one of the 6 txrevise Salmon indices (built in Tier 1).
# Salmon always uses -l A + --validateMappings.
#
# Usage: 03_salmon_alt_tss_polya.sh --config <config.yml> --scripts-dir <dir> <cohort> <sample_id> <group> <position>
#   group:    grp_1 or grp_2
#   position: upstream, contained, or downstream
#
# The driver submits 6 jobs per sample (2 groups x 3 positions).
# Outputs:
#   <cohort_dir>/intermediate/alt_TSS_polyA/<group>.<position>/<sample>/quant.sf
# =============================================================================

#BSUB -q medium
#BSUB -n 16
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 16:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/salmon_alt.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/salmon_alt.%J.err

set -eo pipefail

# ---- Config loading (Tier 2: --config and --scripts-dir as args 1-2) ----
CONFIG="${1:?Usage: --config <config.yml> --scripts-dir <dir> <cohort> <sample> <group> <position>}"
SCRIPTS_DIR="${2:?}"
shift 2
COHORT="${1:?Usage: <cohort> <sample> <group> <position>}"
SAMPLE="${2:?Usage: <cohort> <sample> <group> <position>}"
GROUP="${3:?Usage: group must be grp_1 or grp_2}"
POSITION="${4:?Usage: position must be upstream, contained, or downstream}"

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load all config values as shell vars
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}" --cohort "${COHORT}")"

# ---- Paths from config ----
OUTPUT_BASE="${OUTPUT_BASE}"
REFERENCE_DIR="${REFERENCE_DIR}"
SALMON_INDEX="${REFERENCE_DIR}/txrevise/txrevise.${GROUP}.${POSITION}.salmon_index"
COHORT_DIR="${OUTPUT_BASE}/${COHORT}"
INTERM_DIR="${COHORT_DIR}/intermediate"
FASTQ_MAP="${COHORT_DIR}/fastq_map.txt"
FASTQ_DIR="${COHORT_DIR}/fastq"

ALT_DIR="${INTERM_DIR}/alt_TSS_polyA/${GROUP}.${POSITION}"
OUT_DIR="${ALT_DIR}/${SAMPLE}"
mkdir -p "$OUT_DIR"

if [ ! -d "$SALMON_INDEX" ]; then
    echo "ERROR: Salmon index not found: $SALMON_INDEX" >&2
    echo "  Run Tier 1 reference prep (001/002/003) first." >&2
    exit 1
fi

echo "[$(date)] Salmon quant (alt_TSS_polyA) $SAMPLE $GROUP.$POSITION (cohort $COHORT)"

# ---- Parse FASTQ paths ----
MAP_LINE=$(grep -P "\t${SAMPLE}\s*$" "$FASTQ_MAP" | head -1)
if [ -z "$MAP_LINE" ]; then
    echo "ERROR: Sample $SAMPLE not found in $FASTQ_MAP" >&2
    exit 1
fi

NFIELDS=$(echo "$MAP_LINE" | awk -F'\t' '{print NF}')
if [ "$NFIELDS" -eq 3 ]; then
    R1="${FASTQ_DIR}/$(echo "$MAP_LINE" | awk -F'\t' '{print $1}')"
    R2="${FASTQ_DIR}/$(echo "$MAP_LINE" | awk -F'\t' '{print $2}')"
    SALMON_READS="-1 ${R1} -2 ${R2}"
    echo "  Paired-end: $R1 , $R2"
elif [ "$NFIELDS" -eq 2 ]; then
    R1="${FASTQ_DIR}/$(echo "$MAP_LINE" | awk -F'\t' '{print $1}')"
    SALMON_READS="-r ${R1}"
    echo "  Single-end: $R1"
else
    echo "ERROR: Unexpected field count in fastq_map for $SAMPLE" >&2
    exit 1
fi

# ---- Salmon quant ----
conda activate salmon-1.10.2

salmon quant \
    -i "$SALMON_INDEX" \
    -l A \
    --validateMappings \
    $SALMON_READS \
    -p 16 \
    -o "$OUT_DIR"

conda deactivate 2>/dev/null || true

echo "[$(date)] Done: $SAMPLE $GROUP.$POSITION"
echo "  Output: $OUT_DIR/quant.sf"
