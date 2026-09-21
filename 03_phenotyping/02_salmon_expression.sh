#!/bin/bash
# =============================================================================
# 02_salmon_expression.sh — Per-sample Salmon quant (expression + isoforms)
# =============================================================================
# Uses the pre-built Salmon transcriptome index (HPLRv2).
# Salmon always uses -l A (auto-detect library type) + --validateMappings.
# Emits 20 bootstrap resamples (--numBootstraps 20) so that the isoform
# expression modality can be QU-corrected (quantification-uncertainty / RTA
# overdispersion) downstream via edgeR::catchSalmon in 10_aggregate_expression.sh
# (per Chen et al., NAR 2023, doi:10.1093/nar/gkad1167). Bootstraps land in
# <OUT_DIR>/aux_info/ and are read by catchSalmon from the per-sample dir.
# Outputs quant.sf in the directory expected by assemble_bed.py expression mode.
#
# Usage: 02_salmon_expression.sh --config <config.yml> --scripts-dir <dir> <cohort> <sample_id>
# Outputs:
#   <cohort_dir>/intermediate/expression/<sample>/quant.sf
#   <cohort_dir>/intermediate/expression/<sample>/aux_info/   (20 bootstraps)
# =============================================================================

#BSUB -q medium
#BSUB -n 16
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 16:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/salmon_expr.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/salmon_expr.%J.err

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
SALMON_INDEX="${SALMON_INDEX}"
COHORT_DIR="${OUTPUT_BASE}/${COHORT}"
INTERM_DIR="${COHORT_DIR}/intermediate"
FASTQ_MAP="${COHORT_DIR}/fastq_map.txt"
FASTQ_DIR="${COHORT_DIR}/fastq"

EXPR_DIR="${INTERM_DIR}/expression"
OUT_DIR="${EXPR_DIR}/${SAMPLE}"
mkdir -p "$OUT_DIR"

echo "[$(date)] Salmon quant (expression) for $SAMPLE (cohort $COHORT)"

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
    --numBootstraps 20 \
    $SALMON_READS \
    -p 16 \
    -o "$OUT_DIR"

conda deactivate 2>/dev/null || true

echo "[$(date)] Done: $SAMPLE"
echo "  Output: $OUT_DIR/quant.sf"
head -5 "$OUT_DIR/quant.sf"
