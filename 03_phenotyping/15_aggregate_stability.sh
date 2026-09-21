#!/bin/bash
# =============================================================================
# 15_aggregate_stability.sh — Assemble + normalize RNA stability
# =============================================================================
# Collects per-sample featureCounts outputs (exonic + intronic) and computes
# the exon/intron read ratio = mRNA stability. Assembles into stability.bed
# and normalizes.
#
# Usage: 15_aggregate_stability.sh --config <config.yml> --scripts-dir <dir> <cohort>
# Outputs:
#   <cohort_dir>/output/unnorm/stability.bed
#   <cohort_dir>/output/stability.bed.gz (+ .tbi)
# =============================================================================

#BSUB -q medium
#BSUB -n 8
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 4:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/agg_stab.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/agg_stab.%J.err

set -eo pipefail

# ---- Config loading (Tier 2: --config and --scripts-dir as args 1-2) ----
CONFIG="${1:?Usage: --config <config.yml> --scripts-dir <dir> <cohort>}"
SCRIPTS_DIR="${2:?}"
shift 2
COHORT="${1:?Usage: <cohort>}"

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load all config values as shell vars
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}" --cohort "${COHORT}")"

# ---- Paths from config ----
OUTPUT_BASE="${OUTPUT_BASE}"
# Use the NORMALIZED GTF (has gene_name, gene_biotype, gene features required
# by assemble_bed.py), not the raw SQANTI3 GTF.
REF_ANNO="${NORMALIZED_GTF}"
PANTRY_SCRIPTS="${PANTRY_SCRIPTS}"
COHORT_DIR="${OUTPUT_BASE}/${COHORT}"
INTERM_DIR="${COHORT_DIR}/intermediate"
OUTPUT_DIR="${COHORT_DIR}/output"
UNNORM_DIR="${OUTPUT_DIR}/unnorm"
SAMPLES_FILE="${COHORT_DIR}/samples.txt"

STAB_DIR="${INTERM_DIR}/stability"
mkdir -p "$UNNORM_DIR" "$OUTPUT_DIR"

echo "[$(date)] Aggregating stability (cohort $COHORT)"

conda activate samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

# ---- Assemble stability BED (exon/intron ratio) ----
python3 "${PANTRY_SCRIPTS}/assemble_bed.py" stability \
    --samples "$SAMPLES_FILE" \
    --input-dir "$STAB_DIR" \
    --ref-anno "$REF_ANNO" \
    --output "${UNNORM_DIR}/stability.bed"

# ---- Normalize ----
python3 "${PANTRY_SCRIPTS}/normalize_phenotypes.py" \
    --input "${UNNORM_DIR}/stability.bed" \
    --samples "$SAMPLES_FILE" \
    --output "${OUTPUT_DIR}/stability.bed"

conda deactivate 2>/dev/null || true

# ---- bgzip + tabix ----
# Note: stability has one phenotype per gene, so no phenotype_groups file needed
conda activate samtools-1.16.1

bgzip "${OUTPUT_DIR}/stability.bed"
tabix -p bed "${OUTPUT_DIR}/stability.bed.gz"

conda deactivate 2>/dev/null || true

echo "[$(date)] Done: stability"
ls -lh "${OUTPUT_DIR}/stability.bed.gz"
