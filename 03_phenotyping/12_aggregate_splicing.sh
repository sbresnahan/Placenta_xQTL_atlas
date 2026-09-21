#!/bin/bash
# =============================================================================
# 12_aggregate_splicing.sh — leafCutter cluster + assemble + normalize splicing
# =============================================================================
# Collects per-sample .junc files (from RegTools), clusters them with
# leafCutter, then assembles into splicing.bed and normalizes.
#
# Usage: 12_aggregate_splicing.sh --config <config.yml> --scripts-dir <dir> <cohort>
# Outputs:
#   <cohort_dir>/intermediate/splicing/leafcutter_perind_numers.counts.gz
#   <cohort_dir>/output/unnorm/splicing.bed
#   <cohort_dir>/output/splicing.bed.gz (+ .tbi)
#   <cohort_dir>/output/splicing.phenotype_groups.txt
# =============================================================================

#BSUB -q medium
#BSUB -n 8
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 4:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/agg_splice.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/agg_splice.%J.err

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

SPLICE_DIR="${INTERM_DIR}/splicing"
JUNCFILE_LIST="${SPLICE_DIR}/juncfiles.txt"
LEAFCUTTER_OUT="${SPLICE_DIR}/leafcutter_perind_numers.counts.gz"
mkdir -p "$UNNORM_DIR" "$OUTPUT_DIR"

echo "[$(date)] Aggregating splicing (cohort $COHORT)"

# ---- Build juncfile list from samples ----
python3 -c "
import sys
with open('$SAMPLES_FILE') as f:
    for line in f:
        s = line.strip()
        if s:
            print(f'$SPLICE_DIR/{s}.junc')
" > "$JUNCFILE_LIST"

echo "  $(wc -l < "$JUNCFILE_LIST") junc files"

# ---- leafCutter cluster ----
# leafCutter activate path (infra — hardcoded)
source /rsrch5/home/epi/bhattacharya_lab/software/leafcutter/bin/activate

python3 "${PANTRY_SCRIPTS}/leafcutter_cluster_regtools_py3.py" \
    --juncfiles "$JUNCFILE_LIST" \
    --rundir "$SPLICE_DIR" \
    --maxintronlen 100000 \
    --minclureads 30 \
    --mincluratio 0.001

deactivate 2>/dev/null || true

# ---- Assemble splicing BED ----
conda activate samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

python3 "${PANTRY_SCRIPTS}/assemble_bed.py" splicing \
    --input "$LEAFCUTTER_OUT" \
    --ref-anno "$REF_ANNO" \
    --output "${UNNORM_DIR}/splicing.bed"

# ---- Normalize ----
python3 "${PANTRY_SCRIPTS}/normalize_phenotypes.py" \
    --input "${UNNORM_DIR}/splicing.bed" \
    --samples "$SAMPLES_FILE" \
    --output "${OUTPUT_DIR}/splicing.bed"

conda deactivate 2>/dev/null || true

# ---- bgzip + tabix + phenotype groups ----
conda activate samtools-1.16.1

bgzip "${OUTPUT_DIR}/splicing.bed"
tabix -p bed "${OUTPUT_DIR}/splicing.bed.gz"

zcat < "${OUTPUT_DIR}/splicing.bed.gz" \
    | tail -n +2 \
    | cut -f4 \
    | awk '{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }' \
    > "${OUTPUT_DIR}/splicing.phenotype_groups.txt"

conda deactivate 2>/dev/null || true

echo "[$(date)] Done: splicing"
ls -lh "${OUTPUT_DIR}/splicing.bed.gz" "${OUTPUT_DIR}/splicing.phenotype_groups.txt"
