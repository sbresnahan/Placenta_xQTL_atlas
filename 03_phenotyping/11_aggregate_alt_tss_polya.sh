#!/bin/bash
# =============================================================================
# 11_aggregate_alt_tss_polya.sh — Assemble + normalize alt TSS & polyA
# =============================================================================
# Collects per-sample Salmon quant.sf from the 6 txrevise indices and
# assembles alt_TSS.bed (from upstream positions) and alt_polyA.bed (from
# downstream positions). Each uses grp_1 + grp_2.
#
# Usage: 11_aggregate_alt_tss_polya.sh --config <config.yml> --scripts-dir <dir> <cohort>
# Outputs:
#   <cohort_dir>/output/unnorm/alt_TSS.bed, alt_polyA.bed
#   <cohort_dir>/output/alt_TSS.bed.gz (+ .tbi), alt_polyA.bed.gz (+ .tbi)
#   <cohort_dir>/output/alt_TSS.phenotype_groups.txt, alt_polyA.phenotype_groups.txt
# =============================================================================

#BSUB -q medium
#BSUB -n 8
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 4:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/agg_alt.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/agg_alt.%J.err

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

ALT_DIR="${INTERM_DIR}/alt_TSS_polyA"
mkdir -p "$UNNORM_DIR" "$OUTPUT_DIR"

echo "[$(date)] Aggregating alt TSS + polyA (cohort $COHORT)"

conda activate samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

# ---- alt_TSS: uses upstream position (grp_1.upstream + grp_2.upstream) ----
python3 "${PANTRY_SCRIPTS}/assemble_bed.py" alt-tss-polya \
    --samples "$SAMPLES_FILE" \
    --group1-dir "${ALT_DIR}/grp_1.upstream" \
    --group2-dir "${ALT_DIR}/grp_2.upstream" \
    --ref-anno "$REF_ANNO" \
    --output "${UNNORM_DIR}/alt_TSS.bed"

python3 "${PANTRY_SCRIPTS}/normalize_phenotypes.py" \
    --input "${UNNORM_DIR}/alt_TSS.bed" \
    --samples "$SAMPLES_FILE" \
    --output "${OUTPUT_DIR}/alt_TSS.bed"

# ---- alt_polyA: uses downstream position (grp_1.downstream + grp_2.downstream) ----
python3 "${PANTRY_SCRIPTS}/assemble_bed.py" alt-tss-polya \
    --samples "$SAMPLES_FILE" \
    --group1-dir "${ALT_DIR}/grp_1.downstream" \
    --group2-dir "${ALT_DIR}/grp_2.downstream" \
    --ref-anno "$REF_ANNO" \
    --output "${UNNORM_DIR}/alt_polyA.bed"

python3 "${PANTRY_SCRIPTS}/normalize_phenotypes.py" \
    --input "${UNNORM_DIR}/alt_polyA.bed" \
    --samples "$SAMPLES_FILE" \
    --output "${OUTPUT_DIR}/alt_polyA.bed"

conda deactivate 2>/dev/null || true

# ---- bgzip + tabix + phenotype groups ----
conda activate samtools-1.16.1

for modality in alt_TSS alt_polyA; do
    bgzip "${OUTPUT_DIR}/${modality}.bed"
    tabix -p bed "${OUTPUT_DIR}/${modality}.bed.gz"

    zcat < "${OUTPUT_DIR}/${modality}.bed.gz" \
        | tail -n +2 \
        | cut -f4 \
        | awk '{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }' \
        > "${OUTPUT_DIR}/${modality}.phenotype_groups.txt"
done

conda deactivate 2>/dev/null || true

echo "[$(date)] Done: alt TSS + polyA"
ls -lh "${OUTPUT_DIR}/alt_TSS.bed.gz" "${OUTPUT_DIR}/alt_polyA.bed.gz" \
       "${OUTPUT_DIR}/alt_TSS.phenotype_groups.txt" "${OUTPUT_DIR}/alt_polyA.phenotype_groups.txt"
