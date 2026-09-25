#!/bin/bash
# =============================================================================
# 14_aggregate_rna_editing.sh — Assemble + normalize RNA editing
# =============================================================================
# Collects per-sample edit level files, builds a shared site x sample matrix,
# prepares phenotypes (map to genes, optional clustering), assembles BED,
# and normalizes.
#
# Usage: 14_aggregate_rna_editing.sh --config <config.yml> --scripts-dir <dir> <cohort>
# Outputs:
#   <cohort_dir>/intermediate/RNA_editing/edit_site_matrix.tsv
#   <cohort_dir>/intermediate/RNA_editing/phenotype_matrix.tsv
#   <cohort_dir>/output/unnorm/RNA_editing.bed
#   <cohort_dir>/output/RNA_editing.bed.gz (+ .tbi)
#   <cohort_dir>/output/RNA_editing.phenotype_groups.txt
#   <cohort_dir>/output/RNA_editing.site_to_phenotype.tsv
# =============================================================================

#BSUB -q medium
#BSUB -n 8
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 4:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/agg_rnaedit.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/agg_rnaedit.%J.err

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
REFERENCE_DIR="${REFERENCE_DIR}"
# Use the NORMALIZED GTF (has gene_name, gene_biotype, gene features required
# by assemble_bed.py), not the raw SQANTI3 GTF.
REF_ANNO="${NORMALIZED_GTF}"
PANTRY_SCRIPTS="${PANTRY_SCRIPTS}"
COHORT_DIR="${OUTPUT_BASE}/${COHORT}"
INTERM_DIR="${COHORT_DIR}/intermediate"
OUTPUT_DIR="${COHORT_DIR}/output"
UNNORM_DIR="${OUTPUT_DIR}/unnorm"
SAMPLES_FILE="${COHORT_DIR}/samples.txt"

EDIT_SITES_TO_GENES="${REFERENCE_DIR}/edit_sites_to_genes.tsv"
RNA_EDIT_DIR="${INTERM_DIR}/RNA_editing"
EDIT_LEVELS_DIR="${RNA_EDIT_DIR}/edit_levels"
mkdir -p "$UNNORM_DIR" "$OUTPUT_DIR"

# RNA editing parameters (from config)
MIN_COVERAGE="${RNA_EDITING_EDIT_SITES_MIN_COVERAGE}"
MIN_SAMPLES="${RNA_EDITING_EDIT_SITES_MIN_SAMPLES}"
MIN_SAMPLES_FRACTION="${RNA_EDITING_EDIT_SITES_MIN_SAMPLES_FRACTION}"
CORRELATION_THRESHOLD="${RNA_EDITING_EDIT_SITES_CORRELATION_THRESHOLD}"

# Cluster flag: set to --cluster if edit_sites_cluster is true in config
CLUSTER_FLAG=""
if [ "${RNA_EDITING_EDIT_SITES_CLUSTER}" = "true" ]; then
    CLUSTER_FLAG="--cluster"
fi

# Compute effective min_samples: max(MIN_SAMPLES, ceil(fraction * n_samples))
N_SAMPLES=$(($(wc -l < "$SAMPLES_FILE")))
FRACTION_MIN_SAMPLES=$(python3 -c "import math; print(math.ceil(${MIN_SAMPLES_FRACTION} * ${N_SAMPLES}))")
EFFECTIVE_MIN_SAMPLES=$(python3 -c "print(min(max(${MIN_SAMPLES}, ${FRACTION_MIN_SAMPLES}), ${N_SAMPLES}))")
echo "  Effective min_samples: $EFFECTIVE_MIN_SAMPLES (n_samples=$N_SAMPLES)"

if [ ! -f "$EDIT_SITES_TO_GENES" ]; then
    echo "ERROR: edit_sites_to_genes.tsv not found: $EDIT_SITES_TO_GENES" >&2
    echo "  Run Tier 1 reference prep (001/002/003) first." >&2
    exit 1
fi

echo "[$(date)] Aggregating RNA editing (cohort $COHORT)"

conda activate samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

# ---- Build shared site x sample matrix ----
python3 "${PANTRY_SCRIPTS}/RNA_editing/shared_samples_sites_matrix.py" \
    --path_to_edit_files "$EDIT_LEVELS_DIR" \
    --samples_file "$SAMPLES_FILE" \
    --output_file "${RNA_EDIT_DIR}/edit_site_matrix.tsv" \
    --min_coverage "$MIN_COVERAGE" \
    --min_samples "$EFFECTIVE_MIN_SAMPLES"

# ---- Prepare phenotypes (map to genes, optional clustering) ----
python3 "${PANTRY_SCRIPTS}/RNA_editing/prepare_rna_editing_phenotypes.py" \
    --edit-matrix "${RNA_EDIT_DIR}/edit_site_matrix.tsv" \
    --site-to-gene "$EDIT_SITES_TO_GENES" \
    --output-matrix "${RNA_EDIT_DIR}/phenotype_matrix.tsv" \
    --output-site-map "${OUTPUT_DIR}/RNA_editing.site_to_phenotype.tsv" \
    --correlation-threshold "$CORRELATION_THRESHOLD" \
    $CLUSTER_FLAG

# ---- Assemble BED ----
python3 "${PANTRY_SCRIPTS}/assemble_bed.py" rna-editing \
    --input "${RNA_EDIT_DIR}/phenotype_matrix.tsv" \
    --ref-anno "$REF_ANNO" \
    --output "${UNNORM_DIR}/RNA_editing.bed"

# ---- Canonical BED = unnorm (2026-09 schema: normalization moved to stage 5) ----
cp "${UNNORM_DIR}/RNA_editing.bed" "${OUTPUT_DIR}/RNA_editing.bed"

conda deactivate 2>/dev/null || true

# ---- bgzip + tabix + phenotype groups ----
conda activate samtools-1.16.1

bgzip "${OUTPUT_DIR}/RNA_editing.bed"
tabix -p bed "${OUTPUT_DIR}/RNA_editing.bed.gz"

zcat < "${OUTPUT_DIR}/RNA_editing.bed.gz" \
    | tail -n +2 \
    | cut -f4 \
    | awk '{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }' \
    > "${OUTPUT_DIR}/RNA_editing.phenotype_groups.txt"

conda deactivate 2>/dev/null || true

echo "[$(date)] Done: RNA editing"
ls -lh "${OUTPUT_DIR}/RNA_editing.bed.gz" "${OUTPUT_DIR}/RNA_editing.phenotype_groups.txt" \
       "${OUTPUT_DIR}/RNA_editing.site_to_phenotype.tsv"
