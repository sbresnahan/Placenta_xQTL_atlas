#!/bin/bash
# =============================================================================
# 10_aggregate_expression.sh — Assemble + normalize expression & isoforms
# =============================================================================
# Collects per-sample Salmon quant.sf outputs and assembles them into
# expression.bed (gene-level) and isoforms.bed (isoform ratios), then
# quantile-normalizes both.
#
# Uses the MODIFIED assemble_bed.py that reads Salmon quant.sf format.
#
# QUANTIFICATION-UNCERTAINTY (QU) CORRECTION FOR ISOFORMS:
#   Isoform (transcript) counts are QU-corrected before assembly, per
#   Chen et al., NAR 2023 (doi:10.1093/nar/gkad1167). edgeR::catchSalmon
#   estimates per-transcript RTA overdispersion from the 20 Salmon bootstraps
#   (emitted by 02_salmon_expression.sh), and counts are divided by it. The
#   corrected counts are written as adjusted quant.sf to expression_qu/ and
#   assembled into isoforms.bed from there. Gene-level expression.bed is
#   assembled from the ORIGINAL Salmon dir (gene-level counts are largely
#   unambiguous, so RTA is negligible at gene level).
#
#   The QU step (qu_correct_salmon.R) requires edgeR, tximport, and
#   rtracklayer in the seadragon R library (R_LIBS_USER). rtracklayer is used
#   to build tx2gene from the reference GTF.
#
# Usage: 10_aggregate_expression.sh --config <config.yml> --scripts-dir <dir> <cohort>
# Outputs:
#   <cohort_dir>/intermediate/expression_qu/<sample>/quant.sf   (QU-corrected)
#   <cohort_dir>/output/unnorm/expression.bed
#   <cohort_dir>/output/unnorm/isoforms.bed
#   <cohort_dir>/output/unnorm/isoform_expression.bed
#   <cohort_dir>/output/expression.bed.gz (+ .tbi)
#   <cohort_dir>/output/isoforms.bed.gz (+ .tbi)
#   <cohort_dir>/output/isoforms.phenotype_groups.txt
#   <cohort_dir>/output/isoform_expression.bed.gz (+ .tbi)
#   <cohort_dir>/output/isoform_expression.phenotype_groups.txt
# =============================================================================

#BSUB -q medium
#BSUB -n 8
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 4:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/agg_expr.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/agg_expr.%J.err

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
# Helper scripts (assemble_bed.py, qu_correct_salmon.R) ship next to this
# script in the repo. The config's pantry_scripts / seadragon_scripts keys
# predate the git-clone deployment and may still point at the old
# flat-scripts dir (now the clone root) — self-locate; config is fallback.
_SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -f "${_SELF_DIR}/assemble_bed.py" ]; then
    PANTRY_SCRIPTS="${_SELF_DIR}"
else
    PANTRY_SCRIPTS="${PANTRY_SCRIPTS}"
fi
if [ -f "${_SELF_DIR}/qu_correct_salmon.R" ]; then
    SEADRAGON_SCRIPTS="${_SELF_DIR}"
else
    SEADRAGON_SCRIPTS="${SEADRAGON_SCRIPTS}"
fi
COHORT_DIR="${OUTPUT_BASE}/${COHORT}"
INTERM_DIR="${COHORT_DIR}/intermediate"
OUTPUT_DIR="${COHORT_DIR}/output"
UNNORM_DIR="${OUTPUT_DIR}/unnorm"
SAMPLES_FILE="${COHORT_DIR}/samples.txt"

EXPR_DIR="${INTERM_DIR}/expression"
EXPR_QU_DIR="${INTERM_DIR}/expression_qu"
mkdir -p "$UNNORM_DIR" "$OUTPUT_DIR"

# Singularity R invocation (infra — hardcoded) for the QU step.
SING_R="singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript"
export R_LIBS_USER="/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"

echo "[$(date)] Aggregating expression + isoforms (cohort $COHORT)"

# Use MAJIQ env for Python (pandas, numpy, gtfparse, scikit-learn)
conda activate samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

# ---- Assemble GENE-level BED from ORIGINAL Salmon dir (no QU correction) ----
# Gene-level counts are largely unambiguous (RTA is a transcript-level problem),
# so QU correction is not applied here.
python3 "${PANTRY_SCRIPTS}/assemble_bed.py" expression \
    --samples "$SAMPLES_FILE" \
    --input-dir "$EXPR_DIR" \
    --ref-anno "$REF_ANNO" \
    --output-expression "${UNNORM_DIR}/expression.bed"

# ---- QU correction for isoforms (edgeR::catchSalmon) ----
# Skip only when EVERY sample has an adjusted quant.sf (a first-sample-only
# check can miss a partially cleaned expression_qu; catchSalmon has no
# per-sample resume and would rerun all samples anyway).
N_SAMPLES=$(grep -c . "$SAMPLES_FILE")
N_QU_DONE=$(find "$EXPR_QU_DIR" -mindepth 2 -maxdepth 2 -name quant.sf 2>/dev/null | wc -l)
if [ "$N_SAMPLES" -gt 0 ] && [ "$N_QU_DONE" -eq "$N_SAMPLES" ]; then
    echo "[$(date)] QU correction SKIPPED — adjusted quant.sf for all $N_SAMPLES samples exist in ${EXPR_QU_DIR}"
else
    echo "[$(date)] QU correction (edgeR::catchSalmon) for isoforms ($N_QU_DONE/$N_SAMPLES adjusted quant.sf present)"
    # Run R outside the conda/MAJIQ env stack (singularity provides R + packages).
    conda deactivate 2>/dev/null || true
    $SING_R "${SEADRAGON_SCRIPTS}/qu_correct_salmon.R" \
        --samples "$SAMPLES_FILE" \
        --salmon-dir "$EXPR_DIR" \
        --ref-anno "$REF_ANNO" \
        --out-dir "$EXPR_QU_DIR"
    # Re-activate the Python env for assemble_bed.py.
    conda activate samtools-1.16.1
    source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate
fi

# ---- Assemble ISOFORM-level BEDs from QU-corrected dir ----
# isoforms.bed = within-gene usage ratios; isoform_expression.bed = the same
# QU-corrected transcript TPMs BEFORE the ratio division (abundance modality).
python3 "${PANTRY_SCRIPTS}/assemble_bed.py" expression \
    --samples "$SAMPLES_FILE" \
    --input-dir "$EXPR_QU_DIR" \
    --ref-anno "$REF_ANNO" \
    --output-isoforms "${UNNORM_DIR}/isoforms.bed" \
    --output-isoform-expr "${UNNORM_DIR}/isoform_expression.bed"

# ---- Canonical BEDs = unnorm (2026-09 schema) ----
# Per-cohort QN+INT is discontinued: normalization now happens once, after
# cross-cohort pooling, in stage 5 (19_hcp_factors.sh / 20_combat_modalities.sh).
# The canonical output/<modality>.bed is the unnorm BED so downstream paths
# (16_index_outputs.sh, combine_modalities.sh) are unchanged.
cp "${UNNORM_DIR}/expression.bed" "${OUTPUT_DIR}/expression.bed"
cp "${UNNORM_DIR}/isoforms.bed" "${OUTPUT_DIR}/isoforms.bed"
cp "${UNNORM_DIR}/isoform_expression.bed" "${OUTPUT_DIR}/isoform_expression.bed"

conda deactivate 2>/dev/null || true

# ---- bgzip + tabix ----
conda activate samtools-1.16.1

bgzip -f "${OUTPUT_DIR}/expression.bed"
tabix -p bed "${OUTPUT_DIR}/expression.bed.gz"

bgzip -f "${OUTPUT_DIR}/isoforms.bed"
tabix -p bed "${OUTPUT_DIR}/isoforms.bed.gz"

bgzip -f "${OUTPUT_DIR}/isoform_expression.bed"
tabix -p bed "${OUTPUT_DIR}/isoform_expression.bed.gz"

# ---- Phenotype groups for isoform modalities (gene grouping for tensorQTL) ----
zcat < "${OUTPUT_DIR}/isoforms.bed.gz" \
    | tail -n +2 \
    | cut -f4 \
    | awk '{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }' \
    > "${OUTPUT_DIR}/isoforms.phenotype_groups.txt"

zcat < "${OUTPUT_DIR}/isoform_expression.bed.gz" \
    | tail -n +2 \
    | cut -f4 \
    | awk '{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }' \
    > "${OUTPUT_DIR}/isoform_expression.phenotype_groups.txt"

conda deactivate 2>/dev/null || true

echo "[$(date)] Done: expression + isoforms + isoform_expression"
ls -lh "${OUTPUT_DIR}/expression.bed.gz" "${OUTPUT_DIR}/isoforms.bed.gz" \
       "${OUTPUT_DIR}/isoforms.phenotype_groups.txt" \
       "${OUTPUT_DIR}/isoform_expression.bed.gz" \
       "${OUTPUT_DIR}/isoform_expression.phenotype_groups.txt"
