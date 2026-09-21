#!/bin/bash
# =============================================================================
# 18b_deconvolve_cohort.sh — Build bulk TPM matrix + run MuSiC deconvolution
# =============================================================================
# Stage 2 of the deconvolution pipeline (per-cohort). Builds a gene-symbol x
# sample TPM matrix from Salmon quant.sf files, then runs MuSiC deconvolution
# against the pre-built reference SCE.
#
# This is a shell wrapper that sets up the seadragon environment:
#   - Python: conda activate samtools-1.16.1 + MAJIQ venv (for pandas/numpy)
#   - R:      singularity RStudio 4.3.1 + R_LIBS_USER (for MuSiC/TOAST)
#
# Usage: 18b_deconvolve_cohort.sh --config <config.yml> --scripts-dir <dir> \
#           <cohort> [--reference <path>] [--maternal-threshold <F>]
#
# Outputs (in <deconv_dir>):
#   <cohort>_bulk_tpm_symbols.tsv              (gene-symbol x sample TPM matrix)
#   <cohort>_cell_proportions_full.tsv         (samples x 27 cell types)
#   <cohort>_cell_proportions_collapsed.tsv    (samples x 8 cell types)
#   <cohort>_maternal_flag.tsv                 (sample, maternal_fraction, flag)
# =============================================================================

#BSUB -q medium
#BSUB -n 4
#BSUB -M 16
#BSUB -R "rusage[mem=16]"
#BSUB -W 1:00
#BSUB -o /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs/deconv_%J.out
#BSUB -e /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs/deconv_%J.err

set -eo pipefail

# ---- Parse arguments (PANTRY convention: --config and --scripts-dir first) ----
CONFIG=""
SCRIPTS_DIR=""
COHORT=""
REF_RDS=""
MATERNAL_THRESHOLD="0.10"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)             CONFIG="$2"; shift 2 ;;
        --scripts-dir)        SCRIPTS_DIR="$2"; shift 2 ;;
        --reference)          REF_RDS="$2"; shift 2 ;;
        --maternal-threshold) MATERNAL_THRESHOLD="$2"; shift 2 ;;
        --cohort)             COHORT="$2"; shift 2 ;;
        *)
            # First non-flag argument is the cohort (PANTRY positional convention)
            if [ -z "${COHORT}" ]; then
                COHORT="$1"; shift
            else
                echo "ERROR: Unknown argument: $1" >&2; exit 1
            fi
            ;;
    esac
done

if [ -z "${CONFIG}" ] || [ -z "${SCRIPTS_DIR}" ] || [ -z "${COHORT}" ]; then
    echo "ERROR: --config, --scripts-dir, and <cohort> are required" >&2
    echo "Usage: 18b_deconvolve_cohort.sh --config <cfg> --scripts-dir <dir> <cohort>" >&2
    exit 1
fi

# ---- Global init (seadragon environment) ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load config values as shell vars (config_get.py uses built-in parser)
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}" --cohort "${COHORT}")"

OUTPUT_BASE="${OUTPUT_BASE}"
NORMALIZED_GTF="${NORMALIZED_GTF}"
COHORT_DIR="${OUTPUT_BASE}/${COHORT}"
DECONV_DIR="${OUTPUT_BASE}/deconvolution"
LOG_DIR="${OUTPUT_BASE}/logs"
mkdir -p "${DECONV_DIR}" "${LOG_DIR}"

# Salmon expression dir (matches 02_salmon_expression.sh output convention)
SALMON_DIR="${COHORT_DIR}/intermediate/expression"
# Use config-provided SAMPLES_FILE (source path from cohort config),
# or fall back to the cohort-dir symlink (created by run_pipeline.py staging).
SAMPLES_FILE="${SAMPLES_FILE:-${COHORT_DIR}/samples.txt}"

# Default reference path
REF_RDS="${REF_RDS:-${DECONV_DIR}/placenta_music_reference.rds}"
BULK_TSV="${DECONV_DIR}/${COHORT}_bulk_tpm_symbols.tsv"

# ---- Singularity R invocation (infra — hardcoded) ----
SING_R="singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript"
export R_LIBS_USER="/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"

echo "============================================================"
echo "Stage 2: Deconvolution for cohort ${COHORT}"
echo "  Config:       ${CONFIG}"
echo "  Scripts dir:  ${SCRIPTS_DIR}"
echo "  Cohort:       ${COHORT}"
echo "  Samples:      ${SAMPLES_FILE}"
echo "  Salmon dir:   ${SALMON_DIR}"
echo "  GTF:          ${NORMALIZED_GTF}"
echo "  Reference:    ${REF_RDS}"
echo "  Bulk output:  ${BULK_TSV}"
echo "  Mat threshold:${MATERNAL_THRESHOLD}"
echo "  Started:      $(date)"
echo "  LSF job:      ${LSB_JOBID:-<not under LSF>}"
echo "============================================================"

# ---- Validate inputs ----
if [ ! -f "${REF_RDS}" ]; then
    echo "ERROR: Reference SCE not found: ${REF_RDS}" >&2
    echo "       Run 18a_build_reference.sh first." >&2
    exit 1
fi
if [ ! -f "${NORMALIZED_GTF}" ]; then
    echo "ERROR: Normalized GTF not found: ${NORMALIZED_GTF}" >&2
    exit 1
fi
if [ ! -f "${SAMPLES_FILE}" ]; then
    echo "ERROR: Samples file not found: ${SAMPLES_FILE}" >&2
    exit 1
fi
if [ ! -d "${SALMON_DIR}" ]; then
    echo "ERROR: Salmon expression dir not found: ${SALMON_DIR}" >&2
    exit 1
fi

# ---- Step 1: Build bulk TPM matrix (Python — needs pandas/numpy) ----
echo ""
echo "[$(date)] Step 1: Building bulk TPM matrix (Python)"

# Activate Python env for pandas/numpy
conda activate samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

python3 "${SCRIPTS_DIR}/build_bulk_tpm_matrix.py" \
    --gtf "${NORMALIZED_GTF}" \
    --salmon-dir "${SALMON_DIR}" \
    --samples "${SAMPLES_FILE}" \
    --output "${BULK_TSV}"

if [ $? -ne 0 ]; then
    echo "ERROR: build_bulk_tpm_matrix.py failed" >&2
    exit 1
fi

echo "[$(date)] Bulk TPM matrix: ${BULK_TSV}"
echo "  Size: $(ls -lh "${BULK_TSV}" | awk '{print $5}')"

# Deactivate Python env before R (singularity provides its own R)
conda deactivate 2>/dev/null || true

# ---- Step 2: Run MuSiC deconvolution (R — needs MuSiC/TOAST via singularity) ----
echo ""
echo "[$(date)] Step 2: Running MuSiC deconvolution (R via singularity)"

$SING_R "${SCRIPTS_DIR}/run_music_deconvolution.R" \
    --reference "${REF_RDS}" \
    --bulk "${BULK_TSV}" \
    --cohort "${COHORT}" \
    --output-dir "${DECONV_DIR}" \
    --maternal-threshold "${MATERNAL_THRESHOLD}"

if [ $? -ne 0 ]; then
    echo "ERROR: run_music_deconvolution.R failed" >&2
    exit 1
fi

echo ""
echo "============================================================"
echo "Stage 2 complete for cohort ${COHORT}."
echo "  Outputs in: ${DECONV_DIR}"
echo "    ${COHORT}_cell_proportions_full.tsv"
echo "    ${COHORT}_cell_proportions_collapsed.tsv"
echo "    ${COHORT}_maternal_flag.tsv"
echo "  Ended: $(date)"
echo "============================================================"
