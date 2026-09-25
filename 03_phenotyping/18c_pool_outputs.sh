#!/bin/bash
# =============================================================================
# 18c_pool_outputs.sh — Pool per-cohort deconvolution outputs + QC report
# =============================================================================
# Stage 3 of the deconvolution pipeline. Runs after all per-cohort deconvolution
# jobs complete. Pools the per-cohort proportion tables into all-cohorts tables
# and generates a QC summary report.
#
# This is a shell wrapper that sets up the seadragon Python environment
# (conda + MAJIQ venv for pandas/numpy) before calling
# pool_deconvolution_outputs.py.
#
# Usage: 18c_pool_outputs.sh --config <config.yml> --scripts-dir <dir> \
#           [--ancestry-map <path>]
#
# Outputs (in <deconv_dir>):
#   all_cohorts_cell_proportions_full.tsv       (pooled samples x 27 types)
#   all_cohorts_cell_proportions_collapsed.tsv  (pooled samples x 8 types)
#   all_cohorts_maternal_flag.tsv               (pooled maternal flags)
#   deconvolution_qc_report.md                  (QC summary report)
# =============================================================================

#BSUB -q short
#BSUB -n 2
#BSUB -M 8
#BSUB -R "rusage[mem=8]"
#BSUB -W 0:15
#BSUB -o /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs/deconv_pool.%J.out
#BSUB -e /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs/deconv_pool.%J.err

set -eo pipefail

# ---- Parse arguments (PANTRY convention: --config and --scripts-dir first) ----
CONFIG=""
SCRIPTS_DIR=""
ANCESTRY_MAP=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)        CONFIG="$2"; shift 2 ;;
        --scripts-dir)   SCRIPTS_DIR="$2"; shift 2 ;;
        --ancestry-map)  ANCESTRY_MAP="$2"; shift 2 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# Default --scripts-dir to this script's own directory (repo-clone layout)
if [ -z "${SCRIPTS_DIR}" ]; then
    SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

if [ -z "${CONFIG}" ]; then
    echo "ERROR: --config is required" >&2
    exit 1
fi

# ---- Global init (seadragon environment) ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load config values as shell vars
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}")"

OUTPUT_BASE="${OUTPUT_BASE}"
DECONV_DIR="${OUTPUT_BASE}/deconvolution"
LOG_DIR="${OUTPUT_BASE}/logs"
mkdir -p "${DECONV_DIR}" "${LOG_DIR}"

echo "============================================================"
echo "Stage 3: Pool deconvolution outputs + QC report"
echo "  Config:       ${CONFIG}"
echo "  Scripts dir:  ${SCRIPTS_DIR}"
echo "  Input/Output: ${DECONV_DIR}"
echo "  Ancestry map: ${ANCESTRY_MAP:-<none>}"
echo "  Started:      $(date)"
echo "  LSF job:      ${LSB_JOBID:-<not under LSF>}"
echo "============================================================"

# ---- Activate Python env for pandas/numpy ----
conda activate samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

# ---- Build command ----
POOL_ARGS=(
    --input-dir "${DECONV_DIR}"
    --output-dir "${DECONV_DIR}"
)

if [ -n "${ANCESTRY_MAP}" ]; then
    POOL_ARGS+=(--ancestry-map "${ANCESTRY_MAP}")
fi

# ---- Run pooling ----
python3 "${SCRIPTS_DIR}/pool_deconvolution_outputs.py" "${POOL_ARGS[@]}"

if [ $? -ne 0 ]; then
    echo "ERROR: pool_deconvolution_outputs.py failed" >&2
    exit 1
fi

conda deactivate 2>/dev/null || true

echo ""
echo "============================================================"
echo "Stage 3 complete."
echo "  Pooled outputs in: ${DECONV_DIR}"
echo "    all_cohorts_cell_proportions_full.tsv"
echo "    all_cohorts_cell_proportions_collapsed.tsv"
echo "    all_cohorts_maternal_flag.tsv"
echo "    deconvolution_qc_report.md"
echo "  Ended: $(date)"
echo "============================================================"
