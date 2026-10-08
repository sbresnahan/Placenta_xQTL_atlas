#!/bin/bash
# =============================================================================
# 50a_run_build_gxe_inputs.sh — LSF worker for Module-07 Stage 0
# =============================================================================
# Initializes the Seadragon shell/module environment and the conda environment
# required by 50_build_gxe_inputs.py. Do not rely on the submit shell's active
# conda environment; LSF jobs start with a clean/non-interactive shell context.
#
# Required env: CONFIG, SCRIPTS_DIR, METADATA
# Optional: OUTPUT_BASE, QTL_DIR, RESULTS_DIR, GXE_DIR, ANCESTRIES, MODALITIES,
#   CONDA_EXE, CONDA_ENV.
# =============================================================================
set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
METADATA="${METADATA:?ERROR: METADATA env var required}"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
GXE_DIR="${GXE_DIR:-${RESULTS_DIR}/gxe}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
MODALITIES="${MODALITIES:-expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability}"
CONDA_EXE="${CONDA_EXE:-/risapps/rhel8/miniforge3/24.5.0-0/bin/conda}"
CONDA_ENV="${CONDA_ENV:-tensorqtl}"

source /etc/profile.d/modules.sh
[ -x "$CONDA_EXE" ] || { echo "ERROR: conda executable not found/executable: $CONDA_EXE"; exit 1; }
eval "$("$CONDA_EXE" shell.bash hook)"
conda activate "$CONDA_ENV"

command -v python3 >/dev/null || { echo "ERROR: python3 not found after conda activation"; exit 1; }
command -v bgzip >/dev/null || { echo "ERROR: bgzip not found after conda activation"; exit 1; }
command -v tabix >/dev/null || { echo "ERROR: tabix not found after conda activation"; exit 1; }

python3 - <<'PYENV'
import numpy
import pandas
print("Stage-0 Python environment OK")
print("numpy:", numpy.__version__)
print("pandas:", pandas.__version__)
PYENV

exec python3 "${SCRIPTS_DIR}/50_build_gxe_inputs.py" \
    --qtl-dir "$QTL_DIR" \
    --results-dir "$RESULTS_DIR" \
    --metadata "$METADATA" \
    --ancestries "$ANCESTRIES" \
    --modalities "$MODALITIES" \
    --gxe-dir "$GXE_DIR"
