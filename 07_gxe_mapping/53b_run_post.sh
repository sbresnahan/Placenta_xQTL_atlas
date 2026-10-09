#!/bin/bash
# =============================================================================
# 53b_run_post.sh — LSF worker for the post-scan steps (dispatched by
# POST_MODE): tier1_meta | tier2 | sensitivity | tnt | aggregate. Submitted by
# 53_submit_gxe.sh.
# Optional env: CONDA_EXE, CONDA_ENV.
# =============================================================================
#BSUB -q medium
#BSUB -n 4
#BSUB -M 32G
#BSUB -R "rusage[mem=32G]"
#BSUB -W 8:00

set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR required}"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
GXE_DIR="${GXE_DIR:-${RESULTS_DIR}/gxe}"
POST_MODE="${POST_MODE:?ERROR: POST_MODE required}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
CONDA_EXE="${CONDA_EXE:-/risapps/rhel8/miniforge3/24.5.0-0/bin/conda}"
CONDA_ENV="${CONDA_ENV:-tensorqtl}"
MAF_THRESHOLD="${MAF_THRESHOLD:-0.01}"
MAF_THRESHOLD_INTERACTION="${MAF_THRESHOLD_INTERACTION:-0.05}"

source /etc/profile.d/modules.sh
[ -x "$CONDA_EXE" ] || { echo "ERROR: conda executable not found/executable: $CONDA_EXE"; exit 1; }
eval "$("$CONDA_EXE" shell.bash hook)"
conda activate "$CONDA_ENV"
command -v python3 >/dev/null || { echo "ERROR: python3 not found after conda activation"; exit 1; }

case "$POST_MODE" in
  tier1_meta)
    MODALITY="${MODALITY:?tier1_meta needs MODALITY}"
    EXPOSURE="${EXPOSURE:?tier1_meta needs EXPOSURE}"
    python3 "${SCRIPTS_DIR}/54_tier1_meta.py" \
        --gxe-dir "$GXE_DIR" --ancestries $ANCESTRIES \
        --modality "$MODALITY" --exposure "$EXPOSURE"
    ;;
  tier2)
    MODALITY="${MODALITY:?tier2 needs MODALITY}"
    EXPOSURE="${EXPOSURE:?tier2 needs EXPOSURE}"
    python3 "${SCRIPTS_DIR}/54_tier2_stratified.py" \
        --qtl-dir "$QTL_DIR" --results-dir "$RESULTS_DIR" --gxe-dir "$GXE_DIR" \
        --ancestries $ANCESTRIES --modality "$MODALITY" --exposure "$EXPOSURE" \
        --maf-threshold "$MAF_THRESHOLD" \
        --maf-threshold-interaction "$MAF_THRESHOLD_INTERACTION"
    ;;
  sensitivity)
    python3 "${SCRIPTS_DIR}/55_sensitivity_snpxcov.py" \
        --qtl-dir "$QTL_DIR" --results-dir "$RESULTS_DIR" --gxe-dir "$GXE_DIR" \
        --ancestries $ANCESTRIES ${EXPOSURES:+--exposures $EXPOSURES}
    ;;
  tnt)
    python3 "${SCRIPTS_DIR}/56_transmitted_nontransmitted.py" \
        --qtl-dir "$QTL_DIR" --results-dir "$RESULTS_DIR" --gxe-dir "$GXE_DIR" \
        --ancestries $ANCESTRIES \
        ${MATERNAL_PGEN_DIR:+--maternal-pgen-dir "$MATERNAL_PGEN_DIR"} \
        ${PAIRS:+--pairs "$PAIRS"}
    ;;
  aggregate)
    python3 "${SCRIPTS_DIR}/57_aggregate_gxe.py" \
        --results-dir "$RESULTS_DIR" --gxe-dir "$GXE_DIR" --ancestries $ANCESTRIES
    ;;
  *)
    echo "ERROR: unknown POST_MODE=$POST_MODE"; exit 1;;
esac
