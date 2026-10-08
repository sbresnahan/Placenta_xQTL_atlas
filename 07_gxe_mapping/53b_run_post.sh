#!/bin/bash
# =============================================================================
# 53b_run_post.sh — LSF worker for the post-scan steps (dispatched by
# POST_MODE): tier2 | sensitivity | tnt | aggregate. Submitted by
# 53_submit_gxe.sh.
# =============================================================================
#BSUB -q medium
#BSUB -n 4
#BSUB -M 32G
#BSUB -R "rusage[mem=32G]"
#BSUB -W 8:00

set -euo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR required}"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
GXE_DIR="${GXE_DIR:-${RESULTS_DIR}/gxe}"
POST_MODE="${POST_MODE:?ERROR: POST_MODE required}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"

source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl

case "$POST_MODE" in
  tier2)
    MODALITY="${MODALITY:?tier2 needs MODALITY}"
    EXPOSURE="${EXPOSURE:?tier2 needs EXPOSURE}"
    python3 "${SCRIPTS_DIR}/54_tier2_stratified.py" \
        --qtl-dir "$QTL_DIR" --results-dir "$RESULTS_DIR" --gxe-dir "$GXE_DIR" \
        --ancestries $ANCESTRIES --modality "$MODALITY" --exposure "$EXPOSURE"
    ;;
  sensitivity)
    python3 "${SCRIPTS_DIR}/55_sensitivity_snpxcov.py" \
        --qtl-dir "$QTL_DIR" --results-dir "$RESULTS_DIR" --gxe-dir "$GXE_DIR" \
        --ancestries $ANCESTRIES ${EXPOSURES:+--exposures $EXPOSURES}
    ;;
  tnt)
    python3 "${SCRIPTS_DIR}/56_transmitted_nontransmitted.py" \
        --qtl-dir "$QTL_DIR" --results-dir "$RESULTS_DIR" --gxe-dir "$GXE_DIR" \
        ${MATERNAL_PGEN_DIR:+--maternal-pgen-dir "$MATERNAL_PGEN_DIR"} \
        ${PAIRS:+--pairs "$PAIRS"}
    ;;
  aggregate)
    python3 "${SCRIPTS_DIR}/57_aggregate_gxe.py" \
        --results-dir "$RESULTS_DIR" --gxe-dir "$GXE_DIR"
    ;;
  *)
    echo "ERROR: unknown POST_MODE=$POST_MODE"; exit 1;;
esac
