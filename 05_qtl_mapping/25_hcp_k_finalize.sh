#!/bin/bash
# =============================================================================
# 25_hcp_k_finalize.sh — finalize isolated per-k HCP optimization jobs
# =============================================================================
# Required env:
#   SCRIPTS_DIR, REAL_QTL_DIR, HCP_FINALIZE_MODE, ANCESTRY
#   HCP_FINALIZE_MODE=expression or modality
# For modality mode: MODALITY is also required.
# Optional: K_GRID, FDR, CHR1_MIN
# =============================================================================
#BSUB -q medium
#BSUB -n 2
#BSUB -M 8
#BSUB -R "rusage[mem=8]"
#BSUB -W 4:00

set -euo pipefail

SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR required}"
REAL_QTL_DIR="${REAL_QTL_DIR:?ERROR: REAL_QTL_DIR required}"
HCP_FINALIZE_MODE="${HCP_FINALIZE_MODE:?ERROR: HCP_FINALIZE_MODE required}"
ANCESTRY="${ANCESTRY:?ERROR: ANCESTRY required}"
K_GRID="${K_GRID:-0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100}"
FDR="${FDR:-0.05}"
CHR1_MIN="${CHR1_MIN:-300}"

source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl

ARGS=(
  --mode "$HCP_FINALIZE_MODE"
  --qtl-dir "$REAL_QTL_DIR"
  --scripts-dir "$SCRIPTS_DIR"
  --ancestry "$ANCESTRY"
  --k-grid "$K_GRID"
  --fdr "$FDR"
  --chr1-min "$CHR1_MIN"
)

if [ "$HCP_FINALIZE_MODE" = "modality" ]; then
  MODALITY="${MODALITY:?ERROR: MODALITY required for modality finalization}"
  ARGS+=(--modality "$MODALITY")
fi

python3 "$SCRIPTS_DIR/finalize_hcp_k_grid.py" "${ARGS[@]}"
