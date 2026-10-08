#!/bin/bash
# =============================================================================
# 53a_run_gxe_scan.sh — LSF worker for tier-1 GxE scans (one chromosome) or
# the per-modality merge. Submitted by 53_submit_gxe.sh.
# =============================================================================
# Required env: CONFIG, SCRIPTS_DIR, MODALITY, EXPOSURE, MODE (scan|merge)
# scan:  $LSB_JOBINDEX = chromosome
# merge: concatenates chr parquets + Storey q-values
# Optional: OUTPUT_BASE, QTL_DIR, RESULTS_DIR, GXE_DIR, ANCESTRIES,
#   CIS_WINDOW, MAF_THRESHOLD, PERM_BLOCKS, STOP_P, SEED, DOSAGES=1,
#   R_PACKAGE_LIB, CONDA_EXE, CONDA_ENV.
# =============================================================================
#BSUB -q medium
#BSUB -n 4
#BSUB -M 48G
#BSUB -R "rusage[mem=48G]"
#BSUB -W 12:00

set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR required}"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
GXE_DIR="${GXE_DIR:-${RESULTS_DIR}/gxe}"
INPUTS="${GXE_DIR}/inputs"
TIER1="${GXE_DIR}/tier1"
MODALITY="${MODALITY:?ERROR: MODALITY required}"
EXPOSURE="${EXPOSURE:?ERROR: EXPOSURE required}"
MODE="${MODE:-scan}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
CIS_WINDOW="${CIS_WINDOW:-1000000}"
MAF_THRESHOLD="${MAF_THRESHOLD:-0.01}"
PERM_BLOCKS="${PERM_BLOCKS:-100 400 500 9000}"
STOP_P="${STOP_P:-0.10}"
SEED="${SEED:-12345}"
DOSAGES="${DOSAGES:-0}"
R_PACKAGE_LIB="${R_PACKAGE_LIB:-/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1}"
CONDA_EXE="${CONDA_EXE:-/risapps/rhel8/miniforge3/24.5.0-0/bin/conda}"
CONDA_ENV="${CONDA_ENV:-tensorqtl}"

source /etc/profile.d/modules.sh
[ -x "$CONDA_EXE" ] || { echo "ERROR: conda executable not found/executable: $CONDA_EXE"; exit 1; }
eval "$("$CONDA_EXE" shell.bash hook)"
conda activate "$CONDA_ENV"
command -v python3 >/dev/null || { echo "ERROR: python3 not found after conda activation"; exit 1; }

PGENS=()
for ANC in $ANCESTRIES; do
    PGENS+=("${ANC}=${QTL_DIR}/${ANC}_qtl")
done

mkdir -p "$TIER1"
BASE="${TIER1}/pooled_${MODALITY}_${EXPOSURE}"

if [ "$MODE" = "merge" ]; then
    export QVALUE_RSCRIPT="${QVALUE_RSCRIPT:-$(cd "${SCRIPTS_DIR}/.." && pwd)/bin/Rscript_sif}"
    python3 "${SCRIPTS_DIR}/52_gxe_scan.py" --mode merge \
        --chr-outputs "${BASE}.chr*.gxe_cis.parquet" \
        --rscript "$QVALUE_RSCRIPT" \
        --r-package-lib "$R_PACKAGE_LIB" \
        --out "${BASE}.gxe_cis.parquet"
    exit 0
fi

CHROM="${LSB_JOBINDEX:?scan mode must run under an LSF array}"
OUT="${BASE}.chr${CHROM}.gxe_cis.parquet"
if [ -f "$OUT" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "  $OUT exists — skipping"; exit 0
fi

EXTRA=()
[ "$DOSAGES" = "1" ] && EXTRA+=(--dosages)

python3 "${SCRIPTS_DIR}/52_gxe_scan.py" --mode cis-perm \
    --bed "${INPUTS}/pooled_${MODALITY}.bed.gz" \
    --pgens "${PGENS[@]}" \
    --manifest "${INPUTS}/pooled_sample_manifest.tsv" \
    --covariates "${INPUTS}/pooled_covariates_${MODALITY}.tsv" \
    --exposures "${INPUTS}/exposures.tsv" \
    --exposure "$EXPOSURE" \
    --chrom "$CHROM" \
    --cis-window "$CIS_WINDOW" \
    --maf-threshold "$MAF_THRESHOLD" \
    --perm-blocks $PERM_BLOCKS \
    --stop-p "$STOP_P" \
    --seed "$SEED" \
    "${EXTRA[@]}" \
    --out "$OUT"
