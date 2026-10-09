#!/bin/bash
# =============================================================================
# 53a_run_gxe_scan.sh — LSF worker for ancestry-specific tier-1 GxE scans.
# =============================================================================
# Required env: CONFIG, SCRIPTS_DIR, ANCESTRY, MODALITY, EXPOSURE,
# MODE (scan|merge). scan uses LSB_JOBINDEX as chromosome.
# =============================================================================
#BSUB -q medium
#BSUB -n 4
#BSUB -M 48G
#BSUB -R "rusage[mem=48G]"
#BSUB -W 12:00

set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR required}"
ANCESTRY="${ANCESTRY:?ERROR: ANCESTRY required}"
MODALITY="${MODALITY:?ERROR: MODALITY required}"
EXPOSURE="${EXPOSURE:?ERROR: EXPOSURE required}"
MODE="${MODE:-scan}"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
GXE_DIR="${GXE_DIR:-${RESULTS_DIR}/gxe}"
INPUTS="${GXE_DIR}/inputs"
TIER1="${GXE_DIR}/tier1"
CIS_WINDOW="${CIS_WINDOW:-1000000}"
MAF_THRESHOLD="${MAF_THRESHOLD:-0.01}"
MAF_THRESHOLD_INTERACTION="${MAF_THRESHOLD_INTERACTION:-0.05}"
PERM_BLOCKS="${PERM_BLOCKS:-100 400 500 9000}"
STOP_P="${STOP_P:-0.10}"
SEED="${SEED:-12345}"
CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
DOSAGES="${DOSAGES:-0}"
R_PACKAGE_LIB="${R_PACKAGE_LIB:-/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1}"
CONDA_EXE="${CONDA_EXE:-/risapps/rhel8/miniforge3/24.5.0-0/bin/conda}"
CONDA_ENV="${CONDA_ENV:-tensorqtl}"

source /etc/profile.d/modules.sh
[ -x "$CONDA_EXE" ] || { echo "ERROR: conda executable not found/executable: $CONDA_EXE"; exit 1; }
eval "$("$CONDA_EXE" shell.bash hook)"
conda activate "$CONDA_ENV"
command -v python3 >/dev/null || { echo "ERROR: python3 not found after conda activation"; exit 1; }

mkdir -p "$TIER1"
BASE="${TIER1}/${ANCESTRY}_${MODALITY}_${EXPOSURE}"

if [ "$MODE" = "merge" ]; then
    export QVALUE_RSCRIPT="${QVALUE_RSCRIPT:-$(cd "${SCRIPTS_DIR}/.." && pwd)/bin/Rscript_sif}"
    python3 "${SCRIPTS_DIR}/52_gxe_scan.py" --mode merge \
        --chr-outputs "${BASE}.chr*.gxe_cis.parquet" \
        --rscript "$QVALUE_RSCRIPT" \
        --r-package-lib "$R_PACKAGE_LIB" \
        --out "${BASE}.gxe_cis.parquet"
    exit 0
fi

IDX="${LSB_JOBINDEX:?scan mode must run under an LSF array}"
CHROM="$(echo $CHROMS | awk -v i="$IDX" '{print $i}')"
[ -n "$CHROM" ] || { echo "ERROR: LSB_JOBINDEX=$IDX outside CHROMS='$CHROMS'"; exit 1; }
OUT="${BASE}.chr${CHROM}.gxe_cis.parquet"
if [ -f "$OUT" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "  $OUT exists — skipping"; exit 0
fi

EXTRA=()
[ "$DOSAGES" = "1" ] && EXTRA+=(--dosages)

python3 "${SCRIPTS_DIR}/52_gxe_scan.py" --mode cis-perm \
    --bed "${INPUTS}/${ANCESTRY}_${MODALITY}.bed.gz" \
    --pgens "${ANCESTRY}=${QTL_DIR}/${ANCESTRY}_qtl" \
    --manifest "${INPUTS}/${ANCESTRY}_sample_manifest.tsv" \
    --covariates "${INPUTS}/${ANCESTRY}_covariates_${MODALITY}.tsv" \
    --exposures "${INPUTS}/exposures.tsv" \
    --exposure "$EXPOSURE" \
    --chrom "$CHROM" \
    --cis-window "$CIS_WINDOW" \
    --maf-threshold "$MAF_THRESHOLD" \
    --maf-threshold-interaction "$MAF_THRESHOLD_INTERACTION" \
    --perm-blocks $PERM_BLOCKS \
    --stop-p "$STOP_P" \
    --seed "$SEED" \
    "${EXTRA[@]}" \
    --out "$OUT"
