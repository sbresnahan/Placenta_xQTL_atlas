#!/bin/bash
# =============================================================================
# 18a_build_reference.sh — Build MuSiC reference SCE from GEO GSE182381
# =============================================================================
# Stage 1 of the deconvolution pipeline. Downloads (if needed) the Campbell
# et al. 2023 scRNA-seq reference matrix from GEO and builds a
# SingleCellExperiment object with K pseudo-subjects for MuSiC.
#
# This is a shell wrapper that sets up the seadragon R environment (singularity
# + R_LIBS_USER) before calling build_music_reference.R.
#
# Usage: 18a_build_reference.sh --config <config.yml> --scripts-dir <dir> \
#           [--geo-ref <path>] [--output <path>] [--n-pseudo-subjects <N>] [--seed <N>]
#
# Outputs:
#   <deconv_dir>/placenta_music_reference.rds   (MuSiC SCE reference)
#   <deconv_dir>/reference_gene_symbols.txt     (16,003 gene symbols)
#   <deconv_dir>/GSE182381_reference_sample.txt.gz  (downloaded if absent)
# =============================================================================

#BSUB -q medium
#BSUB -n 4
#BSUB -M 16
#BSUB -R "rusage[mem=16]"
#BSUB -W 0:30
#BSUB -o /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs/deconv_refbuild.%J.out
#BSUB -e /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs/deconv_refbuild.%J.err

set -eo pipefail

# ---- Parse arguments (PANTRY convention: --config and --scripts-dir first) ----
CONFIG=""
SCRIPTS_DIR=""
GEO_REF=""
REF_OUTPUT=""
N_PSEUDO="5"
SEED="42"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)           CONFIG="$2"; shift 2 ;;
        --scripts-dir)      SCRIPTS_DIR="$2"; shift 2 ;;
        --geo-ref)          GEO_REF="$2"; shift 2 ;;
        --output)           REF_OUTPUT="$2"; shift 2 ;;
        --n-pseudo-subjects) N_PSEUDO="$2"; shift 2 ;;
        --seed)             SEED="$2"; shift 2 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [ -z "${CONFIG}" ] || [ -z "${SCRIPTS_DIR}" ]; then
    echo "ERROR: --config and --scripts-dir are required" >&2
    exit 1
fi

# ---- Global init (seadragon environment) ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load config values as shell vars (base conda python3 has no PyYAML, but
# config_get.py uses a built-in parser — no dependencies)
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}")"

OUTPUT_BASE="${OUTPUT_BASE}"
DECONV_DIR="${OUTPUT_BASE}/deconvolution"
LOG_DIR="${OUTPUT_BASE}/logs"
mkdir -p "${DECONV_DIR}" "${LOG_DIR}"

# Default paths
GEO_REF="${GEO_REF:-${DECONV_DIR}/GSE182381_reference_sample.txt.gz}"
REF_OUTPUT="${REF_OUTPUT:-${DECONV_DIR}/placenta_music_reference.rds}"
REF_GENES="${DECONV_DIR}/reference_gene_symbols.txt"

# ---- Singularity R invocation (infra — hardcoded) ----
SING_R="singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript"
export R_LIBS_USER="/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"

echo "============================================================"
echo "Stage 1: Build MuSiC reference SCE"
echo "  Config:       ${CONFIG}"
echo "  Scripts dir:  ${SCRIPTS_DIR}"
echo "  GEO ref:      ${GEO_REF}"
echo "  Output:       ${REF_OUTPUT}"
echo "  Pseudo-subj:  ${N_PSEUDO}"
echo "  Seed:         ${SEED}"
echo "  Started:      $(date)"
echo "  LSF job:      ${LSB_JOBID:-<not under LSF>}"
echo "============================================================"

# ---- Download GEO reference if not present ----
if [ ! -f "${GEO_REF}" ]; then
    echo "[$(date)] Downloading GSE182381_reference_sample.txt.gz from GEO..."
    curl -sL -o "${GEO_REF}" \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182381/suppl/GSE182381_reference_sample.txt.gz"
    if [ $? -ne 0 ] || [ ! -s "${GEO_REF}" ]; then
        echo "ERROR: Failed to download GEO reference" >&2
        exit 1
    fi
    echo "  Downloaded: $(ls -lh "${GEO_REF}" | awk '{print $5}')"
fi

# ---- Extract reference gene symbols for overlap reporting ----
if [ ! -f "${REF_GENES}" ]; then
    echo "[$(date)] Extracting reference gene symbols..."
    zcat "${GEO_REF}" | tail -n +2 | cut -f1 > "${REF_GENES}"
    echo "  $(wc -l < "${REF_GENES}") gene symbols"
fi

# ---- Build MuSiC reference SCE ----
echo "[$(date)] Building MuSiC reference SCE..."
$SING_R "${SCRIPTS_DIR}/build_music_reference.R" \
    --input "${GEO_REF}" \
    --output "${REF_OUTPUT}" \
    --n-pseudo-subjects "${N_PSEUDO}" \
    --seed "${SEED}"

if [ $? -ne 0 ]; then
    echo "ERROR: build_music_reference.R failed" >&2
    exit 1
fi

echo "[$(date)] Reference SCE built: ${REF_OUTPUT}"
echo "  Size: $(ls -lh "${REF_OUTPUT}" | awk '{print $5}')"
echo "============================================================"
echo "Stage 1 complete."
echo "  Ended: $(date)"
echo "============================================================"
