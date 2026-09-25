#!/bin/bash
# =============================================================================
# 25a_optimize_hcp.sh — LSF wrapper for chr1-only HCP-count optimization
# =============================================================================
# Pre-mapping module (devBrain xQTL convention — Wen et al., Science 2024,
# 384:eadh0829, §4.2): selects the number of HCP hidden covariates per
# ancestry by maximizing chr1 cis-eGene discovery (Storey q <= 0.05), with
# HCP re-estimated at each candidate k. Runs AFTER 23/24 and BEFORE the
# canonical 25_build_covariates.py run; installs the winning k* solution as
# ${QTL_DIR}/{ANC}_hcp_factors_harmonized.tsv.
#
# Usage:
#   export CONFIG=/path/to/config.yml
#   export SCRIPTS_DIR=/path/to/scripts
#   export ANCESTRY_MAP=/path/to/pooled_sample_ancestry_RNAseq.tsv
#   export ANCESTRIES="EAS EUR"        # optional; default: all in map
#   bsub < 25a_optimize_hcp.sh
#
# Or interactively (inside: bsub -Is -q medium -n 4 -M 32 -R "rusage[mem=32]" -W 12:00 bash):
#   CONFIG=config.yml SCRIPTS_DIR=/path/to/scripts \
#   ANCESTRY_MAP=/path/to/pooled_sample_ancestry_RNAseq.tsv \
#   bash 25a_optimize_hcp.sh
#
# Optional env vars:
#   ANCESTRIES   — space-separated ancestry labels (default: all in map;
#                  ANCESTRY singular accepted as a fallback alias)
#   K_GRID       — candidate HCP counts (default: "0 5 10 15 20 25 30")
#   QTL_DIR      — canonical QTL inputs dir (default: ${OUTPUT_BASE}/qtl_inputs)
#   HCP_DIR      — script-19 output dir (default: ${OUTPUT_BASE}/hcp)
#   PC_DIR       — genotype PCs dir (default: ${OUTPUT_BASE}/genotype_pcs)
#   WORK_DIR     — staging dir (default: ${QTL_DIR}/hcp_optimization)
#   FDR          — Storey q threshold for eGene counts (default: 0.05)
#   SKIP_EXISTING — set to 1 to reuse existing per-k results (resumable)
#   EXCLUDE_COVARIATES — covariates excluded before correlation pruning in
#                  every per-k model (default: ct_Maternal; "" disables)
# =============================================================================

#BSUB -q long
#BSUB -n 4
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 48:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/hcp_opt.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/hcp_opt.%J.err

set -eo pipefail

# ---- Config / env ----
CONFIG="${CONFIG:?ERROR: CONFIG env var required (path to config.yml)}"
# Default SCRIPTS_DIR to this script's own directory, so the pipeline runs
# directly from the git clone. An explicit SCRIPTS_DIR env var overrides —
# and is REQUIRED when submitting via `bsub < script` (LSF executes a spool
# copy of the script; self-location would resolve to the spool directory).
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# config_get.py lives in ../03_phenotyping in the repo layout; a flat copy in
# SCRIPTS_DIR (legacy deployment) takes precedence.
CONFIG_GET="${SCRIPTS_DIR}/config_get.py"
[ -f "$CONFIG_GET" ] || CONFIG_GET="${SCRIPTS_DIR}/../03_phenotyping/config_get.py"
if [ ! -f "$CONFIG_GET" ]; then
    echo "ERROR: config_get.py not found in $SCRIPTS_DIR or $SCRIPTS_DIR/../03_phenotyping" >&2
    echo "  Submitting via 'bsub <'? Export SCRIPTS_DIR=<repo>/05_qtl_mapping first." >&2
    exit 1
fi
ANCESTRY_MAP="${ANCESTRY_MAP:?ERROR: ANCESTRY_MAP env var required}"
K_GRID="${K_GRID:-0 5 10 15 20 25 30}"
FDR="${FDR:-0.05}"
SKIP_EXISTING="${SKIP_EXISTING:-0}"
# Covariates excluded before correlation pruning in every per-k model
# (coded replacement for the round-3 manual ct_Maternal hand-edit).
# Set to "" to disable.
EXCLUDE_COVARIATES="${EXCLUDE_COVARIATES-ct_Maternal}"

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load config values
eval "$(python3 "$CONFIG_GET" "${CONFIG}")"
OUTPUT_BASE="${OUTPUT_BASE}"

QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
HCP_DIR="${HCP_DIR:-${OUTPUT_BASE}/hcp}"
PC_DIR="${PC_DIR:-${OUTPUT_BASE}/genotype_pcs}"
WORK_DIR="${WORK_DIR:-${QTL_DIR}/hcp_optimization}"
mkdir -p "$WORK_DIR"

# ---- Environments ----
# tensorqtl conda env: python3 + pandas + tensorQTL (+ pysam for bgzip/tabix)
conda activate tensorqtl

# bgzip/tabix fallback if the tensorqtl env lacks pysam/htslib
if ! command -v bgzip >/dev/null 2>&1; then
    echo "  bgzip not found in tensorqtl env; stacking samtools-1.16.1 env"
    conda activate --stack samtools-1.16.1
fi

# Singularity R (sva + Rhcpp) for HCP re-estimation and the Storey q bridge
export R_LIBS_USER="/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
SING_R="singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript"
export QVALUE_RSCRIPT="$SING_R"

# ---- Determine ancestries ----
# Canonical env var is ANCESTRIES (space-separated). ANCESTRY (singular) is
# accepted as a fallback alias — passing ANCESTRY alone previously fell
# through to the all-in-map default and tried to run unprocessed ancestries.
if [ -z "${ANCESTRIES:-}" ] && [ -n "${ANCESTRY:-}" ]; then
    echo "  NOTE: ANCESTRY (singular) set; treating as ANCESTRIES='$ANCESTRY'"
    ANCESTRIES="$ANCESTRY"
fi
if [ -n "${ANCESTRIES:-}" ]; then
    ANCESTRY_LIST="$ANCESTRIES"
else
    echo "  NOTE: ANCESTRIES not set — defaulting to ALL ancestries in the map"
    ANCESTRY_LIST=$(python3 -c "
import pandas as pd
df = pd.read_csv('${ANCESTRY_MAP}', sep='\t')
for a in sorted(df['assigned_ancestry'].unique()):
    print(a)
")
fi

echo "[$(date)] HCP-count optimization (chr1 expression)"
echo "  CONFIG:       $CONFIG"
echo "  SCRIPTS_DIR:  $SCRIPTS_DIR"
echo "  ANCESTRY_MAP: $ANCESTRY_MAP"
echo "  QTL_DIR:      $QTL_DIR"
echo "  HCP_DIR:      $HCP_DIR"
echo "  PC_DIR:       $PC_DIR"
echo "  WORK_DIR:     $WORK_DIR"
echo "  Ancestries:   $ANCESTRY_LIST"
echo "  k grid:       $K_GRID"
echo "  FDR:          $FDR"
echo "  Exclude:      ${EXCLUDE_COVARIATES:-<none>}"

EXTRA_ARGS=""
if [ "$SKIP_EXISTING" = "1" ]; then
    EXTRA_ARGS="--skip-existing"
fi
if [ -n "$EXCLUDE_COVARIATES" ]; then
    EXTRA_ARGS="$EXTRA_ARGS --exclude-covariates $EXCLUDE_COVARIATES"
fi

python3 "${SCRIPTS_DIR}/optimize_hcp_chr1.py" \
    --qtl-dir "$QTL_DIR" \
    --hcp-dir "$HCP_DIR" \
    --pcair-dir "$PC_DIR" \
    --ancestry-map "$ANCESTRY_MAP" \
    --scripts-dir "$SCRIPTS_DIR" \
    --ancestries "$ANCESTRY_LIST" \
    --k-grid "$K_GRID" \
    --work-dir "$WORK_DIR" \
    --fdr "$FDR" \
    --r-cmd "$SING_R" \
    $EXTRA_ARGS

echo ""
echo "[$(date)] HCP optimization complete"
echo "  Per-ancestry results: ${WORK_DIR}/*_optimal_hcp.tsv / .png"
echo "  Canonical HCP files updated: ${QTL_DIR}/*_hcp_factors_harmonized.tsv"
echo "  Excluded covariates (pre-pruning): ${EXCLUDE_COVARIATES:-<none>}"
echo "  Next: run 25_build_covariates.py (canonical) with the same"
echo "        --exclude-covariates setting; no manual edits needed."
