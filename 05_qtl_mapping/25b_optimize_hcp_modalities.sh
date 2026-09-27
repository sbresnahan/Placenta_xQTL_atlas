#!/bin/bash
# =============================================================================
# 25b_optimize_hcp_modalities.sh — LSF wrapper for per-modality HCP optimization
# =============================================================================
# Selects the number of HCP hidden covariates for EACH modality and the
# combined cross-modality arm by maximizing cis-eGene discovery (Storey
# q <= 0.05), with HCPs re-estimated at each candidate k from that
# modality's own harmonized BED (hcp_from_matrix.R — HCP-only; the BEDs are
# already QN+INT+ComBat'd). 25a_optimize_hcp.sh covers the expression-only
# special case. Mapping scope per modality: chr1 subset when the BED has
# >= CHR1_MIN chr1 phenotypes (default 300), else genome-wide.
#
# Intended parallel unit: ONE ancestry x modality per job (20 jobs for
# EAS+EUR x 9 modalities + combined). See docs/runbook_modality_hcp.md for
# the full submission loop.
#
# Usage:
#   bsub -J hcpopt_EAS_splicing -q medium -n 4 -M 32 -R "rusage[mem=32]" -W 24:00 \
#        -o .../logs/hcpopt_EAS_splicing.%J.out -e .../logs/hcpopt_EAS_splicing.%J.err \
#        -env "CONFIG=/path/config.yml,SCRIPTS_DIR=/path/scripts,ANCESTRIES=EAS,MODALITIES=splicing" \
#        < 25b_optimize_hcp_modalities.sh
#
# Required env:
#   CONFIG       — path to config.yml
#   SCRIPTS_DIR  — directory containing optimize_hcp_modalities.py,
#                  hcp_from_matrix.R, 25_build_covariates.py, 27_run_tensorqtl.py
# Optional env:
#   ANCESTRIES   — space-separated ancestry labels (default: EAS EUR)
#   MODALITIES   — space-separated modality labels incl. "combined"
#                  (default: all 9 modalities + combined)
#   K_GRID       — candidate HCP counts (default: "0 5 10 15 20 25 30")
#   QTL_DIR      — canonical QTL inputs dir (default: ${OUTPUT_BASE}/qtl_inputs)
#   QC_METRICS   — pooled Picard QC metrics (default: ${OUTPUT_BASE}/hcp/all_qc_metrics.tsv)
#   PC_DIR       — genotype PCs dir (default: ${OUTPUT_BASE}/genotype_pcs)
#   WORK_DIR     — staging dir (default: ${QTL_DIR}/hcp_optimization_modalities)
#   FDR          — Storey q threshold for eGene counts (default: 0.05)
#   CHR1_MIN     — minimum chr1 phenotypes for chr1-based selection (default: 300)
#   MAX_HCP_PHENOTYPES — phenotype cap for HCP estimation (default: 40000; 0 disables)
#   SKIP_EXISTING — set to 1 to reuse existing per-k results (resumable)
#   EXCLUDE_COVARIATES — covariates excluded before correlation pruning in
#                  every per-k model (default: ct_Maternal; "" disables)
# =============================================================================

#BSUB -q medium
#BSUB -n 4
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 24:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/hcpopt_mod.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/hcpopt_mod.%J.err

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
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
MODALITIES="${MODALITIES:-expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability combined}"
K_GRID="${K_GRID:-0 5 10 15 20 25 30}"
FDR="${FDR:-0.05}"
CHR1_MIN="${CHR1_MIN:-300}"
MAX_HCP_PHENOTYPES="${MAX_HCP_PHENOTYPES:-40000}"
SKIP_EXISTING="${SKIP_EXISTING:-0}"
# Covariates excluded before correlation pruning in every per-k model.
# Set to "" to disable.
EXCLUDE_COVARIATES="${EXCLUDE_COVARIATES-ct_Maternal}"

# Strip any literal quotes that LSF's -env may have preserved in the values
CONFIG="${CONFIG%\"}";           CONFIG="${CONFIG#\"}"
SCRIPTS_DIR="${SCRIPTS_DIR%\"}"; SCRIPTS_DIR="${SCRIPTS_DIR#\"}"
ANCESTRIES="${ANCESTRIES%\"}";   ANCESTRIES="${ANCESTRIES#\"}"
MODALITIES="${MODALITIES%\"}";   MODALITIES="${MODALITIES#\"}"

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load config values
eval "$(python3 "$CONFIG_GET" "${CONFIG}")"
OUTPUT_BASE="${OUTPUT_BASE}"

QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
QC_METRICS="${QC_METRICS:-${OUTPUT_BASE}/hcp/all_qc_metrics.tsv}"
PC_DIR="${PC_DIR:-${OUTPUT_BASE}/genotype_pcs}"
WORK_DIR="${WORK_DIR:-${QTL_DIR}/hcp_optimization_modalities}"
mkdir -p "$WORK_DIR"

# ---- Environments ----
# tensorqtl conda env: python3 + pandas + tensorQTL (+ pysam for bgzip/tabix)
conda activate tensorqtl

# bgzip/tabix fallback if the tensorqtl env lacks pysam/htslib
if ! command -v bgzip >/dev/null 2>&1; then
    echo "  bgzip not found in tensorqtl env; stacking samtools-1.16.1 env"
    conda activate --stack samtools-1.16.1
fi

# Singularity R (Rhcpp) for HCP estimation and the Storey q bridge
export R_LIBS_USER="/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
SING_R="singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript"
export QVALUE_RSCRIPT="$SING_R"

echo "[$(date)] Per-modality HCP-count optimization"
echo "  CONFIG:       $CONFIG"
echo "  SCRIPTS_DIR:  $SCRIPTS_DIR"
echo "  QTL_DIR:      $QTL_DIR"
echo "  QC_METRICS:   $QC_METRICS"
echo "  PC_DIR:       $PC_DIR"
echo "  WORK_DIR:     $WORK_DIR"
echo "  Ancestries:   $ANCESTRIES"
echo "  Modalities:   $MODALITIES"
echo "  k grid:       $K_GRID"
echo "  FDR:          $FDR"
echo "  chr1 minimum: $CHR1_MIN"
echo "  HCP phenotype cap: $MAX_HCP_PHENOTYPES"
echo "  Exclude:      ${EXCLUDE_COVARIATES:-<none>}"

EXTRA_ARGS=""
if [ "$SKIP_EXISTING" = "1" ]; then
    EXTRA_ARGS="--skip-existing"
fi
if [ -n "$EXCLUDE_COVARIATES" ]; then
    EXTRA_ARGS="$EXTRA_ARGS --exclude-covariates $EXCLUDE_COVARIATES"
fi

python3 "${SCRIPTS_DIR}/optimize_hcp_modalities.py" \
    --qtl-dir "$QTL_DIR" \
    --qc-metrics "$QC_METRICS" \
    --pcair-dir "$PC_DIR" \
    --scripts-dir "$SCRIPTS_DIR" \
    --ancestries "$ANCESTRIES" \
    --modalities "$MODALITIES" \
    --k-grid "$K_GRID" \
    --work-dir "$WORK_DIR" \
    --fdr "$FDR" \
    --chr1-min-phenotypes "$CHR1_MIN" \
    --max-hcp-phenotypes "$MAX_HCP_PHENOTYPES" \
    --r-cmd "$SING_R" \
    $EXTRA_ARGS

echo ""
echo "[$(date)] Per-modality HCP optimization complete"
echo "  Results: ${WORK_DIR}/*_optimal_hcp.tsv / .png"
echo "  Installed: ${QTL_DIR}/*_hcp_factors_optimized.tsv"
echo "  Next: canonical per-modality 25_build_covariates.py"
echo "        (see docs/runbook_modality_hcp.md)"
