#!/bin/bash
# =============================================================================
# 19_hcp_factors.sh — LSF wrapper for HCP latent factor extraction pipeline
# =============================================================================
# Orchestrates the full HCP pipeline on seadragon:
#   Stage 1: PicardTools QC metrics (per-sample, per cohort)
#   Stage 2: Pool QC metrics across cohorts
#   Stage 3: Pool expression within ancestry strata
#   Stage 4: ComBat + INT + HCP estimation (per ancestry)
#
# Usage:
#   bsub -env "CONFIG=\"config.yml\",SCRIPTS_DIR=\"/path/to/scripts\",ANCESTRY_MAP=\"/path/to/pooled_sample_ancestry_RNAseq.tsv\"" < 19_hcp_factors.sh
#
# Or run directly (interactive / non-LSF):
#   CONFIG=config.yml SCRIPTS_DIR=/path/to/scripts \
#   ANCESTRY_MAP=/path/to/pooled_sample_ancestry_RNAseq.tsv bash 19_hcp_factors.sh
#
# Required env vars:
#   CONFIG        — path to config.yml
#   SCRIPTS_DIR   — directory containing picard_qc.py, pool_expression_within_ancestry.py,
#                   combat_normalize_hcp.R, and config_get.py
#   ANCESTRY_MAP  — path to pooled_sample_ancestry_RNAseq.tsv
#
# Optional env vars:
#   COHORTS       — space-separated cohort names to process (default: all from config)
#   ANCESTRIES    — space-separated ancestry labels to process (default: all)
#   K             — number of HCP factors (default: 15)
#   OUTPUT_DIR    — output directory (default: ${OUTPUT_BASE}/hcp)
#   REFFLAT       — path to refFlat file (default: ${REFERENCE_DIR}/HPLRv2.refFlat)
#   GENE_ANNOT    — path to gene GC/length annotation TSV (optional)
#   FASTA         — path to genome FASTA (for generating gene annotation if needed)
# =============================================================================

#BSUB -q long
#BSUB -n 4
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 24:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/hcp.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/hcp.%J.err

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
ANCESTRY_MAP="${ANCESTRY_MAP:?ERROR: ANCESTRY_MAP env var required (path to pooled_sample_ancestry_RNAseq.tsv)}"
K="${K:-15}"
OUTPUT_DIR="${OUTPUT_DIR:-}"

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load config values
eval "$(python3 "$CONFIG_GET" "${CONFIG}")"

OUTPUT_BASE="${OUTPUT_BASE}"
REFERENCE_DIR="${REFERENCE_DIR}"
NORMALIZED_GTF="${NORMALIZED_GTF}"
REF_GENOME="${REF_GENOME}"

# Map ancestry-map cohort names -> PANTRY cohort data directories.
# The ancestry map uses real cohort names (NIEHS_RICHS, GUSTO, ...) while the
# data directories use config names (cohort1, ...). The pooler looks up
# directories by the map's cohort values, so translate via --cohort-dirs.
# Override with the COHORT_DIR_MAP env var if the mapping changes.
COHORT_DIR_MAP="${COHORT_DIR_MAP:-NIEHS_RICHS=cohort1 GUSTO=cohort2 SNUH=cohort3 NIGMS=cohort4}"
COHORT_DIRS_ARGS=""
for pair in $COHORT_DIR_MAP; do
    name="${pair%%=*}"
    dir="${pair##*=}"
    COHORT_DIRS_ARGS="${COHORT_DIRS_ARGS} ${name}=${OUTPUT_BASE}/${dir}"
done
echo "  Cohort dir map: $COHORT_DIRS_ARGS"

if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="${OUTPUT_BASE}/hcp"
fi
mkdir -p "$OUTPUT_DIR"

REFFLAT="${REFFLAT:-${REFERENCE_DIR}/HPLRv2.refFlat}"
GENE_ANNOT="${GENE_ANNOT:-${OUTPUT_DIR}/gene_gc_length.tsv}"
FASTA="${FASTA:-${REF_GENOME}}"

echo "[$(date)] HCP Latent Factor Extraction Pipeline"
echo "  CONFIG:       $CONFIG"
echo "  SCRIPTS_DIR:  $SCRIPTS_DIR"
echo "  ANCESTRY_MAP: $ANCESTRY_MAP"
echo "  OUTPUT_DIR:   $OUTPUT_DIR"
echo "  K:            $K"
echo "  REFFLAT:      $REFFLAT"
echo "  GENE_ANNOT:   $GENE_ANNOT"

# ---- Determine cohorts ----
if [ -n "$COHORTS" ]; then
    COHORT_LIST="$COHORTS"
else
    # Parse cohort names from config.yml
    COHORT_LIST=$(python3 -c "
import re
with open('${CONFIG}') as f:
    in_cohorts = False
    for line in f:
        stripped = line.split('#')[0].rstrip()
        if not stripped.strip():
            continue
        if stripped.startswith('cohorts:'):
            in_cohorts = True
            continue
        if in_cohorts:
            if re.match(r'^  \w', stripped) and stripped.strip().endswith(':'):
                print(stripped.strip().rstrip(':'))
            elif not stripped.startswith(' '):
                break
")
fi
echo "  Cohorts: $COHORT_LIST"

# ---- Determine ancestries ----
if [ -n "$ANCESTRIES" ]; then
    ANCESTRY_LIST="$ANCESTRIES"
else
    ANCESTRY_LIST=$(python3 -c "
import pandas as pd
df = pd.read_csv('${ANCESTRY_MAP}', sep='\t')
for a in sorted(df['assigned_ancestry'].unique()):
    print(a)
")
fi
echo "  Ancestries: $ANCESTRY_LIST"

# =============================================================================
# Stage 0: Generate refFlat and gene GC/length annotation (if not present)
# =============================================================================
echo ""
echo "[$(date)] Stage 0: Reference files"

if [ ! -f "$REFFLAT" ]; then
    echo "  Generating refFlat from $NORMALIZED_GTF"
    conda activate samtools-1.16.1
    source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate
    python3 "${SCRIPTS_DIR}/picard_qc.py" \
        --generate-refflat \
        --gtf "$NORMALIZED_GTF" \
        --output "$REFFLAT"
    conda deactivate 2>/dev/null || true
else
    echo "  refFlat exists: $REFFLAT"
fi

if [ ! -f "$GENE_ANNOT" ]; then
    echo "  Generating gene GC/length annotation from $NORMALIZED_GTF + $FASTA"
    conda activate samtools-1.16.1
    source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate
    python3 "${SCRIPTS_DIR}/picard_qc.py" \
        --generate-gene-annot \
        --gtf "$NORMALIZED_GTF" \
        --fasta "$FASTA" \
        --output "$GENE_ANNOT"
    conda deactivate 2>/dev/null || true
else
    echo "  Gene annotation exists: $GENE_ANNOT"
fi

# =============================================================================
# Stage 1: PicardTools QC metrics (per cohort)
# =============================================================================
echo ""
echo "[$(date)] Stage 1: PicardTools QC metrics"

module load picard
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
# picard-2.27.4 env provides the picard binary + java; samtools env + MAJIQ
# venv (stacked on top) provides python3 with pandas/numpy
conda activate picard-2.27.4
conda activate --stack samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

# ---- Verify picard is reachable ---------------------------------------------
# picard_qc.py invokes the literal `picard` executable; if the env stack above
# didn't put it on PATH, every QC call fails per sample (garbage metrics).
# Fail fast with diagnostics instead.
if ! command -v picard >/dev/null 2>&1; then
    echo "ERROR: 'picard' not on PATH after env stack (module load picard +"
    echo "  conda activate picard-2.27.4 + samtools + MAJIQ venv). State:"
    echo "  CONDA_PREFIX=${CONDA_PREFIX:-unset}"
    echo "  picard binary:  $(command -v picard || echo none)"
    echo "  java binary:    $(command -v java || echo none)"
    echo "  env bin/ contents: $(ls "${CONDA_PREFIX:-/nonexistent}"/bin 2>/dev/null | grep -i -m3 picard || echo 'no picard* in $CONDA_PREFIX/bin')"
    echo "Fix the activation (env name/path) — do not shim around it."
    exit 1
fi

QC_DIR="${OUTPUT_DIR}/qc_metrics"
mkdir -p "$QC_DIR"

for COHORT in $COHORT_LIST; do
    echo "  [$(date)] Cohort: $COHORT"
    COHORT_DIR="${OUTPUT_BASE}/${COHORT}"
    SALMON_DIR="${COHORT_DIR}/intermediate/expression"
    SAMPLES_FILE="${COHORT_DIR}/samples.txt"

    QC_OUTPUT="${QC_DIR}/${COHORT}_qc_metrics.tsv"

    if [ -f "$QC_OUTPUT" ]; then
        echo "    Already exists, skipping: $QC_OUTPUT"
        continue
    fi

    python3 "${SCRIPTS_DIR}/picard_qc.py" \
        --config "$CONFIG" \
        --cohort "$COHORT" \
        --samples "$SAMPLES_FILE" \
        --output "$QC_OUTPUT" \
        --refflat "$REFFLAT" \
        --picard-cmd picard \
        --salmon-dir "$SALMON_DIR" \
        --gtf "$NORMALIZED_GTF" \
        --gene-annot "$GENE_ANNOT"
done

conda deactivate 2>/dev/null || true

# =============================================================================
# Stage 2: Pool QC metrics across cohorts
# =============================================================================
echo ""
echo "[$(date)] Stage 2: Pool QC metrics across cohorts"

ALL_QC="${OUTPUT_DIR}/all_qc_metrics.tsv"

python3 -c "
import pandas as pd
import glob
import os

qc_files = sorted(glob.glob('${QC_DIR}/*_qc_metrics.tsv'))
if not qc_files:
    raise SystemExit('ERROR: no QC metric files found in ${QC_DIR}')

dfs = []
for f in qc_files:
    cohort = os.path.basename(f).replace('_qc_metrics.tsv', '')
    df = pd.read_csv(f, sep='\t', index_col=0)
    df['cohort'] = cohort
    dfs.append(df)
    print(f'  {cohort}: {df.shape[0]} samples x {df.shape[1]-1} metrics')

pooled = pd.concat(dfs, axis=0)
# Keep cohort column for reference but don't use it as a metric
pooled_metrics = pooled.drop(columns=['cohort'])
pooled_metrics.to_csv('${ALL_QC}', sep='\t')
print(f'Pooled: {pooled_metrics.shape[0]} samples x {pooled_metrics.shape[1]} metrics')
print(f'Written: ${ALL_QC}')
"

# =============================================================================
# Stage 3: Pool expression within ancestry strata
# =============================================================================
echo ""
echo "[$(date)] Stage 3: Pool expression within ancestry strata"

POOLED_DIR="${OUTPUT_DIR}/pooled_expression"
mkdir -p "$POOLED_DIR"

python3 "${SCRIPTS_DIR}/pool_expression_within_ancestry.py" \
    --ancestry-map "$ANCESTRY_MAP" \
    --config "$CONFIG" \
    --cohort-dirs $COHORT_DIRS_ARGS \
    --modality expression \
    --output-dir "$POOLED_DIR" \
    --ancestries $ANCESTRY_LIST

# =============================================================================
# Stage 4: ComBat + INT + HCP estimation (per ancestry)
# =============================================================================
echo ""
echo "[$(date)] Stage 4: ComBat + INT + HCP estimation"

HCP_DIR="${OUTPUT_DIR}/hcp_factors"
mkdir -p "$HCP_DIR"

export R_LIBS_USER="/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
SING_R="singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript"

for ANCESTRY in $ANCESTRY_LIST; do
    echo "  [$(date)] Ancestry: $ANCESTRY"
    EXPR_FILE="${POOLED_DIR}/${ANCESTRY}_pooled_expression.bed"

    if [ ! -f "$EXPR_FILE" ]; then
        echo "    WARN: pooled expression not found: $EXPR_FILE, skipping"
        continue
    fi

    $SING_R "${SCRIPTS_DIR}/combat_normalize_hcp.R" \
        --expression "$EXPR_FILE" \
        --qc-metrics "$ALL_QC" \
        --ancestry-map "$ANCESTRY_MAP" \
        --ancestry "$ANCESTRY" \
        --k "$K" \
        --output-dir "$HCP_DIR" \
        --qc-cor-threshold 0.9 \
        --lambda1 0.5 \
        --lambda2 1 \
        --lambda3 1
done

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "[$(date)] HCP Pipeline Complete"
echo "  QC metrics:       ${ALL_QC}"
echo "  Pooled expression: ${POOLED_DIR}/"
echo "  HCP factors:       ${HCP_DIR}/"
echo ""
echo "  HCP factor files (for tensorQTL --covariates):"
ls -lh "${HCP_DIR}"/*_hcp_factors.tsv 2>/dev/null || echo "    (none found)"
echo ""
echo "  ComBat+INT expression BEDs:"
ls -lh "${HCP_DIR}"/*_combat_int_expression.bed 2>/dev/null || echo "    (none found)"
echo ""
echo "  Diagnostics:"
ls -lh "${HCP_DIR}"/*_hcp_diagnostics.pdf 2>/dev/null || echo "    (none found)"
