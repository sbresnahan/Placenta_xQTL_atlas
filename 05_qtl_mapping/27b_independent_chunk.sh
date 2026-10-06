#!/bin/bash
# =============================================================================
# 27b_independent_chunk.sh — LSF array task: one chunk of the tensorQTL
#                            cis.map_independent (stepwise) scan
# =============================================================================
# Runs 27_run_tensorqtl.py --independent-only --chunk-index $LSB_JOBINDEX for
# one ancestry × modality, reading the existing map_cis results
# ({ANC}_{MOD}[_ungrouped]_cisqtl.parquet) and writing
# qtl_results/independent_chunks/{ANC}_{MOD}[_ungrouped]/chunk_KKKK.parquet.
#
# Normally submitted by 28b_submit_independent.sh (which sizes the array from
# the number of FDR-significant phenotypes/groups, sets resources, and submits
# the dependent merge job). Direct submission:
#   bsub -J 'inqtl_EAS_splicing[1-42]' \
#        -env "CONFIG=/path/config.yml,SCRIPTS_DIR=/path/scripts,ANCESTRIES=EAS,MODALITY=splicing,MAF_THRESHOLD=0.01" \
#        < 27b_independent_chunk.sh
#
# Required env:
#   CONFIG       — path to config.yml
#   SCRIPTS_DIR  — directory containing 27_run_tensorqtl.py and config_get.py
#   MODALITY     — modality label (expression, splicing, combined, ...)
#   ANCESTRIES   — ONE ancestry label per array (the chunk index must map 1:1
#                  to $LSB_JOBINDEX)
# Optional env:
#   CHUNK_SIZE      — genes/groups (grouped) or phenotypes (ungrouped) per
#                     chunk (default: 100)
#   INDEPENDENT_FDR — FDR entry threshold, must match the merge (default: 0.05)
#   CIS_WINDOW      — cis window in bp (default: 1000000)
#   MAF_THRESHOLD   — minimum MAF (default: 0.01)
#   NO_GROUPS       — 1 = ungrouped layer (default: 0)
#   COVARIATES_FILE — covariates TSV; {ANC} and {MOD} placeholders replaced
#   SEED            — permutation RNG seed (default: unset = unseeded, as
#                     before; set for bit-reproducible chunks)
#
# GPU: no code change needed — tensorQTL auto-uses CUDA when visible. With
# GPU=1, 28b submits this script to GPU_QUEUE with GPU_OPTS; the python logs
# whether CUDA was actually used.
# =============================================================================

#BSUB -q medium
#BSUB -n 8
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 4:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/tensorqtl_ind_chunk.%J.%I.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/tensorqtl_ind_chunk.%J.%I.err

set -eo pipefail

# ---- Config / env ----
CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
# Default SCRIPTS_DIR to this script's own directory (repo-clone layout).
# Explicit SCRIPTS_DIR overrides — 28b_submit_independent.sh always passes it
# via bsub -env (LSF executes a spool copy; self-location would resolve to
# the spool directory).
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
MODALITY="${MODALITY:?ERROR: MODALITY env var required (e.g. expression, splicing)}"
ANCESTRIES="${ANCESTRIES:?ERROR: ANCESTRIES env var required (one ancestry per array)}"
CHUNK_INDEX="${LSB_JOBINDEX:?ERROR: LSB_JOBINDEX not set — submit as an LSF job array, e.g. bsub -J 'inqtl_ANC_MOD[1-N]'}"
CHUNK_SIZE="${CHUNK_SIZE:-100}"
INDEPENDENT_FDR="${INDEPENDENT_FDR:-0.05}"
CIS_WINDOW="${CIS_WINDOW:-1000000}"
MAF_THRESHOLD="${MAF_THRESHOLD:-0.01}"
NO_GROUPS="${NO_GROUPS:-0}"
COVARIATES_FILE="${COVARIATES_FILE:-}"
SEED="${SEED:-}"

# Strip any literal quotes that LSF's -env may have preserved in the values
CONFIG="${CONFIG%\"}";     CONFIG="${CONFIG#\"}"
SCRIPTS_DIR="${SCRIPTS_DIR%\"}"; SCRIPTS_DIR="${SCRIPTS_DIR#\"}"
MODALITY="${MODALITY%\"}"; MODALITY="${MODALITY#\"}"
ANCESTRIES="${ANCESTRIES%\"}"; ANCESTRIES="${ANCESTRIES#\"}"
CHUNK_SIZE="${CHUNK_SIZE%\"}"; CHUNK_SIZE="${CHUNK_SIZE#\"}"
INDEPENDENT_FDR="${INDEPENDENT_FDR%\"}"; INDEPENDENT_FDR="${INDEPENDENT_FDR#\"}"
NO_GROUPS="${NO_GROUPS%\"}"; NO_GROUPS="${NO_GROUPS#\"}"
COVARIATES_FILE="${COVARIATES_FILE%\"}"; COVARIATES_FILE="${COVARIATES_FILE#\"}"
SEED="${SEED%\"}"; SEED="${SEED#\"}"

# ---- Validate inputs before doing any work ----
if [ ! -f "$CONFIG" ]; then
    echo "ERROR: config file not found: '$CONFIG'"
    exit 1
fi
# config_get.py lives in ../03_phenotyping in the repo layout; a flat copy in
# SCRIPTS_DIR (flat deployment) takes precedence.
CONFIG_GET="${SCRIPTS_DIR}/config_get.py"
[ -f "$CONFIG_GET" ] || CONFIG_GET="${SCRIPTS_DIR}/../03_phenotyping/config_get.py"
if [ ! -f "$CONFIG_GET" ]; then
    echo "ERROR: config_get.py not found in $SCRIPTS_DIR or $SCRIPTS_DIR/../03_phenotyping"
    exit 1
fi
if [ ! -f "${SCRIPTS_DIR}/27_run_tensorqtl.py" ]; then
    echo "ERROR: 27_run_tensorqtl.py not found: '${SCRIPTS_DIR}/27_run_tensorqtl.py'"
    exit 1
fi

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load config values
eval "$(python3 "$CONFIG_GET" "${CONFIG}")"
if [ -z "$OUTPUT_BASE" ] || [ ! -d "$OUTPUT_BASE" ]; then
    echo "ERROR: OUTPUT_BASE missing or not a directory: '$OUTPUT_BASE'"
    echo "  (check that config_get.py parsed $CONFIG correctly)"
    exit 1
fi

QTL_DIR="${OUTPUT_BASE}/qtl_inputs"
RESULTS_DIR="${OUTPUT_BASE}/qtl_results"
mkdir -p "$RESULTS_DIR"

# Activate tensorQTL env
conda activate tensorqtl

# Also need bgzip + tabix (from samtools env or system)
conda activate --stack samtools-1.16.1 2>/dev/null || true

echo "[$(date)] tensorQTL chunked independent (stepwise) scan"
echo "  CONFIG:        $CONFIG"
echo "  SCRIPTS_DIR:   $SCRIPTS_DIR"
echo "  QTL_DIR:       $QTL_DIR"
echo "  RESULTS_DIR:   $RESULTS_DIR"
echo "  ANCESTRIES:    $ANCESTRIES"
echo "  MODALITY:      $MODALITY"
echo "  CHUNK_INDEX:   $CHUNK_INDEX (chunk size $CHUNK_SIZE)"
echo "  INDEPENDENT_FDR: $INDEPENDENT_FDR"
echo "  CIS_WINDOW:    ±$((CIS_WINDOW / 1000)) kb"
echo "  MAF_THRESHOLD: $MAF_THRESHOLD"
echo "  NO_GROUPS:     $NO_GROUPS"
echo "  SEED:          ${SEED:-<unset>}"
echo ""

for ANC in $ANCESTRIES; do
    echo "============================================================"
    echo "[$(date)] Ancestry: $ANC / Modality: $MODALITY / chunk $CHUNK_INDEX"
    echo "============================================================"

    # ---- Step 1: sort + bgzip + tabix index phenotype BED ----
    # Should already exist (created by the map_cis job); kept for standalone
    # robustness.
    PHENO_BED="${QTL_DIR}/${ANC}_${MODALITY}_harmonized.bed"
    PHENO_BGZ="${QTL_DIR}/${ANC}_${MODALITY}.bed.gz"

    if [ -f "$PHENO_BGZ" ] && [ -f "${PHENO_BGZ}.tbi" ]; then
        echo "  [Step 1] bgzip + tabix: already exists, skipping"
    elif [ -f "$PHENO_BED" ]; then
        echo "  [Step 1] Sorting + bgzipping + tabix-indexing phenotype BED..."
        # Sort by chrom + start position (tabix requires sorted input).
        # Keep header line first, sort the rest.
        (head -1 "$PHENO_BED" && tail -n +2 "$PHENO_BED" | sort -k1,1 -k2,2n) | bgzip -c > "$PHENO_BGZ"
        tabix -p bed "$PHENO_BGZ"
        echo "    Written: $PHENO_BGZ + .tbi"
    else
        echo "  ERROR: phenotype BED not found: $PHENO_BED"
        echo "  (run 26_harmonize_modalities.py first for non-expression modalities)"
        exit 1
    fi

    # ---- Step 2: run one chunk of cis.map_independent ----
    EXTRA_ARGS="--independent-only --chunk-index $CHUNK_INDEX --chunk-size $CHUNK_SIZE --independent-fdr $INDEPENDENT_FDR"
    if [ "$NO_GROUPS" = "1" ]; then
        EXTRA_ARGS="$EXTRA_ARGS --no-groups"
    fi
    if [ -n "$SEED" ]; then
        EXTRA_ARGS="$EXTRA_ARGS --seed $SEED"
    fi
    # Covariates: {ANC}/{MOD} placeholder substitution
    if [ -n "$COVARIATES_FILE" ]; then
        COV_FILE="${COVARIATES_FILE//\{ANC\}/$ANC}"
        COV_FILE="${COV_FILE//\{MOD\}/$MODALITY}"
        if [ ! -f "$COV_FILE" ]; then
            echo "  ERROR: covariates file not found: $COV_FILE"
            echo "  (COVARIATES_FILE=$COVARIATES_FILE with {ANC} -> $ANC, {MOD} -> $MODALITY)"
            exit 1
        fi
        EXTRA_ARGS="$EXTRA_ARGS --covariates-file $COV_FILE"
        echo "  Covariates: $COV_FILE"
    fi
    python3 "${SCRIPTS_DIR}/27_run_tensorqtl.py" \
        --qtl-dir "$QTL_DIR" \
        --output-dir "$RESULTS_DIR" \
        --ancestry "$ANC" \
        --modality "$MODALITY" \
        --cis-window "$CIS_WINDOW" \
        --maf-threshold "$MAF_THRESHOLD" \
        $EXTRA_ARGS

    echo "  Done: $ANC / $MODALITY / chunk $CHUNK_INDEX"
    echo ""
done

echo "[$(date)] Chunk complete"
