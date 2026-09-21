#!/bin/bash
# =============================================================================
# 27_run_tensorqtl.sh — LSF wrapper for tensorQTL cis-xQTL mapping (any modality)
# =============================================================================
# Steps:
#   1. Sort + bgzip + tabix-index the harmonized phenotype BED (if not done)
#   2. Run tensorQTL cis.map_cis for one ancestry × modality
#
# Usage (NOTE: no inner quotes around -env values — LSF preserves them
# literally and mangles the paths):
#   bsub -env "CONFIG=/path/config.yml,SCRIPTS_DIR=/path/scripts,ANCESTRIES=EAS,MODALITY=splicing" < 27_run_tensorqtl.sh
#
# Required env:
#   CONFIG       — path to config.yml
#   SCRIPTS_DIR  — directory containing 27_run_tensorqtl.py and config_get.py
#   MODALITY     — modality label (expression, alt_polyA, alt_TSS,
#                  intron_retention, isoforms, RNA_editing, splicing,
#                  stability, combined)
# Optional env:
#   ANCESTRIES   — space-separated ancestry labels (default: EAS EUR)
#   CIS_WINDOW   — cis window in bp (default: 1000000)
#   MAF_THRESHOLD — minimum MAF (default: 0.01)
#   NO_GROUPS    — 1 = ignore phenotype_groups.txt and run ungrouped
#                  (per-phenotype lead variants; outputs get '_ungrouped'
#                  suffix). Default: 0 (grouped, one lead variant per gene)
#   INDEPENDENT  — 1 = after map_cis, run cis.map_independent (PANTRY-style
#                  forward-backward stepwise regression) for conditionally
#                  independent xQTLs; outputs get '_cisqtl_independent'
#                  suffix. Default: 0
#   COVARIATES_FILE — covariates TSV (tensorQTL orientation). The literal
#                  placeholder {ANC} is replaced by the ancestry label, so one
#                  value covers both strata. Default: {QTL_DIR}/{ANC}_covariates.tsv
#   QVALUE_METHOD — 'storey' (default; GTEx convention, R qvalue package via
#                  the compute_qvalues.R bridge) or 'bh' (Benjamini-Hochberg
#                  escape hatch, no R dependency)
#   QVALUE_RSCRIPT — Rscript command used for the Storey bridge. Default:
#                  the singularity R 4.3.1 container (same as scripts 19/20),
#                  with qvalue installed into R_LIBS_USER by
#                  21_install_tensorqtl.sh.
#
# Phenotype BED naming: {ANC}_{MODALITY}_harmonized.bed in qtl_inputs
# (expression: {ANC}_expression_harmonized.bed — same pattern).
# =============================================================================

#BSUB -q medium
#BSUB -n 8
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 12:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/tensorqtl.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/tensorqtl.%J.err

set -eo pipefail

# ---- Config / env ----
CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
MODALITY="${MODALITY:?ERROR: MODALITY env var required (e.g. expression, splicing)}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
CIS_WINDOW="${CIS_WINDOW:-1000000}"
MAF_THRESHOLD="${MAF_THRESHOLD:-0.01}"
NO_GROUPS="${NO_GROUPS:-0}"
INDEPENDENT="${INDEPENDENT:-0}"
COVARIATES_FILE="${COVARIATES_FILE:-}"
QVALUE_METHOD="${QVALUE_METHOD:-storey}"

# Strip any literal quotes that LSF's -env may have preserved in the values
CONFIG="${CONFIG%\"}";     CONFIG="${CONFIG#\"}"
SCRIPTS_DIR="${SCRIPTS_DIR%\"}"; SCRIPTS_DIR="${SCRIPTS_DIR#\"}"
MODALITY="${MODALITY%\"}"; MODALITY="${MODALITY#\"}"
ANCESTRIES="${ANCESTRIES%\"}"; ANCESTRIES="${ANCESTRIES#\"}"
NO_GROUPS="${NO_GROUPS%\"}"; NO_GROUPS="${NO_GROUPS#\"}"
INDEPENDENT="${INDEPENDENT%\"}"; INDEPENDENT="${INDEPENDENT#\"}"
COVARIATES_FILE="${COVARIATES_FILE%\"}"; COVARIATES_FILE="${COVARIATES_FILE#\"}"
QVALUE_METHOD="${QVALUE_METHOD%\"}"; QVALUE_METHOD="${QVALUE_METHOD#\"}"

# ---- Validate inputs before doing any work ----
if [ ! -f "$CONFIG" ]; then
    echo "ERROR: config file not found: '$CONFIG'"
    exit 1
fi
if [ ! -f "${SCRIPTS_DIR}/config_get.py" ]; then
    echo "ERROR: config_get.py not found: '${SCRIPTS_DIR}/config_get.py'"
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
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}")"
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

# Storey q-value bridge: R qvalue runs in the singularity R container (same
# setup as scripts 19/20 — proven on this cluster; avoids the conda-R/rpy2
# libstdc++ conflict). qvalue is installed into R_LIBS_USER by
# 21_install_tensorqtl.sh.
export R_LIBS_USER="${R_LIBS_USER:-/rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1}"
export QVALUE_RSCRIPT="${QVALUE_RSCRIPT:-singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript}"

echo "[$(date)] tensorQTL cis-xQTL Mapping"
echo "  CONFIG:        $CONFIG"
echo "  SCRIPTS_DIR:   $SCRIPTS_DIR"
echo "  QTL_DIR:       $QTL_DIR"
echo "  RESULTS_DIR:   $RESULTS_DIR"
echo "  ANCESTRIES:    $ANCESTRIES"
echo "  MODALITY:      $MODALITY"
echo "  CIS_WINDOW:    ±$((CIS_WINDOW / 1000)) kb"
echo "  MAF_THRESHOLD: $MAF_THRESHOLD"
echo "  NO_GROUPS:     $NO_GROUPS"
echo "  INDEPENDENT:   $INDEPENDENT"
echo ""

for ANC in $ANCESTRIES; do
    echo "============================================================"
    echo "[$(date)] Ancestry: $ANC / Modality: $MODALITY"
    echo "============================================================"

    # ---- Step 1: sort + bgzip + tabix index phenotype BED ----
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
        echo "  Skipping $ANC"
        continue
    fi

    # ---- Step 2: Run tensorQTL ----
    EXTRA_ARGS=""
    if [ "$NO_GROUPS" = "1" ]; then
        EXTRA_ARGS="$EXTRA_ARGS --no-groups"
        echo "  [Step 2] Running tensorQTL (UNGROUPED: per-phenotype lead variants)..."
    else
        echo "  [Step 2] Running tensorQTL (grouped if groups file present)..."
    fi
    if [ "$INDEPENDENT" = "1" ]; then
        EXTRA_ARGS="$EXTRA_ARGS --independent"
        echo "           + stepwise regression (cis.map_independent) for conditionally independent xQTLs"
    fi
    # Covariates: per-ancestry {ANC} placeholder substitution
    if [ -n "$COVARIATES_FILE" ]; then
        COV_FILE="${COVARIATES_FILE//\{ANC\}/$ANC}"
        if [ ! -f "$COV_FILE" ]; then
            echo "  ERROR: covariates file not found: $COV_FILE"
            echo "  (COVARIATES_FILE=$COVARIATES_FILE with {ANC} -> $ANC)"
            exit 1
        fi
        EXTRA_ARGS="$EXTRA_ARGS --covariates-file $COV_FILE"
        echo "  Covariates: $COV_FILE"
    fi
    EXTRA_ARGS="${EXTRA_ARGS# }"
    python3 "${SCRIPTS_DIR}/27_run_tensorqtl.py" \
        --qtl-dir "$QTL_DIR" \
        --output-dir "$RESULTS_DIR" \
        --ancestry "$ANC" \
        --modality "$MODALITY" \
        --cis-window "$CIS_WINDOW" \
        --maf-threshold "$MAF_THRESHOLD" \
        --qvalue-method "$QVALUE_METHOD" \
        $EXTRA_ARGS

    echo "  Done: $ANC / $MODALITY"
    echo ""
done

# ---- Summary ----
echo "[$(date)] Pipeline Complete"
echo ""
echo "  cis-xQTL results:"
ls -lh "${RESULTS_DIR}"/*_"${MODALITY}"*_cisqtl.parquet 2>/dev/null || echo "    (none found)"
echo ""
echo "  Top associations:"
ls -lh "${RESULTS_DIR}"/*_"${MODALITY}"*_cisqtl_top.tsv 2>/dev/null || echo "    (none found)"
