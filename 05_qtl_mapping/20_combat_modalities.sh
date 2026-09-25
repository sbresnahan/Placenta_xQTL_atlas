#!/bin/bash
# =============================================================================
# 20_combat_modalities.sh — LSF wrapper for cross-cohort ComBat + pooling
#                            of non-expression RNA modalities (roadmap step 4)
# =============================================================================
# Orchestrates pooling + ComBat + INT for the 6 non-expression modalities
# (isoforms, alt_TSS, alt_polyA, splicing, intron_retention, RNA_editing,
#  stability) across all ancestry strata on seadragon.
#
# Gene-level expression is handled separately by the HCP pipeline
# (19_hcp_factors.sh); HCP factors estimated there are reused as covariates
# for ALL modalities and are NOT re-estimated here.
#
# Pipeline:
#   Stage 1: Pool modalities within ancestry strata
#     - isoforms, alt_TSS, alt_polyA, RNA_editing, stability:
#       per-cohort unnorm BEDs -> phenotype_id intersection (original sample IDs)
#     - splicing, intron_retention:
#       pass through pre-pooled BEDs from harmonize_within_ancestry.py
#       (namespaced sample IDs)
#   Stage 2: QN + INT + ComBat per ancestry x modality (R via singularity)
#     - normalization first, ComBat last (devBrain xQTL schema, 2026-09 revision)
#     - isoforms additionally exclude the expression-outlier samples written
#       by 19_hcp_factors.sh (devBrain §3.3)
#
# Usage:
#   bsub -env "CONFIG=\"config.yml\",SCRIPTS_DIR=\"/path/to/scripts\",ANCESTRY_MAP=\"/path/to/pooled_sample_ancestry_RNAseq.tsv\"" < 20_combat_modalities.sh
#
# Or run directly (interactive / non-LSF):
#   CONFIG=config.yml SCRIPTS_DIR=/path/to/scripts \
#   ANCESTRY_MAP=/path/to/pooled_sample_ancestry_RNAseq.tsv bash 20_combat_modalities.sh
#
# Required env vars:
#   CONFIG        — path to config.yml
#   SCRIPTS_DIR   — directory containing pool_modalities_within_ancestry.py,
#                   combat_normalize_modalities.R, and config_get.py
#   ANCESTRY_MAP  — path to pooled_sample_ancestry_RNAseq.tsv
#
# Optional env vars:
#   MODALITIES      — space-separated modality names to process
#                     (default: isoforms alt_TSS alt_polyA splicing
#                      intron_retention RNA_editing stability)
#   ANCESTRIES      — space-separated ancestry labels to process (default: all)
#   OUTPUT_DIR      — output directory (default: ${OUTPUT_BASE}/combat_modalities)
#   HARMONIZE_DIR   — root dir of harmonize_within_ancestry.py outputs
#                     (default: ${OUTPUT_BASE}; 17 writes
#                      ${OUTPUT_BASE}/<ancestry>/<modality>/unnorm/<modality>.bed)
# =============================================================================

#BSUB -q medium
#BSUB -n 4
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 24:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/combat_mod.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/combat_mod.%J.err

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
OUTPUT_DIR="${OUTPUT_DIR:-}"
HARMONIZE_DIR="${HARMONIZE_DIR:-}"

# Default modality list (all non-expression modalities)
DEFAULT_MODALITIES="isoforms alt_TSS alt_polyA splicing intron_retention RNA_editing stability"
MODALITIES="${MODALITIES:-$DEFAULT_MODALITIES}"

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load config values
eval "$(python3 "$CONFIG_GET" "${CONFIG}")"

OUTPUT_BASE="${OUTPUT_BASE}"

# Map ancestry-map cohort names -> PANTRY cohort data directories.
# The ancestry map uses real cohort names (NIEHS_RICHS, GUSTO, ...) while the
# data directories use config names (cohort1, ...). The pooler looks up
# directories by the map's cohort values, so translate via --cohort-dirs.
# (Only used by intersection modalities; ignored for pre-pooled splicing/IR.)
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
    OUTPUT_DIR="${OUTPUT_BASE}/combat_modalities"
fi
if [ -z "$HARMONIZE_DIR" ]; then
    # 17_harmonize_within_ancestry.sh writes ${OUTPUT_BASE}/<ANC>/<modality>/...
    HARMONIZE_DIR="${OUTPUT_BASE}"
fi
mkdir -p "$OUTPUT_DIR"

echo "[$(date)] Cross-cohort ComBat + Pooling for Non-Expression Modalities"
echo "  CONFIG:        $CONFIG"
echo "  SCRIPTS_DIR:   $SCRIPTS_DIR"
echo "  ANCESTRY_MAP:  $ANCESTRY_MAP"
echo "  OUTPUT_DIR:    $OUTPUT_DIR"
echo "  HARMONIZE_DIR: $HARMONIZE_DIR"
echo "  MODALITIES:    $MODALITIES"

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

# Modalities that use pre-pooled (harmonized) inputs
PREPOOLED_MODALITIES="splicing intron_retention"

# =============================================================================
# Stage 1: Pool modalities within ancestry strata
# =============================================================================
echo ""
echo "[$(date)] Stage 1: Pool modalities within ancestry strata"

POOLED_DIR="${OUTPUT_DIR}/pooled"
mkdir -p "$POOLED_DIR"

for MODALITY in $MODALITIES; do
    echo "  [$(date)] Modality: $MODALITY"

    # Determine if this modality is pre-pooled (splicing/IR) or intersection
    IS_PREPOOLED=false
    for pm in $PREPOOLED_MODALITIES; do
        if [ "$MODALITY" = "$pm" ]; then
            IS_PREPOOLED=true
            break
        fi
    done

    if [ "$IS_PREPOOLED" = "true" ]; then
        # Pre-pooled pass-through (splicing, intron_retention)
        EXTRA_ARGS="--pre-pooled-dir ${HARMONIZE_DIR}"
    else
        # Per-cohort intersection pooling
        EXTRA_ARGS=""
    fi

    python3 "${SCRIPTS_DIR}/pool_modalities_within_ancestry.py" \
        --ancestry-map "$ANCESTRY_MAP" \
        --config "$CONFIG" \
        --cohort-dirs $COHORT_DIRS_ARGS \
        --modality "$MODALITY" \
        --output-dir "$POOLED_DIR" \
        --ancestries $ANCESTRY_LIST \
        $EXTRA_ARGS
done

# =============================================================================
# Stage 2: ComBat + INT per ancestry x modality (R via singularity)
# =============================================================================
echo ""
echo "[$(date)] Stage 2: QN + INT + ComBat per ancestry x modality"

COMBAT_DIR="${OUTPUT_DIR}/combat_int"
mkdir -p "$COMBAT_DIR"

export R_LIBS_USER="/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
SING_R="singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript"

for ANCESTRY in $ANCESTRY_LIST; do
    echo "  [$(date)] Ancestry: $ANCESTRY"

    for MODALITY in $MODALITIES; do
        echo "    [$(date)] Modality: $MODALITY"
        POOLED_BED="${POOLED_DIR}/${ANCESTRY}_${MODALITY}_pooled.bed"

        if [ ! -f "$POOLED_BED" ]; then
            echo "      WARN: pooled BED not found: $POOLED_BED, skipping"
            continue
        fi

        # Determine if this modality is pre-pooled (needs cohort-label sidecar)
        IS_PREPOOLED=false
        for pm in $PREPOOLED_MODALITIES; do
            if [ "$MODALITY" = "$pm" ]; then
                IS_PREPOOLED=true
                break
            fi
        done

        EXTRA_ARGS=""
        if [ "$IS_PREPOOLED" = "true" ]; then
            LABEL_FILE="${POOLED_DIR}/${ANCESTRY}_${MODALITY}_cohort_labels.tsv"
            if [ -f "$LABEL_FILE" ]; then
                EXTRA_ARGS="--cohort-labels ${LABEL_FILE}"
            else
                echo "      WARN: cohort-label sidecar not found: $LABEL_FILE"
                echo "             falling back to ancestry map for batch labels"
            fi
        fi

        # Expression-outlier exclusion for isoforms (devBrain §3.3): apply the
        # connectivity-outlier list written by 19_hcp_factors.sh to isoforms
        # only (not splicing / proportion modalities).
        if [ "$MODALITY" = "isoforms" ]; then
            OUTLIER_FILE="${OUTPUT_BASE}/hcp/hcp_factors/${ANCESTRY}_expression_outliers.tsv"
            if [ -f "$OUTLIER_FILE" ]; then
                EXTRA_ARGS="${EXTRA_ARGS} --exclude-samples ${OUTLIER_FILE}"
            else
                echo "      WARN: expression outlier list not found: $OUTLIER_FILE"
                echo "             isoforms will keep all samples"
            fi
        fi

        $SING_R "${SCRIPTS_DIR}/combat_normalize_modalities.R" \
            --input "$POOLED_BED" \
            --ancestry-map "$ANCESTRY_MAP" \
            --ancestry "$ANCESTRY" \
            --modality "$MODALITY" \
            --output-dir "$COMBAT_DIR" \
            $EXTRA_ARGS
    done
done

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "[$(date)] Pipeline Complete"
echo "  Pooled BEDs:    ${POOLED_DIR}/"
echo "  ComBat+INT BEDs: ${COMBAT_DIR}/"
echo ""
echo "  ComBat+INT phenotype BEDs (for tensorQTL):"
ls -lh "${COMBAT_DIR}"/*_combat_int.bed 2>/dev/null || echo "    (none found)"
echo ""
echo "  Phenotype groups:"
ls -lh "${COMBAT_DIR}"/*.phenotype_groups.txt 2>/dev/null || echo "    (none found)"
echo ""
echo "  Diagnostics:"
ls -lh "${COMBAT_DIR}"/*_combat_diagnostics.pdf 2>/dev/null || echo "    (none found)"
