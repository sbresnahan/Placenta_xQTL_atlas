#!/bin/bash
# =============================================================================
# 53_submit_gxe.sh — Orchestrate the Objective 2.1 GxE pipeline
# =============================================================================
# Stages (STAGES env var selects a subset; default all):
#   0  build pooled inputs        — 50_build_gxe_inputs.py (one job)
#   1  pooled HCPs + covariates   — 51_pooled_hcp.sh array [1..n_modality]
#   2  tier-1 cis-perm scans      — 53a array [1..22] per modality x exposure
#   3  merge + Storey q-values    — 53a MODE=merge per modality x exposure
#   4  tier-2 ancestry-stratified — 53b POST_MODE=tier2 per modality x exposure
#   5  sensitivity + T/NT         — 53b POST_MODE=sensitivity,tnt (status-
#                                   gated: OFF/BLOCKED until the metadata
#                                   columns / maternal pgens exist)
#   6  aggregate                  — 53b POST_MODE=aggregate
#
# Usage:
#   TEST=1 bash 53_submit_gxe.sh    # stage 2-3 pilot: expression x GA chr21
#   bash 53_submit_gxe.sh           # everything
#   STAGES="2 3" MODALITIES="expression" bash 53_submit_gxe.sh
#
# Required env: CONFIG, SCRIPTS_DIR, METADATA (cohort metadata .txt).
# Optional: OUTPUT_BASE, QTL_DIR, RESULTS_DIR, GXE_DIR, ANCESTRIES,
#   MODALITIES, EXPOSURES (default: enabled rows of gxe_config.tsv),
#   QUEUE, WALLTIME, CHROMS, FORCE, R_PACKAGE_LIB.
# =============================================================================
set -euo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
METADATA="${METADATA:?ERROR: METADATA env var required (placenta_QTL_cohort_metadata.txt)}"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
GXE_DIR="${GXE_DIR:-${RESULTS_DIR}/gxe}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
MODALITIES="${MODALITIES:-expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability}"
EXPOSURES="${EXPOSURES:-$(awk -F'\t' 'NR>1 && $6==1 {print $1}' "${SCRIPTS_DIR}/gxe_config.tsv")}"
STAGES="${STAGES:-0 1 2 3 4 5 6}"
QUEUE="${QUEUE:-medium}"
TEST="${TEST:-0}"
FORCE="${FORCE:-0}"
CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
R_PACKAGE_LIB="${R_PACKAGE_LIB:-/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1}"

if [ "$TEST" = "1" ]; then
    MODALITIES="expression"; CHROMS="21"; STAGES="${STAGES:-2 3}"
    echo "TEST=1: expression x first exposure, chr21 only"
fi

N_CHROMS=$(echo $CHROMS | wc -w)
N_MODS=$(echo $MODALITIES | wc -w)
mkdir -p "$LOG_DIR"

echo "=== 53_submit_gxe.sh ==="
echo "  GXE_DIR:    $GXE_DIR"
echo "  MODALITIES: $MODALITIES"
echo "  EXPOSURES:  $EXPOSURES"
echo "  STAGES:     $STAGES"

has_stage() { case " $STAGES " in *" $1 "*) return 0;; *) return 1;; esac; }

# ---- Stage 0: build pooled inputs ----------------------------------------
DEP0=""
if has_stage 0; then
    if [ -f "${GXE_DIR}/inputs/pooled_covariates_base.tsv" ] && [ "$FORCE" != "1" ]; then
        echo "  stage 0: inputs exist — skipping (FORCE=1 to rebuild)"
    else
        J=$(bsub -q "$QUEUE" -n 4 -M 32G -R "rusage[mem=32G]" -W 4:00 \
            -J "gxe_inputs" \
            -o "${LOG_DIR}/gxe_inputs.%J.out" -e "${LOG_DIR}/gxe_inputs.%J.err" \
            -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,ANCESTRIES=$ANCESTRIES,MODALITIES=$MODALITIES" \
            "python3 ${SCRIPTS_DIR}/50_build_gxe_inputs.py --qtl-dir $QTL_DIR --results-dir $RESULTS_DIR --metadata $METADATA --ancestries '$ANCESTRIES' --modalities '$MODALITIES' --gxe-dir $GXE_DIR")
        DEP0=$(echo "$J" | grep -o '[0-9]\+' | head -1)
        echo "  stage 0: gxe_inputs job $DEP0"
    fi
fi

# ---- Stage 1: pooled HCPs -------------------------------------------------
DEP1="$DEP0"
if has_stage 1; then
    DEP=""
    [ -n "$DEP0" ] && DEP="-w done($DEP0)"
    J=$(bsub -q "$QUEUE" -n 4 -M 32G -R "rusage[mem=32G]" -W 8:00 \
        -J "gxe_hcp[1-${N_MODS}]" $DEP \
        -o "${LOG_DIR}/gxe_hcp.%J.%I.out" -e "${LOG_DIR}/gxe_hcp.%J.%I.err" \
        -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,ANCESTRIES=$ANCESTRIES,MODALITIES=$MODALITIES,METADATA=$METADATA,FORCE=$FORCE,R_PACKAGE_LIB=$R_PACKAGE_LIB" \
        "bash ${SCRIPTS_DIR}/51_pooled_hcp.sh")
    DEP1=$(echo "$J" | grep -o '[0-9]\+' | head -1)
    echo "  stage 1: gxe_hcp array $DEP1"
fi

# ---- Stages 2-3: tier-1 scans + merge (per modality x exposure) -----------
declare -a MERGE_JOBS=()
if has_stage 2 || has_stage 3; then
    for MOD in $MODALITIES; do
        for EXP in $EXPOSURES; do
            DEP=""
            [ -n "$DEP1" ] && DEP="-w done($DEP1)"
            SCAN_DEP=""
            if has_stage 2; then
                J=$(bsub -q "$QUEUE" -n 4 -M 48G -R "rusage[mem=48G]" -W 12:00 \
                    -J "gxe_t1_${MOD}_${EXP}[1-${N_CHROMS}]" $DEP \
                    -o "${LOG_DIR}/gxe_t1_${MOD}_${EXP}.%J.%I.out" \
                    -e "${LOG_DIR}/gxe_t1_${MOD}_${EXP}.%J.%I.err" \
                    -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,ANCESTRIES=$ANCESTRIES,MODALITY=$MOD,EXPOSURE=$EXP,MODE=scan,FORCE=$FORCE" \
                    "bash ${SCRIPTS_DIR}/53a_run_gxe_scan.sh")
                SCAN_DEP=$(echo "$J" | grep -o '[0-9]\+' | head -1)
                echo "  stage 2: gxe_t1_${MOD}_${EXP} array $SCAN_DEP"
            fi
            if has_stage 3; then
                MDEP=""
                [ -n "$SCAN_DEP" ] && MDEP="-w done($SCAN_DEP)"
                J=$(bsub -q "$QUEUE" -n 2 -M 16G -R "rusage[mem=16G]" -W 2:00 \
                    -J "gxe_merge_${MOD}_${EXP}" $MDEP \
                    -o "${LOG_DIR}/gxe_merge_${MOD}_${EXP}.%J.out" \
                    -e "${LOG_DIR}/gxe_merge_${MOD}_${EXP}.%J.err" \
                    -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,MODALITY=$MOD,EXPOSURE=$EXP,MODE=merge,R_PACKAGE_LIB=$R_PACKAGE_LIB" \
                    "bash ${SCRIPTS_DIR}/53a_run_gxe_scan.sh")
                MERGE_JOBS+=($(echo "$J" | grep -o '[0-9]\+' | head -1))
                echo "  stage 3: gxe_merge_${MOD}_${EXP} job ${MERGE_JOBS[-1]}"
            fi
        done
    done
fi

# ---- Stage 4: tier-2 ancestry-stratified ----------------------------------
if has_stage 4; then
    for MOD in $MODALITIES; do
        for EXP in $EXPOSURES; do
            J=$(bsub -q "$QUEUE" -n 4 -M 32G -R "rusage[mem=32G]" -W 8:00 \
                -J "gxe_t2_${MOD}_${EXP}" \
                -o "${LOG_DIR}/gxe_t2_${MOD}_${EXP}.%J.out" \
                -e "${LOG_DIR}/gxe_t2_${MOD}_${EXP}.%J.err" \
                -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,ANCESTRIES=$ANCESTRIES,MODALITY=$MOD,EXPOSURE=$EXP,POST_MODE=tier2" \
                "bash ${SCRIPTS_DIR}/53b_run_post.sh")
            echo "  stage 4: gxe_t2_${MOD}_${EXP} job $(echo "$J" | grep -o '[0-9]\+' | head -1)"
        done
    done
fi

# ---- Stage 5: sensitivity + T/NT (status-gated) ---------------------------
if has_stage 5; then
    for PM in sensitivity tnt; do
        J=$(bsub -q "$QUEUE" -n 2 -M 16G -R "rusage[mem=16G]" -W 2:00 \
            -J "gxe_${PM}" \
            -o "${LOG_DIR}/gxe_${PM}.%J.out" -e "${LOG_DIR}/gxe_${PM}.%J.err" \
            -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,ANCESTRIES=$ANCESTRIES,POST_MODE=$PM" \
            "bash ${SCRIPTS_DIR}/53b_run_post.sh")
        echo "  stage 5: gxe_${PM} job $(echo "$J" | grep -o '[0-9]\+' | head -1)"
    done
fi

# ---- Stage 6: aggregate ----------------------------------------------------
if has_stage 6; then
    ADEP=""
    if [ ${#MERGE_JOBS[@]} -gt 0 ]; then
        ADEP="-w $(IFS=,; echo "${MERGE_JOBS[*]/%/}" | sed 's/\([0-9]*\)/done(\1)/g' | tr ',' ' ' | sed 's/ / \&\& /g')"
    fi
    J=$(bsub -q "$QUEUE" -n 2 -M 16G -R "rusage[mem=16G]" -W 2:00 \
        -J "gxe_aggregate" $ADEP \
        -o "${LOG_DIR}/gxe_aggregate.%J.out" -e "${LOG_DIR}/gxe_aggregate.%J.err" \
        -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,POST_MODE=aggregate" \
        "bash ${SCRIPTS_DIR}/53b_run_post.sh")
    echo "  stage 6: gxe_aggregate job $(echo "$J" | grep -o '[0-9]\+' | head -1)"
fi

echo "=== all stages submitted ==="
