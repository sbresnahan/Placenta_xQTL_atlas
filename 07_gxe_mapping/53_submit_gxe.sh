#!/bin/bash
# =============================================================================
# 53_submit_gxe.sh — Orchestrate ancestry-first Objective 2.1 GxE mapping
# =============================================================================
# Stages:
#   0  ancestry-specific inputs       — 50_build_gxe_inputs.py
#   1  ancestry-specific HCPs         — 51_pooled_hcp.sh per ancestry array
#   2  tier-1 ancestry cis-perm scans — chr arrays per ancestry x modality x exposure
#   3  ancestry merge + Storey q      — one merge per ancestry x modality x exposure
#   4  cross-ancestry tier-1 synthesis— ACAT feature evidence + exact-lead IVW
#   5  Aim-1 prioritized tier 2       — ancestry-stratified nominal + IVW
#   6  sensitivity + T/NT             — status-gated
#   7  aggregate atlas tables         — 57_aggregate_gxe.py
# =============================================================================
set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
GXE_DIR="${GXE_DIR:-${RESULTS_DIR}/gxe}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
MODALITIES="${MODALITIES:-expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability}"
COLLAPSED_ANCESTRY_MAP="${COLLAPSED_ANCESTRY_MAP:-${OUTPUT_BASE}/replicate_collapsed/reports/ancestry_map_collapsed.tsv}"
EXPOSURES="${EXPOSURES:-$(awk -F'\t' 'NR>1 && $6==1 {print $1}' "${SCRIPTS_DIR}/gxe_config.tsv")}"
STAGES="${STAGES:-0 1 2 3 4 5 6 7}"
QUEUE="${QUEUE:-medium}"
TEST="${TEST:-0}"
FORCE="${FORCE:-0}"
CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
R_PACKAGE_LIB="${R_PACKAGE_LIB:-/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1}"
CONDA_EXE="${CONDA_EXE:-/risapps/rhel8/miniforge3/24.5.0-0/bin/conda}"
CONDA_ENV="${CONDA_ENV:-tensorqtl}"

if [ "$TEST" = "1" ]; then
    MODALITIES="expression"
    EXPOSURES="$(echo "$EXPOSURES" | awk '{print $1}')"
    CHROMS="21"
    echo "TEST=1: expression x ${EXPOSURES}, chr21"
fi

N_CHROMS=$(echo $CHROMS | wc -w)
N_MODS=$(echo $MODALITIES | wc -w)
mkdir -p "$LOG_DIR"

echo "=== 53_submit_gxe.sh (ancestry-first) ==="
echo "  GXE_DIR:    $GXE_DIR"
echo "  ANCESTRIES: $ANCESTRIES"
echo "  MODALITIES: $MODALITIES"
echo "  EXPOSURES:  $EXPOSURES"
echo "  STAGES:     $STAGES"
echo "  COLLAPSED_ANCESTRY_MAP: $COLLAPSED_ANCESTRY_MAP"

has_stage() { case " $STAGES " in *" $1 "*) return 0;; *) return 1;; esac; }
job_id() { grep -o '[0-9]\+' | head -1; }
dep_expr() {
    local out="" id
    for id in "$@"; do
        [ -z "$id" ] && continue
        [ -n "$out" ] && out+="&&"
        out+="done(${id})"
    done
    printf '%s' "$out"
}

# ---- Stage 0 --------------------------------------------------------------
DEP0=""
if has_stage 0; then
    COMPLETE=1
    for ANC in $ANCESTRIES; do
        [ -f "${GXE_DIR}/inputs/${ANC}_covariates_base.tsv" ] || COMPLETE=0
    done
    if [ "$COMPLETE" = "1" ] && [ "$FORCE" != "1" ]; then
        echo "  stage 0: ancestry inputs exist — skipping (FORCE=1 to rebuild)"
    else
        J=$(bsub -q "$QUEUE" -n 4 -M 32G -R "rusage[mem=32G]" -W 4:00 \
            -J "gxe_inputs" \
            -o "${LOG_DIR}/gxe_inputs.%J.out" -e "${LOG_DIR}/gxe_inputs.%J.err" \
            -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,ANCESTRIES=$ANCESTRIES,MODALITIES=$MODALITIES,COLLAPSED_ANCESTRY_MAP=$COLLAPSED_ANCESTRY_MAP,CONDA_EXE=$CONDA_EXE,CONDA_ENV=$CONDA_ENV" \
            "bash ${SCRIPTS_DIR}/50a_run_build_gxe_inputs.sh")
        DEP0=$(echo "$J" | job_id)
        echo "  stage 0: gxe_inputs job $DEP0"
    fi
fi

# ---- Stage 1 --------------------------------------------------------------
declare -a HCP_JOBS=()
if has_stage 1; then
    for ANC in $ANCESTRIES; do
        DEP_ARGS=()
        [ -n "$DEP0" ] && DEP_ARGS=(-w "done($DEP0)")
        J=$(bsub -q "$QUEUE" -n 4 -M 32G -R "rusage[mem=32G]" -W 8:00 \
            -J "gxe_hcp_${ANC}[1-${N_MODS}]" "${DEP_ARGS[@]}" \
            -o "${LOG_DIR}/gxe_hcp_${ANC}.%J.%I.out" \
            -e "${LOG_DIR}/gxe_hcp_${ANC}.%J.%I.err" \
            -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,ANCESTRY=$ANC,MODALITIES=$MODALITIES,FORCE=$FORCE,R_PACKAGE_LIB=$R_PACKAGE_LIB,CONDA_EXE=$CONDA_EXE,CONDA_ENV=$CONDA_ENV" \
            "bash ${SCRIPTS_DIR}/51_pooled_hcp.sh")
        jid=$(echo "$J" | job_id); HCP_JOBS+=("$jid")
        echo "  stage 1: gxe_hcp_${ANC} array $jid"
    done
fi

# ---- Stages 2-3 -----------------------------------------------------------
declare -A MERGE_BY_KEY=()
declare -a MERGE_JOBS=()
HCP_DEP=$(dep_expr "${HCP_JOBS[@]}")
for ANC in $ANCESTRIES; do
    for MOD in $MODALITIES; do
        for EXP in $EXPOSURES; do
            SCAN_DEP=""
            if has_stage 2; then
                DEP_ARGS=(); [ -n "$HCP_DEP" ] && DEP_ARGS=(-w "$HCP_DEP")
                J=$(bsub -q "$QUEUE" -n 4 -M 48G -R "rusage[mem=48G]" -W 12:00 \
                    -J "gxe_t1_${ANC}_${MOD}_${EXP}[1-${N_CHROMS}]" "${DEP_ARGS[@]}" \
                    -o "${LOG_DIR}/gxe_t1_${ANC}_${MOD}_${EXP}.%J.%I.out" \
                    -e "${LOG_DIR}/gxe_t1_${ANC}_${MOD}_${EXP}.%J.%I.err" \
                    -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,ANCESTRY=$ANC,MODALITY=$MOD,EXPOSURE=$EXP,MODE=scan,CHROMS=$CHROMS,FORCE=$FORCE,CONDA_EXE=$CONDA_EXE,CONDA_ENV=$CONDA_ENV" \
                    "bash ${SCRIPTS_DIR}/53a_run_gxe_scan.sh")
                SCAN_DEP=$(echo "$J" | job_id)
                echo "  stage 2: ${ANC} ${MOD} x ${EXP} array $SCAN_DEP"
            fi
            if has_stage 3; then
                DEP_ARGS=(); [ -n "$SCAN_DEP" ] && DEP_ARGS=(-w "done($SCAN_DEP)")
                J=$(bsub -q "$QUEUE" -n 2 -M 16G -R "rusage[mem=16G]" -W 2:00 \
                    -J "gxe_merge_${ANC}_${MOD}_${EXP}" "${DEP_ARGS[@]}" \
                    -o "${LOG_DIR}/gxe_merge_${ANC}_${MOD}_${EXP}.%J.out" \
                    -e "${LOG_DIR}/gxe_merge_${ANC}_${MOD}_${EXP}.%J.err" \
                    -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,ANCESTRY=$ANC,MODALITY=$MOD,EXPOSURE=$EXP,MODE=merge,R_PACKAGE_LIB=$R_PACKAGE_LIB,CONDA_EXE=$CONDA_EXE,CONDA_ENV=$CONDA_ENV" \
                    "bash ${SCRIPTS_DIR}/53a_run_gxe_scan.sh")
                jid=$(echo "$J" | job_id)
                MERGE_BY_KEY["${ANC}|${MOD}|${EXP}"]="$jid"
                MERGE_JOBS+=("$jid")
                echo "  stage 3: merge ${ANC} ${MOD} x ${EXP} job $jid"
            fi
        done
    done
done

# ---- Stage 4: cross-ancestry tier-1 synthesis -----------------------------
declare -a META_JOBS=()
if has_stage 4; then
    for MOD in $MODALITIES; do
        for EXP in $EXPOSURES; do
            rel=()
            for ANC in $ANCESTRIES; do
                id="${MERGE_BY_KEY[${ANC}|${MOD}|${EXP}]:-}"
                [ -n "$id" ] && rel+=("$id")
            done
            dep=$(dep_expr "${rel[@]}")
            DEP_ARGS=(); [ -n "$dep" ] && DEP_ARGS=(-w "$dep")
            J=$(bsub -q "$QUEUE" -n 2 -M 16G -R "rusage[mem=16G]" -W 2:00 \
                -J "gxe_meta_${MOD}_${EXP}" "${DEP_ARGS[@]}" \
                -o "${LOG_DIR}/gxe_meta_${MOD}_${EXP}.%J.out" \
                -e "${LOG_DIR}/gxe_meta_${MOD}_${EXP}.%J.err" \
                -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,ANCESTRIES=$ANCESTRIES,MODALITY=$MOD,EXPOSURE=$EXP,POST_MODE=tier1_meta,FORCE=$FORCE,CONDA_EXE=$CONDA_EXE,CONDA_ENV=$CONDA_ENV" \
                "bash ${SCRIPTS_DIR}/53b_run_post.sh")
            jid=$(echo "$J" | job_id); META_JOBS+=("$jid")
            echo "  stage 4: meta ${MOD} x ${EXP} job $jid"
        done
    done
fi

# ---- Stage 5: prioritized tier 2 ------------------------------------------
declare -a TIER2_JOBS=()
if has_stage 5; then
    meta_dep=$(dep_expr "${META_JOBS[@]}")
    for MOD in $MODALITIES; do
        for EXP in $EXPOSURES; do
            DEP_ARGS=(); [ -n "$meta_dep" ] && DEP_ARGS=(-w "$meta_dep")
            J=$(bsub -q "$QUEUE" -n 4 -M 32G -R "rusage[mem=32G]" -W 8:00 \
                -J "gxe_t2_${MOD}_${EXP}" "${DEP_ARGS[@]}" \
                -o "${LOG_DIR}/gxe_t2_${MOD}_${EXP}.%J.out" \
                -e "${LOG_DIR}/gxe_t2_${MOD}_${EXP}.%J.err" \
                -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,ANCESTRIES=$ANCESTRIES,MODALITY=$MOD,EXPOSURE=$EXP,POST_MODE=tier2,FORCE=$FORCE,CONDA_EXE=$CONDA_EXE,CONDA_ENV=$CONDA_ENV" \
                "bash ${SCRIPTS_DIR}/53b_run_post.sh")
            jid=$(echo "$J" | job_id); TIER2_JOBS+=("$jid")
            echo "  stage 5: tier2 ${MOD} x ${EXP} job $jid"
        done
    done
fi

# ---- Stage 6: sensitivity + T/NT ------------------------------------------
declare -a FOLLOWUP_JOBS=()
if has_stage 6; then
    for PM in sensitivity tnt; do
        J=$(bsub -q "$QUEUE" -n 2 -M 16G -R "rusage[mem=16G]" -W 2:00 \
            -J "gxe_${PM}" \
            -o "${LOG_DIR}/gxe_${PM}.%J.out" -e "${LOG_DIR}/gxe_${PM}.%J.err" \
            -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,QTL_DIR=$QTL_DIR,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,ANCESTRIES=$ANCESTRIES,POST_MODE=$PM,CONDA_EXE=$CONDA_EXE,CONDA_ENV=$CONDA_ENV" \
            "bash ${SCRIPTS_DIR}/53b_run_post.sh")
        jid=$(echo "$J" | job_id); FOLLOWUP_JOBS+=("$jid")
        echo "  stage 6: gxe_${PM} job $jid"
    done
fi

# ---- Stage 7: aggregate ---------------------------------------------------
if has_stage 7; then
    final_ids=("${MERGE_JOBS[@]}" "${META_JOBS[@]}" "${TIER2_JOBS[@]}" "${FOLLOWUP_JOBS[@]}")
    dep=$(dep_expr "${final_ids[@]}")
    DEP_ARGS=(); [ -n "$dep" ] && DEP_ARGS=(-w "$dep")
    J=$(bsub -q "$QUEUE" -n 2 -M 16G -R "rusage[mem=16G]" -W 2:00 \
        -J "gxe_aggregate" "${DEP_ARGS[@]}" \
        -o "${LOG_DIR}/gxe_aggregate.%J.out" -e "${LOG_DIR}/gxe_aggregate.%J.err" \
        -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,OUTPUT_BASE=$OUTPUT_BASE,RESULTS_DIR=$RESULTS_DIR,GXE_DIR=$GXE_DIR,POST_MODE=aggregate,CONDA_EXE=$CONDA_EXE,CONDA_ENV=$CONDA_ENV" \
        "bash ${SCRIPTS_DIR}/53b_run_post.sh")
    echo "  stage 7: aggregate job $(echo "$J" | job_id)"
fi

echo "=== requested stages submitted ==="
