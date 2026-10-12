#!/bin/bash
# Force-submit FUSION merge jobs for every weight-set x trait pair that has a
# per-shard directory for the current shard plan. Existing merged outputs and
# existing LSF merge jobs are intentionally ignored.
#
# Requires 48c_merge_fusion.sh with ALLOW_MISSING_SHARDS=1 support.
set -euo pipefail

SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR not set}"
ISOTWAS_DIR="${ISOTWAS_DIR:?ERROR: ISOTWAS_DIR not set}"
LOG_DIR="${LOG_DIR:?ERROR: LOG_DIR not set}"
WEIGHT_SETS="${WEIGHT_SETS:-EAS EUR pooled}"
QUEUE="${QUEUE:-medium}"
MERGE_WALLTIME="${MERGE_WALLTIME:-12:00}"
MERGE_MEM_GB="${MERGE_MEM_GB:-4}"
DRY_RUN="${DRY_RUN:-0}"

MERGER="${SCRIPTS_DIR}/48c_merge_fusion.sh"
[ -f "$MERGER" ] || { echo "ERROR: missing merger: $MERGER" >&2; exit 1; }
mkdir -p "$LOG_DIR"

N=0
for WS in $WEIGHT_SETS; do
    POINTER="${ISOTWAS_DIR}/weights/${WS}/fusion_shards/latest_manifest.txt"
    if [ ! -s "$POINTER" ]; then
        echo "SKIP ${WS}: missing manifest pointer: $POINTER" >&2
        continue
    fi

    MANIFEST="$(head -n 1 "$POINTER")"
    if [ ! -s "$MANIFEST" ]; then
        echo "SKIP ${WS}: manifest missing/empty: $MANIFEST" >&2
        continue
    fi

    PLAN_ID="$(basename "$(dirname "$MANIFEST")")"
    TRAIT_ROOT="${ISOTWAS_DIR}/fusion/${WS}/per_shard"
    if [ ! -d "$TRAIT_ROOT" ]; then
        echo "SKIP ${WS}: no per-shard trait directory: $TRAIT_ROOT" >&2
        continue
    fi

    for TRAIT_DIR in "$TRAIT_ROOT"/*; do
        [ -d "$TRAIT_DIR" ] || continue
        TRAIT_ID="$(basename "$TRAIT_DIR")"
        PAIR_DIR="${TRAIT_DIR}/${PLAN_ID}"
        [ -d "$PAIR_DIR" ] || continue

        # A merge can tolerate any number of missing shards, but at least one
        # completed/nonempty shard is needed to establish the FUSION output schema.
        if ! find "$PAIR_DIR" -maxdepth 1 -type f -name '*.twas.tsv' -size +0c -print -quit | grep -q .; then
            echo "SKIP ${WS} x ${TRAIT_ID}: no nonempty shard output exists at all" >&2
            continue
        fi

        MERGE_JOB="fusionmerge_${WS}_${TRAIT_ID}"
        MERGE_ENV="ISOTWAS_DIR=${ISOTWAS_DIR},WEIGHT_SET=${WS},TRAIT_ID=${TRAIT_ID},SHARD_MANIFEST=${MANIFEST},ALLOW_MISSING_SHARDS=1"

        CMD=(
            bsub
            -J "$MERGE_JOB"
            -q "$QUEUE"
            -n 1
            -W "$MERGE_WALLTIME"
            -M "$MERGE_MEM_GB"
            -R "rusage[mem=${MERGE_MEM_GB}]"
            -o "${LOG_DIR}/redo_${MERGE_JOB}.%J.out"
            -e "${LOG_DIR}/redo_${MERGE_JOB}.%J.err"
            -env "$MERGE_ENV"
        )

        echo "SUBMIT ${WS} x ${TRAIT_ID} (ALLOW_MISSING_SHARDS=1)"
        if [ "$DRY_RUN" = "1" ]; then
            printf '  '
            printf '%q ' "${CMD[@]}"
            printf '< %q\n' "$MERGER"
        else
            "${CMD[@]}" < "$MERGER"
        fi
        N=$((N + 1))
    done
done

echo "Submitted ${N} forced FUSION merge jobs."
