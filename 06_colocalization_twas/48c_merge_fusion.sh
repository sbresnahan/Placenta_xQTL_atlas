#!/bin/bash
# =============================================================================
# 48c_merge_fusion.sh — concatenate completed FUSION feature-shard outputs
# =============================================================================
# Submitted by 48_fusion_twas.sh after the entire 48a global shard array ends.
# Refuses to publish a final table unless every row in SHARD_MANIFEST has a
# completed output + success sentinel. This makes partial failures restartable
# instead of leaving a permanently pending merge dependency.
#
# Env: ISOTWAS_DIR, WEIGHT_SET, TRAIT_ID, SHARD_MANIFEST
# =============================================================================
set -eo pipefail

SHARD_MANIFEST="${SHARD_MANIFEST:?ERROR: SHARD_MANIFEST not set}"
if [ ! -s "$SHARD_MANIFEST" ]; then
    echo "ERROR: FUSION shard manifest missing/empty: $SHARD_MANIFEST" >&2
    exit 1
fi

PLAN_ID="$(basename "$(dirname "$SHARD_MANIFEST")")"
PER_SHARD_DIR="${ISOTWAS_DIR}/fusion/${WEIGHT_SET}/per_shard/${TRAIT_ID}/${PLAN_ID}"
OUT_DIR="${ISOTWAS_DIR}/fusion/${WEIGHT_SET}"
FINAL="${OUT_DIR}/${WEIGHT_SET}_${TRAIT_ID}.twas.tsv"
mkdir -p "$OUT_DIR"

TMP_FINAL="${FINAL}.tmp.${LSB_JOBID:-$$}"
trap 'rm -f "$TMP_FINAL"' EXIT
rm -f "$TMP_FINAL"

FIRST=1
N_EXPECTED=0
while IFS=$'\t' read -r ARRAY_INDEX CHR SHARD_INDEX N_FEATURES SHARD_POS; do
    [ "$ARRAY_INDEX" = "array_index" ] && continue
    N_EXPECTED=$((N_EXPECTED + 1))
    printf -v SHARD_STEM 'chr%02d.shard%04d' "$CHR" "$SHARD_INDEX"
    SHARD_OUT="${PER_SHARD_DIR}/${SHARD_STEM}.twas.tsv"
    DONE="${SHARD_OUT}.done"
    if [ ! -f "$DONE" ] || [ ! -s "$SHARD_OUT" ]; then
        echo "ERROR: ${SHARD_STEM} is incomplete: expected $SHARD_OUT and $DONE" >&2
        exit 1
    fi
    if [ "$FIRST" = "1" ]; then
        cat "$SHARD_OUT" > "$TMP_FINAL"
        FIRST=0
    else
        tail -n +2 "$SHARD_OUT" >> "$TMP_FINAL"
    fi
done < "$SHARD_MANIFEST"

if [ "$N_EXPECTED" -eq 0 ] || [ "$FIRST" = "1" ] || [ ! -s "$TMP_FINAL" ]; then
    echo "ERROR: no feature-shard outputs were merged for ${WEIGHT_SET} x ${TRAIT_ID}" >&2
    exit 1
fi

mv -f "$TMP_FINAL" "$FINAL"
trap - EXIT

echo "wrote $FINAL ($(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "$FINAL") model rows from ${N_EXPECTED} chromosome x feature shards)"
