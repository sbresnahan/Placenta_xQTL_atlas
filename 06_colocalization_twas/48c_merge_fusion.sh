#!/bin/bash
# =============================================================================
# 48c_merge_fusion.sh — concatenate completed FUSION chromosome outputs
# =============================================================================
# Submitted by 48_fusion_twas.sh after the entire 48a chromosome array ends.
# Refuses to publish a final table unless all 22 chromosome outputs have
# success sentinels. This makes failed arrays restartable instead of leaving a
# permanently pending merge dependency.
#
# Env: ISOTWAS_DIR, WEIGHT_SET, TRAIT_ID
# =============================================================================
set -eo pipefail

PER_CHR_DIR="${ISOTWAS_DIR}/fusion/${WEIGHT_SET}/per_chr/${TRAIT_ID}"
OUT_DIR="${ISOTWAS_DIR}/fusion/${WEIGHT_SET}"
FINAL="${OUT_DIR}/${WEIGHT_SET}_${TRAIT_ID}.twas.tsv"
mkdir -p "$OUT_DIR"

TMP_FINAL="${FINAL}.tmp.${LSB_JOBID:-$$}"
trap 'rm -f "$TMP_FINAL"' EXIT
rm -f "$TMP_FINAL"

FIRST=1
for CHR in $(seq 1 22); do
    CHR_OUT="${PER_CHR_DIR}/chr${CHR}.twas.tsv"
    DONE="${CHR_OUT}.done"
    if [ ! -f "$DONE" ] || [ ! -s "$CHR_OUT" ]; then
        echo "ERROR: chr${CHR} is incomplete: expected $CHR_OUT and $DONE" >&2
        exit 1
    fi
    if [ "$FIRST" = "1" ]; then
        cat "$CHR_OUT" > "$TMP_FINAL"
        FIRST=0
    else
        tail -n +2 "$CHR_OUT" >> "$TMP_FINAL"
    fi
done

if [ "$FIRST" = "1" ] || [ ! -s "$TMP_FINAL" ]; then
    echo "ERROR: no chromosome outputs were merged for ${WEIGHT_SET} x ${TRAIT_ID}" >&2
    exit 1
fi

mv -f "$TMP_FINAL" "$FINAL"
trap - EXIT

echo "wrote $FINAL ($(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "$FINAL") model rows from 22 chromosomes)"
