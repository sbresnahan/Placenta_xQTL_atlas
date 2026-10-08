#!/bin/bash
# =============================================================================
# 48a_run_fusion.sh — one FUSION.assoc_test.R run per chromosome for a
# weight set x GWAS trait pair; concatenates to the final TWAS table.
# =============================================================================
# Env (from 48_fusion_twas.sh): WEIGHT_SET, TRAIT_ID, SUMSTATS, LDREF_PREFIX
#   (per-chromosome plink1 BED prefix, chromosome number appended),
#   ISOTWAS_DIR, MIN_R2PRED, FUSION_DIR, RSCRIPT, SCRIPTS_DIR
# =============================================================================
set -euo pipefail

POS="${ISOTWAS_DIR}/weights/${WEIGHT_SET}.pos"
WGT_DIR="${ISOTWAS_DIR}/weights/${WEIGHT_SET}/genes"
OUT_DIR="${ISOTWAS_DIR}/fusion/${WEIGHT_SET}"
mkdir -p "$OUT_DIR"
FINAL="${OUT_DIR}/${WEIGHT_SET}_${TRAIT_ID}.twas.tsv"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

# FUSION.assoc_test.R sources utils/ via here(): run from FUSION_DIR
cd "$FUSION_DIR"

for CHR in $(seq 1 22); do
    CHR_OUT="${TMP_DIR}/chr${CHR}.tsv"
    if [ ! -f "${LDREF_PREFIX}${CHR}.bed" ]; then
        echo "  WARNING: no LD ref for chr${CHR}; skipped" >&2
        continue
    fi
    "$RSCRIPT" "${FUSION_DIR}/FUSION.assoc_test.R" \
        --sumstats "$SUMSTATS" \
        --weights "$POS" \
        --weights_dir "$WGT_DIR" \
        --ref_ld_chr "$LDREF_PREFIX" \
        --chr "$CHR" \
        --min_r2pred "$MIN_R2PRED" \
        --out "$CHR_OUT"
done

# concatenate per-chromosome tables (header once)
FIRST=1
for CHR in $(seq 1 22); do
    CHR_OUT="${TMP_DIR}/chr${CHR}.tsv"
    [ -f "$CHR_OUT" ] || continue
    if [ "$FIRST" = "1" ]; then
        cat "$CHR_OUT" > "$FINAL"
        FIRST=0
    else
        tail -n +2 "$CHR_OUT" >> "$FINAL"
    fi
done
echo "wrote $FINAL ($(($(wc -l < "$FINAL") - 1)) model rows)"
