#!/bin/bash
# =============================================================================
# 48a_run_fusion.sh — one chromosome x feature-shard FUSION worker
# =============================================================================
# Submitted as a global LSF shard array by 48_fusion_twas.sh. Each array
# element reads one row of SHARD_MANIFEST, runs exactly one chromosome-specific
# subset of the weight models, and writes a durable per-shard TWAS table.
# 48c_merge_fusion.sh concatenates all completed shard tables for the pair.
#
# Env (from 48_fusion_twas.sh): WEIGHT_SET, TRAIT_ID, SUMSTATS_PREFIX,
#   LDREF_PREFIX (chromosome number appended), SHARD_MANIFEST, ISOTWAS_DIR,
#   MIN_R2PRED, FUSION_DIR, RSCRIPT, SCRIPTS_DIR, FORCE_RUN,
#   LSB_JOBINDEX (1..number of feature shards)
# =============================================================================
set -eo pipefail

ARRAY_INDEX="${LSB_JOBINDEX:?ERROR: LSB_JOBINDEX not set; 48a must run as an LSF array element}"
case "$ARRAY_INDEX" in
    ''|*[!0-9]*) echo "ERROR: invalid array index: $ARRAY_INDEX" >&2; exit 1 ;;
esac

SHARD_MANIFEST="${SHARD_MANIFEST:?ERROR: SHARD_MANIFEST not set}"
if [ ! -s "$SHARD_MANIFEST" ]; then
    echo "ERROR: FUSION shard manifest missing/empty: $SHARD_MANIFEST" >&2
    exit 1
fi

ROW=$(awk -F '\t' -v idx="$ARRAY_INDEX" 'NR == idx + 1 {print $2 "\t" $3 "\t" $4 "\t" $5; exit}' "$SHARD_MANIFEST")
if [ -z "$ROW" ]; then
    echo "ERROR: no manifest row for array index ${ARRAY_INDEX}: $SHARD_MANIFEST" >&2
    exit 1
fi
IFS=$'\t' read -r CHR SHARD_INDEX N_FEATURES SHARD_POS <<< "$ROW"

case "$CHR" in
    ''|*[!0-9]*) echo "ERROR: invalid chromosome in manifest row: $ROW" >&2; exit 1 ;;
esac
if [ "$CHR" -lt 1 ] || [ "$CHR" -gt 22 ]; then
    echo "ERROR: chromosome out of range in manifest row: $ROW" >&2
    exit 1
fi
if [ ! -s "$SHARD_POS" ]; then
    echo "ERROR: shard .pos missing/empty: $SHARD_POS" >&2
    exit 1
fi

PLAN_ID="$(basename "$(dirname "$SHARD_MANIFEST")")"
WGT_DIR="${ISOTWAS_DIR}/weights/${WEIGHT_SET}/genes"
OUT_DIR="${ISOTWAS_DIR}/fusion/${WEIGHT_SET}/per_shard/${TRAIT_ID}/${PLAN_ID}"
mkdir -p "$OUT_DIR"
printf -v SHARD_STEM 'chr%02d.shard%04d' "$CHR" "$SHARD_INDEX"
SHARD_OUT="${OUT_DIR}/${SHARD_STEM}.twas.tsv"
DONE="${SHARD_OUT}.done"
FORCE_RUN="${FORCE_RUN:-0}"
SUMSTATS="${SUMSTATS_PREFIX}${CHR}.fusion.tsv"

if [ -f "$DONE" ] && [ -s "$SHARD_OUT" ] && [ "$FORCE_RUN" != "1" ]; then
    echo "${SHARD_STEM}: already complete; skipping ($SHARD_OUT)"
    exit 0
fi

if [ ! -s "$SUMSTATS" ]; then
    echo "ERROR: chromosome-specific FUSION GWAS missing/empty: $SUMSTATS" >&2
    exit 1
fi
if [ ! -f "${LDREF_PREFIX}${CHR}.bed" ]; then
    echo "ERROR: no FUSION LD reference for chr${CHR}: ${LDREF_PREFIX}${CHR}.bed" >&2
    exit 1
fi

# FUSION.assoc_test.R is an external script without an in-R .libPaths() call,
# and bash-level R_LIBS_* is not reliably propagated here. Run it through the
# repo wrapper 48b_fusion_assoc.R, which sets .libPaths() inside R and then
# sources FUSION.assoc_test.R (optparse still sees the flags below via
# commandArgs(trailingOnly=TRUE)). FUSION_DIR is already in the environment.
cd "$FUSION_DIR"

TMP_OUT="${SHARD_OUT}.tmp.${LSB_JOBID:-$$}.${ARRAY_INDEX}"
trap 'rm -f "$TMP_OUT"' EXIT
rm -f "$DONE" "$TMP_OUT"

echo "${SHARD_STEM}: chr${CHR}, ${N_FEATURES} features, weights=$(basename "$SHARD_POS")"
"$RSCRIPT" "${SCRIPTS_DIR}/48b_fusion_assoc.R" \
    --sumstats "$SUMSTATS" \
    --weights "$SHARD_POS" \
    --weights_dir "$WGT_DIR" \
    --ref_ld_chr "$LDREF_PREFIX" \
    --chr "$CHR" \
    --min_r2pred "$MIN_R2PRED" \
    --out "$TMP_OUT"

if [ ! -s "$TMP_OUT" ]; then
    echo "ERROR: FUSION returned success but ${SHARD_STEM} output is missing/empty: $TMP_OUT" >&2
    exit 1
fi

mv -f "$TMP_OUT" "$SHARD_OUT"
touch "$DONE"
trap - EXIT

N_ROWS=$(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "$SHARD_OUT")
STATUS_COL=$(awk -F '\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="FUSION.STATUS") {print i; exit}}' "$SHARD_OUT")
if [ -n "$STATUS_COL" ]; then
    N_FAILED=$(awk -F '\t' -v c="$STATUS_COL" 'NR>1 && $c!="ok" {n++} END{print n+0}' "$SHARD_OUT")
    N_OK=$((N_ROWS - N_FAILED))
    echo "${SHARD_STEM}: wrote $SHARD_OUT (${N_OK} successful, ${N_FAILED} failed/skipped phenotype rows)"
else
    echo "${SHARD_STEM}: wrote $SHARD_OUT (${N_ROWS} model rows; legacy output without FUSION.STATUS)"
fi
