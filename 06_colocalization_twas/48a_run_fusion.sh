#!/bin/bash
# =============================================================================
# 48a_run_fusion.sh — one FUSION.assoc_test.R chromosome worker
# =============================================================================
# Submitted as a 1-22 LSF array by 48_fusion_twas.sh. Each array element runs
# exactly one chromosome and writes a durable per-chromosome TWAS table. The
# dependent 48c_merge_fusion.sh job concatenates them after the whole array
# completes successfully.
#
# Env (from 48_fusion_twas.sh): WEIGHT_SET, TRAIT_ID, SUMSTATS, LDREF_PREFIX
#   (per-chromosome plink1 BED prefix, chromosome number appended),
#   ISOTWAS_DIR, MIN_R2PRED, FUSION_DIR, RSCRIPT, SCRIPTS_DIR, FORCE_RUN,
#   LSB_JOBINDEX (chromosome 1-22)
# =============================================================================
set -eo pipefail

CHR="${LSB_JOBINDEX:?ERROR: LSB_JOBINDEX not set; 48a must run as an LSF array element}"
case "$CHR" in
    ''|*[!0-9]*) echo "ERROR: invalid chromosome index: $CHR" >&2; exit 1 ;;
esac
if [ "$CHR" -lt 1 ] || [ "$CHR" -gt 22 ]; then
    echo "ERROR: chromosome index out of range: $CHR" >&2
    exit 1
fi

POS="${ISOTWAS_DIR}/weights/${WEIGHT_SET}.pos"
WGT_DIR="${ISOTWAS_DIR}/weights/${WEIGHT_SET}/genes"
OUT_DIR="${ISOTWAS_DIR}/fusion/${WEIGHT_SET}/per_chr/${TRAIT_ID}"
mkdir -p "$OUT_DIR"
CHR_OUT="${OUT_DIR}/chr${CHR}.twas.tsv"
DONE="${CHR_OUT}.done"
FORCE_RUN="${FORCE_RUN:-0}"

if [ -f "$DONE" ] && [ -s "$CHR_OUT" ] && [ "$FORCE_RUN" != "1" ]; then
    echo "chr${CHR}: already complete; skipping ($CHR_OUT)"
    exit 0
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

TMP_OUT="${CHR_OUT}.tmp.${LSB_JOBID:-$$}.${CHR}"
trap 'rm -f "$TMP_OUT"' EXIT
rm -f "$DONE" "$TMP_OUT"

"$RSCRIPT" "${SCRIPTS_DIR}/48b_fusion_assoc.R" \
    --sumstats "$SUMSTATS" \
    --weights "$POS" \
    --weights_dir "$WGT_DIR" \
    --ref_ld_chr "$LDREF_PREFIX" \
    --chr "$CHR" \
    --min_r2pred "$MIN_R2PRED" \
    --out "$TMP_OUT"

if [ ! -s "$TMP_OUT" ]; then
    echo "ERROR: FUSION returned success but chr${CHR} output is missing/empty: $TMP_OUT" >&2
    exit 1
fi

mv -f "$TMP_OUT" "$CHR_OUT"
touch "$DONE"
trap - EXIT

echo "chr${CHR}: wrote $CHR_OUT ($(awk 'END{print (NR > 0 ? NR - 1 : 0)}' "$CHR_OUT") model rows)"
