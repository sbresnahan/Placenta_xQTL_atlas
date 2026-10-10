#!/bin/bash
# =============================================================================
# 40c_run_1kg_ld_reference.sh — LSF worker for one ancestry x chromosome
# =============================================================================
set -eo pipefail

source /etc/profile.d/modules.sh
module load plink
PLINK2_BIN="$(command -v plink2 || true)"
[ -n "$PLINK2_BIN" ] || { echo "ERROR: plink2 not found after module load plink" >&2; exit 1; }

IDX="${LSB_JOBINDEX:?ERROR: LSB_JOBINDEX not set}"
IFS=':' read -r -a ANC_ARR <<< "${ANCESTRIES_CSV:?ERROR: ANCESTRIES_CSV not set}"
N_ANC="${#ANC_ARR[@]}"
TOTAL=$((N_ANC * 22))
if [ "$IDX" -lt 1 ] || [ "$IDX" -gt "$TOTAL" ]; then
    echo "ERROR: array index $IDX outside 1-$TOTAL" >&2
    exit 1
fi
ANC_I=$(( (IDX - 1) / 22 ))
CHR=$(( (IDX - 1) % 22 + 1 ))
ANC="${ANC_ARR[$ANC_I]}"

KEEP="${KG_LD_REF_DIR:?ERROR: KG_LD_REF_DIR not set}/keep/${ANC}.keep"
[ -s "$KEEP" ] || { echo "ERROR: missing/empty keep file: $KEEP" >&2; exit 1; }
OUT_DIR="${KG_LD_REF_DIR:?ERROR: KG_LD_REF_DIR not set}/${ANC}"
PREFIX="${OUT_DIR}/chr${CHR}"
DONE="${PREFIX}.done"
mkdir -p "$OUT_DIR"

if [ -f "$DONE" ] && [ "${FORCE_REF:-0}" != "1" ] \
   && [ -s "${PREFIX}.pgen" ] \
   && { [ -s "${PREFIX}.pvar" ] || [ -s "${PREFIX}.pvar.zst" ]; } \
   && [ -s "${PREFIX}.psam" ]; then
    echo "SKIP: ${ANC} chr${CHR} reference already complete"
    exit 0
fi

TMP="${PREFIX}.tmp.${LSB_JOBID:-manual}_${IDX}"
rm -f "${TMP}.pgen" "${TMP}.pvar" "${TMP}.pvar.zst" "${TMP}.psam" "${TMP}.log"

echo "[$(date)] build 1KG LD reference: ${ANC} chr${CHR}"
echo "  source: $KG_PGEN"
echo "  keep:   $KEEP"
echo "  out:    $PREFIX"

"$PLINK2_BIN" \
    --pfile "${KG_PGEN:?ERROR: KG_PGEN not set}" \
    --keep "$KEEP" \
    --chr "$CHR" \
    --max-alleles 2 \
    --make-pgen \
    --threads 1 \
    --memory "${PLINK_MEMORY_MB:-12000}" \
    --out "$TMP"

[ -s "${TMP}.pgen" ] || { echo "ERROR: missing ${TMP}.pgen" >&2; exit 1; }
[ -s "${TMP}.psam" ] || { echo "ERROR: missing ${TMP}.psam" >&2; exit 1; }
if [ -s "${TMP}.pvar" ]; then
    PVAR_SRC="${TMP}.pvar"
    PVAR_DST="${PREFIX}.pvar"
elif [ -s "${TMP}.pvar.zst" ]; then
    PVAR_SRC="${TMP}.pvar.zst"
    PVAR_DST="${PREFIX}.pvar.zst"
else
    echo "ERROR: missing ${TMP}.pvar(.zst)" >&2
    exit 1
fi

rm -f "${PREFIX}.pgen" "${PREFIX}.pvar" "${PREFIX}.pvar.zst" "${PREFIX}.psam" "$DONE"
mv "${TMP}.pgen" "${PREFIX}.pgen"
mv "$PVAR_SRC" "$PVAR_DST"
mv "${TMP}.psam" "${PREFIX}.psam"
rm -f "${TMP}.log"
{
    echo -e "ancestry\t${ANC}"
    echo -e "chrom\t${CHR}"
    echo -e "source_pgen\t${KG_PGEN}"
    echo -e "keep_file\t${KEEP}"
    echo -e "built_at\t$(date -Is)"
} > "${PREFIX}.metadata.tsv"
touch "$DONE"

echo "[$(date)] DONE: ${ANC} chr${CHR}"
