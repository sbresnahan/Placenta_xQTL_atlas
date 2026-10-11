#!/bin/bash
# =============================================================================
# 48_fusion_twas.sh — Run FUSION.assoc_test.R per weight set x GWAS trait
# =============================================================================
# Objective 1.6: apply the trained isoTWAS/TWAS weights (46/47) to every
# harmonized GWAS (37/38) with ancestry-matched 1KG LD references.
#
# Step 1 converts the already-prebuilt ancestry/chromosome 1KG PGEN references
# from 40b/40c to the PLINK BED format required by FUSION. The full 1KG PGEN
# and ancestry keep files are not touched here.
# Step 2 converts each harmonized GWAS to FUSION format (SNP A1 A2 Z BETA SE).
# Step 3 merges shard .pos fragments into one .pos per weight set.
# Step 4 submits one 22-element LSF chromosome array per weight set x trait.
# Each 48a array element runs one chromosome. A dependent 48c job runs after
# the whole array ends and publishes the final table only if all 22 succeeded.
#
# LD reference ancestry: the trait's primary_ancestry when EAS/EUR; traits
# labeled "both" use the weight set's ancestry (EUR for the pooled set).
#
# Prerequisites: 40b/40c compact 1KG refs, 46/47 weights, FUSION clone (36).
#
# Usage:
#   TEST=1 bash 48_fusion_twas.sh   # prep + ONE pair (22-chr array + merge)
#   bash 48_fusion_twas.sh
#
# Optional overrides:
#   WEIGHT_SETS         — default "EAS EUR pooled"
#   QUEUE               — default medium
#   WALLTIME            — per-chromosome walltime, default 04:00
#   MEM_GB              — memory per chromosome array element, default 128 GB
#   FUSION_MAX_PARALLEL — max concurrent chromosomes per pair, default 4
#   MERGE_WALLTIME      — dependent concatenation walltime, default 00:30
#   MERGE_MEM_GB        — dependent concatenation memory, default 4 GB
#   MIN_R2PRED          — FUSION GWAS Z imputation r2 filter, default 0.7
#   FORCE_RUN=1         — rerun completed chromosomes/final pair outputs
# =============================================================================
set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
REPO_ROOT="$(cd "${SCRIPTS_DIR}/.." && pwd)"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
COLOC_DIR="${COLOC_DIR:-${RESULTS_DIR}/coloc}"
ISOTWAS_DIR="${ISOTWAS_DIR:-${RESULTS_DIR}/isotwas}"
GWAS_DIR="${GWAS_DIR:-${OUTPUT_BASE}/gwas}"
KG_LD_REF_DIR="${KG_LD_REF_DIR:-${COLOC_DIR}/ld_reference/1kg}"
FUSION_DIR="${FUSION_DIR:-/rsrch5/home/epi/bhattacharya_lab/software/fusion_twas}"
LDREF_DIR="${LDREF_DIR:-${GWAS_DIR}/fusion_ldref}"
WEIGHT_SETS="${WEIGHT_SETS:-EAS EUR pooled}"
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-04:00}"
MEM_GB="${MEM_GB:-128}"
FUSION_MAX_PARALLEL="${FUSION_MAX_PARALLEL:-4}"
MERGE_WALLTIME="${MERGE_WALLTIME:-00:30}"
MERGE_MEM_GB="${MERGE_MEM_GB:-4}"
MIN_R2PRED="${MIN_R2PRED:-0.7}"
TEST="${TEST:-0}"
FORCE_RUN="${FORCE_RUN:-0}"
RSCRIPT="${RSCRIPT:-${REPO_ROOT}/bin/Rscript_sif}"
PLINK2="${PLINK2:-plink2}"
N_CHR=22

case "$FUSION_MAX_PARALLEL" in
    ''|*[!0-9]*)
        echo "ERROR: FUSION_MAX_PARALLEL must be an integer in [1, ${N_CHR}]" >&2
        exit 1
        ;;
esac
if [ "$FUSION_MAX_PARALLEL" -lt 1 ] || [ "$FUSION_MAX_PARALLEL" -gt "$N_CHR" ]; then
    echo "ERROR: FUSION_MAX_PARALLEL must be in [1, ${N_CHR}] (got ${FUSION_MAX_PARALLEL})" >&2
    exit 1
fi

# --- Toolchain ---------------------------------------------------------------
# Step 1 builds the 1KG LD references with plink2 on the login node. Load the
# seadragon module when plink2 is not already on PATH.
if ! command -v "$PLINK2" >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh
    module load plink
fi
if ! command -v "$PLINK2" >/dev/null 2>&1; then
    echo "ERROR: '$PLINK2' not on PATH (tried 'module load plink')." >&2
    echo "  Load a module/env providing plink2 or set PLINK2=/path/to/plink2." >&2
    exit 1
fi

mkdir -p "$LDREF_DIR" "${GWAS_DIR}/fusion" "$LOG_DIR"

echo "=== 48_fusion_twas.sh ==="
echo "  WEIGHT_SETS:         $WEIGHT_SETS"
echo "  ISOTWAS_DIR:         $ISOTWAS_DIR"
echo "  FUSION_DIR:          $FUSION_DIR"
echo "  LDREF_DIR:           $LDREF_DIR"
echo "  CHROMOSOME ARRAY:    1-${N_CHR}, max ${FUSION_MAX_PARALLEL} concurrent per pair"
echo "  CHROMOSOME MEMORY:   ${MEM_GB}G per array element"
echo "  MERGE MEMORY:        ${MERGE_MEM_GB}G"

# --- Step 1: FUSION BED views of the prebuilt 1KG references ----------------
for ANC in EAS EUR; do
    for CHR in $(seq 1 "$N_CHR"); do
        SRC="${KG_LD_REF_DIR}/${ANC}/chr${CHR}"
        if [ ! -f "${SRC}.done" ] || [ ! -s "${SRC}.pgen" ]; then
            echo "ERROR: compact 1KG reference missing: ${SRC}; run 40b first" >&2
            exit 1
        fi
        OUT="${LDREF_DIR}/${ANC}.1kg.chr${CHR}"
        if [ -f "${OUT}.bed" ]; then continue; fi
        echo "  building FUSION BED view: ${ANC} chr${CHR}"
        "$PLINK2" --pfile "$SRC" --make-bed --out "$OUT" --silent
    done
done

# --- Step 2: FUSION-format GWAS sumstats --------------------------------------
TRAITS=()
while IFS=$'\t' read -r TRAIT_ID PRIMARY_ANC; do
    HARM="${GWAS_DIR}/${TRAIT_ID}.sumstats.tsv.gz"
    [ -f "$HARM" ] || continue
    FUS="${GWAS_DIR}/fusion/${TRAIT_ID}.fusion.tsv"
    if [ ! -f "$FUS" ]; then
        echo "  converting $TRAIT_ID to FUSION format"
        # harmonized schema: chr pos rsid effect_allele other_allele eaf
        # beta se pval n var_id trait_id  ->  FUSION: SNP A1 A2 Z BETA SE N
        zcat "$HARM" | awk 'NR==1{print "SNP\tA1\tA2\tZ\tBETA\tSE\tN"}
                            NR>1{ if ($7!="NA" && $8!="NA" && $8>0)
                                    printf "%s\t%s\t%s\t%.6g\t%.6g\t%.6g\t%d\n",
                                           $11, $4, $5, $7/$8, $7, $8, $10 }' \
            > "$FUS"
    fi
    TRAITS+=("${TRAIT_ID}:${PRIMARY_ANC}")
done < <(tail -n +2 "${SCRIPTS_DIR}/gwas_catalog.tsv" | cut -f1,12)

if [ "${#TRAITS[@]}" -eq 0 ]; then
    echo "  ERROR: no harmonized GWAS files found in $GWAS_DIR" >&2
    exit 1
fi

# --- Step 3: merge shard .pos fragments ---------------------------------------
for WS in $WEIGHT_SETS; do
    POS_DIR="${ISOTWAS_DIR}/weights/${WS}"
    POS_OUT="${ISOTWAS_DIR}/weights/${WS}.pos"
    if ls "${POS_DIR}"/shard-*.pos >/dev/null 2>&1; then
        head -1 "$(ls "${POS_DIR}"/shard-*.pos | head -1)" > "$POS_OUT"
        for SHARD_POS in "${POS_DIR}"/shard-*.pos; do
            tail -n +2 "$SHARD_POS" >> "$POS_OUT"
        done
        echo "  ${WS}: $(( $(wc -l < "$POS_OUT") - 1 )) models in $(basename "$POS_OUT")"
    else
        echo "  WARNING: no shard .pos files for ${WS} (has 47 run?)"
    fi
done

# --- Step 4: submit chromosome arrays + dependent merges ----------------------
N=0
for WS in $WEIGHT_SETS; do
    POS_OUT="${ISOTWAS_DIR}/weights/${WS}.pos"
    [ -s "$POS_OUT" ] || continue
    for TRAIT_SPEC in "${TRAITS[@]}"; do
        TRAIT_ID="${TRAIT_SPEC%%:*}"
        PRIMARY_ANC="${TRAIT_SPEC##*:}"
        # LD reference ancestry
        if [ "$PRIMARY_ANC" = "EAS" ] || [ "$PRIMARY_ANC" = "EUR" ]; then
            LD_ANC="$PRIMARY_ANC"
        elif [ "$WS" != "pooled" ]; then
            LD_ANC="$WS"
        else
            LD_ANC="EUR"
        fi
        if [ ! -f "${LDREF_DIR}/${LD_ANC}.1kg.chr1.bed" ]; then
            echo "  SKIP ${WS} x ${TRAIT_ID}: no LD ref for ${LD_ANC}"
            continue
        fi
        OUT="${ISOTWAS_DIR}/fusion/${WS}/${WS}_${TRAIT_ID}.twas.tsv"
        if [ -f "$OUT" ] && [ "$FORCE_RUN" != "1" ]; then
            echo "  SKIP ${WS} x ${TRAIT_ID}: output exists"
            continue
        fi

        JOB="fusion_${WS}_${TRAIT_ID}"
        MERGE_JOB="fusionmerge_${WS}_${TRAIT_ID}"
        if bjobs -J "$JOB" 2>/dev/null | grep -q "$JOB" || \
           bjobs -J "$MERGE_JOB" 2>/dev/null | grep -q "$MERGE_JOB"; then
            echo "  SKIP ${WS} x ${TRAIT_ID}: chromosome array or merge running/pending"
            continue
        fi

        ENV_STR="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},REPO_ROOT=${REPO_ROOT},ISOTWAS_DIR=${ISOTWAS_DIR},WEIGHT_SET=${WS},TRAIT_ID=${TRAIT_ID},SUMSTATS=${GWAS_DIR}/fusion/${TRAIT_ID}.fusion.tsv,LDREF_PREFIX=${LDREF_DIR}/${LD_ANC}.1kg.chr,MIN_R2PRED=${MIN_R2PRED},FUSION_DIR=${FUSION_DIR},RSCRIPT=${RSCRIPT},FORCE_RUN=${FORCE_RUN}"
        SUBMIT_OUT=$(bsub -J "${JOB}[1-${N_CHR}]%${FUSION_MAX_PARALLEL}" \
             -q "$QUEUE" -n 1 -W "$WALLTIME" \
             -M "$MEM_GB" -R "rusage[mem=${MEM_GB}]" \
             -o "${LOG_DIR}/${JOB}.%J_%I.out" -e "${LOG_DIR}/${JOB}.%J_%I.err" \
             -env "$ENV_STR" < "${SCRIPTS_DIR}/48a_run_fusion.sh")
        echo "$SUBMIT_OUT"
        ARRAY_JOB_ID=$(printf '%s\n' "$SUBMIT_OUT" | sed -n 's/.*Job <\([0-9][0-9]*\).*/\1/p' | head -1)
        if [ -z "$ARRAY_JOB_ID" ]; then
            echo "ERROR: could not parse LSF job ID for ${JOB}" >&2
            exit 1
        fi

        MERGE_ENV_STR="ISOTWAS_DIR=${ISOTWAS_DIR},WEIGHT_SET=${WS},TRAIT_ID=${TRAIT_ID}"
        bsub -J "$MERGE_JOB" -q "$QUEUE" -n 1 -W "$MERGE_WALLTIME" \
             -M "$MERGE_MEM_GB" -R "rusage[mem=${MERGE_MEM_GB}]" \
             -w "ended(${ARRAY_JOB_ID})" \
             -o "${LOG_DIR}/${MERGE_JOB}.%J.out" -e "${LOG_DIR}/${MERGE_JOB}.%J.err" \
             -env "$MERGE_ENV_STR" < "${SCRIPTS_DIR}/48c_merge_fusion.sh"

        echo "  submitted ${WS} x ${TRAIT_ID}: chr1-${N_CHR} array (max ${FUSION_MAX_PARALLEL}) + dependent merge; LD ref: ${LD_ANC}"
        N=$((N + 1))
        if [ "$TEST" = "1" ]; then
            echo ""
            echo "TEST=1: submitted one pair (${WS} x ${TRAIT_ID}) as a chromosome array + merge."
            exit 0
        fi
    done
done

echo ""
echo "Submitted $N FUSION chromosome arrays (+ dependent merge jobs)."
