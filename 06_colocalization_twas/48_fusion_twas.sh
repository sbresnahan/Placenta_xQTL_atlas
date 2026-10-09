#!/bin/bash
# =============================================================================
# 48_fusion_twas.sh — Run FUSION.assoc_test.R per weight set x GWAS trait
# =============================================================================
# Objective 1.6: apply the trained isoTWAS/TWAS weights (46/47) to every
# harmonized GWAS (37/38) with ancestry-matched 1KG LD references.
#
# Step 1 builds per-chromosome plink1 BED LD references from the 1KG pgen,
# one per superpopulation with a keep file ({COLOC_DIR}/loci/{ANC}.1kg.keep,
# written by 40_prepare_coloc_loci.py).
# Step 2 converts each harmonized GWAS to FUSION format (SNP A1 A2 Z BETA SE).
# Step 3 merges shard .pos fragments into one .pos per weight set.
# Step 4 submits one LSF job per weight set x trait (48a_run_fusion.sh loops
# chromosomes 1-22 and concatenates).
#
# LD reference ancestry: the trait's primary_ancestry when EAS/EUR; traits
# labeled "both" use the weight set's ancestry (EUR for the pooled set).
#
# Prerequisites: 40 (keep files), 46/47 (weights), FUSION clone (36).
#
# Usage:
#   TEST=1 bash 48_fusion_twas.sh   # prep + ONE pair (first WS x first trait)
#   bash 48_fusion_twas.sh
#
# Optional overrides:
#   WEIGHT_SETS  — default "EAS EUR pooled"
#   QUEUE        — default medium;  WALLTIME — default 04:00
#   MIN_R2PRED   — FUSION GWAS Z imputation r2 filter, default 0.7
#   FORCE_RUN=1  — rerun pairs with existing outputs
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
KG_PGEN="${KG_PGEN:-/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/1kGP/1kGP_hg38}"
FUSION_DIR="${FUSION_DIR:-/rsrch5/home/epi/bhattacharya_lab/software/fusion_twas}"
LDREF_DIR="${LDREF_DIR:-${GWAS_DIR}/fusion_ldref}"
WEIGHT_SETS="${WEIGHT_SETS:-EAS EUR pooled}"
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-04:00}"
MIN_R2PRED="${MIN_R2PRED:-0.7}"
TEST="${TEST:-0}"
FORCE_RUN="${FORCE_RUN:-0}"
RSCRIPT="${RSCRIPT:-${REPO_ROOT}/bin/Rscript_sif}"
PLINK2="${PLINK2:-plink2}"

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
echo "  WEIGHT_SETS: $WEIGHT_SETS"
echo "  ISOTWAS_DIR: $ISOTWAS_DIR"
echo "  FUSION_DIR:  $FUSION_DIR"
echo "  LDREF_DIR:   $LDREF_DIR"

# --- Step 1: per-chromosome 1KG LD references --------------------------------
for ANC in EAS EUR; do
    KEEP="${COLOC_DIR}/loci/${ANC}.1kg.keep"
    if [ ! -f "$KEEP" ]; then
        echo "  WARNING: no keep file for ${ANC} ($KEEP); skipping LD ref"
        continue
    fi
    for CHR in $(seq 1 22); do
        OUT="${LDREF_DIR}/${ANC}.1kg.chr${CHR}"
        if [ -f "${OUT}.bed" ]; then continue; fi
        echo "  building LD ref: ${ANC} chr${CHR}"
        "$PLINK2" --pfile "$KG_PGEN" --keep "$KEEP" --chr "$CHR" \
                  --make-bed --out "$OUT" --silent
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

# --- Step 4: submit one job per weight set x trait ----------------------------
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
        if bjobs -J "$JOB" 2>/dev/null | grep -q "$JOB"; then
            echo "  SKIP ${WS} x ${TRAIT_ID}: job running/pending"
            continue
        fi
        ENV_STR="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},REPO_ROOT=${REPO_ROOT},ISOTWAS_DIR=${ISOTWAS_DIR},WEIGHT_SET=${WS},TRAIT_ID=${TRAIT_ID},SUMSTATS=${GWAS_DIR}/fusion/${TRAIT_ID}.fusion.tsv,LDREF_PREFIX=${LDREF_DIR}/${LD_ANC}.1kg.chr,MIN_R2PRED=${MIN_R2PRED},FUSION_DIR=${FUSION_DIR},RSCRIPT=${RSCRIPT}"
        bsub -J "$JOB" -q "$QUEUE" -n 1 -W "$WALLTIME" \
             -o "${LOG_DIR}/${JOB}.%J.out" -e "${LOG_DIR}/${JOB}.%J.err" \
             -env "$ENV_STR" < "${SCRIPTS_DIR}/48a_run_fusion.sh"
        echo "  submitted ${WS} x ${TRAIT_ID} (LD ref: ${LD_ANC})"
        N=$((N + 1))
        if [ "$TEST" = "1" ]; then
            echo ""
            echo "TEST=1: submitted one pair (${WS} x ${TRAIT_ID})."
            exit 0
        fi
    done
done
echo ""
echo "Submitted $N FUSION jobs."
