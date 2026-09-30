#!/bin/bash
# =============================================================================
# 32_submit_sushie.sh — Submit SuSHiE cross-ancestry fine-mapping array jobs
# =============================================================================
# Objective 1.5: joint multi-ancestry fine-mapping of every phenotype with an
# FDR<=5% lead variant in at least one ancestry stratum.
#
# Step 1 builds the per-modality locus list (31_sushie_finemap.py
# prepare-loci: union of significant phenotypes, exact tested windows from
# the mapping parquets, per-locus L from independent-signal counts).
# Step 2 shards each locus list and submits one LSF job per shard; each job
# runs 32a_run_sushie_shard.sh -> 31_sushie_finemap.py run.
#
# Prerequisites:
#   - Scripts 23-28 completed (qtl_inputs has {ANC}_qtl.pgen,
#     {ANC}_{MOD}.bed.gz, {ANC}_covariates_{MOD}.tsv; qtl_results has the
#     cisqtl parquets + top tables from 27/29)
#   - sushie installed into a dedicated conda env (python 3.11):
#       conda create -n sushie python=3.11 && conda activate sushie
#       pip install sushie pandas pyarrow fastparquet
#
# Usage:
#   TEST=1 bash 32_submit_sushie.sh   # prepare loci + submit ONE pilot shard
#   bash 32_submit_sushie.sh          # submit all shards
#
# Optional overrides (environment):
#   ANCESTRIES   — default "EAS EUR" (joint fine-mapping across these strata)
#   MODALITIES   — default all 9 per-modality layers (combined excluded:
#                  it duplicates per-modality phenotypes)
#   SHARD_SIZE   — loci per job, default 50 (~1-2 h per shard expected;
#                  measure with TEST=1 first)
#   QUEUE        — default medium
#   WALLTIME     — default 04:00
#   THREADS      — default 2 (matches -n)
#   FORCE_LOCI=1 — rebuild locus lists even if present
#   FORCE_RUN=1  — re-run loci with existing .done markers
#
# Shards with existing diagnostics or a running/pending job of the same name
# are skipped, so it is safe to rerun this script at any time.
# =============================================================================
set -euo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
MODALITIES="${MODALITIES:-expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability}"
SHARD_SIZE="${SHARD_SIZE:-50}"
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-04:00}"
THREADS="${THREADS:-2}"
TEST="${TEST:-0}"
FORCE_LOCI="${FORCE_LOCI:-0}"
FORCE_RUN="${FORCE_RUN:-0}"

FINEMAP_DIR="${RESULTS_DIR}/finemap"
LOCI_DIR="${FINEMAP_DIR}/loci"
SHARD_DIR="${FINEMAP_DIR}/shards"
mkdir -p "$LOCI_DIR" "$SHARD_DIR" "$LOG_DIR"

echo "=== 32_submit_sushie.sh ==="
echo "  ANCESTRIES:  $ANCESTRIES"
echo "  MODALITIES:  $MODALITIES"
echo "  QTL_DIR:     $QTL_DIR"
echo "  RESULTS_DIR: $RESULTS_DIR"
echo "  SHARD_SIZE:  $SHARD_SIZE   QUEUE: $QUEUE   WALLTIME: $WALLTIME"

# --- Step 1: per-modality locus lists --------------------------------------
for MOD in $MODALITIES; do
    LOCI="${LOCI_DIR}/${MOD}.loci.tsv"
    if [ -f "$LOCI" ] && [ "$FORCE_LOCI" != "1" ]; then
        echo "  loci list exists: $LOCI ($(($(wc -l < "$LOCI") - 1)) loci)"
        continue
    fi
    echo "  building locus list: $MOD"
    python3 "${SCRIPTS_DIR}/31_sushie_finemap.py" prepare-loci \
        --results-dir "$RESULTS_DIR" \
        --ancestries $ANCESTRIES \
        --modality "$MOD" \
        --out "$LOCI"
done

# --- Step 2: shard + submit -------------------------------------------------
N=0
for MOD in $MODALITIES; do
    LOCI="${LOCI_DIR}/${MOD}.loci.tsv"
    if [ ! -s "$LOCI" ]; then
        echo "  SKIP ${MOD}: no loci file"
        continue
    fi
    N_LOCI=$(($(wc -l < "$LOCI") - 1))
    if [ "$N_LOCI" -le 0 ]; then
        echo "  SKIP ${MOD}: 0 significant loci"
        continue
    fi
    # split into shards, preserving the header on each shard file
    MOD_SHARD_DIR="${SHARD_DIR}/${MOD}"
    mkdir -p "$MOD_SHARD_DIR"
    python3 - "$LOCI" "$MOD_SHARD_DIR" "$SHARD_SIZE" <<'PYEOF'
import sys
import pandas as pd
loci_path, out_dir, shard_size = sys.argv[1], sys.argv[2], int(sys.argv[3])
df = pd.read_csv(loci_path, sep="\t")
for i in range(0, len(df), shard_size):
    df.iloc[i:i + shard_size].to_csv(
        f"{out_dir}/shard_{i // shard_size:04d}.tsv", sep="\t", index=False)
print(f"  {loci_path}: {len(df)} loci -> {(len(df) + shard_size - 1) // shard_size} shards")
PYEOF
    for SHARD in "$MOD_SHARD_DIR"/shard_*.tsv; do
        SHARD_NAME="${MOD}_$(basename "$SHARD" .tsv)"
        DIAG="${FINEMAP_DIR}/logs/${SHARD_NAME}.diagnostics.tsv"
        if [ -f "$DIAG" ] && [ "$FORCE_RUN" != "1" ]; then
            echo "  SKIP ${SHARD_NAME}: diagnostics already exist"
            continue
        fi
        if bjobs -J "sushie_${SHARD_NAME}" 2>/dev/null | grep -q "sushie_${SHARD_NAME}"; then
            echo "  SKIP ${SHARD_NAME}: job already running/pending"
            continue
        fi
        ENV_STR="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},QTL_DIR=${QTL_DIR},RESULTS_DIR=${RESULTS_DIR},LOCUS_SHARD=${SHARD},ANCESTRIES='$(echo "$ANCESTRIES" | tr ' ' ',')',THREADS=${THREADS},FORCE_RUN=${FORCE_RUN}"
        bsub -J "sushie_${SHARD_NAME}" -q "$QUEUE" -n "$THREADS" -W "$WALLTIME" \
             -o "${LOG_DIR}/sushie_${SHARD_NAME}.%J.out" \
             -e "${LOG_DIR}/sushie_${SHARD_NAME}.%J.err" \
             -env "$ENV_STR" \
             < "${SCRIPTS_DIR}/32a_run_sushie_shard.sh"
        N=$((N + 1))
        if [ "$TEST" = "1" ]; then
            echo ""
            echo "TEST=1: submitted one pilot shard (sushie_${SHARD_NAME})."
            echo "Check per-locus wall times in its log before submitting the rest:"
            echo "  tail ${LOG_DIR}/sushie_${SHARD_NAME}.<jobid>.out"
            exit 0
        fi
    done
done
echo ""
echo "Submitted $N shard jobs. Aggregate when done:"
echo "  python3 ${SCRIPTS_DIR}/33_aggregate_finemap.py --finemap-dir ${FINEMAP_DIR} --ancestries $ANCESTRIES"
