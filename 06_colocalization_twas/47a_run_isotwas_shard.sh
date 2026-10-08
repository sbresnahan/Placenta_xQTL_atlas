#!/bin/bash
# =============================================================================
# 47a_run_isotwas_shard.sh — LSF array element: train one isoTWAS shard
# =============================================================================
# Invoked by bsub via 47_submit_isotwas.sh with env:
#   WEIGHT_SET, N_SHARDS, QTL_DIR, ISOTWAS_DIR, R2_MIN,
#   FORCE_RUN, RSCRIPT, SCRIPTS_DIR
# $LSB_JOBINDEX selects the shard.
# =============================================================================
set -euo pipefail

SHARD_INDEX="${LSB_JOBINDEX:?ERROR: LSB_JOBINDEX not set (submit as a job array)}"

# FORCE_RUN: remove existing .wgt.RDat for this shard's genes so they retrain
if [ "${FORCE_RUN:-0}" = "1" ]; then
    DIAG="${ISOTWAS_DIR}/diagnostics/isotwas_${WEIGHT_SET}.shard-$(printf '%04d' "$SHARD_INDEX").diagnostics.tsv"
    if [ -f "$DIAG" ]; then
        cut -f1 "$DIAG" | tail -n +2 | while read -r GENE; do
            rm -f "${ISOTWAS_DIR}/weights/${WEIGHT_SET}/genes/${GENE}.wgt.RDat"
        done
    fi
    rm -f "${ISOTWAS_DIR}/weights/${WEIGHT_SET}/shard-$(printf '%04d' "$SHARD_INDEX").pos"
fi

"$RSCRIPT" "${SCRIPTS_DIR}/46_isotwas_train.R" \
    --weight-set "$WEIGHT_SET" \
    --shard-index "$SHARD_INDEX" \
    --n-shards "$N_SHARDS" \
    --qtl-dir "$QTL_DIR" \
    --outdir "$ISOTWAS_DIR" \
    --r2-min "${R2_MIN:-0.01}"
