#!/bin/bash
# =============================================================================
# 47a_run_isotwas_shard.sh — LSF array element: prepare + train isoTWAS shard
# =============================================================================
# Invoked by bsub via 47_submit_isotwas.sh with env:
#   WEIGHT_SET, N_SHARDS, QTL_DIR, ISOTWAS_DIR, R2_MIN,
#   FORCE_RUN, RSCRIPT, SCRIPTS_DIR
# $LSB_JOBINDEX selects the shard.
#
# plink2 runs ONLY in 46_prepare_isotwas_inputs.py on the host compute node.
# 46_isotwas_train.R reads prepared BED subsets/dosages and performs only
# R-native statistical work; it never invokes command-line tools.
# =============================================================================
set -eo pipefail

source /etc/profile.d/modules.sh
module load plink
PLINK2_BIN="$(command -v plink2 || true)"
[ -n "$PLINK2_BIN" ] || { echo "ERROR: plink2 not found after 'module load plink'" >&2; exit 1; }

SHARD_INDEX="${LSB_JOBINDEX:?ERROR: LSB_JOBINDEX not set (submit as a job array)}"

# FORCE_RUN: remove all existing gene + isoform weight files represented by
# this shard's prior diagnostics so stale isoform weights cannot survive reruns.
if [ "${FORCE_RUN:-0}" = "1" ]; then
    DIAG="${ISOTWAS_DIR}/diagnostics/isotwas_${WEIGHT_SET}.shard-$(printf '%04d' "$SHARD_INDEX").diagnostics.tsv"
    if [ -f "$DIAG" ]; then
        cut -f1 "$DIAG" | tail -n +2 | while read -r GENE; do
            rm -f "${ISOTWAS_DIR}/weights/${WEIGHT_SET}/genes/${GENE}.wgt.RDat"
            rm -f "${ISOTWAS_DIR}/weights/${WEIGHT_SET}/genes/${GENE}__"*.wgt.RDat
        done
    fi
    rm -f "${ISOTWAS_DIR}/weights/${WEIGHT_SET}/shard-$(printf '%04d' "$SHARD_INDEX").pos"
fi

PREP_DIR="${ISOTWAS_DIR}/prepared/${WEIGHT_SET}/${LSB_JOBID:-manual}_${SHARD_INDEX}"
MANIFEST="${PREP_DIR}/shard-$(printf '%04d' "$SHARD_INDEX").tsv"
rm -rf "$PREP_DIR"
mkdir -p "$PREP_DIR"

FORCE_ARG=()
if [ "${FORCE_RUN:-0}" = "1" ]; then FORCE_ARG=(--force); fi
python3 "${SCRIPTS_DIR}/46_prepare_isotwas_inputs.py" \
    --weight-set "$WEIGHT_SET" \
    --shard-index "$SHARD_INDEX" \
    --n-shards "$N_SHARDS" \
    --qtl-dir "$QTL_DIR" \
    --outdir "$ISOTWAS_DIR" \
    --work-dir "$PREP_DIR" \
    --manifest "$MANIFEST" \
    --plink2 "$PLINK2_BIN" \
    "${FORCE_ARG[@]}"

"$RSCRIPT" "${SCRIPTS_DIR}/46_isotwas_train.R" \
    --weight-set "$WEIGHT_SET" \
    --shard-index "$SHARD_INDEX" \
    --prepared-manifest "$MANIFEST" \
    --qtl-dir "$QTL_DIR" \
    --outdir "$ISOTWAS_DIR" \
    --r2-min "${R2_MIN:-0.01}"

if [ "${KEEP_PREP:-0}" != "1" ]; then
    rm -rf "$PREP_DIR"
fi
