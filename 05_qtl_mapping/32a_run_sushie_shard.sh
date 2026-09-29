#!/bin/bash
# =============================================================================
# 32a_run_sushie_shard.sh — LSF job wrapper: run one SuSHiE locus shard
# =============================================================================
# Submitted by 32_submit_sushie.sh via bsub -env with:
#   CONFIG, SCRIPTS_DIR, QTL_DIR, RESULTS_DIR, LOCUS_SHARD, ANCESTRIES
#   (comma-separated), THREADS, FORCE_RUN
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
QTL_DIR="${QTL_DIR:?ERROR: QTL_DIR env var required}"
RESULTS_DIR="${RESULTS_DIR:?ERROR: RESULTS_DIR env var required}"
LOCUS_SHARD="${LOCUS_SHARD:?ERROR: LOCUS_SHARD env var required}"
ANCESTRIES="${ANCESTRIES:-EAS,EUR}"
THREADS="${THREADS:-2}"
FORCE_RUN="${FORCE_RUN:-0}"

# tensorqtl conda env (sushie installed there — see runbook); same setup as
# scripts 27/28 — proven on this cluster.
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# conda create -n sushie python=3.11
# conda activate sushie
# pip install sushie pandas pyarrow fastparquet

conda activate sushie

# comma-separated ANCESTRIES -> space-separated for argparse
ANC_ARGS=$(echo "$ANCESTRIES" | tr ',' ' ')

FORCE_ARGS=""
if [ "$FORCE_RUN" = "1" ]; then
    FORCE_ARGS="--force"
fi

echo "Running shard: $LOCUS_SHARD (ancestries: $ANC_ARGS, threads: $THREADS)"
python3 "${SCRIPTS_DIR}/31_sushie_finemap.py" run \
    --loci-file "$LOCUS_SHARD" \
    --qtl-dir "$QTL_DIR" \
    --out-dir "${RESULTS_DIR}/finemap" \
    --ancestries $ANC_ARGS \
    --threads "$THREADS" \
    $FORCE_ARGS
echo "Shard complete: $LOCUS_SHARD"
