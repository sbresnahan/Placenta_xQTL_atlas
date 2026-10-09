#!/bin/bash
# =============================================================================
# 42a_run_coloc_shard.sh — LSF array element: run one coloc shard
# =============================================================================
# Invoked by bsub via 42_submit_coloc.sh with env:
#   TASKS, N_SHARDS, COLOC_DIR, MIN_VARIANTS, PP_H4, FORCE_RUN, RSCRIPT,
#   SCRIPTS_DIR
# $LSB_JOBINDEX selects the shard.
#
# Compute-node environment: seadragon modules 'plink' (plink2) and 'samtools'
# (tabix; 41_susie_coloc.R slices tabix-indexed summary stats and calls plink2
# for LD). R runs via $RSCRIPT (singularity wrapper); the R script sets its
# own .libPaths() — bash-level R_LIBS_* is not relied upon.
# =============================================================================
set -eo pipefail

source /etc/profile.d/modules.sh
module load plink samtools

SHARD_INDEX="${LSB_JOBINDEX:?ERROR: LSB_JOBINDEX not set (submit as a job array)}"

# FORCE_RUN: clear .done sentinels for this shard's tasks so the worker
# re-runs them (worker skips tasks with existing .done otherwise).
# Stdlib python only (csv/re/pathlib) — no conda env needed on compute nodes.
if [ "${FORCE_RUN:-0}" = "1" ]; then
    python3 - "$TASKS" "$SHARD_INDEX" "$N_SHARDS" "${COLOC_DIR}/results" <<'PYEOF'
import csv
import re
import sys
from pathlib import Path

tasks_path, shard_index, n_shards, results_dir = (
    sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), Path(sys.argv[4]))
with open(tasks_path, newline="") as fh:
    tasks = list(csv.DictReader(fh, delimiter="\t"))
rows = [i for i in range(len(tasks)) if (i % n_shards) + 1 == shard_index]
def sanitize(x):
    return re.sub(r"[^A-Za-z0-9._+-]", "_", str(x))
removed = 0
for i in rows:
    t = tasks[i]
    base = f"{t['ancestry']}_{sanitize(t['phenotype_id'])}_{t['trait_id']}"
    done = results_dir / t["modality"] / f"{base}.done"
    if done.exists():
        done.unlink()
        removed += 1
print(f"FORCE_RUN: cleared {removed} .done sentinels for shard {shard_index}")
PYEOF
fi

"$RSCRIPT" "${SCRIPTS_DIR}/41_susie_coloc.R" \
    --tasks "$TASKS" \
    --shard-index "$SHARD_INDEX" \
    --n-shards "$N_SHARDS" \
    --outdir "$COLOC_DIR" \
    --min-variants "${MIN_VARIANTS:-50}" \
    --pp-h4 "${PP_H4:-0.7}"
