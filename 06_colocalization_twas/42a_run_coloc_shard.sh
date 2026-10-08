#!/bin/bash
# =============================================================================
# 42a_run_coloc_shard.sh — LSF array element: run one coloc shard
# =============================================================================
# Invoked by bsub via 42_submit_coloc.sh with env:
#   TASKS, N_SHARDS, COLOC_DIR, MIN_VARIANTS, PP_H4, FORCE_RUN, RSCRIPT,
#   SCRIPTS_DIR
# $LSB_JOBINDEX selects the shard.
# =============================================================================
set -euo pipefail

SHARD_INDEX="${LSB_JOBINDEX:?ERROR: LSB_JOBINDEX not set (submit as a job array)}"

# FORCE_RUN: clear .done sentinels for this shard's tasks so the worker
# re-runs them (worker skips tasks with existing .done otherwise).
if [ "${FORCE_RUN:-0}" = "1" ]; then
    python3 - "$TASKS" "$SHARD_INDEX" "$N_SHARDS" "${COLOC_DIR}/results" <<'PYEOF'
import sys
from pathlib import Path
import pandas as pd
import re

tasks_path, shard_index, n_shards, results_dir = (
    sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), Path(sys.argv[4]))
tasks = pd.read_csv(tasks_path, sep="\t")
rows = [i for i in range(len(tasks)) if (i % n_shards) + 1 == shard_index]
def sanitize(x):
    return re.sub(r"[^A-Za-z0-9._+-]", "_", str(x))
removed = 0
for i in rows:
    t = tasks.iloc[i]
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
