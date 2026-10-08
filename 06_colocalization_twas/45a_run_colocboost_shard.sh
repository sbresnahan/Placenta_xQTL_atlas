#!/bin/bash
# =============================================================================
# 45a_run_colocboost_shard.sh — LSF array element: run one colocBoost shard
# =============================================================================
# Invoked by bsub via 45_submit_colocboost.sh with env:
#   REGIONS, OUTCOMES, N_SHARDS, CB_DIR, LD_XQTL_PGEN, LD_GWAS_PGEN,
#   LD_GWAS_KEEP, FORCE_RUN, RSCRIPT, SCRIPTS_DIR
# $LSB_JOBINDEX selects the shard.
# =============================================================================
set -euo pipefail

SHARD_INDEX="${LSB_JOBINDEX:?ERROR: LSB_JOBINDEX not set (submit as a job array)}"

# FORCE_RUN: clear .done sentinels for this shard's regions
if [ "${FORCE_RUN:-0}" = "1" ]; then
    python3 - "$REGIONS" "$SHARD_INDEX" "$N_SHARDS" "${CB_DIR}/results" <<'PYEOF'
import sys
from pathlib import Path
import pandas as pd

regions_path, shard_index, n_shards, results_dir = (
    sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), Path(sys.argv[4]))
regions = pd.read_csv(regions_path, sep="\t")
rows = [i for i in range(len(regions)) if (i % n_shards) + 1 == shard_index]
removed = 0
for i in rows:
    anc = regions.iloc[i]["ancestry"]
    done = results_dir / anc / f"{regions.iloc[i]['region_id']}.done"
    if done.exists():
        done.unlink()
        removed += 1
print(f"FORCE_RUN: cleared {removed} .done sentinels for shard {shard_index}")
PYEOF
fi

"$RSCRIPT" "${SCRIPTS_DIR}/44_colocboost.R" \
    --regions "$REGIONS" \
    --outcomes "$OUTCOMES" \
    --shard-index "$SHARD_INDEX" \
    --n-shards "$N_SHARDS" \
    --outdir "$CB_DIR" \
    --ld-xqtl-pgen "$LD_XQTL_PGEN" \
    --ld-gwas-pgen "$LD_GWAS_PGEN" \
    --ld-gwas-keep "${LD_GWAS_KEEP:-}"
