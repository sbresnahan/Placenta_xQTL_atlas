#!/bin/bash
# =============================================================================
# 45a_run_colocboost_shard.sh — LSF array element: run one colocBoost shard
# =============================================================================
# Invoked by bsub via 45_submit_colocboost.sh with env:
#   REGIONS, OUTCOMES, N_SHARDS, CB_DIR, LD_XQTL_PGEN, LD_GWAS_PGEN,
#   LD_GWAS_KEEP, FORCE_RUN, RSCRIPT, SCRIPTS_DIR
# $LSB_JOBINDEX selects the shard.
#
# Compute-node environment: seadragon modules 'plink' (plink2) and 'samtools'
# (tabix; 44_colocboost.R slices tabix-indexed summary stats and calls plink2
# for reference dosages). R runs via $RSCRIPT (singularity wrapper); the R
# script sets its own .libPaths() — bash-level R_LIBS_* is not relied upon.
# =============================================================================
set -eo pipefail

source /etc/profile.d/modules.sh
module load plink samtools

SHARD_INDEX="${LSB_JOBINDEX:?ERROR: LSB_JOBINDEX not set (submit as a job array)}"

# FORCE_RUN: clear .done sentinels for this shard's regions.
# Stdlib python only (csv/pathlib) — no conda env needed on compute nodes.
if [ "${FORCE_RUN:-0}" = "1" ]; then
    python3 - "$REGIONS" "$SHARD_INDEX" "$N_SHARDS" "${CB_DIR}/results" <<'PYEOF'
import csv
import sys
from pathlib import Path

regions_path, shard_index, n_shards, results_dir = (
    sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), Path(sys.argv[4]))
with open(regions_path, newline="") as fh:
    regions = list(csv.DictReader(fh, delimiter="\t"))
rows = [i for i in range(len(regions)) if (i % n_shards) + 1 == shard_index]
removed = 0
for i in rows:
    anc = regions[i]["ancestry"]
    done = results_dir / anc / f"{regions[i]['region_id']}.done"
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
