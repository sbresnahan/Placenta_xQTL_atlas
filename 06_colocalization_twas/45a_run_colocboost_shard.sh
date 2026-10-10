#!/bin/bash
# =============================================================================
# 45a_run_colocboost_shard.sh — LSF array element: prepare + run colocBoost
# =============================================================================
# Invoked by bsub via 45_submit_colocboost.sh with env:
#   REGIONS, OUTCOMES, N_SHARDS, CB_DIR, LD_XQTL_PGEN, LD_GWAS_PGEN,
#   LD_GWAS_KEEP, FORCE_RUN, RSCRIPT, SCRIPTS_DIR
# $LSB_JOBINDEX selects the shard.
#
# External command-line tools run ONLY in the host-side Python preparation
# phase. 44_colocboost.R consumes prepared TSV/PLINK .raw files and never
# invokes tabix, plink2, zcat, or any shell command.
# =============================================================================
set -eo pipefail

source /etc/profile.d/modules.sh
module load plink samtools

if ! command -v tabix >/dev/null 2>&1; then
    echo "  tabix not found from samtools module; stacking conda samtools-1.16.1"
    eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
    conda activate --stack samtools-1.16.1
fi
TABIX_BIN="$(command -v tabix || true)"
PLINK2_BIN="$(command -v plink2 || true)"
[ -n "$TABIX_BIN" ] || { echo "ERROR: tabix unavailable after module + conda fallback" >&2; exit 1; }
[ -n "$PLINK2_BIN" ] || { echo "ERROR: plink2 not found after 'module load plink'" >&2; exit 1; }

SHARD_INDEX="${LSB_JOBINDEX:?ERROR: LSB_JOBINDEX not set (submit as a job array)}"

# FORCE_RUN: clear .done sentinels for this shard's regions.
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

PREP_DIR="${CB_DIR}/prepared/${LSB_JOBID:-manual}_${SHARD_INDEX}"
MANIFEST="${PREP_DIR}/shard-$(printf '%04d' "$SHARD_INDEX").tsv"
rm -rf "$PREP_DIR"
mkdir -p "$PREP_DIR"

python3 "${SCRIPTS_DIR}/44_prepare_colocboost_inputs.py" \
    --regions "$REGIONS" \
    --outcomes "$OUTCOMES" \
    --shard-index "$SHARD_INDEX" \
    --n-shards "$N_SHARDS" \
    --outdir "$CB_DIR" \
    --work-dir "$PREP_DIR" \
    --manifest "$MANIFEST" \
    --ld-xqtl-pgen "$LD_XQTL_PGEN" \
    --ld-gwas-pgen "$LD_GWAS_PGEN" \
    --ld-gwas-keep "${LD_GWAS_KEEP:-}" \
    --tabix "$TABIX_BIN" \
    --plink2 "$PLINK2_BIN"

"$RSCRIPT" "${SCRIPTS_DIR}/44_colocboost.R" \
    --prepared-manifest "$MANIFEST" \
    --outdir "$CB_DIR"

if [ "${KEEP_PREP:-0}" != "1" ]; then
    rm -rf "$PREP_DIR"
fi
