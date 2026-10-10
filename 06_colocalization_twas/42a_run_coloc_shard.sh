#!/bin/bash
# =============================================================================
# 42a_run_coloc_shard.sh — LSF array element: prepare + run one coloc shard
# =============================================================================
# Invoked by bsub via 42_submit_coloc.sh with env:
#   TASKS, N_SHARDS, COLOC_DIR, MIN_VARIANTS, PP_H4, FORCE_RUN, RSCRIPT,
#   SCRIPTS_DIR
# $LSB_JOBINDEX selects the shard.
#
# External command-line tools run ONLY in the host-side Python preparation
# phase. 41_susie_coloc.R runs via $RSCRIPT (Singularity) and consumes only
# prepared TSV/LD files; R never invokes tabix, plink/plink2, zcat, or shell commands.
# =============================================================================
set -eo pipefail

source /etc/profile.d/modules.sh
module load plink samtools

# Capture the module-provided PLINK2 before conda changes PATH.  The installed
# PLINK v2.00a3.6LM can export dosages but cannot produce bulk LD matrices.
PLINK2_BIN="$(command -v plink2 || true)"
[ -n "$PLINK2_BIN" ] || { echo "ERROR: plink2 not found after 'module load plink'" >&2; exit 1; }

# The host-side preparer computes signed LD from PLINK2-exported dosages with
# NumPy. Explicitly activate tensorqtl so the batch job uses the same supported
# Python environment as the mapping pipeline instead of the old system Python.
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl
PYTHON_BIN="$(command -v python3 || true)"
[ -n "$PYTHON_BIN" ] || { echo "ERROR: python3 unavailable in tensorqtl env" >&2; exit 1; }
"$PYTHON_BIN" -c 'import numpy' >/dev/null 2>&1 || {
    echo "ERROR: NumPy unavailable in tensorqtl env" >&2
    exit 1
}

# Some samtools module builds do not expose tabix. Match the existing pipeline
# convention and stack the known conda environment only when needed.
if ! command -v tabix >/dev/null 2>&1; then
    echo "  tabix not found from samtools module; stacking conda samtools-1.16.1"
    conda activate --stack samtools-1.16.1
fi
TABIX_BIN="$(command -v tabix || true)"
[ -n "$TABIX_BIN" ] || { echo "ERROR: tabix unavailable after module + conda fallback" >&2; exit 1; }

SHARD_INDEX="${LSB_JOBINDEX:?ERROR: LSB_JOBINDEX not set (submit as a job array)}"

# FORCE_RUN: clear .done sentinels for this shard's tasks so the worker
# re-runs them (worker skips tasks with existing .done otherwise).
# Stdlib python only (csv/re/pathlib) — no conda Python dependency.
if [ "${FORCE_RUN:-0}" = "1" ]; then
    "$PYTHON_BIN" - "$TASKS" "$SHARD_INDEX" "$N_SHARDS" "${COLOC_DIR}/results" <<'PYEOF'
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

# Keep prepared artifacts on the shared /rsrch filesystem so they are visible
# inside the Singularity R container. Remove them after success unless requested.
PREP_DIR="${COLOC_DIR}/prepared/susie/${LSB_JOBID:-manual}_${SHARD_INDEX}"
PREPARED="${PREP_DIR}/shard-$(printf '%04d' "$SHARD_INDEX").tsv"
rm -rf "$PREP_DIR"
mkdir -p "$PREP_DIR"

"$PYTHON_BIN" "${SCRIPTS_DIR}/41_prepare_susie_coloc_inputs.py" \
    --tasks "$TASKS" \
    --shard-index "$SHARD_INDEX" \
    --n-shards "$N_SHARDS" \
    --work-dir "${PREP_DIR}/work" \
    --out "$PREPARED" \
    --tabix "$TABIX_BIN" \
    --plink2 "$PLINK2_BIN" \
    --min-variants "${MIN_VARIANTS:-50}"

"$RSCRIPT" "${SCRIPTS_DIR}/41_susie_coloc.R" \
    --prepared-tasks "$PREPARED" \
    --shard-index "$SHARD_INDEX" \
    --n-shards "$N_SHARDS" \
    --outdir "$COLOC_DIR" \
    --min-variants "${MIN_VARIANTS:-50}" \
    --pp-h4 "${PP_H4:-0.7}"

if [ "${KEEP_PREP:-0}" != "1" ]; then
    rm -rf "$PREP_DIR"
fi
