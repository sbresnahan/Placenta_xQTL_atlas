#!/bin/bash
# =============================================================================
# 51_pooled_hcp.sh — stage Module-05 optimized covariates (one ancestry x modality)
# =============================================================================
# Historical filename retained for compatibility. Module-07 ancestry-specific
# scans must reuse Module-05's finalized/optimized per-modality covariates rather
# than re-estimating HCPs. This worker subsets and reorders
#   qtl_inputs/{ANC}_covariates_{MOD}.tsv
# to the exact sample columns in
#   gxe/inputs/{ANC}_{MOD}.bed.gz
# and writes
#   gxe/inputs/{ANC}_covariates_{MOD}.tsv.
#
# Required env: CONFIG, SCRIPTS_DIR, ANCESTRY. MODALITY may be provided
# directly or resolved from LSB_JOBINDEX against MODALITIES.
# =============================================================================
set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
ANCESTRY="${ANCESTRY:?ERROR: ANCESTRY env var required}"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
GXE_DIR="${GXE_DIR:-${RESULTS_DIR}/gxe}"
INPUTS="${GXE_DIR}/inputs"
MODALITIES="${MODALITIES:-expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability}"
CONDA_EXE="${CONDA_EXE:-/risapps/rhel8/miniforge3/24.5.0-0/bin/conda}"
CONDA_ENV="${CONDA_ENV:-tensorqtl}"

source /etc/profile.d/modules.sh
[ -x "$CONDA_EXE" ] || { echo "ERROR: conda executable not found/executable: $CONDA_EXE"; exit 1; }
eval "$("$CONDA_EXE" shell.bash hook)"
conda activate "$CONDA_ENV"
command -v python3 >/dev/null || { echo "ERROR: python3 not found after conda activation"; exit 1; }

if [ -z "${MODALITY:-}" ]; then
    MODALITY="$(echo $MODALITIES | awk -v i="${LSB_JOBINDEX:?set MODALITY or run under LSF}" '{print $i}')"
fi
[ -n "$MODALITY" ] || { echo "ERROR: unable to resolve modality"; exit 1; }

echo "=== Stage 1: reuse Module-05 covariates: ${ANCESTRY} ${MODALITY} ==="

BED_GZ="${INPUTS}/${ANCESTRY}_${MODALITY}.bed.gz"
SRC_COV="${QTL_DIR}/${ANCESTRY}_covariates_${MODALITY}.tsv"
FINAL_COV="${INPUTS}/${ANCESTRY}_covariates_${MODALITY}.tsv"

for f in "$BED_GZ" "$SRC_COV"; do
    [ -f "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

if [ -f "$FINAL_COV" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "  ${FINAL_COV} exists — skipping (FORCE=1 to restage)"
    exit 0
fi

python3 - "$BED_GZ" "$SRC_COV" "$FINAL_COV" <<'PY'
import gzip
import os
import sys
import tempfile

import pandas as pd

bed_path, src_path, out_path = sys.argv[1:]

with gzip.open(bed_path, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
if len(header) < 5 or header[:4] != ["#chr", "start", "end", "phenotype_id"]:
    raise SystemExit(f"ERROR: unexpected BED header in {bed_path}: {header[:4]}")
samples = header[4:]
if not samples:
    raise SystemExit(f"ERROR: no sample columns in {bed_path}")
if len(samples) != len(set(samples)):
    dup = pd.Index(samples)[pd.Index(samples).duplicated()].unique().tolist()
    raise SystemExit(f"ERROR: duplicated BED sample IDs in {bed_path}: {dup[:5]}")

cov = pd.read_csv(src_path, sep="\t", index_col=0)
cov.columns = cov.columns.astype(str)
cov.index = cov.index.astype(str)
if cov.columns.duplicated().any():
    dup = cov.columns[cov.columns.duplicated()].unique().tolist()
    raise SystemExit(f"ERROR: duplicated sample columns in {src_path}: {dup[:5]}")
if cov.index.duplicated().any():
    dup = cov.index[cov.index.duplicated()].unique().tolist()
    raise SystemExit(f"ERROR: duplicated covariate rows in {src_path}: {dup[:5]}")

missing = [s for s in samples if s not in cov.columns]
if missing:
    raise SystemExit(
        f"ERROR: {len(missing)} Module-07 BED samples are absent from Module-05 "
        f"optimized covariates {src_path}: {missing[:10]}"
    )

extra = [s for s in cov.columns if s not in set(samples)]
staged = cov.loc[:, samples].copy()
if staged.isna().any().any():
    where = staged.isna().stack()
    first = where[where].index[0]
    raise SystemExit(f"ERROR: NaN in staged Module-05 covariates at {first}")

os.makedirs(os.path.dirname(out_path), exist_ok=True)
fd, tmp = tempfile.mkstemp(prefix=os.path.basename(out_path) + ".", suffix=".tmp",
                           dir=os.path.dirname(out_path))
os.close(fd)
try:
    staged.index.name = cov.index.name or "covariate"
    staged.to_csv(tmp, sep="\t", float_format="%.10g")
    os.replace(tmp, out_path)
finally:
    if os.path.exists(tmp):
        os.unlink(tmp)

hcp_n = sum(str(x).startswith("HCP_") for x in staged.index)
print(f"  source: {src_path}")
print(f"  BED samples: {len(samples)}")
print(f"  Module-05 covariates retained: {staged.shape[0]} ({hcp_n} HCP rows)")
print(f"  extra Module-05 sample columns dropped: {len(extra)}")
print(f"  wrote: {out_path}")
PY

echo "=== done: ${ANCESTRY} ${MODALITY} ==="
