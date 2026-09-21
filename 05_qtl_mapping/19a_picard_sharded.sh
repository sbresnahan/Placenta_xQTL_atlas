#!/bin/bash
# =============================================================================
# 19a_picard_sharded.sh — Sharded Stage 1 (PicardTools QC) for 19_hcp_factors.sh
# =============================================================================
# The main wrapper (19_hcp_factors.sh) runs Picard QC serially per cohort with
# no per-sample checkpointing: ~500-700 BAMs x 5 Picard tools (incl.
# MarkDuplicates) will not fit in a 24h LSF job, and a timeout loses ALL of
# Stage 1. This script shards each cohort's samples.txt across an LSF job
# array so chunks run in parallel and each chunk is independently retryable.
#
# Chunk outputs are named  ${COHORT}.chunk${IDX}.qcmetrics.tsv  so they do NOT
# match the main wrapper's Stage-2 glob (*_qc_metrics.tsv). After all chunks
# finish, run this script once with MERGE=1 to write the per-cohort
# ${COHORT}_qc_metrics.tsv files the wrapper expects; the wrapper then skips
# Stage 1 and proceeds to Stages 2-4 (pool QC, pool expression, ComBat+HCP).
#
# Usage:
#   1) Submit one array per cohort (NCHUNKS=16 works for all four):
#        bsub -J "pic_c1[1-16]" -env "CONFIG=...,SCRIPTS_DIR=...,COHORT=cohort1,NCHUNKS=16" < 19a_picard_sharded.sh
#        bsub -J "pic_c2[1-16]" -env "...,COHORT=cohort2,NCHUNKS=16" < 19a_picard_sharded.sh
#        bsub -J "pic_c3[1-16]" -env "...,COHORT=cohort3,NCHUNKS=16" < 19a_picard_sharded.sh
#        bsub -J "pic_c4[1-16]" -env "...,COHORT=cohort4,NCHUNKS=16" < 19a_picard_sharded.sh
#   2) After all arrays finish:
#        CONFIG=... SCRIPTS_DIR=... MERGE=1 bash 19a_picard_sharded.sh
#   3) Run the main wrapper (skips Stage 1):
#        bsub -env "CONFIG=...,SCRIPTS_DIR=...,ANCESTRY_MAP=..." < 19_hcp_factors.sh
#
# Required env vars: CONFIG, SCRIPTS_DIR, COHORT (except MERGE mode: COHORT optional)
# Optional env vars: NCHUNKS (default 16), QC_DIR, REFFLAT, GENE_ANNOT, FASTA
#
# Reference generation (refFlat, gene GC/length) is handled automatically:
# array index 1 generates any missing refs; other indices wait on a sentinel.
# =============================================================================

#BSUB -q long
#BSUB -n 2
#BSUB -M 16
#BSUB -R "rusage[mem=16]"
#BSUB -W 12:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/picard.%J.%I.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/picard.%J.%I.err

set -eo pipefail

# ---- Required env vars ----
CONFIG="${CONFIG:?ERROR: CONFIG env var required (path to config.yml)}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
NCHUNKS="${NCHUNKS:-16}"
IDX="${LSB_JOBINDEX:-1}"

# ---- Global init (same as 19_hcp_factors.sh) ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load config values as shell vars
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}")"

OUTPUT_BASE="${OUTPUT_BASE}"
REFERENCE_DIR="${REFERENCE_DIR}"
NORMALIZED_GTF="${NORMALIZED_GTF}"
REF_GENOME="${REF_GENOME}"

QC_DIR="${QC_DIR:-${OUTPUT_BASE}/hcp/qc_metrics}"
REFFLAT="${REFFLAT:-${REFERENCE_DIR}/HPLRv2.refFlat}"
GENE_ANNOT="${GENE_ANNOT:-${OUTPUT_BASE}/hcp/gene_gc_length.tsv}"
FASTA="${FASTA:-${REF_GENOME}}"
mkdir -p "$QC_DIR"
export QC_DIR

# =============================================================================
# MERGE mode: combine chunk outputs into the per-cohort files the main
# wrapper expects, then exit. Run on a login node after all arrays finish.
# =============================================================================
if [ "${MERGE:-0}" = "1" ]; then
    python3 - <<'PYEOF'
import glob, os
import pandas as pd

qc_dir = os.environ["QC_DIR"]
chunk_files = sorted(glob.glob(os.path.join(qc_dir, "*.chunk*.qcmetrics.tsv")))
if not chunk_files:
    raise SystemExit(f"ERROR: no chunk files found in {qc_dir}")

cohorts = sorted(set(os.path.basename(f).split(".chunk")[0] for f in chunk_files))
for co in cohorts:
    files = sorted(glob.glob(os.path.join(qc_dir, f"{co}.chunk*.qcmetrics.tsv")))
    dfs = [pd.read_csv(f, sep="\t", index_col=0) for f in files]
    merged = pd.concat(dfs)
    merged = merged[~merged.index.duplicated(keep="first")]
    out = os.path.join(qc_dir, f"{co}_qc_metrics.tsv")
    merged.to_csv(out, sep="\t")
    print(f"{co}: {len(files)} chunks -> {merged.shape[0]} samples x "
          f"{merged.shape[1]} metrics -> {out}")
print("Merge complete. Now run 19_hcp_factors.sh (it will skip Stage 1).")
PYEOF
    exit 0
fi

# =============================================================================
# Array mode: process chunk $IDX of $COHORT
# =============================================================================
COHORT="${COHORT:?ERROR: COHORT env var required in array mode (e.g. cohort1)}"

echo "[$(date)] Picard QC shard: cohort=${COHORT} chunk=${IDX}/${NCHUNKS}"
echo "  CONFIG:      $CONFIG"
echo "  QC_DIR:      $QC_DIR"
echo "  REFFLAT:     $REFFLAT"
echo "  GENE_ANNOT:  $GENE_ANNOT"

# ---- Reference files: index 1 generates, others wait on sentinel ----
SENTINEL="${QC_DIR}/.refs_done"

if [ ! -f "$REFFLAT" ] || [ ! -f "$GENE_ANNOT" ]; then
    if [ "$IDX" = "1" ]; then
        echo "[$(date)] Generating reference files (index 1)"
        # samtools env + MAJIQ venv provides python3 with pandas/numpy/pysam
        # (same stack as scripts 10-15 and 17; the picard-2.27.4 conda env
        # has no pandas and must not be used for python)
        conda activate samtools-1.16.1
        source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate
        if [ ! -f "$REFFLAT" ]; then
            python3 "${SCRIPTS_DIR}/picard_qc.py" --generate-refflat \
                --gtf "$NORMALIZED_GTF" --output "$REFFLAT"
        fi
        if [ ! -f "$GENE_ANNOT" ]; then
            python3 "${SCRIPTS_DIR}/picard_qc.py" --generate-gene-annot \
                --gtf "$NORMALIZED_GTF" --fasta "$FASTA" --output "$GENE_ANNOT"
        fi
        conda deactivate 2>/dev/null || true
        touch "$SENTINEL"
        echo "[$(date)] References ready"
    else
        echo "[$(date)] Waiting for index 1 to generate references..."
        for i in $(seq 1 180); do
            [ -f "$SENTINEL" ] && break
            sleep 60
        done
        if [ ! -f "$SENTINEL" ]; then
            echo "ERROR: references not ready after 3h wait" >&2
            exit 1
        fi
    fi
fi

# ---- Extract this chunk's samples (deterministic round-robin) ----
SAMPLES_FILE="${OUTPUT_BASE}/${COHORT}/samples.txt"
CHUNK_FILE="${QC_DIR}/${COHORT}.chunk${IDX}.samples.txt"
awk -v i="$IDX" -v n="$NCHUNKS" '(NR-1) % n == i-1' "$SAMPLES_FILE" > "$CHUNK_FILE"

if [ ! -s "$CHUNK_FILE" ]; then
    echo "[$(date)] No samples in chunk ${IDX} (cohort has fewer than ${NCHUNKS} chunks) — exiting cleanly"
    exit 0
fi
echo "[$(date)] Chunk ${IDX}: $(wc -l < "$CHUNK_FILE") samples"

# ---- Run Picard QC on this chunk ----
# 'module load picard' only prints a pointer to the picard-2.27.4 conda env;
# that env holds the picard binary + java but a pandas-less python3. Stack the
# envs: picard env provides picard/java on PATH; the samtools env + MAJIQ venv
# (prepended last, so first on PATH) provides python3 with pandas/numpy.
module load picard
conda activate picard-2.27.4
conda activate --stack samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

python3 "${SCRIPTS_DIR}/picard_qc.py" \
    --config "$CONFIG" \
    --cohort "$COHORT" \
    --samples "$CHUNK_FILE" \
    --output "${QC_DIR}/${COHORT}.chunk${IDX}.qcmetrics.tsv" \
    --refflat "$REFFLAT" \
    --picard-cmd picard \
    --salmon-dir "${OUTPUT_BASE}/${COHORT}/intermediate/expression" \
    --gtf "$NORMALIZED_GTF" \
    --fasta "$FASTA" \
    --gene-annot "$GENE_ANNOT"

echo "[$(date)] Done: cohort=${COHORT} chunk=${IDX}"
echo "Output: ${QC_DIR}/${COHORT}.chunk${IDX}.qcmetrics.tsv"
