#!/bin/bash
# =============================================================================
# 19a_picard_sharded.sh — Sample-parallel Stage 1 (PicardTools QC) for 19_hcp_factors.sh
# =============================================================================
# Sample-parallel version (2026-10-04 rework). The previous chunk-serial
# version split each cohort's samples.txt across a 16-way array and ran
# Picard SERIALLY within each chunk (~30-45 samples x 5 single-threaded
# Picard tools per 2-core job): chunks needed ~12h+ and a timeout lost the
# whole chunk. Picard tools are single-threaded, so the unit of parallelism
# is now the SAMPLE: one array index = one sample (1 core, ~18 min typical).
# The remaining samples of a cohort finish in ~1-2h of wall time instead of
# another 12+h, and a failed/killed job loses exactly one sample.
#
# BACKWARDS COMPATIBLE with the chunk-serial version: samples already present
# in any ${COHORT}.chunk*.qcmetrics.tsv are treated as done and are never
# resubmitted; their rows are reused VERBATIM at merge time (they carry the
# 19b fixes). New samples produce per-sample outputs in
# ${QC_DIR}/per_sample/${COHORT}/<sample>.qcmetrics.tsv (a subdirectory, so
# the main wrapper's Stage-2 glob '*_qc_metrics.tsv' never sees them).
# Re-running one sample = delete its per-sample file and resubmit.
#
# Modes:
#   SUBMIT (default; run on a login node): scan completed outputs, write a
#     timestamped todo list per cohort, and submit one right-sized LSF array
#     per cohort (indices 1..n_todo). DRYRUN=1 prints bsub commands only.
#   WORKER (LSB_JOBINDEX set; entered via the array submission only): run
#     Picard QC for the single sample at line $LSB_JOBINDEX of TODO_FILE.
#   MERGE=1 (login node): combine legacy chunk rows + new per-sample rows
#     into the per-cohort ${COHORT}_qc_metrics.tsv files the main wrapper
#     expects, then run 19_hcp_factors.sh (it skips Stage 1).
#
# Usage:
#   CONFIG=... SCRIPTS_DIR=... bash 19a_picard_sharded.sh           # submit all cohorts
#   CONFIG=... SCRIPTS_DIR=... COHORTS="cohort3 cohort4" bash 19a_picard_sharded.sh
#   DRYRUN=1 CONFIG=... SCRIPTS_DIR=... bash 19a_picard_sharded.sh  # preview only
#   CONFIG=... SCRIPTS_DIR=... MERGE=1 bash 19a_picard_sharded.sh   # merge after arrays
#
# Required env vars: CONFIG
# Optional env vars: SCRIPTS_DIR (default: this script's directory),
#   COHORTS (default: all cohorts in config), QUEUE (long), MEM_GB (12),
#   WALLTIME (2:00), MAXCONC (unset = no array concurrency cap),
#   LOG_DIR (default: PANTRY logs dir), QC_DIR, REFFLAT, GENE_ANNOT, FASTA,
#   PARSE_ONLY (1 = re-extract from persisted raw outputs, no Picard re-run)
#
# Notes:
#   - NCHUNKS is obsolete and ignored. Direct `bsub -J "pic_c1[1-16]" ...`
#     submissions of the OLD form fail fast in worker mode (TODO_FILE
#     required) — use SUBMIT mode instead.
#   - Todo lists are timestamped and passed to the array via TODO_FILE, so a
#     later resubmission never changes the sample mapping of an in-flight
#     array. Still, let an array finish (or bkill it) before resubmitting.
#   - Reference generation (refFlat, gene GC/length) is handled in worker
#     mode: array index 1 generates any missing refs; other indices wait on a
#     sentinel. No-op once the refs exist.
#   - Median imputation of missing metric values moved to MERGE mode
#     (cohort-level medians over the full cohort). Per-sample jobs run
#     picard_qc.py --no-impute: a 1-sample "column median" is degenerate and
#     would fill failed-tool cells with 0.0.
# =============================================================================

# Fallback resources for direct `bsub < 19a_picard_sharded.sh` submission.
# SUBMIT mode passes -q/-n/-M/-R/-W on the bsub command line (which override
# these), so edit the QUEUE/MEM_GB/WALLTIME env vars, not these lines.
#BSUB -q long
#BSUB -n 1
#BSUB -M 12
#BSUB -R "rusage[mem=12]"
#BSUB -W 2:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/picard.%J.%I.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/picard.%J.%I.err

set -eo pipefail

# ---- Required env vars ----
CONFIG="${CONFIG:?ERROR: CONFIG env var required (path to config.yml)}"
# Default SCRIPTS_DIR to this script's own directory, so the pipeline runs
# directly from the git clone. An explicit SCRIPTS_DIR env var overrides —
# and is REQUIRED when submitting via `bsub < script` (LSF executes a spool
# copy of the script; self-location would resolve to the spool directory).
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# config_get.py lives in ../03_phenotyping in the repo layout; a flat copy in
# SCRIPTS_DIR (flat deployment) takes precedence.
CONFIG_GET="${SCRIPTS_DIR}/config_get.py"
[ -f "$CONFIG_GET" ] || CONFIG_GET="${SCRIPTS_DIR}/../03_phenotyping/config_get.py"
if [ ! -f "$CONFIG_GET" ]; then
    echo "ERROR: config_get.py not found in $SCRIPTS_DIR or $SCRIPTS_DIR/../03_phenotyping" >&2
    echo "  Submitting via 'bsub <'? Export SCRIPTS_DIR=<repo>/05_qtl_mapping first." >&2
    exit 1
fi
PERSAMPLE="${SCRIPTS_DIR}/picard_persample.py"
if [ ! -f "$PERSAMPLE" ]; then
    echo "ERROR: picard_persample.py not found in $SCRIPTS_DIR" >&2
    echo "  It ships with the sample-parallel 19a rework — update your checkout." >&2
    exit 1
fi
# PARSE_ONLY=1: skip Picard execution and re-extract metrics from the
# persisted raw outputs under ${QC_DIR}/raw/${COHORT}/<sample>/ (written by
# earlier runs). Used to widen the metric panel without re-running Picard.
PARSE_ONLY="${PARSE_ONLY:-0}"

# ---- Global init (same as 19_hcp_factors.sh) ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load config values as shell vars
eval "$(python3 "$CONFIG_GET" "${CONFIG}")"

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
# MERGE mode: combine legacy chunk outputs + new per-sample outputs into the
# per-cohort files the main wrapper expects, then exit. Run on a login node
# after all arrays finish.
# =============================================================================
if [ "${MERGE:-0}" = "1" ]; then
    COHORTS_ARG=""
    if [ -n "${COHORTS:-}" ]; then
        COHORTS_ARG="--cohorts ${COHORTS}"
    fi
    # shellcheck disable=SC2086
    python3 "$PERSAMPLE" merge --qc-dir "$QC_DIR" --output-base "$OUTPUT_BASE" $COHORTS_ARG
    exit 0
fi

# =============================================================================
# WORKER mode: process the single sample at line $LSB_JOBINDEX of TODO_FILE
# =============================================================================
if [ -n "${LSB_JOBINDEX:-}" ]; then
    COHORT="${COHORT:?ERROR: COHORT env var required in worker mode}"
    TODO_FILE="${TODO_FILE:?ERROR: TODO_FILE env var required in worker mode.
  Direct 'bsub -J \"pic_c1[1-16]\" ...' submissions of the OLD chunk-serial
  form no longer work. Use SUBMIT mode instead:
    CONFIG=... SCRIPTS_DIR=... bash 19a_picard_sharded.sh}"
    IDX="$LSB_JOBINDEX"

    echo "[$(date)] Picard QC per-sample worker: cohort=${COHORT} index=${IDX}"
    echo "  CONFIG:      $CONFIG"
    echo "  TODO_FILE:   $TODO_FILE"
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

    # ---- This worker's sample (line $IDX of the todo file) ----
    SAMPLE="$(sed -n "${IDX}p" "$TODO_FILE" | tr -d '[:space:]')"
    if [ -z "$SAMPLE" ]; then
        echo "[$(date)] No sample at line ${IDX} of ${TODO_FILE} — exiting cleanly"
        exit 0
    fi

    PS_DIR="${QC_DIR}/per_sample/${COHORT}"
    OUT="${PS_DIR}/${SAMPLE}.qcmetrics.tsv"
    mkdir -p "$PS_DIR"

    # Skip if already done (e.g. overlapping resubmission). A valid output
    # has a header plus at least one data row.
    if [ -f "$OUT" ] && [ "$(wc -l < "$OUT")" -ge 2 ]; then
        echo "[$(date)] ${SAMPLE}: per-sample output already exists, skipping: $OUT"
        exit 0
    fi

    echo "[$(date)] Processing sample: ${SAMPLE}"

    # ---- Environment: picard binary + java from picard env, python3 from
    # samtools env + MAJIQ venv (stacked last, so first on PATH) ----
    module load picard
    conda activate picard-2.27.4
    conda activate --stack samtools-1.16.1
    source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

    # ---- Verify picard is reachable (same check as 19_hcp_factors.sh) ----
    if [ "$PARSE_ONLY" != "1" ] && ! command -v picard >/dev/null 2>&1; then
        echo "ERROR: 'picard' not on PATH after env stack (module load picard +"
        echo "  conda activate picard-2.27.4 + samtools + MAJIQ venv). State:"
        echo "  CONDA_PREFIX=${CONDA_PREFIX:-unset}"
        echo "  picard binary:  $(command -v picard || echo none)"
        echo "  java binary:    $(command -v java || echo none)"
        echo "  env bin/ contents: $(ls "${CONDA_PREFIX:-/nonexistent}"/bin 2>/dev/null | grep -i -m3 picard || echo 'no picard* in $CONDA_PREFIX/bin')"
        echo "Fix the activation (env name/path) — do not shim around it."
        exit 1
    fi

    PARSE_FLAG=""
    if [ "$PARSE_ONLY" = "1" ]; then
        PARSE_FLAG="--parse-only"
    fi

    TMP_SAMPLES="$(mktemp /tmp/picard_sample.XXXXXX)"
    echo "$SAMPLE" > "$TMP_SAMPLES"
    TMP_OUT="${OUT}.tmp.${LSB_JOBID:-$$}"

    # --no-impute: a 1-sample run has a degenerate column median that would
    # fill failed-tool cells with 0.0. Imputation happens at merge time with
    # cohort-level medians (picard_persample.py merge).
    if python3 "${SCRIPTS_DIR}/picard_qc.py" \
        --config "$CONFIG" \
        --cohort "$COHORT" \
        --samples "$TMP_SAMPLES" \
        --output "$TMP_OUT" \
        --refflat "$REFFLAT" \
        --picard-cmd picard \
        --salmon-dir "${OUTPUT_BASE}/${COHORT}/intermediate/expression" \
        --gtf "$NORMALIZED_GTF" \
        --fasta "$FASTA" \
        --gene-annot "$GENE_ANNOT" \
        --raw-dir "${QC_DIR}/raw/${COHORT}" \
        --no-impute \
        $PARSE_FLAG; then
        if [ -s "$TMP_OUT" ] && [ "$(wc -l < "$TMP_OUT")" -ge 2 ]; then
            mv "$TMP_OUT" "$OUT"
        else
            echo "ERROR: picard_qc.py produced no data rows for ${SAMPLE}" >&2
            rm -f "$TMP_OUT" "$TMP_SAMPLES"
            exit 1
        fi
    else
        echo "ERROR: picard_qc.py failed for ${SAMPLE} (see above)" >&2
        rm -f "$TMP_OUT" "$TMP_SAMPLES"
        exit 1
    fi
    rm -f "$TMP_SAMPLES"

    echo "[$(date)] Done: cohort=${COHORT} sample=${SAMPLE}"
    echo "Output: ${OUT}"
    exit 0
fi

# =============================================================================
# SUBMIT mode (default; run on a login node): scan completed outputs, write a
# timestamped todo list per cohort, submit one right-sized array per cohort.
# =============================================================================
if [ -n "${LSB_JOBID:-}" ]; then
    echo "ERROR: running inside an LSF job but LSB_JOBINDEX is unset." >&2
    echo "  19a arrays are submitted by SUBMIT mode; run it from a login node:" >&2
    echo "    CONFIG=... SCRIPTS_DIR=... bash 19a_picard_sharded.sh" >&2
    exit 1
fi
QUEUE="${QUEUE:-long}"
MEM_GB="${MEM_GB:-12}"
WALLTIME="${WALLTIME:-2:00}"
MAXCONC="${MAXCONC:-}"
DRYRUN="${DRYRUN:-0}"
LOG_DIR="${LOG_DIR:-/rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs}"
mkdir -p "$LOG_DIR"

# ---- Determine cohorts (same parsing as 19_hcp_factors.sh) ----
if [ -n "${COHORTS:-}" ]; then
    COHORT_LIST="$COHORTS"
else
    COHORT_LIST=$(python3 -c "
import re
with open('${CONFIG}') as f:
    in_cohorts = False
    for line in f:
        stripped = line.split('#')[0].rstrip()
        if not stripped.strip():
            continue
        if stripped.startswith('cohorts:'):
            in_cohorts = True
            continue
        if in_cohorts:
            if re.match(r'^  \w', stripped) and stripped.strip().endswith(':'):
                print(stripped.strip().rstrip(':'))
            elif not stripped.startswith(' '):
                break
")
fi

SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"

echo "[$(date)] Picard QC sample-parallel submission"
echo "  CONFIG:      $CONFIG"
echo "  SCRIPTS_DIR: $SCRIPTS_DIR"
echo "  QC_DIR:      $QC_DIR"
echo "  Cohorts:     $COHORT_LIST"
echo "  Resources:   queue=${QUEUE} cpus=1 mem=${MEM_GB}GB walltime=${WALLTIME} maxconc=${MAXCONC:-none}"
[ "$DRYRUN" = "1" ] && echo "  DRYRUN=1 — no jobs will be submitted"

for COHORT in $COHORT_LIST; do
    SAMPLES_FILE="${OUTPUT_BASE}/${COHORT}/samples.txt"
    if [ ! -f "$SAMPLES_FILE" ]; then
        echo "WARN: no samples.txt for ${COHORT} at ${SAMPLES_FILE} — skipping" >&2
        continue
    fi

    TS="$(date +%Y%m%d%H%M%S)"
    TODO_FILE="${QC_DIR}/${COHORT}.todo.${TS}.txt"
    N_TODO="$(python3 "$PERSAMPLE" scan \
        --samples "$SAMPLES_FILE" \
        --qc-dir "$QC_DIR" \
        --cohort "$COHORT" \
        --todo-out "$TODO_FILE")"

    if [ "$N_TODO" = "0" ]; then
        echo "${COHORT}: all samples already complete — nothing to submit"
        rm -f "$TODO_FILE"
        continue
    fi

    JOBNAME="pic_${COHORT/cohort/c}"
    CONC=""
    [ -n "$MAXCONC" ] && CONC="%${MAXCONC}"

    BSUB_CMD=(bsub
        -J "${JOBNAME}[1-${N_TODO}]${CONC}"
        -q "$QUEUE"
        -n 1
        -M "$MEM_GB"
        -R "rusage[mem=${MEM_GB}]"
        -W "$WALLTIME"
        -o "${LOG_DIR}/picard.%J.%I.out"
        -e "${LOG_DIR}/picard.%J.%I.err"
        -env "CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},COHORT=${COHORT},TODO_FILE=${TODO_FILE},QC_DIR=${QC_DIR},REFFLAT=${REFFLAT},GENE_ANNOT=${GENE_ANNOT},FASTA=${FASTA},PARSE_ONLY=${PARSE_ONLY}")

    echo "${COHORT}: ${N_TODO} samples to run -> array ${JOBNAME}[1-${N_TODO}]${CONC}"
    if [ "$DRYRUN" = "1" ]; then
        echo "  DRYRUN: ${BSUB_CMD[*]} < ${SCRIPT_PATH}"
    else
        "${BSUB_CMD[@]}" < "$SCRIPT_PATH"
    fi
done

echo "[$(date)] Submission pass complete."
echo "After all arrays finish, merge with:"
echo "  CONFIG=${CONFIG} SCRIPTS_DIR=${SCRIPTS_DIR} MERGE=1 bash ${SCRIPT_PATH}"
