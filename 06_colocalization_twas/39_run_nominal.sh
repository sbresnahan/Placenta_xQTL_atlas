#!/bin/bash
# =============================================================================
# 39_run_nominal.sh — submit genome-wide nominal cis-xQTL scans (module 06)
# =============================================================================
# Objective 1.6: full cis-window summary statistics per ancestry x modality
# for colocalization (SuSiE-coloc / colocBoost) and TWAS QC. One LSF array
# job per ancestry x modality (22 chromosome tasks) plus a dependent merge
# job that writes the bgzipped/tabix'd per-modality nominal store.
#
# Prerequisites: module-05 steps 23-25 completed (qtl_inputs has
#   {ANC}_qtl.pgen, {ANC}_{MOD}.bed.gz, {ANC}_covariates_{MOD}.tsv);
#   tensorqtl conda env (module-05 script 21).
#
# Usage:
#   TEST=1 bash 39_run_nominal.sh    # one pilot chromosome (chr21, expression)
#   bash 39_run_nominal.sh           # all ancestries x modalities
#
# Optional overrides: ANCESTRIES, MODALITIES, CHROMS, QUEUE, WALLTIME,
#   MERGE_WALLTIME, THREADS, MEM, GPU=1 (submit to GPU queue; tensorQTL
#   auto-uses CUDA), FORCE=1, MAF_THRESHOLD, SKIP_COLLAPSE_CHECK=1,
#   MERGE_ONLY=1 (shards already done: submit only the merge jobs, no
#   dependency — the recovery path after a merge-step failure).
#
# Worker environment (set inside the bsub heredocs): conda env 'tensorqtl'
# (module-05 script 21) for 39_run_nominal.py; merge jobs additionally
# 'module load samtools' for bgzip/tabix.
# =============================================================================
set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
MODALITIES="${MODALITIES:-expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability}"
CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-12:00}"
# seadragon queue runlimits: short <3h, medium >3h and <24h, long 24-120h
# (esub rejects out-of-range requests) — the merge default must exceed 3h
MERGE_WALLTIME="${MERGE_WALLTIME:-04:00}"
THREADS="${THREADS:-4}"
MEM="${MEM:-32G}"
GPU="${GPU:-0}"
GPU_QUEUE="${GPU_QUEUE:-gpu}"
MAF_THRESHOLD="${MAF_THRESHOLD:-0.01}"
TEST="${TEST:-0}"
FORCE="${FORCE:-0}"

if [ "$TEST" = "1" ]; then
    ANCESTRIES="EAS"; MODALITIES="expression"; CHROMS="21"
    echo "TEST=1: pilot = EAS expression chr21 only"
fi

# --- Collapsed-replicate input audit ------------------------------------------
# The nominal scans consume sample-level qtl_inputs (BEDs, covariates, pgen).
# Verify those sample sets match the module-05 collapsed-replicate contract
# (the same ancestry_map_collapsed.tsv authority module 07 enforces) before
# submitting. Stdlib-only script: no conda env required on the login node.
if [ "${SKIP_COLLAPSE_CHECK:-0}" != "1" ]; then
    python3 "${SCRIPTS_DIR}/check_collapsed_inputs.py" \
        --qtl-dir "$QTL_DIR" \
        --output-base "$OUTPUT_BASE" \
        --ancestries $ANCESTRIES
fi

mkdir -p "$LOG_DIR"
N_CHROMS=$(echo $CHROMS | wc -w)

# Submit the merge job for the current ANC/MOD; extra bsub args (e.g. the
# -w dependency on the array) are passed through.
submit_merge() {
    bsub -J "$MERGE_JOB" -q "$QUEUE" -n 2 -W "$MERGE_WALLTIME" -M 64G -R "rusage[mem=64G]" \
         "$@" \
         -o "${LOG_DIR}/${MERGE_JOB}.%J.out" -e "${LOG_DIR}/${MERGE_JOB}.%J.err" \
         -env "$ENV_STR" \
         <<'EOF'
#!/bin/bash
set -eo pipefail
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl
module load samtools
python3 "${SCRIPTS_DIR}/39_run_nominal.py" \
    --qtl-dir "$QTL_DIR" --output-dir "$RESULTS_DIR" \
    --ancestry "$ANC" --modality "$MOD" --merge
EOF
}

for ANC in $ANCESTRIES; do
  for MOD in $MODALITIES; do
    JOB="nom_${ANC}_${MOD}"
    MERGE_JOB="nomm_${ANC}_${MOD}"
    SHARD_GLOB="${RESULTS_DIR}/nominal/${ANC}/${ANC}_${MOD}.nominal.chr*.parquet"
    MERGED="${RESULTS_DIR}/nominal/${ANC}/${ANC}_${MOD}.nominal.tsv.gz"
    if [ -f "$MERGED" ] && [ "$FORCE" != "1" ]; then
        echo "  SKIP ${ANC}/${MOD}: merged nominal store exists"
        continue
    fi
    if bjobs -J "$JOB" 2>/dev/null | grep -q "$JOB"; then
        echo "  SKIP ${ANC}/${MOD}: array already running/pending"
        continue
    fi
    # ENV_STR must be set before any submit_merge call (MERGE_ONLY branch
    # below) — an empty -env value is rejected by LSF
    ENV_STR="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},QTL_DIR=${QTL_DIR},RESULTS_DIR=${RESULTS_DIR},ANC=${ANC},MOD=${MOD},MAF_THRESHOLD=${MAF_THRESHOLD},FORCE=${FORCE}"
    if [ "${MERGE_ONLY:-0}" = "1" ]; then
        # recovery mode (e.g. after a merge-step failure): shards are done,
        # submit only the merge job with no dependency
        n_shards=$(ls ${SHARD_GLOB} 2>/dev/null | wc -l)
        if [ "$n_shards" -eq 0 ]; then
            echo "  SKIP ${ANC}/${MOD}: no shards found (run the array first)"
            continue
        fi
        submit_merge
        echo "  submitted merge-only ${MERGE_JOB} (${n_shards} shards)"
        continue
    fi
    EXTRA=""
    if [ "$GPU" = "1" ]; then
        QUEUE_USE="$GPU_QUEUE"
        EXTRA='-gpu "num=1"'
    else
        QUEUE_USE="$QUEUE"
    fi
    # array of chromosome shards
    bsub -J "${JOB}[1-${N_CHROMS}]" -q "$QUEUE_USE" -n "$THREADS" -W "$WALLTIME" \
         -M "$MEM" -R "rusage[mem=${MEM}]" $EXTRA \
         -o "${LOG_DIR}/${JOB}.%J.%I.out" -e "${LOG_DIR}/${JOB}.%J.%I.err" \
         -env "$ENV_STR,CHROMS_LIST=${CHROMS}" \
         <<'EOF'
#!/bin/bash
set -eo pipefail
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl
CHROM=$(echo $CHROMS_LIST | cut -d' ' -f${LSB_JOBINDEX})
FORCE_FLAG=""; [ "$FORCE" = "1" ] && FORCE_FLAG="--force"
python3 "${SCRIPTS_DIR}/39_run_nominal.py" \
    --qtl-dir "$QTL_DIR" --output-dir "$RESULTS_DIR" \
    --ancestry "$ANC" --modality "$MOD" --chrom "$CHROM" \
    --maf-threshold "$MAF_THRESHOLD" $FORCE_FLAG
EOF
    # dependent merge job — skipped in TEST mode: merging a single pilot
    # chromosome would write a partial nominal store that later full runs
    # would mistake for complete (the merged-file check above)
    if [ "$TEST" != "1" ]; then
        submit_merge -w "done(${JOB})"
        echo "  submitted ${JOB}[1-${N_CHROMS}] + merge ${MERGE_JOB}"
    else
        echo "  submitted ${JOB}[1-${N_CHROMS}] (TEST mode: no merge job)"
    fi
    if [ "$TEST" = "1" ]; then
        echo "TEST=1: one pilot array submitted; check ${LOG_DIR}/${JOB}.*.out"
        exit 0
    fi
  done
done
