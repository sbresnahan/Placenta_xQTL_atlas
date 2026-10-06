#!/bin/bash
# =============================================================================
# 28b_submit_independent.sh — Submit CHUNKED tensorQTL cis.map_independent
#                             (stepwise) scans + merge, per ancestry × modality
# =============================================================================
# Splits the conditionally-independent (forward-backward stepwise) scan — the
# >24h bottleneck of the monolithic --independent run — into an LSF job array
# of ~CHUNK_SIZE-gene chunks (27b_independent_chunk.sh), followed by a merge
# job (27c_merge_independent.py) that rebuilds the same final outputs the
# monolithic mode produced:
#   qtl_results/{ANC}_{MOD}[_ungrouped]_cisqtl_independent.parquet / _top.tsv
#
# Chunking is exact: each chunk passes the FULL map_cis results table to
# cis.map_independent (so its internal significance threshold is unchanged)
# and subsets only the phenotype matrix. Merged output is statistically
# identical to the monolithic scan (bit-identical with a fixed SEED).
#
# If the map_cis results ({ANC}_{MOD}[_ungrouped]_cisqtl.parquet) do not exist
# yet, a map_cis-only job is submitted first (or an already-running eqtl_* job
# is reused) and this script re-invokes itself as a chained job once the
# mapping finishes — so one command handles the whole flow. If the map_cis job
# FAILS, the chained planner never runs (done() dependency); fix the failure
# and rerun this script.
#
# Usage:
#   MAF_THRESHOLD=0.01 ANCESTRIES=EAS MODALITIES="combined splicing" \
#       bash 28b_submit_independent.sh
#   MAF_THRESHOLD=0.01 GPU=1 GPU_QUEUE=gpu bash 28b_submit_independent.sh
#   CHUNK_INDICES="3,17" MAF_THRESHOLD=0.01 MODALITIES=combined \
#       bash 28b_submit_independent.sh    # retry failed chunks, then re-merge
#   TEST=1 MAF_THRESHOLD=0.01 MODALITIES=splicing bash 28b_submit_independent.sh
#
# Combos with existing final results or a running/pending job of the same name
# are skipped, so it is safe to rerun this script at any time.
#
# Optional env:
#   ANCESTRIES    — default "EAS EUR"
#   MODALITIES    — default: the 8 non-expression modalities (set explicitly,
#                   e.g. MODALITIES="expression" or "combined")
#   GROUPED       — 1 (default) = grouped layer; 0 = ungrouped (_ungrouped)
#   CHUNK_SIZE    — genes/groups (grouped) or phenotypes (ungrouped) per chunk
#                   (default: 100; consider 500 with GPU=1, 50 for combined if
#                   chunks approach the walltime)
#   INDEPENDENT_FDR — FDR entry threshold (default: 0.05; must match between
#                   array and merge — this script handles that)
#   QUEUE         — chunk queue (default: medium)
#   WALLTIME      — chunk walltime (default: 4:00)
#   NCPU / MEM    — chunk resources (default: 8 cores / 32 GB; the floor is
#                   set by loading the full pgen, same as the map_cis job)
#   GPU           — 1 = submit chunks to a GPU queue (tensorQTL auto-uses CUDA;
#                   typically 10-50x faster per chunk). Default: 0
#   GPU_QUEUE     — default gpu (seadragon also has egpu)
#   GPU_OPTS      — value passed to bsub's -gpu flag (default: num=1; e.g.
#                   "num=1:mode=shared". Space-free so it survives bsub -env
#                   passthrough when chaining)
#   GPU_WALLTIME  — default 2:00
#   MERGE_QUEUE / MERGE_WALLTIME — default short / 1:00
#   MAP_QUEUE / MAP_WALLTIME — for the map_cis job IF one must be submitted
#                   (default: medium / 12:00; use long / 48:00 for combined)
#   CIS_WINDOW    — default 1000000
#   MAF_THRESHOLD — REQUIRED (no default), e.g. 0.01
#   COVARIATES_FILE — covariates TSV with {ANC}/{MOD} placeholders (same
#                   convention as 28_submit_modalities.sh)
#   SEED          — permutation RNG seed for reproducible chunks (default:
#                   unset = unseeded, matching previous behavior)
#   QVALUE_METHOD — only used if a map_cis job is submitted (default: storey)
#   TEST          — 1 = submit a single-chunk pilot array (+merge) for the
#                   first combo only
#
# NOTE: the -env strings deliberately contain NO inner quotes — LSF preserves
# them literally in the variable values, which mangles paths.
# =============================================================================

set -eo pipefail

CONFIG="${CONFIG:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/config.yml}"
# Default SCRIPTS_DIR to this script's own directory (repo-clone layout).
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
OUTPUT_BASE="${OUTPUT_BASE:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY}"
LOG_DIR="${LOG_DIR:-/rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs}"

ANCESTRIES="${ANCESTRIES:-EAS EUR}"
MODALITIES="${MODALITIES:-alt_polyA alt_TSS intron_retention isoforms isoform_expression RNA_editing splicing stability}"
GROUPED="${GROUPED:-1}"
CHUNK_SIZE="${CHUNK_SIZE:-100}"
INDEPENDENT_FDR="${INDEPENDENT_FDR:-0.05}"
CIS_WINDOW="${CIS_WINDOW:-1000000}"
MAF_THRESHOLD="${MAF_THRESHOLD:?ERROR: set MAF_THRESHOLD explicitly (e.g. MAF_THRESHOLD=0.01) — it is never defaulted, to avoid silently running the wrong threshold}"
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-4:00}"
NCPU="${NCPU:-8}"
MEM="${MEM:-32}"
GPU="${GPU:-0}"
GPU_QUEUE="${GPU_QUEUE:-gpu}"
GPU_OPTS="${GPU_OPTS:-num=1}"
GPU_WALLTIME="${GPU_WALLTIME:-2:00}"
MERGE_QUEUE="${MERGE_QUEUE:-short}"
MERGE_WALLTIME="${MERGE_WALLTIME:-1:00}"
MAP_QUEUE="${MAP_QUEUE:-medium}"
MAP_WALLTIME="${MAP_WALLTIME:-12:00}"
COVARIATES_FILE="${COVARIATES_FILE:-}"
SEED="${SEED:-}"
QVALUE_METHOD="${QVALUE_METHOD:-storey}"
TEST="${TEST:-0}"
CHUNK_INDICES="${CHUNK_INDICES:-}"

# Strip any literal quotes that LSF's -env may have preserved in the values
# (this script re-invokes itself via bsub -env when chaining).
for v in CONFIG SCRIPTS_DIR OUTPUT_BASE ANCESTRIES MODALITIES GROUPED \
         CHUNK_SIZE INDEPENDENT_FDR CIS_WINDOW MAF_THRESHOLD QUEUE WALLTIME \
         NCPU MEM GPU GPU_QUEUE GPU_OPTS GPU_WALLTIME MERGE_QUEUE \
         MERGE_WALLTIME MAP_QUEUE MAP_WALLTIME COVARIATES_FILE SEED \
         QVALUE_METHOD TEST CHUNK_INDICES; do
    eval "$v=\"\${$v%\\\"}\"; $v=\"\${$v#\\\"}\""
done

mkdir -p "$LOG_DIR"

# ---- Sanity checks ----
for f in 27b_independent_chunk.sh 27c_merge_independent.py 27_run_tensorqtl.py 27_run_tensorqtl.sh; do
    if [ ! -f "${SCRIPTS_DIR}/$f" ]; then
        echo "ERROR: ${SCRIPTS_DIR}/$f not found"
        exit 1
    fi
done

# The row-counting snippet needs pandas/pyarrow — use the tensorqtl conda env
# (works on submit hosts and inside the chained planner job alike).
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl
python3 -c "import pandas, pyarrow" 2>/dev/null || {
    echo "ERROR: pandas/pyarrow not importable after 'conda activate tensorqtl'."
    echo "  (run 21_install_tensorqtl.sh first)"
    exit 1
}

# Mode-dependent labels (same convention as 28_submit_modalities.sh)
if [ "$GROUPED" != "1" ]; then
    NO_GROUPS=1
    BASE_SUFFIX="_ungrouped"
    MODE_DESC="UNGROUPED"
else
    NO_GROUPS=0
    BASE_SUFFIX=""
    MODE_DESC="grouped"
fi

echo "Submitting CHUNKED independent (stepwise) xQTL scans:"
echo "  CONFIG:      $CONFIG"
echo "  SCRIPTS_DIR: $SCRIPTS_DIR"
echo "  OUTPUT_BASE: $OUTPUT_BASE"
echo "  Ancestries:  $ANCESTRIES"
echo "  Modalities:  $MODALITIES"
echo "  Mode:        $MODE_DESC, chunk size $CHUNK_SIZE, FDR $INDEPENDENT_FDR"
echo "  MAF:         $MAF_THRESHOLD"
if [ "$GPU" = "1" ]; then
    echo "  Chunks:      GPU queue $GPU_QUEUE ($GPU_OPTS), walltime $GPU_WALLTIME"
else
    echo "  Chunks:      CPU queue $QUEUE, $NCPU cores / ${MEM} GB, walltime $WALLTIME"
fi
[ -n "$SEED" ] && echo "  SEED:        $SEED"
[ -n "$CHUNK_INDICES" ] && echo "  Rerun mode:  chunks [$CHUNK_INDICES] only"
echo ""

N=0
for ANC in $ANCESTRIES; do
    for MOD in $MODALITIES; do
        LABEL="${MOD}${BASE_SUFFIX}"
        CIS_PARQUET="${OUTPUT_BASE}/qtl_results/${ANC}_${LABEL}_cisqtl.parquet"
        FINAL_PARQUET="${OUTPUT_BASE}/qtl_results/${ANC}_${LABEL}_cisqtl_independent.parquet"
        ARRAY_NAME="inqtl_${ANC}_${LABEL}"
        MERGE_NAME="inqtlm_${ANC}_${LABEL}"

        # ---- Skip logic (same convention as 28_submit_modalities.sh) ----
        if [ -f "$FINAL_PARQUET" ]; then
            echo "  SKIP ${ANC}/${LABEL}: independent results already exist"
            continue
        fi
        if bjobs -J "$ARRAY_NAME" 2>/dev/null | grep -q "$ARRAY_NAME"; then
            echo "  SKIP ${ANC}/${LABEL}: chunk array already running/pending"
            continue
        fi
        if bjobs -J "$MERGE_NAME" 2>/dev/null | grep -q "$MERGE_NAME"; then
            echo "  SKIP ${ANC}/${LABEL}: merge job already running/pending"
            continue
        fi

        # ---- Ensure map_cis results exist; chain via self-reinvocation ----
        if [ ! -f "$CIS_PARQUET" ]; then
            OLD_ID=$(bjobs -J "eqtl_${ANC}_${LABEL}_ind" -noheader -o jobid 2>/dev/null | head -1 || true)
            if [ -n "$OLD_ID" ]; then
                echo "  WARNING ${ANC}/${LABEL}: old-style MONOLITHIC job eqtl_${ANC}_${LABEL}_ind ($OLD_ID) is running."
                echo "    It runs the slow unchunked stepwise scan. Its map_cis parquet is written ~1h in,"
                echo "    so:  bkill $OLD_ID   then rerun this script — the chunked scan will reuse the parquet."
                continue
            fi
            MAP_JOBID=$(bjobs -J "eqtl_${ANC}_${LABEL}" -noheader -o jobid 2>/dev/null | head -1 || true)
            if [ -z "$MAP_JOBID" ]; then
                # Submit the map_cis-only job ourselves
                if [ ! -f "${OUTPUT_BASE}/qtl_inputs/${ANC}_${MOD}_harmonized.bed" ] && \
                   [ ! -f "${OUTPUT_BASE}/qtl_inputs/${ANC}_${MOD}.bed.gz" ]; then
                    echo "  SKIP ${ANC}/${LABEL}: harmonized BED not found (run 26 first)"
                    continue
                fi
                MAP_ENV="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},ANCESTRIES=${ANC},MODALITY=${MOD},NO_GROUPS=${NO_GROUPS},INDEPENDENT=0,MAF_THRESHOLD=${MAF_THRESHOLD},CIS_WINDOW=${CIS_WINDOW},QVALUE_METHOD=${QVALUE_METHOD}"
                if [ -n "$COVARIATES_FILE" ]; then
                    MAP_ENV="${MAP_ENV},COVARIATES_FILE=${COVARIATES_FILE}"
                fi
                MAP_JOBID=$(bsub -J "eqtl_${ANC}_${LABEL}" -q "$MAP_QUEUE" -n 8 -M 32 -R "rusage[mem=32]" -W "$MAP_WALLTIME" \
                    -o "${LOG_DIR}/eqtl_${ANC}_${LABEL}.%J.out" \
                    -e "${LOG_DIR}/eqtl_${ANC}_${LABEL}.%J.err" \
                    -env "$MAP_ENV" \
                    < "${SCRIPTS_DIR}/27_run_tensorqtl.sh" | sed -n 's/.*<\([0-9][0-9]*\)>.*/\1/p')
                if [ -z "$MAP_JOBID" ]; then
                    echo "  ERROR ${ANC}/${LABEL}: failed to submit map_cis job"
                    continue
                fi
                echo "  ${ANC}/${LABEL}: submitted map_cis job $MAP_JOBID (eqtl_${ANC}_${LABEL})"
            else
                echo "  ${ANC}/${LABEL}: map_cis job $MAP_JOBID already running/pending"
            fi
            # Chain: re-invoke this script once map_cis completes successfully.
            PLAN_ENV="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},OUTPUT_BASE=${OUTPUT_BASE},ANCESTRIES=${ANC},MODALITIES=${MOD},GROUPED=${GROUPED},CHUNK_SIZE=${CHUNK_SIZE},INDEPENDENT_FDR=${INDEPENDENT_FDR},CIS_WINDOW=${CIS_WINDOW},MAF_THRESHOLD=${MAF_THRESHOLD},QUEUE=${QUEUE},WALLTIME=${WALLTIME},NCPU=${NCPU},MEM=${MEM},GPU=${GPU},GPU_QUEUE=${GPU_QUEUE},GPU_OPTS=${GPU_OPTS},GPU_WALLTIME=${GPU_WALLTIME},MERGE_QUEUE=${MERGE_QUEUE},MERGE_WALLTIME=${MERGE_WALLTIME},MAP_QUEUE=${MAP_QUEUE},MAP_WALLTIME=${MAP_WALLTIME},QVALUE_METHOD=${QVALUE_METHOD},TEST=${TEST}"
            if [ -n "$COVARIATES_FILE" ]; then
                PLAN_ENV="${PLAN_ENV},COVARIATES_FILE=${COVARIATES_FILE}"
            fi
            if [ -n "$SEED" ]; then
                PLAN_ENV="${PLAN_ENV},SEED=${SEED}"
            fi
            bsub -J "inqtl_plan_${ANC}_${LABEL}" -q short -n 1 -M 4 -R "rusage[mem=4]" -W 0:30 \
                -w "done(${MAP_JOBID})" \
                -o "${LOG_DIR}/inqtl_plan_${ANC}_${LABEL}.%J.out" \
                -e "${LOG_DIR}/inqtl_plan_${ANC}_${LABEL}.%J.err" \
                -env "$PLAN_ENV" \
                < "${SCRIPTS_DIR}/28b_submit_independent.sh" > /dev/null
            echo "    chained: chunk submission will run when job $MAP_JOBID finishes (inqtl_plan_${ANC}_${LABEL})"
            N=$((N + 1))
            if [ "$TEST" = "1" ]; then
                echo ""
                echo "TEST=1: chained one pilot. After the map_cis job finishes, the planner submits a 1-chunk array + merge."
                exit 0
            fi
            continue
        fi

        # ---- Count FDR-significant rows -> array size ----
        N_SIG=$(python3 - "$CIS_PARQUET" "$INDEPENDENT_FDR" <<'PYEOF'
import sys
import pandas as pd
try:
    df = pd.read_parquet(sys.argv[1], columns=["qval"])
except Exception as e:
    sys.exit(f"ERROR: could not read a 'qval' column from {sys.argv[1]}: {e}\n"
             "  Re-run the map_cis step with the current 27_run_tensorqtl.py.")
print(int((pd.to_numeric(df["qval"], errors="coerce") <= float(sys.argv[2])).sum()))
PYEOF
)
        N_CHUNKS=$(( (N_SIG + CHUNK_SIZE - 1) / CHUNK_SIZE ))

        # ---- Merge-job submission (inline script; values via -env) ----
        MERGE_ENV="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},OUTPUT_BASE=${OUTPUT_BASE},ANC=${ANC},MOD=${MOD},NO_GROUPS=${NO_GROUPS},CHUNK_SIZE=${CHUNK_SIZE},INDEPENDENT_FDR=${INDEPENDENT_FDR}"
        submit_merge() {
            local dep="$1"
            local dep_arg=""
            if [ -n "$dep" ]; then
                dep_arg="-w ended($dep)"
            fi
            bsub -J "$MERGE_NAME" -q "$MERGE_QUEUE" -n 1 -M 8 -R "rusage[mem=8]" -W "$MERGE_WALLTIME" \
                $dep_arg \
                -o "${LOG_DIR}/inqtlm_${ANC}_${LABEL}.%J.out" \
                -e "${LOG_DIR}/inqtlm_${ANC}_${LABEL}.%J.err" \
                -env "$MERGE_ENV" <<'MERGE_EOF' > /dev/null
#!/bin/bash
set -eo pipefail
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl
EXTRA=""
if [ "$NO_GROUPS" = "1" ]; then EXTRA="--no-groups"; fi
python3 "${SCRIPTS_DIR}/27c_merge_independent.py" \
    --output-dir "${OUTPUT_BASE}/qtl_results" \
    --ancestry "$ANC" \
    --modality "$MOD" \
    --chunk-size "$CHUNK_SIZE" \
    --independent-fdr "$INDEPENDENT_FDR" $EXTRA
MERGE_EOF
        }

        if [ "$N_SIG" -eq 0 ] && [ -z "$CHUNK_INDICES" ]; then
            submit_merge ""
            echo "  ${ANC}/${LABEL}: 0 significant at FDR <= ${INDEPENDENT_FDR} — submitted merge only (writes empty outputs)"
            N=$((N + 1))
            continue
        fi

        # ---- Submit the chunk array ----
        if [ -n "$CHUNK_INDICES" ]; then
            IDX_SPEC="$CHUNK_INDICES"
        else
            IDX_SPEC="1-${N_CHUNKS}"
        fi
        if [ "$TEST" = "1" ]; then
            IDX_SPEC="1"
        fi

        if [ "$GPU" = "1" ]; then
            Q="$GPU_QUEUE"; W="$GPU_WALLTIME"
            GPU_ARGS=(-gpu "$GPU_OPTS")
        else
            Q="$QUEUE"; W="$WALLTIME"
            GPU_ARGS=()
        fi

        ENV_STR="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},ANCESTRIES=${ANC},MODALITY=${MOD},NO_GROUPS=${NO_GROUPS},CHUNK_SIZE=${CHUNK_SIZE},INDEPENDENT_FDR=${INDEPENDENT_FDR},CIS_WINDOW=${CIS_WINDOW},MAF_THRESHOLD=${MAF_THRESHOLD}"
        if [ -n "$COVARIATES_FILE" ]; then
            ENV_STR="${ENV_STR},COVARIATES_FILE=${COVARIATES_FILE}"
        fi
        if [ -n "$SEED" ]; then
            ENV_STR="${ENV_STR},SEED=${SEED}"
        fi

        ARRAY_JOBID=$(bsub -J "${ARRAY_NAME}[${IDX_SPEC}]" -q "$Q" -n "$NCPU" -M "$MEM" -R "rusage[mem=$MEM]" -W "$W" "${GPU_ARGS[@]}" \
            -o "${LOG_DIR}/inqtl_${ANC}_${LABEL}.%J.%I.out" \
            -e "${LOG_DIR}/inqtl_${ANC}_${LABEL}.%J.%I.err" \
            -env "$ENV_STR" \
            < "${SCRIPTS_DIR}/27b_independent_chunk.sh" | sed -n 's/.*<\([0-9][0-9]*\)>.*/\1/p')
        if [ -z "$ARRAY_JOBID" ]; then
            echo "  ERROR ${ANC}/${LABEL}: failed to submit chunk array"
            continue
        fi

        # Merge depends on the WHOLE array (ended, not done: the merge itself
        # verifies chunk completeness and prints a resubmission command if any
        # chunk failed, instead of pending forever).
        submit_merge "$ARRAY_JOBID"

        echo "  ${ANC}/${LABEL}: ${N_SIG} significant -> array ${ARRAY_NAME}[${IDX_SPEC}] (job $ARRAY_JOBID) + merge $MERGE_NAME"
        N=$((N + 1))

        if [ "$TEST" = "1" ]; then
            echo ""
            echo "TEST=1: submitted a single-chunk pilot array + merge for ${ANC}/${LABEL}."
            echo "Check the chunk log for 'CUDA available' (GPU mode) and 'Written: .../chunk_0001.parquet',"
            echo "and the merge log for the final independent parquets. Then rerun without TEST=1."
            exit 0
        fi
    done
done

echo ""
echo "Submitted/chained $N combos ($MODE_DESC independent layer, chunk size $CHUNK_SIZE)."
echo "Monitor with: bjobs; logs: ${LOG_DIR}/inqtl_*.<jobid>.<index>.out"
echo "If chunks fail, rerun only those: CHUNK_INDICES=\"3,17\" ... bash 28b_submit_independent.sh"
