#!/bin/bash
# =============================================================================
# 03_refprep_post.sh — Reference prep STAGE 3 of 3 (post-array)
# =============================================================================
# Runs after the STAGE 2 txrevise array completes (LSF job-name dependency on
# the array "txrevise"). It:
#   - 5e: merges the per-batch GFF3s into the 6 final txrevise GFF3 files
#   - 6:  extracts transcript sequences (gffread) and builds the 6 Salmon indices
#
# Depends on STAGE 2. `done(txrevise)` is satisfied only when ALL array elements
# have finished.
#
# Config loading: CONFIG and SCRIPTS_DIR env vars must be set via LSF -env at
# submission time:
#   bsub -env "CONFIG=\"config.yml\",SCRIPTS_DIR=\"/path/to/scripts\"" < 03_refprep_post.sh
# =============================================================================

# ---- LSF directives ----
#BSUB -J refprep_post
#BSUB -q long
#BSUB -n 12
#BSUB -M 120
#BSUB -R "rusage[mem=120]"
#BSUB -W 25:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/refprep_post.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/refprep_post.%J.err

set -eo pipefail

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# ---- Config loading (Tier 1: CONFIG and SCRIPTS_DIR via LSF -env) ----
CONFIG="${CONFIG:?Set CONFIG via bsub -env}"
SCRIPTS_DIR="${SCRIPTS_DIR:?Set SCRIPTS_DIR via bsub -env}"

# Load all config values as shell vars
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}")"

# ---- Paths from config ----
REF_GENOME="${REF_GENOME}"
REFERENCE_DIR="${REFERENCE_DIR}"
PANTRY_SCRIPTS="${PANTRY_SCRIPTS}"

# Scratch log dir (infra — hardcoded, used by #BSUB directives above)
LOG_DIR="/rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs"

# Batch count from config (must match 01_refprep_pre.sh and 002_txrevise_array.sh)
TXREVISE_N_BATCHES="${TXREVISE_N_BATCHES}"
TXREVISE_DIR="$REFERENCE_DIR/txrevise"

mkdir -p "$LOG_DIR"
cd "$REFERENCE_DIR"

echo "[$(date)] STAGE 3 (post-array) reference prep starting"

# =============================================================================
# 5e. Merge batch outputs into 6 GFF3 files (grp_1/grp_2 x upstream/contained/downstream)
# =============================================================================
# Section-level check: skip merge entirely if all 6 merged GFF3s already exist.
txrevise_all_exist=true
for group in grp_1 grp_2; do
    for position in upstream contained downstream; do
        [ -s "$TXREVISE_DIR/txrevise.${group}.${position}.gff3" ] || txrevise_all_exist=false
    done
done

if [ "$txrevise_all_exist" = true ]; then
    echo "[$(date)] Section 5e SKIPPED — all 6 txrevise GFF3s already exist"
else
    if ! ls "$TXREVISE_DIR/batch/"*.gff3 >/dev/null 2>&1; then
        echo "ERROR: no batch GFF3s found in $TXREVISE_DIR/batch."
        echo "  The STAGE 2 array (02_txrevise_array.sh) must complete before this runs."
        exit 1
    fi
    echo "[$(date)] Merging txrevise batch outputs"
    for group in grp_1 grp_2; do
        for position in upstream contained downstream; do
            out_file="$TXREVISE_DIR/txrevise.${group}.${position}.gff3"
            if [ -s "$out_file" ]; then
                echo "  [$group.$position] SKIPPED — merged GFF3 already exists: $out_file"
                continue
            fi
            : > "$out_file"
            for batch in $(seq 1 $TXREVISE_N_BATCHES); do
                batch_file="$TXREVISE_DIR/batch/txrevise.${group}.${position}.${batch}_${TXREVISE_N_BATCHES}.gff3"
                if [ -f "$batch_file" ]; then
                    grep -v "^#" "$batch_file" >> "$out_file"
                fi
            done
            echo "  Merged $out_file ($(wc -l < "$out_file") lines)"
        done
    done
    echo "[$(date)] txrevise merge done"
fi

# =============================================================================
# 6. Transcript sequences from txrevise GFF3 (gffread) + Salmon indices
# =============================================================================
# Section-level check: skip entirely if all 6 Salmon indices already exist.
salmon_all_exist=true
for group in grp_1 grp_2; do
    for position in upstream contained downstream; do
        idx_dir="$TXREVISE_DIR/txrevise.${group}.${position}.salmon_index"
        if [ ! -d "$idx_dir" ] || [ -z "$(ls -A "$idx_dir" 2>/dev/null)" ]; then
            salmon_all_exist=false
        fi
    done
done

if [ "$salmon_all_exist" = true ]; then
    echo "[$(date)] Section 6 SKIPPED — all 6 txrevise Salmon indices already exist"
else
    echo "[$(date)] Building txrevise transcript sequences + Salmon indices"
    # gffread runs from its own sourced venv (stays active throughout). bgzip
    # lives in the samtools env and salmon in its own env, and the two conda
    # envs can't be co-active, so extraction (gffread | bgzip) and indexing
    # (salmon) switch conda env per command.
    # gffread activate path (infra — hardcoded)
    source /rsrch5/home/epi/bhattacharya_lab/software/gffread/bin/activate

    for group in grp_1 grp_2; do
        for position in upstream contained downstream; do
            gff3="$TXREVISE_DIR/txrevise.${group}.${position}.gff3"
            fa="$TXREVISE_DIR/txrevise.${group}.${position}.fa.gz"
            idx_dir="$TXREVISE_DIR/txrevise.${group}.${position}.salmon_index"

            if [ -d "$idx_dir" ] && [ -n "$(ls -A "$idx_dir" 2>/dev/null)" ]; then
                echo "  [$group.$position] SKIPPED — Salmon index already exists: $idx_dir"
                continue
            fi

            echo "  [$group.$position] Extracting transcript sequences"
            conda activate samtools-1.16.1   # bgzip
            gffread "$gff3" -g "$REF_GENOME" -w - | bgzip -c > "$fa"
            conda deactivate 2>/dev/null || true

            echo "  [$group.$position] Building Salmon index"
            conda activate salmon-1.10.2
            salmon index -t "$fa" -i "$idx_dir"
            conda deactivate 2>/dev/null || true

            echo "  [$group.$position] done"
        done
    done
    echo "[$(date)] txrevise Salmon indices done"
fi

# =============================================================================
# Done
# =============================================================================
echo "[$(date)] STAGE 3 (post-array) COMPLETE — Tier 1 reference prep finished"
echo "Outputs in: $REFERENCE_DIR"
ls -la "$REFERENCE_DIR"
echo "---"
ls -la "$TXREVISE_DIR"
