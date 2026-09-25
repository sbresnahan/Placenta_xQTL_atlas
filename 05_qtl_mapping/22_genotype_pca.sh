#!/bin/bash
# =============================================================================
# 22_genotype_pca.sh — Genotype PCs via cohort-only PCA (GTEx convention)
# =============================================================================
# REPLACES 22_genotype_pcair.sh. Per GTEx convention, genotype PCs are
# computed by PCA on the cohort's OWN LD-pruned variants — no 1KG merge,
# no KING kinship, no PC-AiR, no reference-panel projection. This is much
# lighter than the old script and needs only plink2 + python.
#
# Per ancestry:
#   1. LD-prune the pooled per-ancestry pgen (--indep-pairwise 200 50 0.2,
#      same parameters as the retired PC-AiR script).
#   2. plink2 --pca 20 exact on the pruned variants (cohort samples only).
#      If this cluster's plink2 alpha rejects 'exact', drop the modifier —
#      the default approximate PCA is fine at this sample size.
#   3. genotype_pca_format.py reshapes .eigenvec/.eigenval into
#      {ANC}_genotype_pcs.tsv (sample_id, PC1..PC20 — the schema
#      24_outlier_exclusion.py already reads) + a scree plot.
#
# Usage:
#   bsub -env "CONFIG=\"config.yml\",SCRIPTS_DIR=\"/path/to/scripts\",ANCESTRIES=\"EAS EUR\"" < 22_genotype_pca.sh
#
# Required env:
#   CONFIG       — path to config.yml
#   SCRIPTS_DIR  — directory containing genotype_pca_format.py and config_get.py
# Optional env:
#   ANCESTRIES   — space-separated ancestry labels (default: EAS EUR)
#   GENO_DIR     — directory with {ANC}_pooled.{pgen,pvar,psam}
#   OUTPUT_DIR   — output directory (default: ${OUTPUT_BASE}/genotype_pcs)
#   N_PCS        — number of PCs to compute (default: 20)
# =============================================================================

#BSUB -q long
#BSUB -n 12
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 25:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/pca.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/pca.%J.err

set -eo pipefail

# ---- Config / env ----
CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
# Default SCRIPTS_DIR to this script's own directory, so the pipeline runs
# directly from the git clone. An explicit SCRIPTS_DIR env var overrides —
# and is REQUIRED when submitting via `bsub < script` (LSF executes a spool
# copy of the script; self-location would resolve to the spool directory).
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# config_get.py lives in ../03_phenotyping in the repo layout; a flat copy in
# SCRIPTS_DIR (legacy deployment) takes precedence.
CONFIG_GET="${SCRIPTS_DIR}/config_get.py"
[ -f "$CONFIG_GET" ] || CONFIG_GET="${SCRIPTS_DIR}/../03_phenotyping/config_get.py"
if [ ! -f "$CONFIG_GET" ]; then
    echo "ERROR: config_get.py not found in $SCRIPTS_DIR or $SCRIPTS_DIR/../03_phenotyping" >&2
    echo "  Submitting via 'bsub <'? Export SCRIPTS_DIR=<repo>/05_qtl_mapping first." >&2
    exit 1
fi
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
N_PCS="${N_PCS:-20}"

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
module load plink

# Load config values
eval "$(python3 "$CONFIG_GET" "${CONFIG}")"
OUTPUT_BASE="${OUTPUT_BASE}"

# ---- Paths ----
GENO_DIR="${GENO_DIR:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/genotypes}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_BASE}/genotype_pcs}"
WORK_DIR="${OUTPUT_DIR}/intermediate"

mkdir -p "$OUTPUT_DIR" "$WORK_DIR"

echo "[$(date)] Genotype PCA Pipeline (cohort-only, GTEx convention)"
echo "  CONFIG:       $CONFIG"
echo "  SCRIPTS_DIR:  $SCRIPTS_DIR"
echo "  GENO_DIR:     $GENO_DIR"
echo "  OUTPUT_DIR:   $OUTPUT_DIR"
echo "  ANCESTRIES:   $ANCESTRIES"
echo "  N_PCS:        $N_PCS"
echo ""

# tensorQTL env provides pandas + matplotlib for the format/scree step
conda activate tensorqtl

# ---- Per-ancestry loop ----
for ANC in $ANCESTRIES; do
    echo "============================================================"
    echo "[$(date)] Ancestry: $ANC"
    echo "============================================================"

    POOLED_PREFIX="${GENO_DIR}/${ANC}_pooled"
    if [ ! -f "${POOLED_PREFIX}.pgen" ]; then
        echo "  ERROR: pooled pgen not found: ${POOLED_PREFIX}.pgen — skipping $ANC"
        continue
    fi

    OUT_TSV="${OUTPUT_DIR}/${ANC}_genotype_pcs.tsv"
    if [ -f "$OUT_TSV" ]; then
        echo "  PCs already exist: $OUT_TSV — skipping (delete to regenerate)"
        continue
    fi

    PRUNE_PREFIX="${WORK_DIR}/${ANC}_prune_stats"
    PCA_PREFIX="${WORK_DIR}/${ANC}_pca"

    # ---- Step 1: LD pruning (cohort-only) ----
    if [ -f "${PRUNE_PREFIX}.prune.in" ]; then
        echo "  [Step 1] LD pruning: already exists, skipping"
    else
        echo "  [Step 1] LD pruning (200 kb / 50 kb / r2 0.2)..."
        plink2 --pfile "${POOLED_PREFIX}" \
            --indep-pairwise 200 50 0.2 \
            --threads 12 \
            --out "${PRUNE_PREFIX}"
        if [ ! -f "${PRUNE_PREFIX}.prune.in" ]; then
            echo "    ERROR: LD pruning failed for $ANC"
            exit 1
        fi
    fi
    echo "    Pruned variants: $(wc -l < "${PRUNE_PREFIX}.prune.in")"

    # ---- Step 2: PCA on pruned cohort variants ----
    if [ -f "${PCA_PREFIX}.eigenvec" ]; then
        echo "  [Step 2] PCA: already exists, skipping"
    else
        echo "  [Step 2] plink2 --pca ${N_PCS} exact (cohort-only)..."
        if ! plink2 --pfile "${POOLED_PREFIX}" \
                --extract "${PRUNE_PREFIX}.prune.in" \
                --pca ${N_PCS} exact \
                --threads 12 \
                --out "${PCA_PREFIX}"; then
            echo "    WARN: '--pca exact' failed on this plink2 build; retrying with default (approximate) PCA"
            plink2 --pfile "${POOLED_PREFIX}" \
                --extract "${PRUNE_PREFIX}.prune.in" \
                --pca ${N_PCS} \
                --threads 12 \
                --out "${PCA_PREFIX}"
        fi
        if [ ! -f "${PCA_PREFIX}.eigenvec" ]; then
            echo "    ERROR: PCA failed for $ANC"
            exit 1
        fi
    fi

    # ---- Step 3: Format PCs + scree plot ----
    echo "  [Step 3] Formatting PCs + scree plot..."
    python3 "${SCRIPTS_DIR}/genotype_pca_format.py" \
        --eigenvec "${PCA_PREFIX}.eigenvec" \
        --eigenval "${PCA_PREFIX}.eigenval" \
        --ancestry "${ANC}" \
        --n-pcs "${N_PCS}" \
        --output-dir "${OUTPUT_DIR}"

    echo "  Done: $ANC"
    echo ""
done

# ---- Summary ----
echo "[$(date)] Pipeline Complete"
echo ""
echo "  Genotype PCs:"
ls -lh "${OUTPUT_DIR}"/*_genotype_pcs.tsv 2>/dev/null || echo "    (none found)"
echo ""
echo "  Scree plots:"
ls -lh "${OUTPUT_DIR}"/*_genotype_pcs_scree.png 2>/dev/null || echo "    (none found)"
