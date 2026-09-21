#!/bin/bash
# =============================================================================
# 01_refprep_pre.sh — Reference prep STAGE 1 of 3 (pre-array)
# =============================================================================
# This is the first of three stage scripts that replace the monolithic
# 00_reference_prep.sh. The txrevise event-construction step (formerly 5d) is
# now an LSF job array, so the pipeline is split into three jobs chained by
# LSF job-name dependencies:
#
#   01_refprep_pre.sh    -J refprep_pre      (this script)
#   02_txrevise_array.sh -J "txrevise[1-N]"  -w "done(refprep_pre)"
#   03_refprep_post.sh   -J refprep_post     -w "done(txrevise)"
#
# STAGE 1 produces everything needed BEFORE event construction:
#   - HPLRv2.0.annotated.PANTRY.gtf   (section 0)
#   - rediportal_hg38.bed             (section 1)
#   - exonic.gtf, intronic.gtf        (section 2)
#   - majiq/annotation.gff3           (section 3)
#   - edit_sites_to_genes.tsv         (section 4)
#   - txrevise/gene_annotations.gtf   (5a)
#   - txrevise/transcript_tags.txt    (5b)
#   - txrevise/txrevise_annotations.rds (5c)  <-- input to the STAGE 2 array
#
# Config loading: CONFIG and SCRIPTS_DIR env vars must be set via LSF -env at
# submission time:
#   bsub -env "CONFIG=\"config.yml\",SCRIPTS_DIR=\"/path/to/scripts\"" < 01_refprep_pre.sh
#
# DEPENDENCY CAVEAT: 02 and 03 depend on their predecessors BY JOB NAME
# (done(refprep_pre) / done(txrevise)). LSF resolves a name dependency only if a
# job with that name exists in the system at submission time. Submit all three
# together, in order; do not wait for 01 to finish and clear before submitting
# 02, or the name lookup can fail and the job may never leave PEND.
# =============================================================================

# ---- LSF directives ----
#BSUB -J refprep_pre
#BSUB -q long
#BSUB -n 12
#BSUB -M 120
#BSUB -R "rusage[mem=120]"
#BSUB -W 25:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/refprep_pre.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/refprep_pre.%J.err

set -eo pipefail

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# ---- Config loading (Tier 1: CONFIG and SCRIPTS_DIR via LSF -env) ----
CONFIG="${CONFIG:?Set CONFIG via bsub -env, e.g. -env \"CONFIG=\\\"config.yml\\\"\"}"
SCRIPTS_DIR="${SCRIPTS_DIR:?Set SCRIPTS_DIR via bsub -env}"

# Load all config values as shell vars (OUTPUT_BASE, REF_GENOME, REF_ANNO,
# NORMALIZED_GTF, REFERENCE_DIR, PANTRY_SCRIPTS, SEADRAGON_SCRIPTS,
# CLASSIFICATION, REFERENCE_GTF, REDIPORTAL_INPUT, TXREVISE_N_BATCHES,
# SOURCE_PREFIX, etc.)
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}")"

# ---- Paths from config ----
REF_GENOME="${REF_GENOME}"
RAW_GTF="${REF_ANNO}"
CLASSIFICATION="${CLASSIFICATION}"
REFERENCE_GTF="${REFERENCE_GTF}"
REFERENCE_DIR="${REFERENCE_DIR}"
PANTRY_SCRIPTS="${PANTRY_SCRIPTS}"
SEADRAGON_SCRIPTS="${SEADRAGON_SCRIPTS}"

# Source prefix for prep_PANTRY_gtf.R (GTF column 2 default + ##source header).
# Defaults to HPLR if not set in config.
SOURCE_PREFIX="${SOURCE_PREFIX:-HPLR}"

# Scratch log dir (infra — hardcoded, used by #BSUB directives above)
LOG_DIR="/rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs"

# Normalized GTF — produced by prep_PANTRY_gtf.R (step 0) or pre-existing
NORMALIZED_GTF="${NORMALIZED_GTF}"

# REDIportal raw table (from config; env var override still supported)
REDIPORTAL_INPUT="${REDIPORTAL_INPUT:-${REDIPORTAL_INPUT}}"

# Singularity R invocation (bind both /rsrch5 and /rsrch9) — infra, hardcoded
SING_R="singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript"
export R_LIBS_USER="/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"

# txrevise batch count (from config; must match #BSUB -J array range in 002)
TXREVISE_N_BATCHES="${TXREVISE_N_BATCHES}"
TXREVISE_DIR="$REFERENCE_DIR/txrevise"

# All downstream steps use the NORMALIZED GTF, not the raw GTF.
REF_ANNO="$NORMALIZED_GTF"

mkdir -p "$REFERENCE_DIR" "$LOG_DIR"
cd "$REFERENCE_DIR"

echo "[$(date)] STAGE 1 (pre-array) reference prep starting"
echo "  REF_GENOME:      $REF_GENOME"
echo "  RAW_GTF:         $RAW_GTF"
echo "  NORMALIZED_GTF:  $NORMALIZED_GTF"
echo "  REFERENCE_DIR:   $REFERENCE_DIR"
echo "  REDIPORTAL_INPUT: $REDIPORTAL_INPUT"
echo "  SEADRAGON_SCRIPTS: $SEADRAGON_SCRIPTS"
echo "  SOURCE_PREFIX:   $SOURCE_PREFIX"

# =============================================================================
# 0. GTF normalization (prep_PANTRY_gtf.R)
# =============================================================================
# Normalizes the SQANTI3-annotated HPLRv2 GTF to be PANTRY-compatible by adding:
#   - gene feature lines (from transcript boundaries, overridden by reference
#     GENCODE GTF for standard Ensembl gene IDs)
#   - gene_biotype / transcript_biotype attributes (from SQANTI3 coding status)
#   - gene_name attribute (from SQANTI3 associated_gene / reference GTF)
#   - CDS feature lines (from SQANTI3 CDS_genomic_start/end intersected with exons)
#   - tag "basic" on transcript lines
#
# The --source-prefix arg controls the GTF column 2 (source) default for gene
# feature lines whose input line has no source, and the ##source header comment.
# Default: HPLR. Set source_prefix in config.yml to override.
#
# This step is OPTIONAL: if you have already run prep_PANTRY_gtf.R interactively
# in RStudio (recommended — it gives you visibility into the SQANTI3 stats and
# any Ensembl gene validation errors), the normalized GTF will already exist and
# this step is skipped.
if [ -f "$NORMALIZED_GTF" ]; then
    echo "[$(date)] Step 0 SKIPPED — normalized GTF already exists: $NORMALIZED_GTF"
    echo "  (To regenerate, delete it and re-run, or run prep_PANTRY_gtf.R in RStudio)"
else
    echo "[$(date)] Step 0: Normalizing GTF with prep_PANTRY_gtf.R"
    $SING_R "$SEADRAGON_SCRIPTS/prep_PANTRY_gtf.R" \
        --gtf "$RAW_GTF" \
        --classification "$CLASSIFICATION" \
        --reference "$REFERENCE_GTF" \
        --out "$NORMALIZED_GTF" \
        --source-prefix "$SOURCE_PREFIX"
    if [ ! -f "$NORMALIZED_GTF" ]; then
        echo "ERROR: prep_PANTRY_gtf.R did not produce $NORMALIZED_GTF"
        echo "  Run it interactively in RStudio to see the full error output."
        exit 1
    fi
    echo "[$(date)] Step 0 done — normalized GTF: $NORMALIZED_GTF"
fi

# =============================================================================
# 1. REDIportal preprocessing
# =============================================================================
if [ -f "$REFERENCE_DIR/rediportal_hg38.bed" ]; then
    echo "[$(date)] Section 1 SKIPPED — REDIportal BED already exists: $REFERENCE_DIR/rediportal_hg38.bed"
elif [ -n "$REDIPORTAL_INPUT" ] && [ -f "$REDIPORTAL_INPUT" ]; then
    echo "[$(date)] Preprocessing REDIportal edit sites (keeping chr prefix)"
    # Use MAJIQ env for Python (has pandas etc.)
    conda activate samtools-1.16.1
    source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate
    python3 "$SEADRAGON_SCRIPTS/rediportal_preprocess.py" \
        --input "$REDIPORTAL_INPUT" \
        --output "$REFERENCE_DIR/rediportal_hg38.bed" \
        --verbose
    conda deactivate 2>/dev/null || true
    echo "[$(date)] REDIportal preprocessing done"
else
    echo "WARNING: REDIPORTAL_INPUT not set or not found: '$REDIPORTAL_INPUT'"
    echo "  RNA editing modality will not be available."
    echo "  Set rediportal_input in config.yml or REDIPORTAL_INPUT env var."
fi

EDIT_SITES_BED="$REFERENCE_DIR/rediportal_hg38.bed"
if [ ! -f "$EDIT_SITES_BED" ]; then
    echo "WARNING: $EDIT_SITES_BED not found. RNA editing modality will not be available."
    echo "  Set rediportal_input in config.yml to your REDIportal hg38 table."
fi

# =============================================================================
# 2. Exonic + intronic GTF (for featureCounts / stability)
# =============================================================================
if [ -f "$REFERENCE_DIR/exonic.gtf" ] && [ -f "$REFERENCE_DIR/intronic.gtf" ]; then
    echo "[$(date)] Section 2 SKIPPED — exonic/intronic GTFs already exist"
else
    echo "[$(date)] Building exonic/intronic GTFs"
    conda activate samtools-1.16.1
    source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate
    python3 "$PANTRY_SCRIPTS/exonic_intronic_from_gtf.py" \
        "$REF_ANNO" \
        "$REFERENCE_DIR/exonic.gtf" \
        "$REFERENCE_DIR/intronic.gtf" \
        --min-exon-fraction 0.8 \
        --verbose
    conda deactivate 2>/dev/null || true
    echo "[$(date)] Exonic/intronic GTFs done"
fi

# =============================================================================
# 3. MAJIQ GFF3 (for intron retention)
# =============================================================================
if [ -f "$REFERENCE_DIR/majiq/annotation.gff3" ]; then
    echo "[$(date)] Section 3 SKIPPED — MAJIQ GFF3 already exists: $REFERENCE_DIR/majiq/annotation.gff3"
else
    echo "[$(date)] Building MAJIQ GFF3"
    mkdir -p "$REFERENCE_DIR/majiq"
    conda activate samtools-1.16.1
    source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate
    python3 "$PANTRY_SCRIPTS/gtf_to_majiq_gff3.py" \
        --input "$REF_ANNO" \
        --output "$REFERENCE_DIR/majiq/annotation.gff3"
    conda deactivate 2>/dev/null || true
    echo "[$(date)] MAJIQ GFF3 done"
fi

# =============================================================================
# 4. Edit sites -> genes map (for RNA editing)
# =============================================================================
# Uses the ORIGINAL PANTRY map_edit_sites_to_genes.py. The normalized GTF now
# contains `gene` feature lines, so the original script works without modification.
if [ -f "$REFERENCE_DIR/edit_sites_to_genes.tsv" ]; then
    echo "[$(date)] Section 4 SKIPPED — edit sites -> genes map already exists: $REFERENCE_DIR/edit_sites_to_genes.tsv"
elif [ -f "$EDIT_SITES_BED" ]; then
    echo "[$(date)] Mapping edit sites to genes"
    conda activate samtools-1.16.1
    source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate
    python3 "$PANTRY_SCRIPTS/RNA_editing/map_edit_sites_to_genes.py" \
        "$EDIT_SITES_BED" \
        "$REF_ANNO" \
        "$REFERENCE_DIR/edit_sites_to_genes.tsv"
    conda deactivate 2>/dev/null || true
    echo "[$(date)] Edit sites -> genes map done"
fi

# =============================================================================
# 5 (pre). txrevise annotation prep: 5a preprocess, 5b tags, 5c prepareAnnotations
# =============================================================================
# Event construction (5d) runs as the STAGE 2 LSF job array
# (02_txrevise_array.sh); the merge (5e) + Salmon indices run in STAGE 3
# (03_refprep_post.sh).
mkdir -p "$TXREVISE_DIR"

# Skip the whole prep if the 6 final merged GFF3s already exist.
txrevise_all_exist=true
for group in grp_1 grp_2; do
    for position in upstream contained downstream; do
        [ -s "$TXREVISE_DIR/txrevise.${group}.${position}.gff3" ] || txrevise_all_exist=false
    done
done

if [ "$txrevise_all_exist" = true ]; then
    echo "[$(date)] txrevise prep SKIPPED — all 6 txrevise GFF3s already exist"
else
    mkdir -p "$TXREVISE_DIR/batch"

    # 5a. Preprocess GTF for txrevise (Python, needs gtfparse)
    # Uses the ORIGINAL PANTRY preprocess_gtf.py. The normalized GTF now contains
    # gene_biotype attributes and gene feature lines, so the original script works
    # without modification.
    if [ -f "$TXREVISE_DIR/gene_annotations.gtf" ]; then
        echo "  [5a] SKIPPED — gene_annotations.gtf already exists"
    else
        conda activate samtools-1.16.1
        source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate
        python3 "$PANTRY_SCRIPTS/txrevise/preprocess_gtf.py" \
            --input "$REF_ANNO" \
            --output "$TXREVISE_DIR/gene_annotations.gtf"
        conda deactivate 2>/dev/null || true
    fi

    # 5b. Extract transcript tags (Python)
    if [ -f "$TXREVISE_DIR/transcript_tags.txt" ]; then
        echo "  [5b] SKIPPED — transcript_tags.txt already exists"
    else
        conda activate samtools-1.16.1
        source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate
        python3 "$PANTRY_SCRIPTS/txrevise/extractTranscriptTags.py" \
            --gtf <(gzip -c "$TXREVISE_DIR/gene_annotations.gtf") \
            > "$TXREVISE_DIR/transcript_tags.txt"
        conda deactivate 2>/dev/null || true
    fi

    # 5c. Prepare annotations (R, via singularity) -> txrevise_annotations.rds
    if [ -f "$TXREVISE_DIR/txrevise_annotations.rds" ]; then
        echo "  [5c] SKIPPED — txrevise_annotations.rds already exists"
    else
        echo "[$(date)] txrevise prepareAnnotations.R"
        $SING_R "$PANTRY_SCRIPTS/txrevise/prepareAnnotations.R" \
            --gtf "$TXREVISE_DIR/gene_annotations.gtf" \
            --tags "$TXREVISE_DIR/transcript_tags.txt" \
            --out "$TXREVISE_DIR/txrevise_annotations.rds"
    fi
    echo "[$(date)] txrevise prep done"
fi

echo "[$(date)] STAGE 1 (pre-array) COMPLETE"
echo "  Next: 02_txrevise_array.sh (event construction), then 03_refprep_post.sh"
