# Module 06 Runbook — GWAS Colocalization and isoTWAS

Module: `06_colocalization_twas`
Repository: `Placenta_xQTL_atlas`

This module integrates the Module 05 xQTL/fine-mapping outputs with external GWAS using:

1. GWAS acquisition and harmonization.
2. Genome-wide nominal xQTL statistics.
3. Pairwise SuSiE-coloc.
4. Multi-trait colocBoost.
5. isoTWAS/TWAS weight training.
6. FUSION association testing.
7. Aggregation and cell-type annotation.

The executable sequence is:

```text
36  install R/FUSION environment
37  fetch GWAS
38  harmonize GWAS
39  genome-wide nominal xQTL statistics
40  prepare SuSiE-coloc tasks
41  SuSiE-coloc worker
42  submit SuSiE-coloc arrays
43  prepare colocBoost regions/outcomes
44  colocBoost worker
45  submit colocBoost arrays
46  isoTWAS/TWAS training worker
47  submit isoTWAS/TWAS training arrays
48  FUSION TWAS association testing
49  aggregate results + cell-type annotation
```

---

## Revision notes (2026-10-09 patch)

This revision fixes four operational issues:

1. **R package library is always set inside R.** Every R entry point in this
   module now guarantees `.libPaths()` includes the lab R 4.3.1 library
   (`/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1`)
   from within the R process. Bash-level `R_LIBS`/`R_LIBS_USER` are **not**
   relied upon — they are not reliably propagated to R on this cluster.
   `41_susie_coloc.R`, `44_colocboost.R`, and `46_isotwas_train.R` already
   contained the in-R `.libPaths()` call; the external `FUSION.assoc_test.R`
   does not, so step 48 now runs it through the repo wrapper
   `48b_fusion_assoc.R`, which sets `.libPaths()` and then sources FUSION
   (FUSION's `optparse` flags pass through unchanged).

2. **Every script specifies its required conda environment / seadragon
   module.** Scripts that invoke `python3` with non-stdlib imports activate
   the pipeline conda env themselves (`PYENV`, default `tensorqtl` from
   module-05 script 21) and fail fast with a clear message if pandas is still
   missing. Compute-node workers load `module load plink samtools` (plink2 +
   tabix) themselves instead of relying on the submission shell's PATH. See
   the per-step table in Section 0.

3. **`TEST=1` is now a safe pilot for steps 42, 45, and 47.** Previously
   these submitters replaced `N_SHARDS=1` when `TEST=1`, and because the R
   workers assign tasks round-robin by `N_SHARDS`, the "pilot" shard received
   **every** task/region/gene. The submitters now keep the true `N_SHARDS`
   in the worker environment and submit only array index `[1]`, so a pilot
   processes approximately `SHARD_SIZE` items. The manual pilot commands in
   the older runbook are no longer needed (they remain valid).

4. **Module 06 outputs are compatible with module 07.** Two contract fixes:
   (a) module 07's tier-2 resolver (`54_tier2_stratified.py`) now looks in
   `coloc/aggregated/` — the directory 49 actually writes — before its older
   fallback locations, and accepts the `GENE` column (FUSION `.pos`
   convention) in `twas_gene_results.tsv.gz`; (b) a new preflight audit,
   `check_collapsed_inputs.py`, verifies that module 06's sample-level inputs
   (`qtl_inputs` BEDs, covariates, pgen) match the collapsed-replicate
   retained-run authority (`replicate_collapsed/reports/
   ancestry_map_collapsed.tsv`) that module 07 enforces. See Section 10.

Also fixed: `37_fetch_gwas.sh` no longer prints a stale next-step hint
pointing at `${GWAS_DIR}/harmonized`, and `42_submit_coloc.sh` now points to
the correct aggregator (`49_aggregate_coloc_twas.py`, run with `python3`).

### GWAS directory layout (unchanged contract)

```text
$GWAS_DIR/
├── raw/
│   └── {trait_id}.txt.gz
├── {trait_id}.sumstats.tsv.gz
├── {trait_id}.sumstats.tsv.gz.tbi
├── harmonization_qc.tsv
├── 1kg_sample_superpop.tsv
├── fusion/
└── fusion_ldref/
```

Harmonized GWAS files are written **directly into `$GWAS_DIR`**; steps
40/42/43/45/48 expect `${GWAS_DIR}/{trait_id}.sumstats.tsv.gz`.

---

# 0. Base environment

Start from the repository root.

```bash
cd /path/to/Placenta_xQTL_atlas

export REPO_ROOT="$PWD"
export SCRIPTS_DIR="${REPO_ROOT}/06_colocalization_twas"

export CONFIG=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/config.yml
export OUTPUT_BASE=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY

export QTL_DIR="${OUTPUT_BASE}/qtl_inputs"
export RESULTS_DIR="${OUTPUT_BASE}/qtl_results"
export LOG_DIR="${OUTPUT_BASE}/logs"

export GWAS_DIR="${OUTPUT_BASE}/gwas"
export COLOC_DIR="${RESULTS_DIR}/coloc"
export CB_DIR="${COLOC_DIR}/colocboost"
export ISOTWAS_DIR="${RESULTS_DIR}/isotwas"

# Module-01 / pooled genotype resources used for GWAS harmonization.
export PGEN_DIR=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/genotypes
export CHAIN=/home/stbresnahan/bhattacharya_lab/data/GenomicReferences/liftover/hg19ToHg38.over.chain

# 1000 Genomes GRCh38 reference.
export KG_PGEN=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/1kGP/1kGP_hg38
export KG_SAMPLE_MAP="${GWAS_DIR}/1kg_sample_superpop.tsv"

# Module-06 R/FUSION software.
export RSCRIPT="${REPO_ROOT}/bin/Rscript_sif"
export FUSION_DIR=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/fusion_twas

# Collapsed-replicate retained-run authority (module 05; shared with module 07).
export COLLAPSED_ANCESTRY_MAP="${OUTPUT_BASE}/replicate_collapsed/reports/ancestry_map_collapsed.tsv"

mkdir -p "$GWAS_DIR" "$LOG_DIR"
```

`REPO_ROOT`, `PGEN_DIR`, and `CHAIN` above are runbook convenience variables.
Individual module drivers do not necessarily consume them directly.

## 0.1 Environments at a glance

Every step's required conda environment or seadragon module. "In-script"
means the script activates/loads it itself; "manual" means you must activate
it before running (the step has no shell driver).

| Step | Runs on | Required environment | How it is provided |
|---|---|---|---|
| 36 install | login | R 4.3.1 singularity container; `git`; internet | `$RSCRIPT` wrapper; no conda env |
| 37 fetch | login | conda env `tensorqtl` (pandas); internet | **in-script** (`PYENV` override; `PYENV=none` skips) |
| 38 harmonize | login | conda env `tensorqtl` (pandas, numpy) + `bgzip`/`tabix` | **manual**: `conda activate tensorqtl` + `module load samtools` |
| 39 nominal | login driver | driver: `python3` stdlib only (collapse audit) + LSF | — |
| 39 workers | compute | conda env `tensorqtl` (tensorQTL, torch); merge job adds `module load samtools` | in-script (bsub heredocs) |
| 40/42 coloc | login driver | conda env `tensorqtl` (pandas for step 40) | **in-script** (`PYENV` override) |
| 41/42a workers | compute | `module load plink samtools` (plink2, tabix); R via `$RSCRIPT` | in-script (42a) |
| 43/45 colocBoost | login driver | conda env `tensorqtl` (pandas for step 43) | **in-script** (`PYENV` override) |
| 44/45a workers | compute | `module load plink samtools`; R via `$RSCRIPT` | in-script (45a) |
| 46/47 isoTWAS | login driver | `python3` stdlib only (gene counts, collapse audit) | — |
| 46/47a workers | compute | `module load plink samtools`; R via `$RSCRIPT` | in-script (47a) |
| 48 FUSION | login driver | `plink2` (LD reference build) | in-script: `module load plink` attempted if missing |
| 48a workers | compute | R via `$RSCRIPT` through `48b_fusion_assoc.R` (sets `.libPaths()` in R) | in-script |
| 49 aggregate | login | conda env `tensorqtl` (pandas, numpy, scipy) | **manual**: `conda activate tensorqtl` |

The `tensorqtl` conda env is created by module-05 `21_install_tensorqtl.sh`
(python 3.10; provides pandas/numpy/scipy). Any env with those packages works;
override with `PYENV=<name>` for the scripts that self-activate.

**R package library rule for this module:** any R script or `Rscript -e`
snippet must begin with

```r
.libPaths(c("/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1", .libPaths()))
```

Do not rely on exported `R_LIBS`/`R_LIBS_USER` — they do not reliably reach R
here. All module R entry points (36 heredocs, 41, 44, 46, 48b) already follow
this rule.

## 0.2 Preflight

```bash
test -f "$CONFIG"
test -d "$SCRIPTS_DIR"
test -d "$QTL_DIR"
test -d "$RESULTS_DIR"

for ANC in EAS EUR; do
    test -f "${QTL_DIR}/${ANC}_qtl.pgen"
    test -f "${QTL_DIR}/${ANC}_qtl.pvar"
    test -f "${QTL_DIR}/${ANC}_qtl.psam"

    test -f "${PGEN_DIR}/${ANC}_pooled.pvar"

    test -f "${QTL_DIR}/${ANC}_expression.bed.gz"
    test -f "${QTL_DIR}/${ANC}_isoform_expression.bed.gz"
    test -f "${QTL_DIR}/${ANC}_covariates_expression.tsv"
    test -f "${QTL_DIR}/${ANC}_covariates_isoform_expression.tsv"
done

test -f "${KG_PGEN}.pgen"
test -f "${KG_PGEN}.pvar"
test -f "${KG_PGEN}.psam"
test -f "$KG_SAMPLE_MAP"
test -f "$CHAIN"

command -v bsub
command -v bjobs
command -v singularity   # required by $RSCRIPT on login and compute nodes
command -v plink2 || module load plink
command -v bgzip || module load samtools
command -v tabix
```

### Collapsed-replicate audit (required)

Module 06 must operate on the same collapsed-replicate sample set as module
07. Run the audit once before the first compute step (it also runs
automatically inside `39_run_nominal.sh` and `47_submit_isotwas.sh`):

```bash
python3 "${SCRIPTS_DIR}/check_collapsed_inputs.py" \
    --qtl-dir "$QTL_DIR" \
    --output-base "$OUTPUT_BASE" \
    --ancestries EAS EUR
```

The script is stdlib-only (no conda env needed). It fails closed if any
qtl_inputs sample column / psam IID is absent from the collapsed ancestry
map, or if duplicate sample columns (uncollapsed replicates) are present.
Escape hatch for cohorts with no technical replicates where
`collapse_replicates.py` was intentionally not run: `SKIP_COLLAPSE_CHECK=1`.
See Section 10 for the full contract.

Fine-mapping locus lists from Module 05 are required by coloc and colocBoost:

```bash
ls -lh "${RESULTS_DIR}/finemap/loci/"*.loci.tsv
```

Expected modalities are:

```text
expression
isoforms
isoform_expression
splicing
intron_retention
alt_TSS
alt_polyA
RNA_editing
stability
```

---

# 1. Step 36 — install Module 06 software

Run this **once on a login node**, because package installation and the
FUSION clone require internet access. No conda environment is needed: R runs
through the `$RSCRIPT` singularity wrapper and every R heredoc in the script
sets `.libPaths()` internally.

```bash
cd "$SCRIPTS_DIR"

bash 36_install_coloc_env.sh
```

## Environment controls

| Variable | Default | Purpose |
|---|---|---|
| `SCRIPTS_DIR` | directory containing script | Locate repository/module |
| `RSCRIPT` | `$REPO_ROOT/bin/Rscript_sif` | R executable/wrapper |
| `FUSION_DIR` | `/rsrch5/.../software/fusion_twas` | FUSION checkout |
| `SKIP_R` | `0` | Set `1` to skip R package installation |

Examples:

```bash
# Normal install
bash 36_install_coloc_env.sh

# R packages already installed; only ensure FUSION exists
SKIP_R=1 bash 36_install_coloc_env.sh

# Override R launcher
RSCRIPT=/path/to/Rscript bash 36_install_coloc_env.sh
```

The script installs/checks:

```text
coloc
susieR
glmnet
optparse
data.table
R.utils
remotes
matrixStats
irlba
colocboost
plink2R
```

and clones `gusevlab/fusion_twas`.

Verification — note the mandatory in-R `.libPaths()` line:

```bash
"$RSCRIPT" -e '
.libPaths(c("/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1", .libPaths()))
library(coloc)
library(susieR)
library(glmnet)
library(colocboost)
cat("coloc:", as.character(packageVersion("coloc")), "\n")
cat("susieR:", as.character(packageVersion("susieR")), "\n")
cat("glmnet:", as.character(packageVersion("glmnet")), "\n")
cat("colocboost:", as.character(packageVersion("colocboost")), "\n")
'

test -f "${FUSION_DIR}/FUSION.assoc_test.R"
```

---

# 2. Step 37 — fetch GWAS summary statistics

Run on the **login node** because compute nodes do not have internet access.

**Required environment:** `37_fetch_gwas.py` imports pandas, which the
login-node system `python3` does not provide. The driver now activates the
`tensorqtl` conda env itself (override with `PYENV=<env>`; `PYENV=none`
skips activation) and fails fast with a clear message if pandas is still not
importable.

```bash
cd "$SCRIPTS_DIR"

bash 37_fetch_gwas.sh
```

Raw files are written to:

```text
$GWAS_DIR/raw/{trait_id}.txt.gz
```

## Environment controls

| Variable | Default | Purpose |
|---|---|---|
| `SCRIPTS_DIR` | script directory | Contains `gwas_catalog.tsv` |
| `OUTPUT_BASE` | PANTRY output root | Used to derive `GWAS_DIR` |
| `GWAS_DIR` | `$OUTPUT_BASE/gwas` | GWAS working directory |
| `TRAITS` | empty | Space-separated subset of trait IDs |
| `FORCE` | `0` | `1` = re-download existing files |
| `PYENV` | `tensorqtl` | Conda env providing pandas (`none` skips activation) |

Examples:

```bash
# All catalog traits
bash 37_fetch_gwas.sh

# One trait
TRAITS="egg_bw_fetal_2019" bash 37_fetch_gwas.sh

# Multiple selected traits
TRAITS="egg_bw_fetal_2019 egg_ga_maternal_2023" \
    bash 37_fetch_gwas.sh

# Re-download
FORCE=1 TRAITS="egg_bw_fetal_2019" \
    bash 37_fetch_gwas.sh
```

Check the manifest:

```bash
column -t -s $'\t' "${GWAS_DIR}/raw/fetch_manifest.tsv" | less -S
```

Manual-access traits must be placed at the filename expected by the
catalog/step 37 contract under:

```text
$GWAS_DIR/raw/
```

Before proceeding, confirm the raw files you intend to analyze exist:

```bash
ls -lh "${GWAS_DIR}/raw/"*.txt.gz
```

---

# 3. Step 38 — harmonize GWAS to GRCh38 and atlas variants

There is no shell driver for Step 38. Run the Python script directly.

**Required environment (manual):**

```bash
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl     # pandas, numpy
module load samtools         # bgzip, tabix — required: 41/44 slice the
                             # tabix-indexed outputs
```

The runbook deliberately writes harmonized GWAS files **directly under
`$GWAS_DIR`** so that Steps 40–48 find them without changing `GWAS_DIR`.

```bash
python3 "${SCRIPTS_DIR}/38_harmonize_gwas.py" \
    --catalog "${SCRIPTS_DIR}/gwas_catalog.tsv" \
    --raw-dir "${GWAS_DIR}/raw" \
    --out-dir "${GWAS_DIR}" \
    --chain "${CHAIN}" \
    --pgen-dir "${PGEN_DIR}" \
    --ancestries EAS EUR
```

Expected output:

```text
$GWAS_DIR/{trait_id}.sumstats.tsv.gz
$GWAS_DIR/{trait_id}.sumstats.tsv.gz.tbi
$GWAS_DIR/harmonization_qc.tsv
```

## Step 38 command-line controls

Step 38 does not consume analysis controls through environment variables. Its
flags are:

| Flag | Required/default | Purpose |
|---|---|---|
| `--catalog` | required | `gwas_catalog.tsv` |
| `--raw-dir` | required | Raw GWAS directory |
| `--out-dir` | required | Harmonized output directory |
| `--chain` | required when any build-37 trait is included | hg19 → GRCh38 chain |
| `--pgen-dir` | required | Directory containing `{ANC}_pooled.pvar` |
| `--ancestries` | `EAS EUR` | Pooled variant panels to index |
| `--traits` | all | Restrict harmonization to selected traits |
| `--column-map` | none | Per-trait column overrides |

Subset example:

```bash
python3 "${SCRIPTS_DIR}/38_harmonize_gwas.py" \
    --catalog "${SCRIPTS_DIR}/gwas_catalog.tsv" \
    --raw-dir "${GWAS_DIR}/raw" \
    --out-dir "${GWAS_DIR}" \
    --chain "${CHAIN}" \
    --pgen-dir "${PGEN_DIR}" \
    --ancestries EAS EUR \
    --traits egg_bw_fetal_2019
```

Column override example:

```bash
--column-map \
    my_trait:beta=BETA_COLUMN \
    my_trait:se=SE_COLUMN
```

QC:

```bash
column -t -s $'\t' "${GWAS_DIR}/harmonization_qc.tsv" | less -S

for F in "${GWAS_DIR}"/*.sumstats.tsv.gz; do
    echo "$F"
    zcat "$F" | head -2
    test -f "${F}.tbi" || echo "WARNING: missing tabix index"
done
```

A `WARNING: missing tabix index` here means bgzip/tabix were not available
during harmonization (`harmonization_qc.tsv` will show `tabix = False`);
rerun with `module load samtools` — the coloc/colocBoost workers cannot slice
unindexed files.

---

# 4. Step 39 — genome-wide nominal xQTL statistics

This reruns tensorQTL in nominal mode for every ancestry × modality ×
chromosome and merges the chromosome shards into tabix-indexed files.

**Required environment:** the driver itself needs only LSF and stdlib
`python3` (for the collapsed-replicate audit). The compute-node workers
activate the `tensorqtl` conda env inside the bsub heredocs, and the merge
jobs additionally `module load samtools` (bgzip/tabix). No manual activation
is needed.

Before any submission, the driver runs `check_collapsed_inputs.py` (Section
0.2) and aborts if the qtl_inputs sample sets do not match the
collapsed-replicate contract. Override with `SKIP_COLLAPSE_CHECK=1`.

## Environment controls

| Variable | Default | Purpose |
|---|---|---|
| `CONFIG` | required | Pipeline config |
| `SCRIPTS_DIR` | required | Module 06 directory |
| `OUTPUT_BASE` | `dirname $CONFIG` | Output root |
| `QTL_DIR` | `$OUTPUT_BASE/qtl_inputs` | Module-05 QTL inputs |
| `RESULTS_DIR` | `$OUTPUT_BASE/qtl_results` | Result root |
| `LOG_DIR` | `$OUTPUT_BASE/logs` | LSF logs |
| `ANCESTRIES` | `EAS EUR` | Ancestry subset |
| `MODALITIES` | all 9 | Modality subset |
| `CHROMS` | `1 ... 22` | Chromosomes |
| `QUEUE` | `medium` | LSF queue |
| `WALLTIME` | `12:00` | Per chromosome |
| `THREADS` | `4` | LSF CPU request |
| `MEM` | `32G` | LSF memory request |
| `GPU` | `0` | Set `1` for GPU submission |
| `GPU_QUEUE` | `gpu` | GPU queue |
| `MAF_THRESHOLD` | `0.01` | tensorQTL MAF threshold |
| `TEST` | `0` | `1` = EAS/expression/chr21 |
| `FORCE` | `0` | `1` = rerun existing chromosome shards |
| `SKIP_COLLAPSE_CHECK` | `0` | `1` = skip the collapsed-replicate audit |

### Safe pilot

```bash
cd "$SCRIPTS_DIR"

TEST=1 bash 39_run_nominal.sh
```

`TEST=1` overrides:

```text
ANCESTRIES=EAS
MODALITIES=expression
CHROMS=21
```

Check:

```bash
bjobs -J 'nom_EAS_expression*'

ls -lh \
    "${RESULTS_DIR}/nominal/EAS/EAS_expression.nominal.chr21.parquet"
```

### Full run

```bash
TEST=0 bash 39_run_nominal.sh
```

Examples of selective runs:

```bash
ANCESTRIES="EAS" \
MODALITIES="expression splicing" \
bash 39_run_nominal.sh
```

```bash
ANCESTRIES="EUR" \
MODALITIES="isoform_expression" \
CHROMS="1 2 3 4 5" \
bash 39_run_nominal.sh
```

GPU example:

```bash
GPU=1 \
GPU_QUEUE=gpu \
THREADS=4 \
MEM=32G \
bash 39_run_nominal.sh
```

Force rerun:

```bash
FORCE=1 \
ANCESTRIES="EAS" \
MODALITIES="expression" \
CHROMS="21" \
bash 39_run_nominal.sh
```

### Worker options not exposed by `39_run_nominal.sh`

`39_run_nominal.py` additionally supports:

```text
--cis-window
--covariates-file
--allow-missing-chroms
```

The shell driver does **not** currently provide environment-variable mappings
for those options.

If one must be changed, invoke the Python worker directly or patch the shell
driver.

### Completion check

```bash
for ANC in EAS EUR; do
    for MOD in \
        expression isoforms isoform_expression splicing \
        intron_retention alt_TSS alt_polyA RNA_editing stability
    do
        F="${RESULTS_DIR}/nominal/${ANC}/${ANC}_${MOD}.nominal.tsv.gz"
        printf "%-4s %-20s " "$ANC" "$MOD"
        if [ -f "$F" ] && [ -f "${F}.tbi" ]; then
            echo OK
        else
            echo MISSING
        fi
    done
done
```

---

# 5. Steps 40–42 — pairwise SuSiE-coloc

Step 40 constructs:

```text
$COLOC_DIR/loci/{MOD}.tasks.tsv
$COLOC_DIR/loci/EAS.1kg.keep
$COLOC_DIR/loci/EUR.1kg.keep
```

Step 42 normally runs Step 40 automatically and then submits Step 41 workers.

**Required environment:** the driver activates the `tensorqtl` conda env
itself (pandas for step 40; `PYENV` override). The compute-node workers
(`42a`) load `module load plink samtools` (plink2 for LD, tabix for summary
-statistic slicing) and run R through `$RSCRIPT`; `41_susie_coloc.R` sets
`.libPaths()` internally.

## Environment controls for `42_submit_coloc.sh`

| Variable | Default | Purpose |
|---|---|---|
| `CONFIG` | required | Pipeline config |
| `SCRIPTS_DIR` | required | Module 06 directory |
| `OUTPUT_BASE` | `dirname $CONFIG` | Output root |
| `QTL_DIR` | `$OUTPUT_BASE/qtl_inputs` | xQTL inputs |
| `RESULTS_DIR` | `$OUTPUT_BASE/qtl_results` | Results |
| `LOG_DIR` | `$OUTPUT_BASE/logs` | LSF logs |
| `COLOC_DIR` | `$RESULTS_DIR/coloc` | Coloc output |
| `GWAS_DIR` | `$OUTPUT_BASE/gwas` | **Must contain harmonized GWAS directly** |
| `KG_PGEN` | hard-coded 1KG default | 1KG pgen prefix |
| `KG_SAMPLE_MAP` | `$GWAS_DIR/1kg_sample_superpop.tsv` | Sample → superpopulation map |
| `MODALITIES` | all 9 | Modalities to process |
| `SHARD_SIZE` | `25` | Target tasks per worker |
| `QUEUE` | `medium` | LSF queue |
| `WALLTIME` | `04:00` | Per shard |
| `THREADS` | `2` | LSF CPU request |
| `MIN_VARIANTS` | `50` | Minimum common variants |
| `PP_H4` | `0.7` | Colocalization threshold |
| `TEST` | `0` | `1` = submit only shard `[1]` of the first modality (safe pilot) |
| `FORCE_TASKS` | `0` | `1` = rebuild task tables |
| `FORCE_RUN` | `0` | `1` = clear `.done` and rerun |
| `RSCRIPT` | `$REPO_ROOT/bin/Rscript_sif` | R launcher |
| `PYENV` | `tensorqtl` | Conda env providing pandas for step 40 |

### Prepare task lists explicitly

Doing Step 40 explicitly makes the workflow easier to inspect before
submission.

```bash
python3 "${SCRIPTS_DIR}/40_prepare_coloc_loci.py" \
    --results-dir "$RESULTS_DIR" \
    --qtl-dir "$QTL_DIR" \
    --gwas-dir "$GWAS_DIR" \
    --catalog "${SCRIPTS_DIR}/gwas_catalog.tsv" \
    --kg-pgen "$KG_PGEN" \
    --kg-sample-map "$KG_SAMPLE_MAP" \
    --ancestries EAS EUR \
    --out-dir "$COLOC_DIR"
```

(If run standalone as above, activate the pandas env first:
`conda activate tensorqtl`. Via `42_submit_coloc.sh` this is automatic.)

Check:

```bash
ls -lh "${COLOC_DIR}/loci/"*.tasks.tsv
ls -lh "${COLOC_DIR}/loci/"*.1kg.keep

for F in "${COLOC_DIR}/loci/"*.tasks.tsv; do
    echo "$(basename "$F"): $(( $(wc -l < "$F") - 1 )) tasks"
done
```

### Safe pilot

```bash
TEST=1 bash 42_submit_coloc.sh
```

`TEST=1` submits only array index `[1]` of the first modality while keeping
the true `N_SHARDS` in the worker environment, so the pilot processes
approximately `SHARD_SIZE` tasks (25 by default) — not the whole modality.
(This was **not** true in earlier revisions; see the revision notes.)

### Full submission

After the pilot is satisfactory:

```bash
TEST=0 \
SHARD_SIZE=25 \
MIN_VARIANTS=50 \
PP_H4=0.7 \
bash 42_submit_coloc.sh
```

Subset:

```bash
MODALITIES="expression isoform_expression splicing" \
SHARD_SIZE=25 \
PP_H4=0.7 \
bash 42_submit_coloc.sh
```

Rebuild task lists:

```bash
FORCE_TASKS=1 bash 42_submit_coloc.sh
```

Rerun completed tasks:

```bash
FORCE_RUN=1 bash 42_submit_coloc.sh
```

Be cautious with `FORCE_RUN=1`: the worker deletes existing `.done` sentinels
for its assigned tasks before rerunning them.

### Step-41 flags not exposed through the Step-42 environment

`41_susie_coloc.R` additionally supports:

```text
--plink2
--gwas-L
--max-L
--tmp-dir
```

Current `42_submit_coloc.sh` does not map environment variables to those
flags.

Defaults are:

```text
--plink2       plink2
--gwas-L       10
--max-L        20
--tmp-dir      R tempdir()
```

---

# 6. Steps 43–45 — multi-trait colocBoost

Step 43 constructs per-ancestry region and outcome manifests.

Step 45 submits Step 44 workers.

**Required environment:** same pattern as steps 40–42 — the driver activates
the `tensorqtl` conda env itself (pandas for step 43; `PYENV` override), and
the workers (`45a`) load `module load plink samtools` and run R through
`$RSCRIPT` (`44_colocboost.R` sets `.libPaths()` internally).

Prerequisite: the 1KG ancestry keep files from Step 40 must already exist:

```bash
test -f "${COLOC_DIR}/loci/EAS.1kg.keep"
test -f "${COLOC_DIR}/loci/EUR.1kg.keep"
```

## Environment controls for `45_submit_colocboost.sh`

| Variable | Default | Purpose |
|---|---|---|
| `CONFIG` | required | Pipeline config |
| `SCRIPTS_DIR` | required | Module directory |
| `OUTPUT_BASE` | `dirname $CONFIG` | Output root |
| `QTL_DIR` | `$OUTPUT_BASE/qtl_inputs` | QTL inputs |
| `RESULTS_DIR` | `$OUTPUT_BASE/qtl_results` | Results |
| `LOG_DIR` | `$OUTPUT_BASE/logs` | Logs |
| `COLOC_DIR` | `$RESULTS_DIR/coloc` | Coloc directory |
| `CB_DIR` | `$COLOC_DIR/colocboost` | colocBoost directory |
| `GWAS_DIR` | `$OUTPUT_BASE/gwas` | Harmonized GWAS |
| `KG_PGEN` | hard-coded 1KG default | GWAS-side LD |
| `ANCESTRIES` | `EAS EUR` | Ancestries |
| `MODALITIES` | all 9 | xQTL modalities |
| `SHARD_SIZE` | `10` | Target regions per worker |
| `QUEUE` | `medium` | LSF queue |
| `WALLTIME` | `06:00` | Per worker |
| `THREADS` | `2` | CPUs |
| `MERGE_GAP` | `100000` | Locus merge distance in bp |
| `TEST` | `0` | `1` = submit only shard `[1]` of the first ancestry (safe pilot) |
| `FORCE_PREP` | `0` | `1` = rebuild manifests |
| `FORCE_RUN` | `0` | `1` = rerun completed regions |
| `RSCRIPT` | repository R wrapper | R launcher |
| `PYENV` | `tensorqtl` | Conda env providing pandas for step 43 |
| `KG_SAMPLE_MAP` | defined by script | Currently unused by Step 45 |

### Prepare manifests explicitly

```bash
python3 "${SCRIPTS_DIR}/43_prepare_colocboost.py" \
    --results-dir "$RESULTS_DIR" \
    --qtl-dir "$QTL_DIR" \
    --gwas-dir "$GWAS_DIR" \
    --catalog "${SCRIPTS_DIR}/gwas_catalog.tsv" \
    --ancestries EAS EUR \
    --modalities \
        expression isoforms isoform_expression splicing \
        intron_retention alt_TSS alt_polyA RNA_editing stability \
    --merge-gap 100000 \
    --out-dir "$CB_DIR"
```

(Standalone: activate `tensorqtl` first. Via `45_submit_colocboost.sh` this is
automatic.)

Check:

```bash
for ANC in EAS EUR; do
    echo "$ANC"
    wc -l "${CB_DIR}/${ANC}.regions.tsv"
    wc -l "${CB_DIR}/${ANC}.outcomes.tsv"
done
```

If manually preparing the manifests, create the sentinel expected by Step 45:

```bash
touch "${CB_DIR}/manifest.done"
```

### Safe pilot

```bash
TEST=1 bash 45_submit_colocboost.sh
```

Submits only array index `[1]` of the first ancestry with the true
`N_SHARDS`, so the pilot processes approximately `SHARD_SIZE` regions (10 by
default).

### Full submission

```bash
TEST=0 \
ANCESTRIES="EAS EUR" \
SHARD_SIZE=10 \
MERGE_GAP=100000 \
bash 45_submit_colocboost.sh
```

Examples:

```bash
ANCESTRIES="EAS" \
MODALITIES="expression splicing" \
bash 45_submit_colocboost.sh
```

```bash
FORCE_PREP=1 bash 45_submit_colocboost.sh
```

```bash
FORCE_RUN=1 bash 45_submit_colocboost.sh
```

### Step-43/44 flags not exposed through the Step-45 environment

`43_prepare_colocboost.py` additionally has:

```text
--max-outcomes     default 500
```

`44_colocboost.R` additionally has:

```text
--plink2           default plink2
--min-variants     default 50
--min-outcomes     default 2
--M                default 500
--tmp-dir          default R tempdir()
```

These currently have no environment-variable mapping in
`45_submit_colocboost.sh`.

---

# 7. Steps 46–47 — train isoTWAS/TWAS weights

Weight sets are:

```text
EAS
EUR
pooled
```

The pooled set stacks EAS and EUR ancestry blocks.

**Required environment:** the driver needs only stdlib `python3` (gene counts
and the collapsed-replicate audit). The workers (`47a`) load `module load
plink samtools` (plink2 for genotype dosages, tabix for BED access) and run R
through `$RSCRIPT`; `46_isotwas_train.R` sets `.libPaths()` internally.

Before submitting, the driver runs `check_collapsed_inputs.py` (Section 0.2):
weight training fits sample-level expression/isoform data, so training on a
non-collapsed sample set would silently produce weights inconsistent with
module 07. Override with `SKIP_COLLAPSE_CHECK=1`.

## Environment controls for `47_submit_isotwas.sh`

| Variable | Default | Purpose |
|---|---|---|
| `CONFIG` | required | Pipeline config |
| `SCRIPTS_DIR` | required | Module directory |
| `OUTPUT_BASE` | `dirname $CONFIG` | Root |
| `QTL_DIR` | `$OUTPUT_BASE/qtl_inputs` | QTL data |
| `RESULTS_DIR` | `$OUTPUT_BASE/qtl_results` | Results |
| `LOG_DIR` | `$OUTPUT_BASE/logs` | Logs |
| `ISOTWAS_DIR` | `$RESULTS_DIR/isotwas` | Weight/output directory |
| `WEIGHT_SETS` | `EAS EUR pooled` | Weight sets |
| `SHARD_SIZE` | `25` | Target genes per worker |
| `QUEUE` | `medium` | LSF queue |
| `WALLTIME` | `04:00` | Worker walltime |
| `THREADS` | `2` | CPUs |
| `R2_MIN` | `0.01` | CV R² retention threshold |
| `TEST` | `0` | `1` = submit only shard `[1]` of the first weight set (safe pilot) |
| `FORCE_RUN` | `0` | Retrain existing gene models |
| `SKIP_COLLAPSE_CHECK` | `0` | `1` = skip the collapsed-replicate audit |
| `RSCRIPT` | repository R wrapper | R launcher |
| `PLINK2` | `plink2` | Defined by submitter but currently **not passed to Step 46** |

### Safe pilot

```bash
TEST=1 bash 47_submit_isotwas.sh
```

Submits only array index `[1]` of the first weight set with the true
`N_SHARDS`, so the pilot trains approximately `SHARD_SIZE` genes (25 by
default) — not the entire gene universe as in earlier revisions.

### Full training

```bash
TEST=0 \
WEIGHT_SETS="EAS EUR pooled" \
SHARD_SIZE=25 \
R2_MIN=0.01 \
bash 47_submit_isotwas.sh
```

Individual set:

```bash
WEIGHT_SETS="EAS" \
SHARD_SIZE=25 \
bash 47_submit_isotwas.sh
```

Rerun completed models:

```bash
FORCE_RUN=1 \
WEIGHT_SETS="EAS" \
bash 47_submit_isotwas.sh
```

### Step-46 options not exposed through the Step-47 environment

`46_isotwas_train.R` additionally supports:

| Flag | Default |
|---|---:|
| `--cis-window` | `1e6` |
| `--alpha` | `0.5` |
| `--nfolds` | `5` |
| `--min-variants` | `10` |
| `--max-isoforms` | `50` |
| `--gene-list` | none |
| `--plink2` | `plink2` |
| `--tmp-dir` | R `tempdir()` |

The current submitter does not expose environment variables for these
controls.

If these parameters need to be configurable during production, the preferred
solution is to add explicit environment-to-CLI mappings to
`47_submit_isotwas.sh` rather than manually editing `46_isotwas_train.R`.

### Completion check

```bash
for WS in EAS EUR pooled; do
    echo "=== $WS ==="
    find "${ISOTWAS_DIR}/weights/${WS}/genes" \
        -maxdepth 1 -name '*.wgt.RDat' | wc -l

    ls "${ISOTWAS_DIR}/diagnostics/isotwas_${WS}".shard-*.diagnostics.tsv \
        2>/dev/null | wc -l
done
```

---

# 8. Step 48 — FUSION TWAS association testing

Step 48 performs four operations:

```text
1. Build ancestry-specific chromosome-level 1KG PLINK1 LD references.
2. Convert harmonized GWAS to FUSION summary-statistic format.
3. Merge isoTWAS shard .pos files.
4. Submit one FUSION job per weight-set × GWAS trait.
```

**Required environment:** the driver needs `plink2` on the login node for the
LD-reference build — it attempts `module load plink` automatically when
`plink2` is not already on `PATH` (override with `PLINK2=/path/to/plink2`).
The workers (`48a`) run R through `$RSCRIPT`, invoking FUSION via the repo
wrapper `48b_fusion_assoc.R`, which sets `.libPaths()` **inside R** before
sourcing `FUSION.assoc_test.R` (the external FUSION script has no in-R
library setup, and bash-level `R_LIBS_*` is not reliably propagated here).
FUSION's `optparse` flags are unaffected by the wrapper.

## Environment controls

| Variable | Default | Purpose |
|---|---|---|
| `CONFIG` | required | Pipeline config |
| `SCRIPTS_DIR` | required | Module directory |
| `OUTPUT_BASE` | `dirname $CONFIG` | Root |
| `QTL_DIR` | `$OUTPUT_BASE/qtl_inputs` | Defined, but not materially used by Step 48 |
| `RESULTS_DIR` | `$OUTPUT_BASE/qtl_results` | Results |
| `LOG_DIR` | `$OUTPUT_BASE/logs` | Logs |
| `COLOC_DIR` | `$RESULTS_DIR/coloc` | Contains 1KG keep files |
| `ISOTWAS_DIR` | `$RESULTS_DIR/isotwas` | Weights and TWAS outputs |
| `GWAS_DIR` | `$OUTPUT_BASE/gwas` | Harmonized GWAS |
| `KG_PGEN` | hard-coded 1KG default | Source 1KG pgen |
| `FUSION_DIR` | hard-coded FUSION default | FUSION scripts |
| `LDREF_DIR` | `$GWAS_DIR/fusion_ldref` | PLINK1 LD references |
| `WEIGHT_SETS` | `EAS EUR pooled` | Weight sets |
| `QUEUE` | `medium` | LSF queue |
| `WALLTIME` | `04:00` | Per pair |
| `MIN_R2PRED` | `0.7` | FUSION Z-imputation R² threshold |
| `TEST` | `0` | `1` = submit one weight-set × trait pair |
| `FORCE_RUN` | `0` | Rerun existing TWAS output |
| `RSCRIPT` | repository R wrapper | R launcher |
| `PLINK2` | `plink2` | Used to construct LD reference files |

### Safe pilot

```bash
TEST=1 \
WEIGHT_SETS="EAS EUR pooled" \
MIN_R2PRED=0.7 \
bash 48_fusion_twas.sh
```

`TEST=1` here genuinely limits submission to one weight-set × trait pair.

Monitor:

```bash
bjobs -J 'fusion_*'
```

Check output:

```bash
find "${ISOTWAS_DIR}/fusion" -name '*.twas.tsv' -ls
```

### Full run

```bash
TEST=0 \
WEIGHT_SETS="EAS EUR pooled" \
MIN_R2PRED=0.7 \
bash 48_fusion_twas.sh
```

Subset:

```bash
WEIGHT_SETS="EAS" \
bash 48_fusion_twas.sh
```

Force rerun:

```bash
FORCE_RUN=1 \
WEIGHT_SETS="EAS" \
bash 48_fusion_twas.sh
```

Expected LD references:

```text
$LDREF_DIR/EAS.1kg.chr1.bed
...
$LDREF_DIR/EAS.1kg.chr22.bed

$LDREF_DIR/EUR.1kg.chr1.bed
...
$LDREF_DIR/EUR.1kg.chr22.bed
```

Expected association output:

```text
$ISOTWAS_DIR/fusion/{WEIGHT_SET}/{WEIGHT_SET}_{trait_id}.twas.tsv
```

---

# 9. Step 49 — aggregate coloc/TWAS results

Step 49 aggregates:

```text
SuSiE-coloc
colocBoost
FUSION TWAS
gene-level ACAT
cell-type annotation
```

There is no shell wrapper; analysis controls are direct command-line flags.

**Required environment (manual):**

```bash
conda activate tensorqtl     # pandas, numpy, scipy
```

Default run:

```bash
python3 "${SCRIPTS_DIR}/49_aggregate_coloc_twas.py" \
    --results-dir "$RESULTS_DIR" \
    --qtl-dir "$QTL_DIR"
```

## Step 49 flags

| Flag | Default | Purpose |
|---|---:|---|
| `--results-dir` | required | `$RESULTS_DIR` |
| `--qtl-dir` | required | `$QTL_DIR` |
| `--ancestries` | `EAS EUR` | Cell-type annotation strata |
| `--pp-h4` | `0.7` | Coloc call threshold |
| `--twas-q` | `0.05` | Gene-level TWAS threshold for annotation |
| `--skip-celltype` | false | Skip cell-type annotation |

Explicit thresholds:

```bash
python3 "${SCRIPTS_DIR}/49_aggregate_coloc_twas.py" \
    --results-dir "$RESULTS_DIR" \
    --qtl-dir "$QTL_DIR" \
    --ancestries EAS EUR \
    --pp-h4 0.7 \
    --twas-q 0.05
```

Skip cell-type annotation:

```bash
python3 "${SCRIPTS_DIR}/49_aggregate_coloc_twas.py" \
    --results-dir "$RESULTS_DIR" \
    --qtl-dir "$QTL_DIR" \
    --skip-celltype
```

Cell-type annotation additionally expects, for each ancestry:

```text
$QTL_DIR/{ANC}_deconvolution_harmonized.tsv
$QTL_DIR/{ANC}_expression.bed.gz
$QTL_DIR/{ANC}_covariates_expression.tsv
```

If those are missing, that ancestry is skipped.

Expected aggregated outputs:

```text
$COLOC_DIR/aggregated/
├── coloc_results.tsv.gz
├── coloc_best.tsv.gz
├── colocboost_clusters.tsv.gz
├── coloc_diagnostics_summary.tsv
├── twas_results.tsv.gz
├── twas_gene_results.tsv.gz
└── gene_celltype_annotation.tsv
```

Check:

```bash
ls -lh "${COLOC_DIR}/aggregated/"
```

---

# 10. Module 07 compatibility

Module 07 (`07_gxe_mapping`) consumes module 06 in two ways. Both are now
covered by explicit contracts.

## 10.1 Aggregated result tables (tier-2 prioritized loci)

`54_tier2_stratified.py` builds its tier-2 prioritized phenotype set from:

| File (under `$RESULTS_DIR`) | Columns module 07 uses |
|---|---|
| `coloc/aggregated/coloc_results.tsv.gz` | `modality`, `phenotype_id`, `coloc_call` (fallback: `PP.H4.abf >= --pp-h4`) |
| `coloc/aggregated/twas_gene_results.tsv.gz` | `acat_q`, `GENE` |

Contract notes:

- The directory is `coloc/aggregated/` (with an "s"), matching module 05's
  `finemap/aggregated/` convention. Module 07's resolver now checks
  `coloc/aggregated/` first, then its older `coloc/aggregate/` and `coloc/`
  fallbacks — do not rename the directory.
- `coloc_call` is written by R as `TRUE`/`FALSE`; pandas reads these as
  booleans automatically, so module 07's `astype(bool)` filter is correct.
- The TWAS gene column is `GENE` (uppercase, FUSION `.pos` convention);
  module 07's resolver now accepts it alongside `gene`/`gene_id`/`id`.
- `twas_gene_results.tsv.gz` has one row per `weight_set` × `trait_id` ×
  `GENE`; module 07 takes the union of genes with `acat_q <= --twas_q`
  across all rows.

## 10.2 Collapsed-replicate sample set

Module 07 subsets every module-05 product to the retained RNA run per
individual using
`replicate_collapsed/reports/ancestry_map_collapsed.tsv` as the authority,
and fails closed when that contract is violated. Module 06 reads the same
`qtl_inputs` directly, so it must be held to the same contract — otherwise
coloc LD, nominal scans, and isoTWAS weights would be computed on a different
(larger or unaveraged) sample set than module 07's GxE scans.

`check_collapsed_inputs.py` enforces this for module 06. Per ancestry it
verifies:

1. the collapsed ancestry map exists with `sample_id` /
   `assigned_ancestry` / `cohort` columns and no duplicate retained runs;
2. exactly one retained run per `array_id` (joining the map to module-05
   `{ANC}_metadata.tsv` on `rnaseq_id`, the same join module 07 performs);
3. every sample column in `{ANC}_*.bed.gz` and `{ANC}_covariates_*.tsv`, and
   every IID in `{ANC}_qtl.psam`, is a retained `array_id` — extras are
   **errors** (module 07 excludes those samples);
4. no duplicated sample columns (a signature of uncollapsed technical
   replicates);
5. retained `array_id`s absent from qtl_inputs are **warnings** only (the
   module-05 sample intersection legitimately drops individuals).

The audit runs automatically before submission in `39_run_nominal.sh` and
`47_submit_isotwas.sh`, and can be run any time:

```bash
python3 "${SCRIPTS_DIR}/check_collapsed_inputs.py" \
    --qtl-dir "$QTL_DIR" \
    --output-base "$OUTPUT_BASE" \
    --ancestries EAS EUR
```

It is stdlib-only (no conda env). Escape hatch for cohorts without technical
replicates (collapse step intentionally not run): `SKIP_COLLAPSE_CHECK=1`.

**Scope of the audit.** The check enforces sample-set identity with the
collapsed map. The count-exact value collapse (summing per-run counts and
recomputing ratios) happens upstream in module 05 — steps 17–23 must have
been run against the collapse staging tree per the module-05 README
("Technical replicates"). If module 05 was run without the collapse
workflow, same-cohort replicate pairs are row-wise averaged by the
step-23 guard instead; the audit cannot distinguish those values, but it
does catch the cross-protocol individuals and duplicate columns that the
guard leaves behind.

---

# 11. Worker-only environment variables

These are normally set by the submission scripts. They should not need to be
exported globally.

## Step 39 workers

```text
ANC
MOD
CHROMS_LIST
LSB_JOBINDEX
```

## Step 42 / 42a workers

```text
TASKS
N_SHARDS
COLOC_DIR
MIN_VARIANTS
PP_H4
FORCE_RUN
RSCRIPT
SCRIPTS_DIR
LSB_JOBINDEX
```

## Step 45 / 45a workers

```text
REGIONS
OUTCOMES
N_SHARDS
CB_DIR
LD_XQTL_PGEN
LD_GWAS_PGEN
LD_GWAS_KEEP
FORCE_RUN
RSCRIPT
SCRIPTS_DIR
LSB_JOBINDEX
```

## Step 47 / 47a workers

```text
WEIGHT_SET
N_SHARDS
QTL_DIR
ISOTWAS_DIR
R2_MIN
FORCE_RUN
RSCRIPT
SCRIPTS_DIR
LSB_JOBINDEX
```

## Step 48 / 48a workers

```text
WEIGHT_SET
TRAIT_ID
SUMSTATS
LDREF_PREFIX
ISOTWAS_DIR
MIN_R2PRED
FUSION_DIR
RSCRIPT
SCRIPTS_DIR
```

The 42a/45a/47a workers additionally load the seadragon modules `plink` and
`samtools` themselves; the 48a worker relies on `FUSION_DIR` from the
submitter environment (consumed by `48b_fusion_assoc.R`).

---

# 12. Recommended production sequence

Run the module in this order:

```bash
cd "$SCRIPTS_DIR"

# One-time software installation on login node (no conda env needed).
bash 36_install_coloc_env.sh

# Collapsed-replicate audit (stdlib python3; also runs inside 39/47).
python3 check_collapsed_inputs.py \
    --qtl-dir "$QTL_DIR" --output-base "$OUTPUT_BASE" --ancestries EAS EUR

# Login node: network access required; script self-activates the pandas env.
bash 37_fetch_gwas.sh

# Harmonize (manual env: conda activate tensorqtl; module load samtools).
python3 38_harmonize_gwas.py \
    --catalog gwas_catalog.tsv \
    --raw-dir "${GWAS_DIR}/raw" \
    --out-dir "${GWAS_DIR}" \
    --chain "${CHAIN}" \
    --pgen-dir "${PGEN_DIR}" \
    --ancestries EAS EUR

# Genome-wide nominal QTLs.
TEST=1 bash 39_run_nominal.sh
TEST=0 bash 39_run_nominal.sh

# Pairwise coloc: pilot shard, then full arrays.
TEST=1 bash 42_submit_coloc.sh
TEST=0 bash 42_submit_coloc.sh

# colocBoost: pilot shard, then full arrays.
TEST=1 bash 45_submit_colocboost.sh
TEST=0 bash 45_submit_colocboost.sh

# isoTWAS training: pilot shard, then full arrays.
TEST=1 bash 47_submit_isotwas.sh
TEST=0 bash 47_submit_isotwas.sh

# FUSION pilot followed by full production.
TEST=1 bash 48_fusion_twas.sh
TEST=0 bash 48_fusion_twas.sh

# Final aggregation (manual env: conda activate tensorqtl).
python3 49_aggregate_coloc_twas.py \
    --results-dir "$RESULTS_DIR" \
    --qtl-dir "$QTL_DIR"
```

---

# 13. Module tests

The repository provides:

```bash
pytest \
    test_susie_coloc.py \
    test_colocboost.py \
    test_isotwas_train.py \
    test_aggregate.py
```

The tests require the relevant Python/R/tool dependencies, including
`plink2`, `bgzip`, `tabix`, `coloc`, `colocboost`, and `glmnet`. Run them
inside the `tensorqtl` conda env with `module load plink samtools` and the
module-06 R library available to R (the R workers under test set `.libPaths()`
themselves).

---

# 14. Restart / force semantics

The module is largely incremental.

Use the narrowest force flag necessary:

```text
37   FORCE=1
      Re-download selected raw GWAS.

39   FORCE=1
      Rerun nominal chromosome shards.

42   FORCE_TASKS=1
      Rebuild coloc task tables.

42   FORCE_RUN=1
      Delete matching .done markers and rerun pairwise coloc tasks.

45   FORCE_PREP=1
      Rebuild colocBoost region/outcome manifests.

45   FORCE_RUN=1
      Delete matching region .done markers and rerun colocBoost.

47   FORCE_RUN=1
      Remove existing weight output for assigned genes and retrain.

48   FORCE_RUN=1
      Rerun existing weight-set × GWAS TWAS output.
```

Avoid exporting a force flag globally for the entire module. Prefer
command-scoped assignments such as:

```bash
FORCE_RUN=1 WEIGHT_SETS="EAS" bash 47_submit_isotwas.sh
```

rather than:

```bash
export FORCE_RUN=1
```

This reduces the chance of accidentally rerunning a later step.

---

# 15. Parameters currently not exposed by the submitters

Several scientifically meaningful worker parameters cannot presently be
changed through a submitter environment variable.

These include:

```text
39:
  cis window
  covariate-file override
  allow-missing-chromosomes

41/42:
  plink2 path
  GWAS SuSiE L
  xQTL SuSiE maximum L
  temporary directory

43/45:
  maximum outcomes per colocBoost region

44/45:
  plink2 path
  minimum variants
  minimum outcomes
  colocBoost M
  temporary directory

46/47:
  cis window
  elastic-net alpha
  number of CV folds
  minimum variants
  maximum isoforms
  gene list
  plink2 path
  temporary directory
```

For reproducible production use, these should be exposed by the corresponding
submitter scripts as explicit environment variables rather than changed
inside the worker scripts.
