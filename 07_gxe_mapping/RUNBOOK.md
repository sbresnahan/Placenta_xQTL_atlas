# Module 07 G×E Mapping Runbook — ancestry-first design

## 1. Purpose and analysis contract

Module 07 performs placental cis genotype × exposure mapping with the primary model

```text
Y = b_G G + b_E E + b_GxE (G×E) + Z_cov γ + ε
```

Primary discovery is performed **within ancestry on the full autosomal feature set**.
No EAS∩EUR phenotype intersection is imposed before scanning. This is particularly
important for splicing and intron retention, whose feature spaces differ substantially
between ancestries and whose upstream Module-03 IDs contain ancestry-local numeric
suffixes.

The workflow is:

```text
Module-05 final qtl_inputs
        │
        ├── EAS full autosomal features ──> EAS G×E scan ──> EAS q-values
        │
        └── EUR full autosomal features ──> EUR G×E scan ──> EUR q-values
                                                       │
                                                       ▼
                                      post-hoc canonical feature matching
                                                       │
                          ┌────────────────────────────┴──────────────────────┐
                          │                                                   │
                 ancestry-only features                           shared features
                                                                      │
                                               ACAT feature evidence + exact-lead IVW
```

Ancestry-specific tier-1 scans are the discovery results. Cross-ancestry synthesis is
performed after those scans; it does not shrink the discovery feature space.

## 2. Production paths

On Seadragon:

```bash
export REPO_ROOT=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/Pantry/phenotyping/scripts
export SCRIPTS_DIR=${REPO_ROOT}/07_gxe_mapping

export OUTPUT_BASE=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY
export CONFIG=${OUTPUT_BASE}/config.yml
export QTL_DIR=${OUTPUT_BASE}/qtl_inputs
export RESULTS_DIR=${OUTPUT_BASE}/qtl_results
export GXE_DIR=${RESULTS_DIR}/gxe
export LOG_DIR=${OUTPUT_BASE}/logs

export COLLAPSED_ANCESTRY_MAP=${OUTPUT_BASE}/replicate_collapsed/reports/ancestry_map_collapsed.tsv
export R_PACKAGE_LIB=/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1
export CONDA_EXE=/risapps/rhel8/miniforge3/24.5.0-0/bin/conda
export CONDA_ENV=tensorqtl

cd "$SCRIPTS_DIR"
```

Default ancestries are `EAS EUR`. Default chromosomes are autosomes 1–22.

## 3. Preflight

### 3.1 Repository/script checks

```bash
set -euo pipefail

for f in \
  50_build_gxe_inputs.py \
  50a_run_build_gxe_inputs.sh \
  51_pooled_hcp.sh \
  52_gxe_scan.py \
  53_submit_gxe.sh \
  53a_run_gxe_scan.sh \
  53b_run_post.sh \
  54_tier1_meta.py \
  54_tier2_stratified.py \
  55_sensitivity_snpxcov.py \
  56_transmitted_nontransmitted.py \
  57_aggregate_gxe.py; do
  test -f "$SCRIPTS_DIR/$f" || { echo "MISSING: $f"; exit 1; }
done

bash -n 50a_run_build_gxe_inputs.sh
bash -n 51_pooled_hcp.sh
bash -n 53_submit_gxe.sh
bash -n 53a_run_gxe_scan.sh
bash -n 53b_run_post.sh

python3 -m py_compile \
  50_build_gxe_inputs.py 52_gxe_scan.py 54_tier1_meta.py \
  54_tier2_stratified.py 55_sensitivity_snpxcov.py \
  56_transmitted_nontransmitted.py 57_aggregate_gxe.py
```

Module 07 launches Module-05 R helpers through `${REPO_ROOT}/bin/Rscript_sif`.
Verify it is executable:

```bash
test -x "${REPO_ROOT}/bin/Rscript_sif" || chmod u+x "${REPO_ROOT}/bin/Rscript_sif"
```

The R package library is prepended **inside the R scripts** by the Module-07 wrappers;
it is not sufficient to rely only on `R_LIBS_USER` in the shell.

### 3.2 Module-05 provenance checks

The ancestry-first G×E inputs must derive from the final Module-05 QTL inputs and the
same retained-run map used by the replicate-collapse workflow.

```bash
test -s "$COLLAPSED_ANCESTRY_MAP"
head -3 "$COLLAPSED_ANCESTRY_MAP"

for ANC in EAS EUR; do
  test -s "$QTL_DIR/${ANC}_metadata.tsv"
  test -s "$QTL_DIR/${ANC}_qtl.pgen"
  test -s "$QTL_DIR/${ANC}_qtl.pvar"
  test -s "$QTL_DIR/${ANC}_qtl.psam"
  test -s "$QTL_DIR/${ANC}_selected_pcs.txt"
done
```

### 3.3 Old pooled outputs

Older Module-07 runs may have produced `inputs/pooled_*` and
`tier1/pooled_*.gxe_cis.parquet`. The redesigned code does not use those files, and
the aggregator only accepts configured ancestry prefixes. They may remain for audit
purposes, but do not use them as inputs to the ancestry-first workflow.

## 4. Stage 0 — build ancestry-specific inputs

Submit:

```bash
STAGES="0" FORCE=1 bash 53_submit_gxe.sh
```

Monitor:

```bash
bjobs -w
ls -lt "$LOG_DIR"/gxe_inputs.*.out "$LOG_DIR"/gxe_inputs.*.err | head
```

Stage 0 performs the following independently for EAS and EUR:

1. reads `${QTL_DIR}/{ANC}_{MOD}.bed.gz`;
2. retains only chr1–22;
3. canonicalizes splicing/IR ancestry-local suffixes;
4. **does not intersect phenotypes across ancestries**;
5. z-scores each phenotype within ancestry;
6. writes `${GXE_DIR}/inputs/{ANC}_{MOD}.bed.gz` + `.tbi`;
7. builds `{ANC}_sample_manifest.tsv` from the retained Module-05 metadata;
8. writes `{ANC}_metadata.tsv` for HCP QC aggregation;
9. builds `{ANC}_covariates_base.tsv`;
10. writes one common `exposures.tsv` whose continuous exposures are scaled across
    all retained ancestries.

Expected output examples:

```text
${GXE_DIR}/inputs/EAS_expression.bed.gz
${GXE_DIR}/inputs/EUR_expression.bed.gz
${GXE_DIR}/inputs/EAS_splicing.bed.gz
${GXE_DIR}/inputs/EUR_splicing.bed.gz
${GXE_DIR}/inputs/EAS_intron_retention.bed.gz
${GXE_DIR}/inputs/EUR_intron_retention.bed.gz
${GXE_DIR}/inputs/EAS_sample_manifest.tsv
${GXE_DIR}/inputs/EUR_sample_manifest.tsv
${GXE_DIR}/inputs/EAS_covariates_base.tsv
${GXE_DIR}/inputs/EUR_covariates_base.tsv
${GXE_DIR}/inputs/exposures.tsv
```

### 4.1 Critical Stage-0 QC

Confirm splicing and IR retain ancestry-specific feature counts rather than collapsing
to an intersection:

```bash
for ANC in EAS EUR; do
  for MOD in expression splicing intron_retention; do
    F="$GXE_DIR/inputs/${ANC}_${MOD}.bed.gz"
    printf '%s\t%s\t' "$ANC" "$MOD"
    printf 'features='; zcat "$F" | awk 'END{print NR-1}'
  done
done
```

A valid run should show nonzero, potentially very different EAS and EUR feature counts.
No log line should report a `cross-ancestry phenotype intersection` during Stage 0.

Confirm BED layout and tabix indexing:

```bash
for F in "$GXE_DIR"/inputs/{EAS,EUR}_*.bed.gz; do
  zcat "$F" | head -1 | cut -f1-5
  tabix -l "$F" | head
  echo
 done
```

Column order must begin:

```text
#chr    start    end    phenotype_id    <sample1> ...
```

### 4.2 Metadata/collapse QC

```bash
for ANC in EAS EUR; do
  echo "=== $ANC ==="
  cut -f1 "$GXE_DIR/inputs/${ANC}_sample_manifest.tsv" | tail -n +2 | sort | uniq -d | head
  wc -l "$GXE_DIR/inputs/${ANC}_sample_manifest.tsv"
  wc -l "$GXE_DIR/inputs/${ANC}_metadata.tsv"
done
```

The manifest should have one row per `array_id`. `{ANC}_metadata.tsv` may have more
RNA-seq rows than the manifest because same-cohort technical replicate rows are kept
for HCP QC averaging.

### 4.3 Cell-fraction QC

Stage 0 should log, separately by ancestry, the dominant cell type dropped as the
compositional reference and the exclusion of the maternal fraction, for example:

```text
EAS cell types: 6 (dropped dominant: Syncytiotrophoblast; excluded maternal fraction: Maternal)
```

The dominant reference may differ by ancestry; that is acceptable because covariates
are fit within ancestry.

## 5. Stage 1 — ancestry/modality-specific HCP factors

Submit:

```bash
STAGES="1" FORCE=1 bash 53_submit_gxe.sh
```

The historical filename `51_pooled_hcp.sh` is retained for compatibility, but the
worker now runs **one ancestry × modality**. Expected outputs include:

```text
inputs/EAS_hcp_expression.tsv
inputs/EUR_hcp_expression.tsv
inputs/EAS_covariates_expression.tsv
inputs/EUR_covariates_expression.tsv
...
```

For each modality, finalized covariates are restricted to that modality's actual
sample columns. This is required because modalities can differ by one or more samples.

QC:

```bash
for ANC in EAS EUR; do
  for MOD in expression splicing intron_retention; do
    C="$GXE_DIR/inputs/${ANC}_covariates_${MOD}.tsv"
    test -s "$C" || echo "MISSING $C"
    printf '%s %s: ' "$ANC" "$MOD"
    awk 'END{print NR-1 " covariates"}' "$C"
  done
done
```

## 6. Pilot before the full scan

A recommended pilot is expression × the first enabled exposure on chr21:

```bash
TEST=1 STAGES="2 3 4" FORCE=1 bash 53_submit_gxe.sh
```

`TEST=1` restricts the chromosome list to `21`, and the worker maps the LSF array index
back to the requested chromosome list. It therefore really scans chr21 rather than
array index 1.

Expected primary outputs:

```text
tier1/EAS_expression_GA.chr21.gxe_cis.parquet
tier1/EUR_expression_GA.chr21.gxe_cis.parquet
tier1/EAS_expression_GA.gxe_cis.parquet
tier1/EUR_expression_GA.gxe_cis.parquet
tier1_meta/expression_GA.meta.tsv.gz
```

Inspect logs before scaling up:

```bash
grep -H -E 'ERROR|Traceback|wrote|analysis samples|qval' \
  "$LOG_DIR"/gxe_t1_* "$LOG_DIR"/gxe_merge_* "$LOG_DIR"/gxe_meta_* 2>/dev/null | tail -100
```

## 7. Stage 2 — full ancestry-specific tier-1 scans

Submit chromosome arrays:

```bash
STAGES="2" bash 53_submit_gxe.sh
```

The job unit is:

```text
ANCESTRY × MODALITY × EXPOSURE × CHROMOSOME
```

For example:

```text
EAS × splicing × GA × chr1
...
EUR × splicing × GA × chr22
```

The scanner reads only that ancestry's pgen. MAF/monomorphic filtering is therefore
performed within ancestry rather than after cross-ancestry variant intersection.

Per-chromosome output:

```text
${GXE_DIR}/tier1/{ANC}_{MOD}_{EXP}.chr{CHR}.gxe_cis.parquet
```

Each phenotype row contains the ancestry-specific top cis interaction variant and
adaptive permutation/Beta-approximation statistics.

## 8. Stage 3 — merge chromosomes and control FDR within ancestry

```bash
STAGES="3" bash 53_submit_gxe.sh
```

Outputs:

```text
${GXE_DIR}/tier1/{ANC}_{MOD}_{EXP}.gxe_cis.parquet
${GXE_DIR}/tier1/{ANC}_{MOD}_{EXP}.gxe_cis_top.tsv
```

Storey q-values are computed on `pval_beta` separately for each
ancestry × modality × exposure family.

QC example:

```bash
python3 - <<'PY'
import glob, pandas as pd, os
for p in sorted(glob.glob(os.environ['GXE_DIR'] + '/tier1/*.gxe_cis.parquet')):
    if '.chr' in p:
        continue
    d = pd.read_parquet(p)
    print(os.path.basename(p), len(d), 'q<=0.05=', int((d.qval <= .05).sum()))
PY
```

## 9. Stage 4 — cross-ancestry tier-1 synthesis

```bash
STAGES="4" bash 53_submit_gxe.sh
```

Output:

```text
${GXE_DIR}/tier1_meta/{MOD}_{EXP}.meta.tsv.gz
```

Interpret the `scope` column carefully:

- `ancestry_only`: phenotype was testable in only one ancestry;
- `shared_concordant_lead`: both ancestries tested the phenotype and selected the
  same lead variant; exact fixed-effect IVW beta/se and Q/I² are reported;
- `shared_discordant_lead`: phenotype was shared but ancestry-specific lead variants
  differed; no exact effect meta-analysis is claimed.

`p_feature_acat` combines the ancestry-specific `pval_beta` values at the **feature
level**, after each ancestry's cis search/permutation correction. `q_feature_acat`
is BH across shared features. This score tests cross-ancestry feature evidence; it
must not be interpreted as proof that the same causal variant operates in both
ancestries.

For `shared_concordant_lead`, inspect:

```text
meta_variant_id
beta_meta
se_meta
p_meta
p_het
i2
```

The ancestry-specific scans remain the primary discovery calls.

## 10. Stage 5 — Aim-1-prioritized tier 2

```bash
STAGES="5" bash 53_submit_gxe.sh
```

Tier 2 retains the existing prioritized-locus strategy (Module-06 coloc / TWAS
phenotypes) but now uses the same ancestry-specific Module-07 BEDs and covariates as
tier 1.

Outputs:

```text
${GXE_DIR}/tier2/{MOD}_{EXP}.tier2.tsv
${GXE_DIR}/tier2/{MOD}_{EXP}.tier2_variants.tsv.gz
```

## 11. Stage 6 — sensitivity and transmitted/non-transmitted analyses

```bash
STAGES="6" bash 53_submit_gxe.sh
```

These are status-gated:

- SNP × covariate sensitivity remains `OFF` until sensitivity covariates are enabled;
- T/NT remains `BLOCKED` until maternal pgens and mother-child pairs are supplied.

When enabled, both follow-ups now operate on ancestry-specific tier-1 significant
hits and ancestry-specific Module-07 BED/covariate files.

## 12. Stage 7 — aggregate

```bash
STAGES="7" bash 53_submit_gxe.sh
```

Outputs:

```text
${GXE_DIR}/aggregate/gxe_tier1_results.tsv.gz
${GXE_DIR}/aggregate/gxe_tier1_significant.tsv
${GXE_DIR}/aggregate/gxe_tier1_meta.tsv.gz
${GXE_DIR}/aggregate/gxe_tier2_results.tsv.gz
${GXE_DIR}/aggregate/gxe_summary.tsv
${GXE_DIR}/aggregate/gxe_run_status.tsv
```

`gxe_tier1_results.tsv.gz` and `gxe_tier1_significant.tsv` contain an explicit
`ancestry` column. An ancestry-specific discovery must not be labeled heterogeneous
simply because the feature or variant was not testable in the other ancestry.

## 13. Full production submission

After the pilot is clean:

```bash
STAGES="0 1 2 3 4 5 6 7" FORCE=0 bash 53_submit_gxe.sh
```

If Stage 0 or Stage 1 was already generated under the old pooled design, rebuild them:

```bash
STAGES="0 1" FORCE=1 bash 53_submit_gxe.sh
```

Then submit scans/merges/synthesis:

```bash
STAGES="2 3 4" bash 53_submit_gxe.sh
```

Tier 2 and aggregation can then be submitted separately:

```bash
STAGES="5 6 7" bash 53_submit_gxe.sh
```

## 14. Restart and FORCE behavior

- Stage-2 chromosome workers skip an existing parquet unless `FORCE=1`.
- Stage-1 HCP workers skip an existing finalized ancestry/modality covariate file
  unless `FORCE=1`.
- Stage-4 synthesis skips an existing `.meta.tsv.gz` unless `FORCE=1`.
- Use `FORCE=1` after any input/covariate/metadata logic change.

Examples:

```bash
# Redo one ancestry/modality/exposure scan manually
ANCESTRY=EAS MODALITY=splicing EXPOSURE=GA MODE=scan \
CHROMS="21" LSB_JOBINDEX=1 FORCE=1 \
bash 53a_run_gxe_scan.sh

# Redo only the ancestry merges through the submitter
STAGES="3" FORCE=1 MODALITIES="splicing" EXPOSURES="GA" bash 53_submit_gxe.sh

# Recompute cross-ancestry synthesis
STAGES="4" FORCE=1 MODALITIES="splicing" EXPOSURES="GA" bash 53_submit_gxe.sh
```

## 15. Validation tests

Static checks:

```bash
python3 -m py_compile \
  50_build_gxe_inputs.py 52_gxe_scan.py 54_tier1_meta.py \
  54_tier2_stratified.py 55_sensitivity_snpxcov.py \
  56_transmitted_nontransmitted.py 57_aggregate_gxe.py

for f in 50a_run_build_gxe_inputs.sh 51_pooled_hcp.sh \
         53_submit_gxe.sh 53a_run_gxe_scan.sh 53b_run_post.sh; do
  bash -n "$f"
done
```

The repository's `test_gxe_scan.py` imports a fixture helper from another module
(`test_susie_coloc`). Run the full pytest suite from the repository checkout where
that helper is available, rather than from an isolated `07_gxe_mapping` archive:

```bash
cd "$REPO_ROOT"
pytest 07_gxe_mapping/test_gxe_scan.py -v
```

## 16. Interpretation rules for the atlas

1. **Ancestry-specific significant, feature absent in other ancestry:** report as an
   ancestry-specific discovery/testability pattern; do not call effect heterogeneity.
2. **Feature present in both, significant in one only:** inspect effect estimates,
   standard errors and power before calling ancestry specificity.
3. **Shared phenotype + same lead variant:** use the IVW effect and heterogeneity
   fields from Stage 4.
4. **Shared phenotype + different leads:** retain both ancestry-specific discoveries;
   Stage 4 deliberately does not force an effect meta-analysis.
5. **Tier 2:** use the prioritized all-variant ancestry-stratified analysis for deeper
   locus-level follow-up at Aim-1 loci.
