# Module 06 — cell-type attribution diagnostics (stage 50)

Stage 50 is a **QC/power submodule** for the genotype × cell-proportion
attribution produced by stage 49. It does not change the primary cell-type
assignments.

The key premise is fixed: **EUR contains 136 samples in total.** Therefore an
EUR model N of ~136 is expected and is *not* an unexplained sample-intersection
loss. Sample-overlap tables are retained only to detect any additional loss
below the cohort size. The primary question is whether the much lower EUR
cell-type assignment rate is the expected consequence of testing interaction
terms at N=136, or whether genotype/cell-composition support, effect
heterogeneity, or model calibration also contribute.

## What it tests

1. **EAS → EUR-sized repeated downsampling (primary diagnostic).** For every EAS
   ancestry × modality G×cell test, stage 50 repeatedly samples a common set of
   EAS analysis-eligible subjects down to the EUR cohort size and refits the
   exact classical model `Y ~ G + C + G×C + Z`. BH FDR is recomputed within
   modality for every replicate. By default, the target N is read from
   `EUR_qtl.psam` (expected to be 136) and 100 replicates are run. Outputs report
   the resulting nominal-positive rate, FDR-positive rate, and **feature-level
   cell-type assignment rate**. This directly tests how much of the EAS/EUR
   contrast is attributable to sample size alone.

2. **Genotype × cell-fraction information.** For every tested xQTL feature ×
   cell type, the diagnostic reports MAF, genotype counts, minor-allele carrier
   counts in low/middle/high cell-fraction tertiles, genotype–cell correlation,
   cell-fraction SD/IQR, and the variance/information in `G×C` after projecting
   it on `G`, `C`, and the Module-05 covariates. This asks whether EUR has less
   interaction information even after acknowledging its smaller N.

3. **Matched EAS/EUR feature comparison.** Tests present in both ancestries are
   matched on modality + phenotype + cell type. Because Module 05 can select a
   different lead xQTL in each ancestry, both lead variants are retained and a
   `same_variant` flag is written. Direct beta concordance is summarized only
   for same-variant pairs.

4. **HC3 robust standard errors.** Every valid model is refit with HC3
   heteroskedasticity-robust SEs unless disabled. BH q-values are recomputed
   within ancestry × modality. If the large EAS signal excess disappears under
   HC3, the ancestry contrast is not simply a power issue.

5. **Residual-permutation calibration.** By default, 100 reduced-model residual
   permutations are run for the 20 strongest and 40 random tests in each
   ancestry × modality stratum. Residuals are shuffled within cohort when a
   cohort label is available. This provides an empirical calibration check,
   especially for the unusually strong EAS signal distribution.

6. **Sample overlap QC.** Genotype, phenotype, covariate, and deconvolution IDs
   are still compared, but only as a sanity check. For EUR, the expected upper
   bound is the known cohort N=136.

## Run

After stage 49 is complete:

```bash
export CONFIG=/path/to/config.yml
export SCRIPTS_DIR=/path/to/Placenta_xQTL_atlas/06_colocalization_twas
bash 50_submit_celltype_diagnostics.sh
```

The default run performs 100 EAS→EUR-sized refits per modality, HC3 fits for all
valid tests, and the selected-test permutation diagnostic. Useful overrides:

```bash
# Faster first pass
DOWNSAMPLE_REPS=20 PERMUTATIONS=0 THREADS=8 \
  bash 50_submit_celltype_diagnostics.sh

# Full 100-replicate sample-size diagnostic (default)
DOWNSAMPLE_REPS=100 THREADS=12 bash 50_submit_celltype_diagnostics.sh

# Explicitly force N=136 rather than deriving it from EUR_qtl.psam
DOWNSAMPLE_N=136 DOWNSAMPLE_REPS=100 THREADS=12 \
  bash 50_submit_celltype_diagnostics.sh

# More precise permutation calibration
PERMUTATIONS=500 PERMUTE_TOP=20 PERMUTE_RANDOM=40 THREADS=12 \
  bash 50_submit_celltype_diagnostics.sh

# Rebuild an existing bundle
FORCE=1 bash 50_submit_celltype_diagnostics.sh
```

A smoke test is available:

```bash
TEST=1 bash 50_submit_celltype_diagnostics.sh
```

## Outputs

The default directory is:

```text
$RESULTS_DIR/coloc/aggregated/celltype_diagnostics/
```

| File | Purpose |
|---|---|
| `gxc_downsample_replicates.tsv` | **Primary output:** per-replicate EAS→EUR-sized nominal/FDR/feature-assignment rates |
| `gxc_downsample_tests.tsv.gz` | per-test significance/effect-direction retention across downsample replicates |
| `gxc_downsample_features.tsv.gz` | per-feature probability of receiving a cell-type assignment at EUR-sized N |
| `gxc_matched_ancestry.tsv.gz` | matched EAS/EUR modality + phenotype + cell-type tests with ancestry-specific variants/support |
| `gxc_matched_ancestry_summary.tsv` | compact matched-feature ancestry comparison |
| `gxc_support.tsv.gz` | per-test MAF, cell-range, genotype-by-cell support, interaction information, leverage/conditioning |
| `gxc_hc3.tsv.gz` | classical vs HC3 SE/P/FDR comparison |
| `gxc_permutation.tsv.gz` | empirical P values for selected top/random tests |
| `celltype_distribution.tsv` | ancestry-specific raw/transformed cell-fraction distributions |
| `sample_overlap_summary.tsv` | sample-overlap sanity check; not an explanation for EUR N=136 |
| `sample_overlap_by_sample.tsv.gz` | exact sample membership flags |
| `diagnostic_summary.tsv` | compact ancestry × modality and ancestry-wide summary |
| `run_metadata.json` | command/provenance/settings |
| `celltype_diagnostics_bundle.zip` | all reportable outputs bundled for review |

## Interpretation order

Start with `gxc_downsample_replicates.tsv`. If EAS assignment rates collapse
close to the observed EUR rates when EAS is reduced to N=136, then sample size
alone explains most of the ancestry difference. If downsampled EAS remains far
above EUR, examine `gxc_matched_ancestry_summary.tsv` and `gxc_support.tsv.gz`
for lower EUR cell-fraction variance, sparse minor-allele carriers across the
cell-composition range, or lower residualized G×C information. Then use
`gxc_hc3.tsv.gz` and `gxc_permutation.tsv.gz` to determine whether the EAS excess
is robust to heteroskedasticity and empirical null calibration.

**Upload `celltype_diagnostics_bundle.zip` for summarization.**
