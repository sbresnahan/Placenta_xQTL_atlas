# 07_gxe_mapping — ancestry-first SNP × exposure interaction mapping

Module 07 maps placental cis-xQTL genotype × exposure effects with

```text
Y = b_G G + b_E E + b_GxE (G×E) + Z_cov γ + ε
```

The primary discovery design is **ancestry-specific**, not pooled-intersection:

1. build full autosomal GxE inputs separately for each ancestry;
2. run the adaptive cis-permutation scan on every ancestry's full feature set;
3. merge chromosomes and control FDR within ancestry × modality × exposure;
4. synthesize results across ancestries after discovery;
5. run the existing Aim-1-prioritized ancestry-stratified/IVW tier-2 analysis;
6. run sensitivity/T/NT follow-ups when their inputs are available.

This prevents splicing/IR and other ancestry-specific features from being discarded
before testing. Splicing and intron-retention ancestry-local suffixes are
canonicalized only so results can be matched after scanning.

## Metadata and replicate provenance

Stage 0 reads Module-05 `{ANC}_metadata.tsv` files and uses
`replicate_collapsed/reports/ancestry_map_collapsed.tsv` as the retained-run
authority. Individual-level exposures/demographics therefore follow the same
representative RNA run used by Module 05. `{ANC}_metadata.tsv` is retained in the
Module-07 input directory for provenance/QC only; it is not used to re-estimate HCPs.

## Covariates

Module 07 reuses Module-05's finalized, optimized per-modality covariates
`${QTL_DIR}/{ANC}_covariates_{MOD}.tsv`, including the HCP count selected in Module 5
and Module-05's fixed-covariate/correlation-pruning decisions. Stage 1 only subsets
and reorders those covariates to the exact samples in each Module-07 ancestry BED;
it does not re-estimate HCPs or re-prune covariates. A covariate row with the same
name as the active exposure (for example `GA`) is removed by the scanner because
the exposure enters the interaction model as a main effect.

Continuous exposures are range-filtered/transformed and then standardized once
across all retained ancestries so ancestry-specific interaction coefficients remain
on the same exposure-SD scale. Phenotypes are standardized within ancestry.


## Variant-frequency filters for interaction scans

Tier 1 and tier 2 retain the existing overall in-sample MAF prefilter
(`MAF_THRESHOLD`, default `0.01`) and additionally apply tensorQTL-style
interaction filtering with `MAF_THRESHOLD_INTERACTION` (default `0.05`).
For a single continuous exposure, samples are stable-sorted by the exposure
and split into lower and upper halves; a variant must have MAF >= 0.05 in
**both** halves to enter the interaction model. This mirrors tensorQTL's
`maf_threshold_interaction` behavior and prevents sparse minor alleles
concentrated at one end of the exposure distribution from driving GxE tests.

The production default can be overridden when submitting Module 07, e.g.
`MAF_THRESHOLD_INTERACTION=0.05`, but the 0.05 default is intentional and
should be relaxed only for a documented sensitivity analysis.

## Tier-1 cross-ancestry synthesis

`54_tier1_meta.py` outer-joins ancestry-specific discovery results by canonical
`phenotype_id`:

- ancestry-only phenotypes are retained;
- shared phenotypes receive an ACAT combination of their ancestry-specific
  `pval_beta` values;
- exact fixed-effect IVW beta/se, Cochran Q and I² are reported only when all
  contributing ancestries selected the same lead `variant_id`.

Different ancestry-specific lead variants are explicitly labeled
`shared_discordant_lead`; Module 07 does **not** claim they are the same
feature–variant test. The ancestry-specific scans remain the primary discovery
results.

## Stages

| Stage | Action |
|---|---|
| 0 | Build ancestry-specific BEDs, manifests, metadata provenance, exposures |
| 1 | Stage/validate Module-05 optimized covariates per ancestry × modality |
| 2 | Tier-1 adaptive cis-permutation scans per ancestry × modality × exposure × chromosome |
| 3 | Merge chromosomes and compute Storey q-values within ancestry |
| 4 | Cross-ancestry tier-1 synthesis |
| 5 | Aim-1 prioritized tier-2 ancestry-stratified scan + IVW |
| 6 | Sensitivity and transmitted/non-transmitted follow-ups |
| 7 | Aggregate atlas tables/status |

See `RUNBOOK.md` for the production commands, preflight checks, expected outputs,
and restart procedures.
