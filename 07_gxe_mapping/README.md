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
representative RNA run used by Module 05. Same-cohort technical replicate rows may
remain in `{ANC}_metadata.tsv` only for HCP QC aggregation.

## Covariates

Covariates are built within ancestry: genotype PCs, sex, cohort dummies, GA, cell
fractions, and ancestry/modality-specific HCP factors. The maternal cell fraction
is excluded; the dominant placental cell type is also dropped as the compositional
reference. A covariate row with the same name as the active exposure (for example
`GA`) is removed by the scanner because the exposure enters the interaction model
as a main effect.

Continuous exposures are range-filtered/transformed and then standardized once
across all retained ancestries so ancestry-specific interaction coefficients remain
on the same exposure-SD scale. Phenotypes are standardized within ancestry.

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
| 0 | Build ancestry-specific BEDs, manifests, metadata, exposures, base covariates |
| 1 | Estimate HCPs and finalize covariates per ancestry × modality |
| 2 | Tier-1 adaptive cis-permutation scans per ancestry × modality × exposure × chromosome |
| 3 | Merge chromosomes and compute Storey q-values within ancestry |
| 4 | Cross-ancestry tier-1 synthesis |
| 5 | Aim-1 prioritized tier-2 ancestry-stratified scan + IVW |
| 6 | Sensitivity and transmitted/non-transmitted follow-ups |
| 7 | Aggregate atlas tables/status |

See `RUNBOOK.md` for the production commands, preflight checks, expected outputs,
and restart procedures.
