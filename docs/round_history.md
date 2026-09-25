# Round history — design evolution of the xQTL mapping

The mapping design went through three rounds (plus one intermediate) for the
EAS + EUR strata, followed by a normalization-schema revision (round 4). Round 4
("devBrain schema") is the version shipped in this repository; round-3 results
stand until the round-4 rerun completes. This document records what changed,
why, and which approaches were retired — so that superseded results are never
mistaken for current ones.

**Update (2026-09-25) — round 4: normalization-schema revision (devBrain xQTL
conventions).** Coauthor review determined the phenotype
normalization/transformation schema was wrong. Round 4 reorders it to match the
devBrain xQTL atlas (Wen et al., Science 2024, 384:eadh0829): **pool (unnorm)
→ QN → INT → ComBat(batch=cohort, last)**, with devBrain detection filters
(TPM > 0.1 in > 25% for expression/isoforms; detected in ≥ 40% for proportion
modalities), connectivity-based expression-outlier removal (signed bicor
network, z < −3; extended to isoforms), no log2/logit pre-transforms (no-ops
under rank-based QN+INT), and a new pre-mapping module (25a) that optimizes the
HCP count per ancestry by maximizing chr1 eGene discovery (k ∈ {0,5,…,30},
HCP re-estimated at each k). Per-cohort QN+INT in stage 3 is retired — stage-3
`output/<modality>.bed.gz` files are now unnormalized. Round-3 results (below)
are superseded once the round-4 rerun completes.

**Update (2026-09-21) — EAS stratum re-mapped after a genotype-pooling fix.**
A post-round-3 audit found that GUSTO EAS/SAS samples had been silently dropped
from the pooled EAS genotypes (a missing ancestry-assignment table). After
rebuilding the pooled VCFs and re-running the mapping chain (scripts 22→29),
the EAS stratum grew from n ≈ 108 to n ≈ 272 and the round-1-style
low-carrier artifact class is fully cured in EAS (minimum MAC now 6). The
design below is unchanged; the "current" numbers in this document reflect the
refreshed EAS stratum.

## Rounds at a glance

| Round | Design | Result |
|---|---|---|
| 1 | MAF ≥ 0.01, PC-AiR/1KG PCs, original covariates, BH q-values | 13 EAS splicing "hits" — **all 2–3-carrier artifacts** (median 3 carriers) |
| 1.5 | MAF ≥ 0.05, original covariates | Null (0 FDR hits) |
| 2 | MAF ≥ 0.05, optimized covariates (24; genetically correlated HCPs/cell types removed after an HCP–genotype diagnostic) | 1 locus: EAS splicing at *LARP4B* (combined q = 0.007; lead 10:125824:C:G, 27 carriers) |
| 3 | **GTEx conventions** (below) | **2 FDR hits**: EAS isoforms at *TSPAN3* (ungrouped q = 0.021, grouped q = 0.022) and EUR RNA editing at *NCOA4* (ungrouped q = 0.029); *LARP4B* collapsed (p 5.8e-7 → 0.91) — treated as retracted |
| 4 | **devBrain schema** (below): pool → QN → INT → ComBat last; devBrain filters + outlier removal; HCP k optimized on chr1 | rerun pending |

## Round-4 design (current shipped version)

- **Normalization order**: pool unnorm BEDs within ancestry → quantile
  normalization (samples) → rank-based INT (features) → ComBat (batch=cohort)
  **last** (Wen et al., Science 2024, 384:eadh0829, §3.4 splicing convention,
  applied to all modalities). Retired: per-cohort QN+INT (stage 3), log2/logit
  pre-transforms, ComBat-before-INT.
- **Feature filters** (devBrain §3.3/§3.4): expression/isoforms TPM > 0.1 in
  > 25% of stratum samples; proportion modalities + stability detected
  (non-NA) in ≥ 40%; zeros are real PSI values (the >50%-zeros rule is
  retired). No-variance and per-batch ≥ 2 non-NA guards retained.
- **Sample outliers**: connectivity outliers (signed biweight midcorrelation
  network, z < −3) removed from expression and isoforms (devBrain §3.3);
  splicing/proportion modalities keep all samples.
- **HCP count**: optimized per ancestry by maximizing chr1 cis-eGenes
  (Storey q ≤ 0.05) over k ∈ {0,5,10,15,20,25,30}, HCP re-estimated at each k
  (devBrain §4.2; module 25a). Script 19 runs with provisional k=15 to feed
  the script-23 intersection; 25a installs the k\* solution before canonical
  covariate assembly (25).
- **ComBat-last consequence**: final phenotypes are approximately but not
  exactly N(0,1) per feature — accepted; QQ calibration is checked in the
  rerun report.
- All other round-3 design elements unchanged (MAC ≥ 5 floor, cohort-only
  genotype PCA + 5 PCs, 6-SD PC outliers, covariate pruning/cap, Storey q,
  grouped/ungrouped + independent mapping).

## Round-3 design

- **Variant floor**: MAF ≥ 0.01 **+ MAC ≥ 5** at the genotype intersection
  (script 23) — blocks the round-1 low-carrier artifact class while keeping the
  variant pool large (effective floor ≈ MAF 0.018–0.022 at current n).
- **Genotype PCs**: cohort-only PCA on LD-pruned variants (plink2
  `--indep-pairwise 200 50 0.2` → `--pca 20 exact`; script 22); first 5 PCs used;
  6-SD outlier exclusion on PC1–5 removed 2 EAS + 5 EUR samples (script 24).
- **Covariates**: 24 per ancestry after correlation pruning (|r| > 0.9), a 25-cap,
  and manual removal of the maternal cell-fraction covariate (script 25 + manual
  edit). Both ancestries: PC1–5, HCP_1–11, sex, GA, 6 cell types (no
  correlation pruning triggered in either stratum). No ppBMI.
- **Multiple testing**: Storey q-values on `pval_beta` (R `qvalue`, default lambda)
  via the `compute_qvalues.R` bridge (script 27).
- **Isoforms**: log2(x + 1) isoform *expression* (the earlier logit transform, which
  is for ratios, was retired).
- **Mapping n**: EAS ≈ 272, EUR ≈ 136 (varies slightly by modality layer).

## Current results (verified)

- **2 FDR ≤ 5% hits of 32 ancestry × modality × layer cells**:
  - **EAS isoforms, *TSPAN3*** (ENSG00000140391), isoform ENST00000267970.
    Lead 15:77800023:G:A, AF = 0.36, 202 carriers, slope = −0.54 (r² ≈ 0.13),
    nominal p = 6.8e-9, q = 0.021 ungrouped / 0.022 grouped, 716 kb from the
    isoform TSS. EAS-specific (EUR isoforms q ≈ 0.91) and isoform-specific
    (EAS expression q ≈ 0.99); combined q ≈ 0.15.
  - **EUR ungrouped RNA editing**, site chr10:46,006,565 inside *NCOA4*
    (ENSG00000266412). Lead 10:45322442:G:A, AF = 0.05, 14 carriers,
    slope = −1.78 (r² ≈ 0.3), nominal p = 1.4e-8, q = 0.0287, 684 kb from the
    site. Gene-level tests do **not** pass (EUR grouped q = 0.17, combined
    q = 0.97); no EAS counterpart.

  Both are candidates requiring visual QC + replication, not established
  discoveries.
- **Storey π1 = 0–4% per cell** (trace signal fraction, concentrated in the
  layers carrying the hits); Storey q-values are very close to BH on this data
  (max |q − BH| = 0.006, combined layer).
- **QQ plots clean** (no inflation/deflation); cross-ancestry rank concordance ≈ 0
  (Spearman ρ = −0.012, combined layer).
- **Stepwise (`cis.map_independent`) found no secondary signals** — the only
  FDR-qualifying group (EAS isoforms, *TSPAN3*) entered the scan and yielded
  its primary lead variant alone; all other cells had no entrant at FDR ≤ 5%
  (verified from the 2026-09-21 script-27 logs; see watch-list item 1).
- **Power floor**: at n ≈ 272/136 only r² ≳ 0.10/0.19 is detectable (vs ≈ 0.06
  at PANTRY n = 445). The residual null is a power statement, not a pipeline
  failure.

## Retired approaches (do not resurrect without a reason)

| Retired | Replaced by | Why |
|---|---|---|
| Mapping-stage PC-AiR on 1KG-projected genotypes (`22_genotype_pcair.sh`, `genotype_pcair.R`) | Cohort-only plink2 PCA (`22_genotype_pca.sh`) | GTEx convention; the 1KG-projected PC-AiR solution made eigenvalue order meaningless and complicated PC selection |
| Logit transform of isoform ratios | log2(x + 1) isoform expression (script 20) | Logit is for ratios; isoform phenotypes are expression values |
| BH-only q-values | Storey q-values (GTEx convention), BH retained as escape hatch (`--qvalue-method bh`) | Convention; on current data the two are very close (π1 ≤ 4% per cell, max |q − BH| = 0.006 on the combined layer) |
| rpy2/conda-R for Storey q-values | File-based singularity-R bridge (`compute_qvalues.R`) | rpy2 hits `GLIBCXX_3.4.30 not found` on seadragon (see MANIFEST landmines) |
| ppBMI covariate | Removed | Study decision |
| Within-stratum variance-based PC selection | First-5-PCs rule | Eigenvalue order is meaningful under cohort-only PCA |

Note: **imputation-stage PC-AiR** (`01_genotype_imputation/PC_AiR.R`) is *not*
retired — it remains the ancestry-assignment method for genotype pooling.

## Open issues / watch list

1. **Stepwise layer verified (2026-09-21 logs)**: the EAS isoforms *TSPAN3*
   group (grouped q = 0.022) correctly entered `cis.map_independent` and
   yielded exactly one association — the primary lead variant; no group has a
   conditionally independent secondary signal. All other cells had no
   FDR ≤ 5% entrant, so their `_independent` outputs are empty by design.
   (The `_independent` tables in the earlier results upload were all empty
   because they predate this run; the current cluster outputs are correct.)
2. **Carrier-floor propagation (EUR only)**: MAC ≥ 5 was applied at the
   intersection *before* outlier exclusion, so a small number of tested pairs
   retain 3–4 carriers (EUR: 2,052 combined / 8,722 grouped / 20,809 ungrouped
   tests). The refreshed EAS stratum is clean (minimum MAC = 6; zero tests
   below 5 carriers). No FDR result is carrier-limited. Post-filter EUR hit
   lists to ma_count ≥ 5 for any deliverable.
3. ***LARP4B* interpretation**: the round-2 hit collapsed under GTEx conventions
   despite its lead variant remaining well-powered (23 carriers). Most parsimonious:
   it tagged residual stratification that round-2 covariate tuning did not absorb.
   Retracted unless it re-emerges at larger n.
4. **Modality-level n differs slightly by layer** within an ancestry
   (harmonization differences); harmless but noted.

## Next steps (priority order)

1. Extend mapping to **AFR/AMR/SAS** strata (pooled genotypes exist; run scripts
   19→20→22→23→24→25→26→30→28 per ancestry; watch small-n strata).
2. Add remaining cohorts when phenotyping/genotyping land (power is the binding
   constraint: doubling n roughly halves detectable r²).
3. Hit QC: plot *TSPAN3* isoform levels and *NCOA4* editing levels by genotype;
   seek replication in larger ancestry-matched placental cohorts or
   meta-analysis.
4. Aim-1 downstream once hits exist: cross-ancestry fine-mapping,
   colocalization/isoTWAS vs the R21 GWAS panel (birth weight, gestational age,
   childhood anthropometrics, youth-onset T2D). Aim-2 GxE scans (SNP × GDM,
   pre-pregnancy BMI, OGTT glucose) are underpowered at current n but can be
   piloted — see [`gxe_power_analysis.html`](gxe_power_analysis.html) for
   full-cohort (N = 2,126) power calculations.
