# 07_gxe_mapping — SNP × exposure (GxE) interaction mapping (Objective 2.1)

Maps genotype × exposure interaction effects on placental cis-xQTL
phenotypes, per the aims' model

```
Y = b_SNP * X + b_E * E + b_int * (X * E) + Z_cov * gamma + eps
```

with `Z_cov` = cell-type proportions, cohort, sex, gestational age, ancestry
PCs, and HCP latent factors. Two tiers:

1. **Tier 1 — genome-wide discovery on pooled multi-ancestry samples**
   (maximal power; >80% power for b_int >= 0.3 at common variants per the
   aims' Fig. 4A). Adaptive phenotype-permutation scan with the GTEx /
   tensorQTL multiple-testing scheme (empirical `pval_perm` + Beta-approx
   `pval_beta` + Storey `qval`).
2. **Tier 2 — ancestry-stratified scans at Aim-1 prioritized loci**
   (phenotypes colocalizing with GWAS at PP.H4 >= 0.7 or with significant
   isoTWAS gene associations from module 06), meta-analyzed across
   ancestries by fixed-effect inverse-variance weighting with Cochran's Q /
   I² heterogeneity stats.

Plus two status-gated follow-ups: SNP × covariate sensitivity refits
(`55`, OFF until maternal age / parity / smoking columns exist) and
transmitted/non-transmitted allele decomposition (`56`, BLOCKED until
maternal pgens + mother-child pairs exist; logic implemented and
fixture-tested).

## Statistical notes

* **Interaction math mirrors tensorQTL 1.0.10 exactly**
  (`tensorqtl.core.calculate_interaction_nominal`): g, E, g:E and the
  phenotype are centered and residualized on the covariates (QR
  projection), the per-variant design is X = [g0, e0, ge0], and
  dof = n − k − 4. `test_gxe_scan.py` verifies b / b_se / t match
  tensorQTL to float32 tolerance on identical inputs.
* **Permutations permute the phenotype** (GTEx convention), re-residualizing
  each permuted phenotype against the fixed covariate projection — the
  genotypes/exposure design stays fixed, so each permutation costs one
  batched matvec (validated against brute-force full refits in the tests).
  Caveat: phenotype permutation is slightly conservative for the
  interaction null when a main G effect exists (the permuted null carries
  neither signal). This is the same scheme tensorQTL's interaction mode
  uses.
* **Adaptive schedule**: permutation blocks of 100 / 400 / 500 / 9000
  (cumulative 100 / 500 / 1000 / 10000); a phenotype stops early when the
  running empirical p̂ > 0.10. Per-phenotype RNG seeded by
  `--seed + <global BED row>` — deterministic and chromosome-shard
  invariant.
* **Pooled genotypes without a pooled pgen**: per chromosome, variant IDs
  (which encode chr:pos:ref:alt) are intersected across the per-ancestry
  pgens and ancestry blocks are stacked; missing hardcalls are imputed to
  the within-ancestry-block variant mean; the MAF filter is pooled.
  Covariate ancestry dummies make the genotype main effect within-ancestry.
* **Pooled phenotypes** are z-scored *within* each ancestry before stacking;
  pooled HCPs are re-estimated on the pooled matrix (k = max of the
  per-ancestry optimized k, default 20).
* **The exposure is a model main effect**, so a covariate row named like the
  exposure (e.g. `GA`) is dropped automatically by the scanner. Continuous
  exposures are z-scored pooled (b_int per exposure SD); binary exposures
  stay 0/1.
* Tier 2 uses BH q-values across the small prioritized set (Storey's pi0 is
  unstable on small/discrete sets); tier 1 uses Storey q-values on
  `pval_beta` via the module-05 `compute_qvalues.R` bridge.

## Exposures

Configured in `gxe_config.tsv` (`enabled` flag + valid range). The pilot
exposure is **GA** (gestational age, already in the cohort metadata).
`gdm` and `ogtt` are pre-configured but disabled — enable them after the
columns are added to the cohort metadata (see docs/data_availability.md).

## Order of operations

```bash
export CONFIG=/path/to/config.yml
export SCRIPTS_DIR=$PWD
export METADATA=/path/to/placenta_QTL_cohort_metadata.txt

# everything (stages 0-6):
bash 53_submit_gxe.sh

# or a pilot (expression x GA, chr21 scan + merge only):
TEST=1 STAGES="2 3" bash 53_submit_gxe.sh
```

Stages (see `53_submit_gxe.sh` header): 0 = pooled inputs (`50`), 1 =
pooled HCPs + finalized covariates (`51`), 2 = tier-1 chromosome arrays
(`53a` → `52 --mode cis-perm`), 3 = merge + Storey q (`52 --mode merge`),
4 = tier 2 (`54`), 5 = sensitivity + T/NT (`55`, `56`), 6 = aggregate
(`57`). Per-chromosome parquets and `.done`-style skip checks make reruns
incremental.

## Outputs (`{qtl_results}/gxe/`)

| Path | Contents |
|---|---|
| `inputs/` | pooled BEDs, covariates, exposures, manifest, pooled metadata |
| `tier1/pooled_{MOD}_{EXP}.chr{C}.gxe_cis.parquet` | per-chromosome scan rows |
| `tier1/pooled_{MOD}_{EXP}.gxe_cis.parquet` / `_top.tsv` | merged, with `qval` |
| `tier2/{MOD}_{EXP}.tier2.tsv` | top meta variant per prioritized phenotype (+ per-ancestry stats, Q/I², BH q) |
| `tier2/{MOD}_{EXP}.tier2_variants.tsv.gz` | all tested variant × phenotype rows |
| `sensitivity/`, `tnt/` | status.tsv (OFF/BLOCKED until inputs exist) + results when enabled |
| `aggregate/` | `gxe_tier1_results.tsv.gz`, `gxe_tier1_significant.tsv`, `gxe_tier2_results.tsv.gz`, `gxe_summary.tsv`, `gxe_summary.png`, `gxe_run_status.tsv` |

Tier-1 parquet columns mirror `cis.map_cis` where applicable
(`num_var`, `variant_id`, `start_distance`, `end_distance`, `af`,
`ma_samples`, `ma_count`, `pval_nominal`, `pval_perm`, `pval_beta`,
`beta_shape1/2`, `true_df`, `qval`) plus the interaction fit at the top
variant (`b_g`, `b_e`, `b_gi` + standard errors, `tstat_gi`, `dof`,
`nperm`). `pval_nominal`/`pval_beta`/`qval` all refer to the **interaction
term**; the reported variant is the top *interaction* variant, which need
not coincide with the top marginal xQTL variant.

## Tests

```bash
cd 07_gxe_mapping && pytest test_gxe_scan.py -v
```

8 fixture tests: tensorQTL equivalence, permutation fast-path vs brute
force, interaction recovery / null separation, adaptive-schedule behavior,
Beta-approx wrapper, IVW meta-analysis, the full 3×3 T/NT table, and an
end-to-end CLI run on tiny two-ancestry pgens.
