# 06_colocalization_twas — GWAS colocalization and isoTWAS (Objective 1.6)

> **Operator runbook:** [RUNBOOK.md](RUNBOOK.md) — step-by-step commands,
> required conda envs / seadragon modules per step, safe pilots, the
> collapsed-replicate input contract shared with module 07, and the
> module-07 output contract.

Integrates the placental xQTL atlas (module 05) with external GWAS of birth
weight, gestational duration, glycemic traits, and childhood adiposity:

1. **GWAS curation** (`gwas_catalog.tsv`, 37, 38): fetch and harmonize GWAS
   summary statistics to GRCh38 and the atlas variant-ID convention
   (`chr:pos:ref:alt`).
2. **Genome-wide nominal xQTL stats** (39): `tensorQTL cis.map_nominal` per
   ancestry x modality x chromosome, merged + tabix-indexed.
3. **Pairwise SuSiE-coloc** (40-42): per fine-mapped locus x ancestry x
   trait, `coloc::runsusie` + `coloc::coloc.susie` with ancestry-matched LD
   on both sides (in-sample pgen for xQTL; prebuilt ancestry/chromosome 1KG
   reference for GWAS). GWAS-side LD is cached by exact variant set.
   Colocalization call: PP.H4 >= 0.7. Convergence failures are flagged, not
   rescued (no eCAVIAR fallback).
4. **Multi-trait colocBoost** (43-45): per merged region x ancestry, joint
   colocalization across all overlapping xQTL phenotypes (9 modalities) and
   all ancestry-matched GWAS, using reference genotype matrices (X_ref mode).
5. **isoTWAS** (46-48): multivariate elastic net (`glmnet` mgaussian,
   alpha = 0.5) per gene across its isoforms plus per-gene expression
   models; 5-fold CV, retention at CV R^2 > 0.01; weight sets EAS / EUR /
   pooled (ancestry-stacked, ancestry-indicator residualized); FUSION-format
   per-model `.wgt.RDat`; association via `FUSION.assoc_test.R` with
   ancestry-matched 1KG LD references.
6. **Aggregation + cell-type attribution** (49): coloc/colocBoost/TWAS
   result tables, gene-level ACAT combination across isoform models, and the
   R21 Objective 1.4 genotype-by-cell-proportion test at ancestry-specific
   FDR-significant xQTL leads. For each significant phenotype, the exact
   Module-05 discovery lead variant is fit as `Y ~ G + C + G:C + covariates`;
   the target cell-type
   main-effect covariate is replaced by the explicit `C` term, other optimized
   Module-05 covariates remain, and primary cell type is assigned only when
   the `G:C` term survives BH FDR. Maternal fraction is not a target cell type;
   a >10% maternal-fraction exclusion is reported as a sensitivity analysis.

## Script map

| Script | Purpose |
|---|---|
| `36_install_coloc_env.sh` | Install R packages (coloc, susieR, glmnet, colocboost, plink2R) into the pipeline R library; clone FUSION |
| `37_fetch_gwas.py` / `.sh` | Catalog-driven GWAS fetch (direct / page_scrape / GWAS Catalog GCST / manual). JECS and ProDiGY are manual-placement (see `gwas_catalog.tsv` notes) |
| `38_harmonize_gwas.py` | Column standardization, hg19->GRCh38 liftover (pure-Python chain mapper), allele alignment to the pooled pgen, bgzip+tabix output `{trait_id}.sumstats.tsv.gz` |
| `39_run_nominal.py` / `.sh` | Genome-wide `cis.map_nominal` shards + merge to `{ANC}_{MOD}.nominal.tsv.gz` (+ tabix) |
| `40_prepare_coloc_loci.py` | Per-modality coloc task lists (`loci/{MOD}.tasks.tsv`) + optional legacy 1KG superpopulation keep files (`--skip-kg-keep` disables them; active downstream stages do not consume them) |
| `40b_prepare_1kg_reference.py` + `40b_submit_1kg_ld_reference.sh` / `40c_run_1kg_ld_reference.sh` | One-time validated ancestry keep generation + LSF build of ancestry-specific chromosome PGENs (`ld_reference/1kg/{ANC}/chr{1..22}`); runtime stages never `--keep` the full 1KG PGEN |
| `41_prepare_susie_coloc_inputs.py` | Host-side tabix slicing + in-sample xQTL LD + cached GWAS LD from compact ancestry/chromosome 1KG references |
| `41_susie_coloc.R` | Pure-R pairwise SuSiE-coloc statistical worker; reads prepared summary-statistic/LD files only |
| `42_submit_coloc.sh` / `42a_run_coloc_shard.sh` | LSF array driver (one array per modality); host preparation runs before Singularity R |
| `43_prepare_colocboost.py` | Region builder: per-ancestry merged union of fine-mapped loci + outcome manifests (`colocboost/{ANC}.regions.tsv`, `{ANC}.outcomes.tsv`) |
| `44_prepare_colocboost_inputs.py` | Host-side tabix slicing + plink2 dosage export using the same compact ancestry/chromosome 1KG references |
| `44_colocboost.R` | Pure-R colocBoost statistical worker (X_ref dosage mode) |
| `45_submit_colocboost.sh` / `45a_run_colocboost_shard.sh` | LSF array driver (one array per ancestry); host preparation runs before Singularity R |
| `46_prepare_isotwas_inputs.py` | Host-side BED extraction + plink2 cis-dosage export for one isoTWAS shard |
| `46_isotwas_train.R` | Pure-R isoTWAS/TWAS weight trainer; FUSION-format per-model `.wgt.RDat` + shard `.pos` |
| `47_submit_isotwas.sh` / `47a_run_isotwas_shard.sh` | LSF array driver (one array per weight set); host preparation runs before Singularity R |
| `48_fusion_twas.sh` / `48a_run_fusion.sh` | FUSION association testing per weight set x trait; converts the prebuilt 40b/40c chromosome PGENs to FUSION BED views |
| `49_aggregate_coloc_twas.py` | Aggregation + gene-level ACAT + R21 genotype x cell-proportion attribution at significant xQTL leads |
| `coloc_common.R` | Pure-R shared helpers only; no command execution |
| `gwas_catalog.tsv` | GWAS source catalog (see below) |

## GWAS catalog

`gwas_catalog.tsv` columns include `trait_id`, `access`
(`direct`/`page_scrape`/`gcst`/`manual`), `primary_ancestry` (EUR / EAS /
both), `trait_type` (`quantitative`/`cc`), `prop_cases` (for cc), and
`notes`. Manual traits:

- **JECS** (birth weight + gestational age, EAS): request via the form URL in
  `notes`; place files under `$GWAS_DIR/raw/jecs_*` and rerun 37/38.
- **ProDiGY** (T2D, multi-ancestry): T2D Knowledge Portal download
  (registration); place under `$GWAS_DIR/raw/prodigy_*`.

Known discrepancies vs the aims document (flagged in `notes`): the EGG
SEM-partitioned birth-weight files are EUR-only; the EGG gestational-duration
N differs from the aims table; the "childhood anthropometrics" trait is the
UKB recalled-body-size GWAS (GCST010988, 100% EUR), not multi-ancestry.

## Running

All drivers follow the module-05 conventions: `CONFIG`, `SCRIPTS_DIR`,
`OUTPUT_BASE` env vars; `TEST=1` pilots; skip-if-done guards; reruns are
incremental (`.done` sentinels / existing diagnostics).

```bash
export CONFIG=... SCRIPTS_DIR=$PWD/06_colocalization_twas
bash 36_install_coloc_env.sh
bash 37_fetch_gwas.sh                 # login node (network)
python3 38_harmonize_gwas.py ...      # after raw GWAS are placed
TEST=1 bash 39_run_nominal.sh         # pilot: EAS expression chr21
bash 39_run_nominal.sh
TEST=1 bash 40b_submit_1kg_ld_reference.sh  # validate one compact reference
bash 40b_submit_1kg_ld_reference.sh         # one-time EAS/EUR chr1-22 build
TEST=1 bash 42_submit_coloc.sh
bash 42_submit_coloc.sh
TEST=1 bash 45_submit_colocboost.sh
bash 45_submit_colocboost.sh
TEST=1 bash 47_submit_isotwas.sh
bash 47_submit_isotwas.sh
TEST=1 bash 48_fusion_twas.sh
bash 48_fusion_twas.sh
python3 49_aggregate_coloc_twas.py --results-dir $RESULTS_DIR --qtl-dir $QTL_DIR
```

## Outputs

```
$RESULTS_DIR/
  nominal/{ANC}/{ANC}_{MOD}.nominal.tsv.gz(+ .tbi)   # genome-wide nominal stats
  coloc/
    loci/{MOD}.tasks.tsv, {ANC}.1kg.keep   # legacy compatibility only
    ld_reference/1kg/{ANC}/chr{1..22}.{pgen,pvar,psam} + .done
    ld_cache/1kg/{ANC}/chr*/.../*.ld(.vars)
    results/{MOD}/{ANC}_{phenotype}_{trait}.coloc.tsv / .variants.tsv / .done
    diagnostics/{MOD}.shard-*.diagnostics.tsv
    colocboost/{ANC}.regions.tsv, {ANC}.outcomes.tsv
    colocboost/results/{ANC}/{region}.clusters.tsv / .vcp.tsv / .done
    aggregated/coloc_results.tsv.gz, coloc_best.tsv.gz,
               colocboost_clusters.tsv.gz, coloc_diagnostics_summary.tsv,
               xqtl_celltype_interactions.tsv.gz,
               xqtl_celltype_annotation.tsv.gz,
               gene_celltype_annotation.tsv
  isotwas/
    weights/{WS}/genes/{model}.wgt.RDat, shard-*.pos, {WS}.pos
    fusion/{WS}/{WS}_{trait}.twas.tsv
    diagnostics/isotwas_{WS}.shard-*.diagnostics.tsv
  (aggregated TWAS tables land in coloc/aggregated/: twas_results.tsv.gz,
   twas_gene_results.tsv.gz)
```

## Tests

`pytest test_1kg_ld_reference.py test_susie_coloc.py test_colocboost.py test_isotwas_train.py
test_no_r_cli_downstream.py test_aggregate.py` — fixture-based: shared vs distinct causal variants for
coloc/colocBoost; weight recovery + R^2 gate + FUSION round-trip for isoTWAS;
aggregation smoke test. Fixtures need `plink2`, `bgzip`, `tabix`, and R with
`coloc`, `colocboost`, `glmnet` (sandbox: run 36 or use the conda R).
