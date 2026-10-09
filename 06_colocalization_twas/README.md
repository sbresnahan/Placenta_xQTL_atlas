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
   on both sides (in-sample pgen for xQTL; 1KG superpopulation for GWAS).
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
6. **Aggregation** (49): coloc/colocBoost/TWAS result tables, gene-level
   ACAT combination across isoform models, and primary cell-type annotation
   (max |Spearman rho| between residualized expression and MuSiC
   deconvolution proportions; "unassigned" when |rho| < 0.2).

## Script map

| Script | Purpose |
|---|---|
| `36_install_coloc_env.sh` | Install R packages (coloc, susieR, glmnet, colocboost, plink2R) into the pipeline R library; clone FUSION |
| `37_fetch_gwas.py` / `.sh` | Catalog-driven GWAS fetch (direct / page_scrape / GWAS Catalog GCST / manual). JECS and ProDiGY are manual-placement (see `gwas_catalog.tsv` notes) |
| `38_harmonize_gwas.py` | Column standardization, hg19->GRCh38 liftover (pure-Python chain mapper), allele alignment to the pooled pgen, bgzip+tabix output `{trait_id}.sumstats.tsv.gz` |
| `39_run_nominal.py` / `.sh` | Genome-wide `cis.map_nominal` shards + merge to `{ANC}_{MOD}.nominal.tsv.gz` (+ tabix) |
| `40_prepare_coloc_loci.py` | Per-modality coloc task lists (`loci/{MOD}.tasks.tsv`) + 1KG superpopulation keep files (`loci/{ANC}.1kg.keep`) |
| `41_susie_coloc.R` | Pairwise SuSiE-coloc worker (one shard of a task list) |
| `42_submit_coloc.sh` / `42a_run_coloc_shard.sh` | LSF array driver (one array per modality) |
| `43_prepare_colocboost.py` | Region builder: per-ancestry merged union of fine-mapped loci + outcome manifests (`colocboost/{ANC}.regions.tsv`, `{ANC}.outcomes.tsv`) |
| `44_colocboost.R` | colocBoost worker (one shard of regions; X_ref dosage mode) |
| `45_submit_colocboost.sh` / `45a_run_colocboost_shard.sh` | LSF array driver (one array per ancestry) |
| `46_isotwas_train.R` | isoTWAS/TWAS weight trainer (one shard of genes; FUSION-format per-model `.wgt.RDat` + shard `.pos`) |
| `47_submit_isotwas.sh` / `47a_run_isotwas_shard.sh` | LSF array driver (one array per weight set) |
| `48_fusion_twas.sh` / `48a_run_fusion.sh` | FUSION association testing per weight set x trait (per-chromosome 1KG LD refs built here) |
| `49_aggregate_coloc_twas.py` | Aggregation + gene-level ACAT + cell-type annotation |
| `coloc_common.R` | Shared helpers (tabix slicing, plink2 LD/dosage export, psam N) |
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
    loci/{MOD}.tasks.tsv, {ANC}.1kg.keep
    results/{MOD}/{ANC}_{phenotype}_{trait}.coloc.tsv / .variants.tsv / .done
    diagnostics/{MOD}.shard-*.diagnostics.tsv
    colocboost/{ANC}.regions.tsv, {ANC}.outcomes.tsv
    colocboost/results/{ANC}/{region}.clusters.tsv / .vcp.tsv / .done
    aggregated/coloc_results.tsv.gz, coloc_best.tsv.gz,
               colocboost_clusters.tsv.gz, coloc_diagnostics_summary.tsv,
               gene_celltype_annotation.tsv
  isotwas/
    weights/{WS}/genes/{model}.wgt.RDat, shard-*.pos, {WS}.pos
    fusion/{WS}/{WS}_{trait}.twas.tsv
    diagnostics/isotwas_{WS}.shard-*.diagnostics.tsv
  (aggregated TWAS tables land in coloc/aggregated/: twas_results.tsv.gz,
   twas_gene_results.tsv.gz)
```

## Tests

`pytest test_susie_coloc.py test_colocboost.py test_isotwas_train.py
test_aggregate.py` — fixture-based: shared vs distinct causal variants for
coloc/colocBoost; weight recovery + R^2 gate + FUSION round-trip for isoTWAS;
aggregation smoke test. Fixtures need `plink2`, `bgzip`, `tabix`, and R with
`coloc`, `colocboost`, `glmnet` (sandbox: run 36 or use the conda R).
