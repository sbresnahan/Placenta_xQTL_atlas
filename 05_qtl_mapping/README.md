# 05_qtl_mapping — ancestry-stratified cis-xQTL mapping (tensorQTL)

GTEx-conventions cis-xQTL mapping per ancestry stratum: HCP latent factors,
cross-cohort ComBat + INT of all modalities, cohort-only genotype PCA,
genotype×phenotype intersection with a minor-allele-count floor, PC-outlier
exclusion, covariate optimization, modality harmonization, and tensorQTL mapping
with Storey q-values. Current baseline = the GTEx-conventions design (see
[`../docs/round_history.md`](../docs/round_history.md)).

Full stage documentation: [`../MANIFEST.md`](../MANIFEST.md) stage 5. These
scripts share `config.yml` / `config_get.py` with `../03_phenotyping/` (one
scripts directory on seadragon).

## Run order

| # | Script(s) | Purpose |
|---|---|---|
| 19 | `19_hcp_factors.sh`, `19a_picard_sharded.sh`, `19b_fix_missing_metrics.sh` | HCP pipeline: Picard QC (`picard_qc.py`) → pool metrics (`fix_missing_metrics.py`) → pool expression within ancestry (`pool_expression_within_ancestry.py`) → ComBat + INT + HCP k=15 (`combat_normalize_hcp.R`, `peer_factors.py`) |
| 20 | `20_combat_modalities.sh` | Pool + ComBat + INT, 7 non-expression modalities (`pool_modalities_within_ancestry.py`, `combat_normalize_modalities.R`); isoforms log2(x+1); splicing/IR via stage-17 pre-pooled BEDs |
| 21 | `21_install_tensorqtl.sh` | One-time setup: tensorqtl conda env + R `qvalue` + Storey-bridge smoke test |
| 22 | `22_genotype_pca.sh` | Cohort-only genotype PCA: LD-prune → `plink2 --pca 20 exact` → `genotype_pca_format.py` → `{ANC}_genotype_pcs.tsv` + scree (`PCA_scree.R`) |
| 23 | `23_prepare_intersection.py` | Genotype×phenotype intersection; pgen → intersection samples with MAC ≥ 5; IDs rnaseq_id → array_id |
| 24 | `24_outlier_exclusion.py` | First-5-PC selection; 6-SD PC outliers removed from all inputs (**edits in place — rerun 23 first**) |
| 25 | `25_build_covariates.py` | Covariates: PC1–5 + HCP_1–15 + sex + GA + cell types; \|r\| > 0.9 pruning; cap ≤ 25 (**+ manual maternal-fraction removal afterward — see MANIFEST landmine 4**) |
| 26 | `26_harmonize_modalities.py` | Harmonize 7 non-expression BEDs to the final array_id sample set |
| 30 | `30_combine_modalities.py` | Combined cross-modality BED (`{modality}__{id}` namespacing; cross-modality gene groups) |
| 27 | `27_run_tensorqtl.sh` + `27_run_tensorqtl.py` | tensorQTL `cis.map_cis` per ancestry × modality (grouped when `phenotype_groups.txt` exists); `--independent` stepwise conditional mode; Storey q-values via `compute_qvalues.R` bridge (`QVALUE_RSCRIPT`), `--qvalue-method bh` escape hatch |
| 28 | `28_submit_modalities.sh` | Submission driver: one LSF job per ancestry × modality (`TEST=1` pilot) |
| 29 | `29_make_top_tables.py` | Rebuild sorted `*_cisqtl_top.tsv` from parquets |

Diagnostics/utilities: `hcp_diagnostic.R`, `hcp_diagnostic2.R`,
`check_hcp_chunks.sh`, `test_hcp.py`, `test_combat_modalities.py`.

## Inputs / outputs

**Inputs** (from earlier stages): `{ANC}_pooled.pgen` (01), modality BEDs +
deconvolution proportions (03), `rnaseq_to_array_id_map.csv` (04), cohort
metadata (`rnaseq_id, array_id, ancestry, cohort, sex, GA, ppBMI`).

**Outputs**: `$QTL_DIR` (intersected pgens, genotype PCs, covariates, harmonized
BEDs) and `$RESULTS_DIR` (`{ANC}_{modality}[_ungrouped]_cisqtl.parquet` +
`_top.tsv` per layer — grouped / ungrouped / combined — plus `_independent`
stepwise outputs).

## Environment

`tensorqtl` conda env (python-only: tensorqtl 1.0.10, pandas, pyarrow,
genotypeio) + R 4.3.1 singularity image with `qvalue` for the Storey bridge.
Do **not** attempt rpy2/conda-R on seadragon (GLIBCXX conflict — MANIFEST
landmine 1).
