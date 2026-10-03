# MANIFEST

File inventory. See the top-level README for the workflow overview and the
per-module READMEs for walkthroughs.

## Top level

| File | Purpose |
|------|---------|
| README.md | Project overview, module map, quick start |
| MANIFEST.md | This file inventory |
| LICENSE | GPL v3 |
| docs/data_availability.md | Cohorts, accessions, ancestry x cohort sample counts, reference data |
| docs/software_environments.md | Tool, environment, and package versions |
| docs/porting.md | Adapting the pipeline to a different cluster |
| docs/gxe_power_analysis.html | Power calculations for SNP x exposure (GxE) scans |

## 01_genotype_imputation/

Per-cohort genotype retrieval, TOPMed imputation, PC-AiR ancestry
assignment, and cross-cohort pooling into per-ancestry plink2 files.

| File | Purpose |
|------|---------|
| pull_RICHS_genotypes.lsf | dbGaP array genotype retrieval (NIEHS_RICHS) |
| pull_Yale_genotypes.lsf | Array genotype retrieval (GUSTO) |
| retrieve_SNUH.lsf | Genotype retrieval (SNUH) |
| pull_NIGMS_exomes.lsf | SRA WXS retrieval (NIGMS) |
| NIGMS_preprocessing.lsf | NIGMS WXS: SRA -> FASTQ -> alignment -> gVCF |
| loop_NIGMS_preprocessing.sh | Submit one NIGMS_preprocessing job per sample |
| merge_NIGMS_gvcfs_setup.lsf | Index NIGMS gVCFs for joint genotyping |
| merge_NIGMS_gvcfs_chr.lsf | Per-chromosome joint genotyping of NIGMS gVCFs |
| liftover_RICHS_genotypes.lsf | Lift array genotypes to GRCh38 |
| BWA_ref.lsf | Build BWA reference index |
| prep_dbSNP.lsf | Prepare dbSNP reference |
| extract_topmed.lsf | Download TOPMed imputation results (DIR/PASSWORD via bsub -env) |
| extract_G3A.lsf | Download G3A imputation results (DIR/PASSWORD via bsub -env) |
| filter_vcf_for_imputation.lsf | Pre-imputation filtering (SNPs only, biallelic, missingness/MAF/HWE) |
| filter_vcf_for_imputation_{NIEHS,NIGMS,SNUH,MALI}.lsf | Per-cohort variants of the above |
| pull_imputed_{RICHS,NIGMS,MALI,MALI_G3A}_genotypes.lsf | Pull TOPMed imputation server outputs |
| prep_QC_{GUSTO,MALI,MALI_G3A,NIEHS_RICHS,NIGMS,SNUH}.lsf | Post-imputation QC per cohort (Rsq >= 0.4 -> plink2 pgen) |
| prep_GUSTO_redo.lsf | GUSTO QC rerun |
| cohort_imputationQC.R | Per-cohort imputation QC summaries |
| mega_imputationQC.R | Pooled imputation QC summaries |
| check_vcf_quality.R | Pre-imputation VCF checks |
| prep_pcAiR.lsf | Build 1000 Genomes reference pgen for PC-AiR |
| PC_AiR.R | GENESIS PC-AiR ancestry assignment against 1KG |
| AFR_check.R | AFR-stratum assignment checks |
| AFR_LD.bed | LD regions for AFR checks |
| mega_merge_and_filter.lsf | Cross-cohort pooling: per-ancestry subsetting, multiallelic splitting, bcftools merge, GRCh38 normalization, pgen conversion, pooled filters |
| mega_merge_and_filter_check.lsf | Pooling QC checks |
| cleanup.sh | Housekeeping |

## 02_rnaseq_quantification/

SRA retrieval, trimming, and Salmon quantification for cohorts processed
outside the phenotyping driver.

| File | Purpose |
|------|---------|
| pull_sra_SNUH.lsf | prefetch + fasterq-dump + gzip for one accession (bsub -env ACC=SRR...) |
| loop_salmon_NIEHS_RICHS.sh | Submit one Salmon job per FASTQ (NIEHS_RICHS) |
| run_salmon_NIEHS_RICHS.lsf | Trim Galore -> salmon quant (single-end, 50 bootstraps) |
| loop_salmon_SNUH.sh | Submit one Salmon job per accession (SNUH) |
| run_salmon_SNUH.lsf | Trim Galore -> salmon quant (paired-end) |

## 03_phenotyping/

Per-cohort RNA-seq processing: alignment, quantification, per-modality
aggregation into unnormalized BED matrices, cross-cohort harmonization,
cell-type deconvolution.

| File | Purpose |
|------|---------|
| config.yml | Shared reference paths + per-cohort blocks (also read by 05_qtl_mapping) |
| config_get.py | config.yml parser; emits shell exports |
| run_pipeline.py | Driver: stages cohort inputs, submits per-cohort jobs with LSF dependencies, monitors |
| 00_run_pipeline.sh | Shell wrapper around run_pipeline.py |
| 001_refprep_pre.sh | Reference prep stage 1: normalized GTF, REDIportal BED, exonic/intronic GTFs, MAJIQ GFF3, edit-site map, txrevise inputs |
| 002_txrevise_array.sh | Reference prep stage 2: txrevise event construction (100-batch job array) |
| 003_refprep_post.sh | Reference prep stage 3: merge per-batch GFF3s; build 6 txrevise Salmon indices |
| 01_star_align.sh | STAR alignment per sample |
| 02_salmon_expression.sh | Salmon quantification (expression index, 20 bootstraps) |
| 03_salmon_alt_tss_polya.sh | Salmon quantification against 6 txrevise indices |
| 04_regtools_junctions.sh | RegTools splice-junction extraction per sample |
| 05_featureCounts.sh | featureCounts exon/intron counts (stability) |
| 06_rna_editing_pileup.sh | Per-sample mpileup at REDIportal editing sites |
| 07_majiq_build.sh | MAJIQ splice-graph build |
| 10_aggregate_expression.sh | Aggregate expression + isoform BEDs (qu_correct_salmon.R QU correction for isoforms) |
| 11_aggregate_alt_tss_polya.sh | Aggregate alt_TSS + alt_polyA BEDs |
| 12_aggregate_splicing.sh | leafCutter clustering and splicing BED |
| 13_aggregate_intron_retention.sh | MAJIQ intron-retention PSI BED |
| 14_aggregate_rna_editing.sh | RNA-editing BED |
| 15_aggregate_stability.sh | Stability (exon/intron ratio) BED |
| 16_index_outputs.sh | Verify, bgzip, and tabix-index modality BEDs |
| 17_harmonize_within_ancestry.sh | Cross-cohort coordinate harmonization for splicing and intron retention (per ancestry) |
| harmonize_within_ancestry.py | Implementation of 17 |
| 18_deconvolution.sh | Deconvolution orchestrator |
| 18a_build_reference.sh | Build MuSiC single-cell reference (build_music_reference.R) |
| 18b_deconvolve_cohort.sh | MuSiC deconvolution per cohort (build_bulk_tpm_matrix.py + run_music_deconvolution.R) |
| 18c_pool_outputs.sh | Pool deconvolution outputs (pool_deconvolution_outputs.py) |
| assemble_bed.py | Shared BED assembly helpers used by 10-15 (modified PANTRY script: reads Salmon quant.sf) |
| normalize_phenotypes.py | Per-cohort normalization helpers |
| prep_PANTRY_gtf.R | GTF normalization (step 001) |
| qu_correct_salmon.R | edgeR catchSalmon quantification-uncertainty correction |
| exonic_intronic_from_gtf.py | Exonic/intronic GTF extraction |
| gtf_to_majiq_gff3.py | GTF to MAJIQ GFF3 conversion |
| rediportal_preprocess.py | REDIportal table to BED |
| leafcutter_cluster_regtools_py3.py | leafCutter clustering (python 3 port) |
| build_bulk_tpm_matrix.py | Bulk TPM matrix for deconvolution |
| build_music_reference.R | MuSiC reference construction |
| run_music_deconvolution.R | MuSiC deconvolution |
| pool_deconvolution_outputs.py | Pool per-cohort cell-type fractions |
| check_cell_proportions.R | Deconvolution QC |
| reference_gene_symbols.txt | Gene symbols for the deconvolution reference |
| combine_modalities.sh | Concatenate modality BEDs into a cross-modality BED |
| regenerate_refflat.sh | Regenerate the 11-field refFlat for Picard |
| write_laddr_config.py | Optional LaDDR config writer (not part of the production flow) |
| RNA_editing/ | Editing-site helpers (map to genes, phenotype prep, level query, matrix) |
| intron_retention/extract_ir_psi.py | MAJIQ PSI extraction |
| txrevise/ | Bundled txrevise event construction (constructEvents.R, prepareAnnotations.R, extractTranscriptTags.py, preprocess_gtf.py) |
| test_harmonize.py | Test suite for 17 |
| test_assemble_bed_isoform_expr.py | Test for assemble_bed.py isoform expression |
| test_prepare_data.sh, test_subset_fastq.py | Pilot-run helpers |

## 04_sample_linkage/

| File | Purpose |
|------|---------|
| build_RNA_to_DNA_map.R | Build rnaseq_to_array_id_map.csv linking PANTRY BED sample IDs to pooled-genotype sample IDs (base R, run interactively) |

## 05_qtl_mapping/

Cross-cohort pooling, normalization, covariates, tensorQTL mapping, SuSHiE
fine-mapping, and annotation enrichment.

| File | Purpose |
|------|---------|
| collapse_replicates.py | Collapse technical-replicate runs to one column per individual |
| 19_hcp_factors.sh | Expression pooling, Picard QC, QN + within-cohort INT, HCP factors |
| 19a_picard_sharded.sh | Sharded Picard CollectRnaSeqMetrics |
| 19b_fix_missing_metrics.sh | Retry failed Picard shards |
| picard_qc.py, fix_missing_metrics.py | Picard helpers (also refFlat generation) |
| check_hcp_chunks.sh | QC for Picard shard completeness |
| pool_expression_within_ancestry.py | Expression pooling implementation |
| pool_modalities_within_ancestry.py | Modality pooling implementation |
| normalize_expression_hcp.R | Pooled QN + within-cohort INT + HCP estimation (expression) |
| normalize_modalities.R | Pooled QN + within-cohort INT for non-expression modalities |
| hcp_from_matrix.R | HCP estimation from a matrix (used by 25b) |
| hcp_diagnostic.R, hcp_diagnostic2.R | HCP diagnostic plots |
| peer_factors.py | PEER-style factor helper |
| 20_normalize_modalities.sh | Pool + QN + within-cohort INT for non-expression modalities |
| 21_install_tensorqtl.sh | Build the tensorqtl conda env and install R qvalue |
| 22_genotype_pca.sh | Per-ancestry genotype PCA (plink2) |
| genotype_pca_format.py | Reshape plink2 eigenvec/eigenval to covariate TSV + scree input |
| PCA_scree.R | Scree plots per ancestry |
| 23_prepare_intersection.py | Genotype x phenotype intersection; pgen filtered to intersection samples (MAC >= 5); IDs renamed rnaseq_id -> array_id |
| 24_outlier_exclusion.py | 6-SD outlier exclusion on genotype PC1-5 (edits intersection files in place) |
| 25_build_covariates.py | Assemble covariate table (PCs + HCPs + sex + cell types), prune correlated covariates |
| 25a_optimize_hcp.sh | Expression-only HCP k grid search on chr1 |
| optimize_hcp_chr1.py | Implementation of 25a |
| 25b_optimize_hcp_modalities.sh | Per-modality HCP k grid search |
| optimize_hcp_modalities.py | Implementation of 25b |
| 26_harmonize_modalities.py | Harmonize modality BEDs to the final array_id sample set |
| 27_run_tensorqtl.py | tensorQTL cis mapping for one ancestry x modality |
| 27_run_tensorqtl.sh | Job wrapper for 27 |
| 28_submit_modalities.sh | Submit all ancestry x modality tensorQTL jobs (grouped / ungrouped / combined layers; stepwise independent scans via INDEPENDENT=1) |
| compute_qvalues.R | Storey q-values via the R qvalue package (file bridge) |
| 29_make_top_tables.py | Top-association tables from parquet outputs (all layers: grouped, ungrouped, combined, plus independent scans) |
| 29b_expression_diagnostics.py | Expression diagnostic inputs: sample PCs, per-cohort residual variances, k=15 vs k=45 lead-stability table |
| 29c_choi_comparison.py | Choi 2024 comparison: gene-top / exact-lead tables, significant-set retention, common-test catalog |
| 29d_cohort_heterogeneity.py | Per-cohort cis scans for pooled-significant pairs; inverse-variance meta-analysis (Cochran Q, I²); genotype x cohort interaction scan |
| 29e_plot_diagnostics.R | Diagnostic figures + baseline-vs-rerun validation summary |
| 29f_run_diagnostics.sh | LSF driver chaining 29b-29e per ancestry (gate before fine-mapping) |
| 30_combine_modalities.py | Combined cross-modality BED (namespaced phenotype IDs, gene groups) |
| 31_sushie_finemap.py | SuSHiE locus preparation (ungrouped phenotype-level discovery layer; cis windows reconstructed from phenotype BEDs) and per-shard fine-mapping |
| 32_submit_sushie.sh | SuSHiE LSF driver (locus prep, sharding, array submission) |
| 32a_run_sushie_shard.sh | Per-shard SuSHiE job |
| 33_aggregate_finemap.py | Aggregate per-locus SuSHiE outputs (PIPs, credible sets, locus summary) |
| 34_pip_annotation_enrichment.py | fastVEP consequence + cCRE/OCR enrichment of high-PIP variants |
| 35_extract_report_extras.py | Targeted extractions for DevBrain-style report panels (cross-ancestry/cross-modality lead lookups, showcase-locus scans, lead-variant annotations, gnomAD constraint, PC-AiR staging) and QC summaries (before/after phenotype values/PCAs, genotype stage counts, Picard metrics) |
| 25a_submit_hcp_k_jobs.sh, 25b_submit_hcp_k_jobs.sh | Restart-safe one-job-per-k HCP grid submitters (25a expression, 25b modalities + combined); `LAMBDA1` selects the HCP prior strength (non-default values use a lambda-suffixed sandbox tree and skip finalization unless `FINALIZE=1`) |
| 25_hcp_k_finalize.sh, finalize_hcp_k_grid.py | Per-ancestry/arm grid finalizers: combine completed k points, select k\\*, install canonical HCP factors |
| test_build_covariates.py, test_collapse_replicates.py, test_normalize_modalities.py, test_hcp.py, test_optimize_hcp.py, test_optimize_hcp_modalities.py, test_extract_report_extras.py, test_prepare_intersection.py, test_picard_qc_fullpanel.py, test_union_pooling.py | Unit tests |
| test_sushie_finemap.py | End-to-end fine-mapping fixture (locus prep, SuSHiE recovery, aggregation, enrichment, report render) |

## reports/

| File | Purpose |
|------|---------|
| report_placenta_xqtl.Rmd | Main report source |
| report_placenta_xqtl.html | Rendered self-contained main report |
| xqtl_report_functions.R | Helper module (loaders, QQ/power/concordance helpers, report theme) |
| gene_map.tsv | gene_id -> symbol/biotype/description map |
| report_finemap.Rmd | Fine-mapping diagnostic report source |
| make_report_archive.sh | Build the report-input tarball from mapping outputs |
| runbook_report_extras_seadragon.md | Step-by-step seadragon runbook for the report-extras extraction (incl. QC section) and re-archiving |
