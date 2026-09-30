# 02_rnaseq_quantification

SRA retrieval, adapter trimming, and Salmon quantification for the cohorts
whose RNA-seq FASTQs were not already local (SNUH) or were quantified
outside the phenotyping driver (NIEHS_RICHS). LSF jobs for the seadragon
cluster.

## Prerequisites

- Software: sratoolkit, Trim Galore 0.6.10, Salmon 1.10.2 (see
  [../docs/software_environments.md](../docs/software_environments.md))
- A Salmon index for the quantification reference

## Steps

1. **SNUH retrieval** — one job per accession:

   ```bash
   bsub -env "ACC=SRR..." < pull_sra_SNUH.lsf   # prefetch + fasterq-dump --split-files + gzip
   ```

2. **Quantification** — submit one job per FASTQ (NIEHS_RICHS, single-end)
   or per accession (SNUH, paired-end):

   ```bash
   bash loop_salmon_NIEHS_RICHS.sh   # submits run_salmon_NIEHS_RICHS.lsf per FASTQ
   bash loop_salmon_SNUH.sh          # submits run_salmon_SNUH.lsf per accession
   ```

   Each run script: Trim Galore, then
   `salmon quant -l A --validateMappings --seqBias --numBootstraps 50`.

## Outputs

Per-sample Salmon `quant.sf` directories alongside the FASTQs.

## Operational notes

- GUSTO FASTQs were already local; NIGMS RNA-seq runs come from SRA
  PRJNA671171 (see
  [../docs/data_availability.md](../docs/data_availability.md)).
- Production phenotyping re-quantifies every sample inside the phenotyping
  pipeline (`../03_phenotyping/02_salmon_expression.sh`, 20 bootstraps,
  HPLRv2 index from `config.yml`). The wrappers here document the initial
  fetch/trim/quant steps.

Next: [../03_phenotyping/README.md](../03_phenotyping/README.md).
