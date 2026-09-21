# 02_rnaseq_quantification — SRA retrieval, trimming, Salmon quantification

Fetch/trim/quant wrappers for the cohorts whose RNA-seq FASTQs were not already
local (SNUH) or needed (re)quantification outside the PANTRY driver
(NIEHS_RICHS). LSF jobs for the seadragon cluster.

| Script | Purpose |
|---|---|
| `pull_sra_SNUH.lsf` | `prefetch` + `fasterq-dump --split-files` + gzip for one accession (`bsub -env ACC=SRR...`) |
| `loop_salmon_NIEHS_RICHS.sh` | Submit one `run_salmon_NIEHS_RICHS.lsf` job per FASTQ |
| `run_salmon_NIEHS_RICHS.lsf` | Trim Galore → `salmon quant -l A --validateMappings --seqBias --numBootstraps 50` (single-end) |
| `loop_salmon_SNUH.sh` | Submit one `run_salmon_SNUH.lsf` job per accession |
| `run_salmon_SNUH.lsf` | Trim Galore → Salmon quant (paired-end) |

## Notes

- GUSTO FASTQs were already local; NIGMS RNA-seq runs come from SRA
  PRJNA671171 (see [`../docs/data_availability.md`](../docs/data_availability.md)).
- Production phenotyping re-quantifies every sample inside the PANTRY pipeline
  (`03_phenotyping/02_salmon_expression.sh`, 20 bootstraps, HPLRv2 index from
  `config.yml`). The wrappers here document the initial fetch/trim/quant steps.
- Software: sratoolkit, Trim Galore 0.6.10, Salmon 1.10.2.
