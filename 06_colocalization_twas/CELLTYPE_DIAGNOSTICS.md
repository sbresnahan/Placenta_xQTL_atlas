# Module 06 cell-type attribution diagnostics — sharded LSF implementation

This stage diagnoses the ancestry contrast in the Module-06 genotype × cell-proportion (`G×C`) attribution model. **EUR N=136 is the complete EUR cohort size and is expected.** The main power diagnostic repeatedly downsamples EAS to N=136 and asks how much of the EAS assignment yield survives under the EUR-sized sample.

## Why the original diagnostic job was slow

The first implementation parallelized over features *inside one LSF job*, but it still processed modalities serially. More importantly, the downsampling loop performed a separate regression for every `feature × cell type × replicate`. With 3,606 EAS features, seven cell types, and 100 downsampling replicates, that is roughly **2.5 million regression refits** before HC3 and permutation diagnostics. A 12-core job therefore spent most of its wall time inside one modality at a time and could not exploit the cluster across modalities.

The revised implementation changes both levels of parallelism:

1. **One LSF job array per modality.** Each array contains phenotype shards from both ancestries. Different modalities can run concurrently.
2. **Phenotype sharding.** By default an array element handles 100 xQTL features and uses a small number of local worker processes.
3. **Batched downsampling.** Within a replicate, all features for a given cell type are fitted together using a Frisch–Waugh–Lovell matrix formulation instead of thousands of independent `lstsq` calls.
4. **Replicate parallelism.** Downsample replicates are distributed across the CPUs assigned to that shard.
5. **Reusable input cache.** PLINK lead-variant dosages and selected phenotype BED rows are extracted once, not independently by every array task.
6. **Modality-wide aggregation.** Raw downsample P values from every phenotype shard are combined before BH correction, preserving the original ancestry × modality multiple-testing unit.
7. **Resume support.** Completed shards and modality aggregations are skipped unless `FORCE=1`.
8. **BLAS threads are pinned to one** inside each worker to prevent worker × BLAS oversubscription.

## Files

The sharded implementation adds:

```text
50a_prepare_celltype_diagnostics.py
50b_run_celltype_diag_shard.py
50c_aggregate_celltype_diag_modality.py
50d_finalize_celltype_diagnostics.py
```

`50_celltype_diagnostics.py` remains the shared statistical core and contains the exact full-N diagnostics, HC3, permutation, batched FWL downsampling, ancestry matching, and final summary helpers.

## Submit

Recommended production settings:

```bash
DOWNSAMPLE_REPS=100 \
SHARD_FEATURES=100 \
CPUS_PER_SHARD=2 \
ARRAY_CONCURRENCY=8 \
bash 50_submit_celltype_diagnostics.sh
```

`THREADS` is no longer the principal parallelism setting. The useful controls are:

| Variable | Default | Meaning |
|---|---:|---|
| `SHARD_FEATURES` | 100 | unique xQTL features per array element |
| `CPUS_PER_SHARD` | 2 | local workers within an array element |
| `ARRAY_CONCURRENCY` | 8 | maximum simultaneously running elements **per modality** |
| `SHARD_MEM_GB` | 8 | memory per element |
| `SHARD_WALLTIME` | 04:00 | walltime per element |
| `DOWNSAMPLE_REPS` | 100 | EAS→EUR-sized replicate fits |
| `PERMUTATIONS` | 100 | residual permutations for selected tests |
| `HC3` | 1 | robust-SE diagnostic |
| `FORCE` | 0 | resume existing work |

Because every modality has its own array, `ARRAY_CONCURRENCY=8` is an upper bound per modality, not globally. LSF still controls how many elements can actually run based on available resources and fair-share policy. If this submits more concurrency than desired, lower it to 3–4.

## Execution graph

```text
plan generation (submission host)
          |
          v
ctd_prepare  -- cache all required inputs once
          |
          +------------------+-------------------+ ...
          v                  v
ctd_expression[]       ctd_isoform_expression[]   ...  one array/modality
          |                  |
          v                  v
ctd_expression_agg     ctd_isoform_expression_agg ...
          \                  /
           \                /
              ctd_finalize
                   |
                   v
celltype_diagnostics_bundle.zip
```

Each modality plan is stored under:

```text
${OUT_DIR}/work/plans/<modality>.tsv
```

Intermediate shard outputs are under:

```text
${OUT_DIR}/work/shards/<modality>/
```

and modality aggregates under:

```text
${OUT_DIR}/work/modality/<modality>/
```

These intermediates make interrupted runs resumable.

## Final outputs

The final public/reportable outputs remain:

```text
sample_overlap_summary.tsv
sample_overlap_by_sample.tsv.gz
celltype_distribution.tsv
gxc_support.tsv.gz
gxc_hc3.tsv.gz
gxc_permutation.tsv.gz
gxc_downsample_replicates.tsv
gxc_downsample_tests.tsv.gz
gxc_downsample_features.tsv.gz
gxc_matched_ancestry.tsv.gz
gxc_matched_ancestry_summary.tsv
diagnostic_summary.tsv
run_metadata.json
celltype_diagnostics_bundle.zip
```

Upload `celltype_diagnostics_bundle.zip` for interpretation.

## Monitoring

```bash
bjobs -w | grep ctd_
```

Inspect one array element with:

```bash
bpeek 'JOBID[INDEX]'
```

The modality aggregator does not start until every element in that modality array has completed successfully. The finalizer waits for all pending modality aggregators.
