#!/bin/bash
# =============================================================================
# 07_majiq_build.sh — Per-sample MAJIQ build (splice graph + .sj coverage)
# =============================================================================
# MAJIQ build produces the splice graph (splicegraph.zarr) and per-sample
# junction coverage (.sj files). The splice graph is shared across all samples
# in a cohort, so only the FIRST sample's build creates it; subsequent samples
# append .sj files. To keep jobs independent (no cross-job file locking), each
# sample runs its own build into a per-sample directory, and the aggregation
# step (13_aggregate_intron_retention.sh) merges them.
#
# Actually, MAJIQ build takes an experiments.tsv listing ALL BAMs and produces
# the splice graph + all .sj files in one run. So this per-sample script is
# NOT used for build. Instead, the driver writes the experiments.tsv and the
# aggregation script runs `majiq build` once for the whole cohort.
#
# This script is kept as a placeholder for potential per-sample .sj pre-compute
# but the standard PANTRY flow runs majiq build at the cohort level.
# See 13_aggregate_intron_retention.sh for the actual MAJIQ build + psi flow.
# =============================================================================

echo "07_majiq_build.sh: MAJIQ build runs at the cohort level (see 13_aggregate_intron_retention.sh)"
echo "This per-sample script is a no-op placeholder."
echo "MAJIQ build requires all BAMs listed in one experiments.tsv, so it cannot be split per-sample."
exit 0
