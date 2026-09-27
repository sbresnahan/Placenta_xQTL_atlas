#!/usr/bin/env python3
"""
optimize_hcp_modalities.py — per-modality HCP-count optimization (module 25b)

Generalizes optimize_hcp_chr1.py (25a, expression-only) to every modality and
the combined cross-modality arm: HCPs are estimated from, and optimized for,
exactly the matrix each tensorQTL run maps. One ancestry x modality per
invocation is the intended parallel unit (the 25b LSF wrapper submits one job
per pair); multiple pairs per invocation are supported for serial runs.

Design (per ancestry x modality, per k in --k-grid):
  1. HCP-only estimation at k on the FULL harmonized BED
     ({ANC}_{MOD}_harmonized.bed — QN+INT+ComBat already applied) via
     hcp_from_matrix.R. k=0 writes a header-only file (no HCP covariates).
     --max-hcp-phenotypes bounds the combined arm's HCP cost with a
     deterministic phenotype subsample (estimation only; mapping always
     uses the full-scope staging BED below).
  2. Covariates: 25_build_covariates.py --hcp-file <per-k> --hcp-k k
     --out-suffix _k{k} against a staging dir (symlinked genotypes/metadata).
     --exclude-covariates is forwarded; no covariate cap.
  3. Mapping-scope staging BED: chr1 subset when the harmonized BED has
     >= --chr1-min-phenotypes chr1 rows, else genome-wide (small modalities
     such as RNA_editing have ~30 chr1 phenotypes — too few for a stable
     k* signal). The groups file, when present, is subset to the staging
     phenotypes so mapping stays grouped.
  4. Map: 27_run_tensorqtl.py --modality <MOD> with per-k covariates.
  5. Count eGenes: unique phenotype_id with Storey qval <= --fdr.

k* = argmax(eGenes); ties break toward the smaller k. Writes
{work_dir}/{ANC}_{MOD}_optimal_hcp.tsv + .png and installs the k* HCP set as
{qtl_dir}/{ANC}_{MOD}_hcp_factors_optimized.tsv (input to the canonical
per-modality 25_build_covariates.py run; see docs/runbook_modality_hcp.md).

Usage (typically via 25b_optimize_hcp_modalities.sh):
  python3 optimize_hcp_modalities.py \
      --qtl-dir <qtl_inputs dir> \
      --qc-metrics <hcp dir>/all_qc_metrics.tsv \
      --pcair-dir <genotype_pcs dir> \
      --scripts-dir <scripts dir> \
      --ancestries "EAS" --modalities "splicing" \
      --k-grid "0 5 10 15 20 25 30" \
      --r-cmd "Rscript"
"""

import argparse
import os
import shlex
import shutil
import sys

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import optimize_hcp_chr1 as base  # shared helpers: run, bgzip_tabix, write_empty_hcp, count_egenes

plt.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']

ALL_MODALITIES = ("expression isoforms isoform_expression splicing "
                  "intron_retention alt_TSS alt_polyA RNA_editing stability combined")


def count_chr1_phenotypes(bed_path):
    """Count chr1 rows in a harmonized BED (header-aware, chr-prefixed)."""
    bed = pd.read_csv(bed_path, sep='\t', usecols=['#chr'])
    return int((bed['#chr'] == 'chr1').sum()), int(bed.shape[0])


def write_staging_bed(bed_path, staging_bed, scope):
    """Write the mapping-scope staging BED: 'chr1' subsets to chr1 rows,
    'genome-wide' copies all rows. Returns the phenotype_id list."""
    bed = pd.read_csv(bed_path, sep='\t')
    if scope == 'chr1':
        bed = bed[bed['#chr'] == 'chr1']
    bed.to_csv(staging_bed, sep='\t', index=False)
    return bed['phenotype_id'].tolist()


def subset_groups_file(groups_path, phenotype_ids, out_path):
    """Filter a phenotype_groups.txt (no header: phenotype_id, group_id) to
    the staging-BED phenotypes so grouped mapping covers the subset."""
    groups = pd.read_csv(groups_path, sep='\t', header=None,
                         names=['phenotype_id', 'group_id'])
    keep = set(phenotype_ids)
    sub = groups[groups['phenotype_id'].isin(keep)]
    sub.to_csv(out_path, sep='\t', index=False, header=False)
    return sub.shape[0]


def plot_optimization(res_df, anc, mod, k_star, fdr, scope_label, out_path):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(res_df['k'], res_df['n_egenes'], 'o-', color='#0279EE',
            markersize=6, linewidth=1.5)
    ax.axvline(k_star, color='#FF9400', linestyle='--', linewidth=1.2,
               label=f'k* = {k_star}')
    for _, row in res_df.iterrows():
        ax.annotate(str(int(row['n_egenes'])),
                    (row['k'], row['n_egenes']),
                    textcoords='offset points', xytext=(0, 7),
                    ha='center', fontsize=8)
        if 'n_hcp_dropped' in res_df.columns and row.get('n_hcp_dropped', 0) > 0:
            ax.annotate(f"({int(row['n_hcp_used'])} kept)",
                        (row['k'], row['n_egenes']),
                        textcoords='offset points', xytext=(0, -13),
                        ha='center', fontsize=7, color='#FF9400')
    ax.set_xlabel('Number of HCP factors (k)')
    ax.set_ylabel(f'eGenes (q <= {fdr})')
    ax.set_title(f'{anc} {mod}: HCP count optimization ({scope_label})')
    ax.legend(frameon=False)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def optimize_one(anc, mod, args, k_grid, work_root, r_cmd):
    """Run the full k grid for one ancestry x modality. Returns the
    per-k results DataFrame with the chosen k flagged."""
    print(f"\n{'='*60}\nAncestry: {anc} / Modality: {mod}\n{'='*60}")
    staging = os.path.join(work_root, anc, mod)
    os.makedirs(staging, exist_ok=True)
    log_path = os.path.join(staging, f'{anc}_{mod}_optimize_hcp.log')

    # ---- Required inputs ----
    meta_path = os.path.join(args.qtl_dir, f'{anc}_metadata.tsv')
    harm_bed = os.path.join(args.qtl_dir, f'{anc}_{mod}_harmonized.bed')
    for f in [meta_path, harm_bed]:
        if not os.path.exists(f):
            sys.exit(f"ERROR: required input not found: {f} "
                     f"(run 23/24/26/30 first)")
    meta = pd.read_csv(meta_path, sep='\t')
    print(f"  Samples (post-outlier, array_id): {meta.shape[0]}")

    # Symlink genotype + covariate-input files into staging
    for fname in [f'{anc}_qtl.pgen', f'{anc}_qtl.pvar', f'{anc}_qtl.psam',
                  f'{anc}_selected_pcs.txt', f'{anc}_metadata.tsv',
                  f'{anc}_deconvolution_harmonized.tsv']:
        src = os.path.join(args.qtl_dir, fname)
        dst = os.path.join(staging, fname)
        if os.path.exists(src) and not os.path.exists(dst):
            os.symlink(os.path.abspath(src), dst)

    # ---- Mapping scope: chr1 subset vs genome-wide fallback ----
    n_chr1, n_total = count_chr1_phenotypes(harm_bed)
    scope = 'chr1' if n_chr1 >= args.chr1_min_phenotypes else 'genome-wide'
    print(f"  Phenotypes: {n_total} total, {n_chr1} chr1 -> mapping scope: {scope} "
          f"(chr1 minimum: {args.chr1_min_phenotypes})")

    staging_gz = os.path.join(staging, f'{anc}_{mod}.bed.gz')
    if not (args.skip_existing and os.path.exists(staging_gz + '.tbi')):
        staging_bed = os.path.join(staging, f'{anc}_{mod}.bed')
        pheno_ids = write_staging_bed(harm_bed, staging_bed, scope)
        print(f"    Staging BED: {len(pheno_ids)} phenotypes ({scope})")
        base.bgzip_tabix(staging_bed, staging_gz, log_path)
        os.remove(staging_bed)
    else:
        # Recover the staging phenotype list for groups subsetting
        import gzip
        with gzip.open(staging_gz, 'rt') as fh:
            pheno_ids = [ln.split('\t')[3] for ln in fh if not ln.startswith('#')]

    # Groups file: subset to staging phenotypes (grouped mapping)
    src_groups = os.path.join(args.qtl_dir, f'{anc}_{mod}.phenotype_groups.txt')
    dst_groups = os.path.join(staging, f'{anc}_{mod}.phenotype_groups.txt')
    if os.path.exists(src_groups):
        if not (args.skip_existing and os.path.exists(dst_groups)):
            n_kept = subset_groups_file(src_groups, pheno_ids, dst_groups)
            print(f"    Groups file: {n_kept} phenotypes covered ({scope})")

    # ---- k grid ----
    results = []
    for k in k_grid:
        print(f"\n  -- k = {k} --")
        per_k_hcp = os.path.join(staging, f'{anc}_{mod}_hcp_k{k}.tsv')
        results_dir = os.path.join(staging, f'results_k{k}')
        parquet = os.path.join(results_dir, f'{anc}_{mod}_cisqtl.parquet')

        # 1. HCP estimation at k (full harmonized BED, capped subsample)
        if not (args.skip_existing and os.path.exists(per_k_hcp)):
            if k == 0:
                base.write_empty_hcp(meta, per_k_hcp)
                print("    k=0: no HCP covariates (sample set only)")
            else:
                base.run(r_cmd + [os.path.join(args.scripts_dir, 'hcp_from_matrix.R'),
                                  '--bed', harm_bed,
                                  '--qc-metrics', args.qc_metrics,
                                  '--metadata', meta_path,
                                  '--k', str(k),
                                  '--output', per_k_hcp,
                                  '--lambda1', str(args.lambda1),
                                  '--lambda2', str(args.lambda2),
                                  '--lambda3', str(args.lambda3),
                                  '--qc-cor-threshold', str(args.qc_cor_threshold),
                                  '--max-phenotypes', str(args.max_hcp_phenotypes),
                                  '--seed', str(args.seed)],
                         f"HCP estimation (k={k})", log_path)

        # 2. Covariates for this k
        cov_file = os.path.join(staging, f'{anc}_covariates_k{k}.tsv')
        if not (args.skip_existing and os.path.exists(cov_file)):
            cov_cmd = [sys.executable,
                       os.path.join(args.scripts_dir, '25_build_covariates.py'),
                       '--qtl-dir', staging,
                       '--pcair-dir', args.pcair_dir,
                       '--ancestries', anc,
                       '--hcp-file', per_k_hcp,
                       '--hcp-k', str(k),
                       '--out-suffix', f'_k{k}']
            if args.exclude_covariates:
                cov_cmd += ['--exclude-covariates'] + list(args.exclude_covariates)
            base.run(cov_cmd, f"covariate assembly (k={k})", log_path)

        # 2b. Effective HCP count after correlation pruning
        n_hcp_dropped = 0
        pruning_file = os.path.join(staging, f'{anc}_covariate_pruning_k{k}.tsv')
        if k > 0 and os.path.exists(pruning_file):
            pruned = pd.read_csv(pruning_file, sep='\t')
            n_hcp_dropped = int(pruned['covariate'].astype(str)
                                .str.startswith('HCP_').sum())
        n_hcp_used = k - n_hcp_dropped
        if n_hcp_dropped:
            print(f"    k={k}: {n_hcp_dropped} HCP(s) dropped by pruning "
                  f"-> {n_hcp_used} enter the model")

        # 3-4. Mapping (grouped when a groups file exists in staging)
        if not (args.skip_existing and os.path.exists(parquet)):
            base.run([sys.executable,
                      os.path.join(args.scripts_dir, '27_run_tensorqtl.py'),
                      '--qtl-dir', staging,
                      '--output-dir', results_dir,
                      '--ancestry', anc,
                      '--modality', mod,
                      '--covariates-file', cov_file,
                      '--cis-window', str(args.cis_window),
                      '--maf-threshold', str(args.maf_threshold)],
                     f"tensorQTL cis mapping (k={k})", log_path)

        # 5. eGene count
        n_eg, n_tested = base.count_egenes(parquet, args.fdr)
        print(f"    k={k}: {n_eg} eGenes (of {n_tested} tested, {scope})")
        results.append({'ancestry': anc, 'modality': mod, 'k': k,
                        'n_egenes': n_eg, 'n_tested': n_tested,
                        'mapping_scope': scope,
                        'n_hcp_used': n_hcp_used,
                        'n_hcp_dropped': n_hcp_dropped})

    # ---- Select k* (max eGenes; ties -> smaller k) ----
    res_df = pd.DataFrame(results).sort_values('k').reset_index(drop=True)
    best = res_df['n_egenes'].max()
    k_star = int(res_df.loc[res_df['n_egenes'] == best, 'k'].min())
    res_df['chosen'] = res_df['k'] == k_star
    print(f"\n  Optimal k for {anc} {mod}: {k_star} ({best} eGenes)")

    tsv_path = os.path.join(work_root, f'{anc}_{mod}_optimal_hcp.tsv')
    res_df.to_csv(tsv_path, sep='\t', index=False)
    print(f"  Written: {tsv_path}")

    png_path = os.path.join(work_root, f'{anc}_{mod}_optimal_hcp.png')
    plot_optimization(res_df, anc, mod, k_star, args.fdr, scope, png_path)
    print(f"  Written: {png_path}")

    # ---- Install k* as the per-modality optimized HCP file ----
    installed = os.path.join(args.qtl_dir, f'{anc}_{mod}_hcp_factors_optimized.tsv')
    kstar_file = os.path.join(staging, f'{anc}_{mod}_hcp_k{k_star}.tsv')
    shutil.copy2(kstar_file, installed)
    print(f"  Installed k*={k_star} HCP solution: {installed}")
    return res_df


def main():
    parser = argparse.ArgumentParser(
        description="Per-modality HCP-count optimization (module 25b)")
    parser.add_argument("--qtl-dir", required=True,
                        help="Canonical QTL inputs dir (post-23/24/26/30 outputs)")
    parser.add_argument("--qc-metrics", required=True,
                        help="Pooled Picard QC metrics TSV "
                             "(e.g. {OUTPUT_BASE}/hcp/all_qc_metrics.tsv)")
    parser.add_argument("--pcair-dir", required=True,
                        help="Genotype PCs directory")
    parser.add_argument("--scripts-dir", required=True,
                        help="Directory with hcp_from_matrix.R, "
                             "25_build_covariates.py, 27_run_tensorqtl.py")
    parser.add_argument("--ancestries", default="EAS EUR",
                        help="Space-separated ancestry labels")
    parser.add_argument("--modalities", default=ALL_MODALITIES,
                        help="Space-separated modality labels, including "
                             "'combined' (default: all 9 modalities + combined)")
    parser.add_argument("--k-grid", default="0 5 10 15 20 25 30",
                        help="Space-separated candidate HCP counts")
    parser.add_argument("--work-dir", default=None,
                        help="Staging dir (default: {qtl-dir}/hcp_optimization_modalities)")
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="Storey q-value threshold for eGene counts (default: 0.05)")
    parser.add_argument("--r-cmd", default="Rscript",
                        help="R command for HCP estimation (may be a "
                             "singularity prefix; default: Rscript)")
    parser.add_argument("--lambda1", type=float, default=0.5,
                        help="HCP prior strength (default: 0.5, as in 19_hcp_factors.sh)")
    parser.add_argument("--lambda2", type=float, default=1.0)
    parser.add_argument("--lambda3", type=float, default=1.0)
    parser.add_argument("--qc-cor-threshold", type=float, default=0.9)
    parser.add_argument("--chr1-min-phenotypes", type=int, default=300,
                        help="Minimum chr1 phenotypes for chr1-based k* selection; "
                             "smaller modalities map genome-wide (default: 300)")
    parser.add_argument("--max-hcp-phenotypes", type=int, default=40000,
                        help="Phenotype cap for HCP ESTIMATION (deterministic "
                             "subsample; 0 disables). Bounds the combined arm's "
                             "HCP cost; mapping is unaffected (default: 40000)")
    parser.add_argument("--seed", type=int, default=1,
                        help="Seed for the HCP phenotype subsample (default: 1)")
    parser.add_argument("--exclude-covariates", nargs='+', default=None,
                        metavar='NAME',
                        help="Covariates excluded before correlation pruning in "
                             "every per-k model (e.g. ct_Maternal). Forwarded to "
                             "25_build_covariates.py --exclude-covariates.")
    parser.add_argument("--cis-window", type=int, default=1000000)
    parser.add_argument("--maf-threshold", type=float, default=0.01)
    parser.add_argument("--skip-existing", action="store_true",
                        help="Reuse existing per-k outputs (resumable reruns)")
    args = parser.parse_args()

    k_grid = [int(k) for k in args.k_grid.split()]
    work_root = args.work_dir or os.path.join(args.qtl_dir, 'hcp_optimization_modalities')
    os.makedirs(work_root, exist_ok=True)
    r_cmd = shlex.split(args.r_cmd)

    if not os.path.exists(args.qc_metrics):
        sys.exit(f"ERROR: pooled QC metrics not found: {args.qc_metrics} (run 19 first)")

    all_results = []
    for anc in args.ancestries.split():
        for mod in args.modalities.split():
            res_df = optimize_one(anc, mod, args, k_grid, work_root, r_cmd)
            all_results.append(res_df)

    if len(all_results) > 1:
        summary = pd.concat(all_results, ignore_index=True)
        summary_path = os.path.join(work_root, 'optimal_hcp_summary.tsv')
        summary.to_csv(summary_path, sep='\t', index=False)
        print(f"\nWritten: {summary_path}")

    print(f"\nNext: run the canonical per-modality 25_build_covariates.py "
          f"(see docs/runbook_modality_hcp.md).")


if __name__ == '__main__':
    main()
