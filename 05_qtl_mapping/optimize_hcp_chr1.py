#!/usr/bin/env python3
"""
optimize_hcp_chr1.py — chr1-only HCP-count optimization (module 25a driver)

Pre-mapping module that selects the number of HCP hidden covariates per
ancestry by maximizing cis-eGene discovery on chromosome 1, following the
devBrain xQTL atlas convention (Wen et al., Science 2024, 384:eadh0829,
§4.2: HCP count chosen by maximizing eGenes at FDR < 0.05, with HCP
re-estimated at each candidate k).

Run AFTER 23_prepare_intersection.py and 24_outlier_exclusion.py (needs
{ANC}_metadata.tsv, {ANC}_expression_harmonized.bed, {ANC}_qtl.pgen,
{ANC}_selected_pcs.txt, {ANC}_deconvolution_harmonized.tsv in --qtl-dir),
and BEFORE the canonical 25_build_covariates.py run.

Per ancestry, per k in --k-grid (default: 0 5 10 15 20 25 30):
  1. Re-estimate HCP at k: combat_normalize_hcp.R --k k into a per-k dir
     (full QN -> INT -> outlier removal -> ComBat -> HCP chain; k=0 skips
     HCP estimation via --skip-hcp).
  2. Harmonize HCP factors to array_id space, restricted to the
     post-outlier sample set in {qtl_dir}/{ANC}_metadata.tsv.
  3. Build covariates: 25_build_covariates.py --hcp-file <per-k> --hcp-k k
     --out-suffix _k{k} against a staging dir (symlinked genotypes/metadata).
  4. Subset the harmonized expression BED to chr1 -> bgzip + tabix.
  5. Map: 27_run_tensorqtl.py --modality expression with per-k covariates
     (cis-window 1 Mb, MAF >= 0.01, 5 genotype PCs via the covariate table).
  6. Count eGenes: unique phenotype_id with Storey qval <= --fdr.

Select k* = argmax(eGenes); ties broken toward the smaller k.
Writes {work_dir}/{ANC}_optimal_hcp.tsv + {ANC}_optimal_hcp.png, and
installs the k* solution as {qtl_dir}/{ANC}_hcp_factors_harmonized.tsv
(the pre-existing file is backed up to *.pre25a_backup.tsv).

Usage (typically via 25a_optimize_hcp.sh):
  python3 optimize_hcp_chr1.py \
      --qtl-dir <qtl_inputs dir> \
      --hcp-dir <hcp dir> \
      --pcair-dir <genotype_pcs dir> \
      --ancestry-map <pooled_sample_ancestry_RNAseq.tsv> \
      --scripts-dir <scripts dir> \
      --ancestries "EAS EUR" \
      --k-grid "0 5 10 15 20 25 30" \
      --r-cmd "Rscript"
"""

import argparse
import os
import shlex
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']


def run(cmd, desc, log_path=None):
    """Run a command (list or shell string), teeing output to an optional log."""
    printable = cmd if isinstance(cmd, str) else ' '.join(cmd)
    print(f"    [{pd.Timestamp.now()}] {desc}")
    print(f"      $ {printable}")
    with open(log_path, 'a') if log_path else open(os.devnull, 'w') as logf:
        proc = subprocess.run(cmd, shell=isinstance(cmd, str),
                              stdout=logf, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        sys.exit(f"ERROR: {desc} failed (exit {proc.returncode}). "
                 f"See log: {log_path or '(none)'}")


def bgzip_tabix(bed_path, out_gz):
    """bgzip + tabix a BED file. Tries pysam, falls back to htslib binaries."""
    try:
        import pysam
        pysam.tabix_compress(bed_path, out_gz, force=True)
        pysam.tabix_index(out_gz, preset='bed', force=True)
        return
    except ImportError:
        pass
    # htslib binaries fallback (bgzip -c writes to stdout -> shell redirect)
    run(f"bgzip -f -c {bed_path} > {out_gz}", f"bgzip {os.path.basename(bed_path)}")
    run(['tabix', '-f', '-p', 'bed', out_gz], f"tabix {os.path.basename(out_gz)}")


def harmonize_hcp(hcp_path, meta, out_path):
    """Rename HCP factor columns rnaseq_id -> array_id and restrict to the
    post-outlier sample set in meta (the {ANC}_metadata.tsv table)."""
    hcp = pd.read_csv(hcp_path, sep='\t', index_col=0)
    id_map = dict(zip(meta['rnaseq_id'], meta['array_id']))
    keep_array = set(meta['array_id'])
    new_cols, seen = [], set()
    for c in hcp.columns:
        mapped = id_map.get(c, c)  # tolerate already-array_id columns
        if mapped in keep_array and mapped not in seen:
            new_cols.append(mapped)
            seen.add(mapped)
        else:
            new_cols.append(None)
    hcp = hcp.loc[:, [c is not None for c in new_cols]]
    hcp.columns = [c for c in new_cols if c is not None]
    hcp.index.name = 'covariate'
    hcp.to_csv(out_path, sep='\t')
    return hcp


def write_empty_hcp(meta, out_path):
    """Header-only HCP file (0 factors) carrying the sample set for k=0."""
    cols = ['covariate'] + sorted(meta['array_id'].tolist())
    with open(out_path, 'w') as f:
        f.write('\t'.join(cols) + '\n')


def count_egenes(parquet_path, fdr):
    """eGenes = unique phenotypes with Storey qval <= fdr."""
    df = pd.read_parquet(parquet_path)
    if 'qval' not in df.columns:
        sys.exit(f"ERROR: no qval column in {parquet_path} "
                 f"(columns: {list(df.columns)})")
    hits = df.loc[df['qval'] <= fdr, 'phenotype_id']
    return int(hits.nunique()), int(df.shape[0])


def plot_optimization(res_df, anc, k_star, fdr, out_path):
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
    ax.set_xlabel('Number of HCP factors (k)')
    ax.set_ylabel(f'chr1 eGenes (q <= {fdr})')
    ax.set_title(f'{anc}: HCP count optimization (chr1 expression)')
    ax.legend(frameon=False)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="chr1-only HCP-count optimization (module 25a)")
    parser.add_argument("--qtl-dir", required=True,
                        help="Canonical QTL inputs dir (post-23/24 outputs)")
    parser.add_argument("--hcp-dir", required=True,
                        help="HCP pipeline dir (script 19 outputs: "
                             "pooled_expression/, all_qc_metrics.tsv)")
    parser.add_argument("--pcair-dir", required=True,
                        help="Genotype PCs directory")
    parser.add_argument("--ancestry-map", required=True,
                        help="pooled_sample_ancestry_RNAseq.tsv")
    parser.add_argument("--scripts-dir", required=True,
                        help="Directory with combat_normalize_hcp.R, "
                             "25_build_covariates.py, 27_run_tensorqtl.py")
    parser.add_argument("--ancestries", default="EAS EUR",
                        help="Space-separated ancestry labels")
    parser.add_argument("--k-grid", default="0 5 10 15 20 25 30",
                        help="Space-separated candidate HCP counts")
    parser.add_argument("--work-dir", default=None,
                        help="Staging dir (default: {qtl-dir}/hcp_optimization)")
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="Storey q-value threshold for eGene counts (default: 0.05)")
    parser.add_argument("--r-cmd", default="Rscript",
                        help="R command for HCP re-estimation (may be a "
                             "singularity prefix; default: Rscript)")
    parser.add_argument("--lambda1", type=float, default=0.5,
                        help="HCP prior strength (default: 0.5, as in 19_hcp_factors.sh)")
    parser.add_argument("--lambda2", type=float, default=1.0)
    parser.add_argument("--lambda3", type=float, default=1.0)
    parser.add_argument("--qc-cor-threshold", type=float, default=0.9)
    parser.add_argument("--cis-window", type=int, default=1000000)
    parser.add_argument("--maf-threshold", type=float, default=0.01)
    parser.add_argument("--skip-existing", action="store_true",
                        help="Reuse existing per-k mapping results (resumable reruns)")
    args = parser.parse_args()

    k_grid = [int(k) for k in args.k_grid.split()]
    work_root = args.work_dir or os.path.join(args.qtl_dir, 'hcp_optimization')
    os.makedirs(work_root, exist_ok=True)
    r_cmd = shlex.split(args.r_cmd)

    pooled_expr = os.path.join(args.hcp_dir, 'pooled_expression')
    all_qc = os.path.join(args.hcp_dir, 'all_qc_metrics.tsv')
    if not os.path.exists(all_qc):
        sys.exit(f"ERROR: pooled QC metrics not found: {all_qc} (run 19 first)")

    for anc in args.ancestries.split():
        print(f"\n{'='*60}\nAncestry: {anc}\n{'='*60}")
        staging = os.path.join(work_root, anc)
        os.makedirs(staging, exist_ok=True)
        log_path = os.path.join(staging, f'{anc}_optimize_hcp.log')

        # ---- Required inputs from the canonical qtl dir ----
        meta_path = os.path.join(args.qtl_dir, f'{anc}_metadata.tsv')
        expr_harm = os.path.join(args.qtl_dir, f'{anc}_expression_harmonized.bed')
        if not os.path.exists(meta_path) or not os.path.exists(expr_harm):
            sys.exit(f"ERROR: run 23/24 first — missing {meta_path} or {expr_harm}")
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

        # ---- chr1 expression BED (phenotypes are k-independent) ----
        chr1_gz = os.path.join(staging, f'{anc}_expression.bed.gz')
        if not (args.skip_existing and os.path.exists(chr1_gz + '.tbi')):
            print("  Subsetting harmonized expression to chr1")
            bed = pd.read_csv(expr_harm, sep='\t')
            chr1 = bed[bed['#chr'] == 'chr1']
            print(f"    chr1 phenotypes: {chr1.shape[0]} of {bed.shape[0]}")
            chr1_bed = os.path.join(staging, f'{anc}_expression.bed')
            chr1.to_csv(chr1_bed, sep='\t', index=False)
            bgzip_tabix(chr1_bed, chr1_gz)
            os.remove(chr1_bed)

        # ---- k grid ----
        results = []
        for k in k_grid:
            print(f"\n  -- k = {k} --")
            per_k_hcp = os.path.join(staging, f'{anc}_hcp_k{k}_harmonized.tsv')
            results_dir = os.path.join(staging, f'results_k{k}')
            parquet = os.path.join(results_dir, f'{anc}_expression_cisqtl.parquet')

            # 1-2. HCP re-estimation + harmonization
            if not (args.skip_existing and os.path.exists(per_k_hcp)):
                if k == 0:
                    write_empty_hcp(meta, per_k_hcp)
                    print(f"    k=0: no HCP covariates (sample set only)")
                else:
                    hcp_out_dir = os.path.join(staging, f'hcp_k{k}')
                    os.makedirs(hcp_out_dir, exist_ok=True)
                    expr_file = os.path.join(pooled_expr, f'{anc}_pooled_expression.bed')
                    if not os.path.exists(expr_file):
                        sys.exit(f"ERROR: pooled expression not found: {expr_file}")
                    run(r_cmd + [os.path.join(args.scripts_dir, 'combat_normalize_hcp.R'),
                                 '--expression', expr_file,
                                 '--qc-metrics', all_qc,
                                 '--ancestry-map', args.ancestry_map,
                                 '--ancestry', anc,
                                 '--k', str(k),
                                 '--output-dir', hcp_out_dir,
                                 '--qc-cor-threshold', str(args.qc_cor_threshold),
                                 '--lambda1', str(args.lambda1),
                                 '--lambda2', str(args.lambda2),
                                 '--lambda3', str(args.lambda3)],
                        f"HCP re-estimation (k={k})", log_path)
                    raw_hcp = os.path.join(hcp_out_dir, f'{anc}_hcp_factors.tsv')
                    hcp_h = harmonize_hcp(raw_hcp, meta, per_k_hcp)
                    print(f"    Harmonized HCP: {hcp_h.shape[0]} factors x "
                          f"{hcp_h.shape[1]} samples")

            # 3. Covariates for this k
            cov_file = os.path.join(staging, f'{anc}_covariates_k{k}.tsv')
            if not (args.skip_existing and os.path.exists(cov_file)):
                run([sys.executable,
                     os.path.join(args.scripts_dir, '25_build_covariates.py'),
                     '--qtl-dir', staging,
                     '--pcair-dir', args.pcair_dir,
                     '--ancestries', anc,
                     '--hcp-file', per_k_hcp,
                     '--hcp-k', str(k),
                     '--out-suffix', f'_k{k}'],
                    f"covariate assembly (k={k})", log_path)

            # 4-5. chr1 cis mapping
            if not (args.skip_existing and os.path.exists(parquet)):
                run([sys.executable,
                     os.path.join(args.scripts_dir, '27_run_tensorqtl.py'),
                     '--qtl-dir', staging,
                     '--output-dir', results_dir,
                     '--ancestry', anc,
                     '--modality', 'expression',
                     '--covariates-file', cov_file,
                     '--cis-window', str(args.cis_window),
                     '--maf-threshold', str(args.maf_threshold)],
                    f"tensorQTL cis mapping (k={k})", log_path)

            # 6. eGene count
            n_eg, n_tested = count_egenes(parquet, args.fdr)
            print(f"    k={k}: {n_eg} eGenes (of {n_tested} chr1 genes tested)")
            results.append({'ancestry': anc, 'k': k, 'n_egenes': n_eg,
                            'n_tested': n_tested})

        # ---- Select k* (max eGenes; ties -> smaller k) ----
        res_df = pd.DataFrame(results).sort_values('k').reset_index(drop=True)
        best = res_df['n_egenes'].max()
        k_star = int(res_df.loc[res_df['n_egenes'] == best, 'k'].min())
        res_df['chosen'] = res_df['k'] == k_star
        print(f"\n  Optimal k for {anc}: {k_star} ({best} eGenes)")

        tsv_path = os.path.join(work_root, f'{anc}_optimal_hcp.tsv')
        res_df.to_csv(tsv_path, sep='\t', index=False)
        print(f"  Written: {tsv_path}")

        png_path = os.path.join(work_root, f'{anc}_optimal_hcp.png')
        plot_optimization(res_df, anc, k_star, args.fdr, png_path)
        print(f"  Written: {png_path}")

        # ---- Install k* as the canonical harmonized HCP file ----
        canonical = os.path.join(args.qtl_dir, f'{anc}_hcp_factors_harmonized.tsv')
        kstar_file = os.path.join(staging, f'{anc}_hcp_k{k_star}_harmonized.tsv')
        if os.path.exists(canonical):
            backup = canonical.replace('.tsv', '.pre25a_backup.tsv')
            if not os.path.exists(backup):
                shutil.copy2(canonical, backup)
                print(f"  Backed up provisional HCP file: {backup}")
        shutil.copy2(kstar_file, canonical)
        print(f"  Installed k*={k_star} HCP solution: {canonical}")
        print(f"  Next: run 25_build_covariates.py (canonical) for {anc}, "
              f"then re-apply the manual maternal-fraction covariate removal.")


if __name__ == '__main__':
    main()
