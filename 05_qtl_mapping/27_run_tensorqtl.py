#!/usr/bin/env python3
"""
27_run_tensorqtl.py — cis-xQTL mapping with tensorQTL (any modality)

Per ancestry × modality:
  1. Loads phenotype BED (bgzipped + tabix-indexed, harmonized to array_id).
  2. Loads genotype (plink2 pgen via genotypeio.load_genotypes).
  3. Loads covariates (TSV in tensorQTL format; shared across modalities).
  4. If {ANC}_{modality}.phenotype_groups.txt exists, loads it and passes
     group_s to cis.map_cis → grouped testing (one lead variant per gene),
     following the PANTRY/Pheast convention. Without a groups file, runs
     ungrouped (one lead variant per phenotype; e.g. expression).
  5. Writes summary statistics (parquet) + top variants (TSV).
  6. Q-values on pval_beta per --qvalue-method: 'storey' (default, GTEx
     convention — R qvalue package via the compute_qvalues.R Rscript
     bridge; QVALUE_RSCRIPT env var sets the R command) or 'bh'
     (Benjamini-Hochberg escape hatch, no R dependency).
     With --independent (PANTRY-style stepwise regression): runs
     cis.map_independent → conditionally independent cis-xQTLs per
     phenotype/group with a 'rank' column. Writes
     {ANC}_{label}_cisqtl_independent.parquet + _top.tsv.
     If no phenotype/group passes the FDR threshold, empty independent
     outputs are written (exit 0).

Usage:
  python3 27_run_tensorqtl.py \
      --qtl-dir <qtl_inputs dir> \
      --output-dir <qtl_results dir> \
      --ancestry EAS \
      --modality splicing \
      --cis-window 1000000 \
      --maf-threshold 0.05 \
      [--no-groups] [--independent]
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd


def bh_qvalues(pvals):
    """Benjamini-Hochberg adjusted p-values (q-values), preserving input
    order. NaN p-values are treated as 1."""
    p = np.asarray(pvals, dtype=float)
    p = np.where(np.isnan(p), 1.0, p)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p, kind='stable')
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]  # enforce monotonicity
    q = np.clip(q, 0.0, 1.0)
    out = np.empty(n)
    out[order] = q
    return out


def load_groups(groups_path, phenotype_ids):
    """Load a PANTRY phenotype_groups.txt file (no header; phenotype_id,
    gene_id) and return a Series mapping phenotype_id -> group_id, indexed
    by phenotype_id. Verifies coverage of all phenotypes in the BED."""
    groups_df = pd.read_csv(groups_path, sep='\t', header=None,
                            names=['phenotype_id', 'group_id'])
    group_s = pd.Series(groups_df['group_id'].values,
                        index=groups_df['phenotype_id'])
    # Drop duplicate phenotype entries if any (keep first)
    group_s = group_s[~group_s.index.duplicated(keep='first')]
    missing = set(phenotype_ids) - set(group_s.index)
    if missing:
        sys.exit(f"ERROR: groups file {groups_path} does not cover "
                 f"{len(missing)} phenotypes in the BED "
                 f"(e.g. {sorted(missing)[:3]}). Regenerate or delete the "
                 f"groups file to run ungrouped.")
    n_groups = group_s.loc[list(phenotype_ids)].nunique()
    print(f"    Groups: {len(group_s)} phenotypes -> {n_groups} genes "
          f"(grouped testing via group_s)")
    return group_s


def main():
    parser = argparse.ArgumentParser(description="Run tensorQTL cis-xQTL mapping")
    parser.add_argument("--qtl-dir", required=True, help="QTL inputs directory")
    parser.add_argument("--output-dir", required=True, help="Output directory for results")
    parser.add_argument("--ancestry", required=True, help="Ancestry label (e.g. EAS)")
    parser.add_argument("--modality", default="expression",
                        help="Modality label (default: expression). Determines "
                             "phenotype BED ({ANC}_{modality}.bed.gz), optional groups "
                             "file ({ANC}_{modality}.phenotype_groups.txt), and output names.")
    parser.add_argument("--cis-window", type=int, default=1000000,
                        help="cis window in bp (default: 1000000 = ±1 Mb)")
    parser.add_argument("--maf-threshold", type=float, default=0.01,
                        help="Minimum MAF (default: 0.01)")
    parser.add_argument("--no-groups", action="store_true",
                        help="Ignore phenotype_groups.txt and run ungrouped "
                             "(per-phenotype lead variants, e.g. transcript-level "
                             "for isoforms); outputs get an '_ungrouped' suffix")
    parser.add_argument("--independent", action="store_true",
                        help="After map_cis, run cis.map_independent (PANTRY-style "
                             "forward-backward stepwise regression) for conditionally "
                             "independent cis-xQTLs; writes "
                             "{ANC}_{label}_cisqtl_independent.parquet/_top.tsv")
    parser.add_argument("--independent-fdr", type=float, default=0.05,
                        help="FDR threshold (q-value on pval_beta) for a "
                             "phenotype/group to enter stepwise regression "
                             "(default: 0.05)")
    parser.add_argument("--qvalue-method", choices=["storey", "bh"], default="storey",
                        help="Q-value method on pval_beta. 'storey' (default, GTEx "
                             "convention): R qvalue package via compute_qvalues.R "
                             "(Rscript bridge; QVALUE_RSCRIPT env var sets the R "
                             "command). 'bh': Benjamini-Hochberg escape hatch, "
                             "no R dependency.")
    parser.add_argument("--covariates-file", default=None,
                        help="Covariates TSV (tensorQTL orientation). "
                             "Default: {qtl-dir}/{ancestry}_covariates.tsv")
    args = parser.parse_args()

    anc = args.ancestry
    mod = args.modality
    qtl_dir = args.qtl_dir
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # ---- File paths ----
    phenotype_bed = os.path.join(qtl_dir, f"{anc}_{mod}.bed.gz")
    groups_path = os.path.join(qtl_dir, f"{anc}_{mod}.phenotype_groups.txt")
    plink_prefix = os.path.join(qtl_dir, f"{anc}_qtl")
    covariates_path = args.covariates_file or os.path.join(qtl_dir, f"{anc}_covariates.tsv")

    for f in [phenotype_bed, f"{plink_prefix}.pgen", covariates_path]:
        if not os.path.exists(f):
            sys.exit(f"ERROR: required file not found: {f}")

    print(f"[{anc} / {mod}] tensorQTL cis-xQTL mapping")
    print(f"  Phenotype: {phenotype_bed}")
    print(f"  Genotype:  {plink_prefix}.pgen")
    print(f"  Covariates: {covariates_path}")
    print(f"  cis window: ±{args.cis_window // 1000} kb")
    print(f"  MAF threshold: {args.maf_threshold}")

    # ---- Import tensorQTL ----
    import tensorqtl
    from tensorqtl import genotypeio, cis

    # ---- Load phenotype ----
    print(f"\n  Loading phenotype BED...")
    phenotypes, phenotypes_pos = tensorqtl.read_phenotype_bed(phenotype_bed)
    print(f"    {phenotypes.shape[0]} phenotypes x {phenotypes.shape[1]} samples")
    print(f"    Phenotype positions: {phenotypes_pos.shape[0]}")

    # ---- Load phenotype groups (unless --no-groups) ----
    group_s = None
    if args.no_groups:
        print(f"  --no-groups: running ungrouped (per-phenotype lead variants)")
    elif os.path.exists(groups_path):
        print(f"  Loading phenotype groups: {groups_path}")
        group_s = load_groups(groups_path, phenotypes.index)
    else:
        print(f"  No groups file ({groups_path}) — running ungrouped")
    out_label = f"{mod}_ungrouped" if args.no_groups else mod

    # ---- Load genotype ----
    print(f"  Loading genotype...")
    genotype_df, variant_df = genotypeio.load_genotypes(plink_prefix)
    print(f"    {genotype_df.shape[0]} variants x {genotype_df.shape[1]} samples")

    # ---- Reconcile chromosome naming (phenotype BED may use 'chr1', genotype .pvar may use '1') ----
    pheno_has_chr_prefix = phenotypes_pos['chr'].astype(str).str.startswith('chr').any()
    variant_has_chr_prefix = variant_df['chrom'].astype(str).str.startswith('chr').any()
    if pheno_has_chr_prefix and not variant_has_chr_prefix:
        phenotypes_pos = phenotypes_pos.copy()
        phenotypes_pos['chr'] = phenotypes_pos['chr'].str.replace('^chr', '', regex=True)
        print("  Normalized phenotype chromosome naming: stripped 'chr' prefix to match genotype convention")
    elif variant_has_chr_prefix and not pheno_has_chr_prefix:
        variant_df = variant_df.copy()
        variant_df['chrom'] = 'chr' + variant_df['chrom'].astype(str)
        print("  Normalized genotype chromosome naming: added 'chr' prefix to match phenotype convention")

    # ---- Load covariates ----
    print(f"  Loading covariates...")
    covariates = pd.read_csv(covariates_path, sep='\t', index_col=0).T
    # covariates is now samples x covariates
    print(f"    {covariates.shape[0]} samples x {covariates.shape[1]} covariates")

    # ---- Verify sample overlap ----
    pheno_samples = set(phenotypes.columns)
    geno_samples = set(genotype_df.columns)
    cov_samples = set(covariates.index)
    common = pheno_samples & geno_samples & cov_samples
    print(f"\n  Sample overlap:")
    print(f"    Phenotype: {len(pheno_samples)}")
    print(f"    Genotype:  {len(geno_samples)}")
    print(f"    Covariates: {len(cov_samples)}")
    print(f"    Common:    {len(common)}")

    if len(common) < 10:
        sys.exit(f"ERROR: Too few common samples ({len(common)}) for QTL mapping")

    # ---- Run cis-xQTL mapping ----
    print(f"\n  Running cis.map_cis...")

    # cis.map_cis requires covariates_df.index to exactly equal phenotype_df.columns
    # (same samples, same order), so align all three inputs to the common sample set.
    common_ordered = [s for s in phenotypes.columns if s in common]
    phenotypes = phenotypes[common_ordered]
    genotype_df = genotype_df[common_ordered]
    covariates = covariates.loc[common_ordered]

    result = cis.map_cis(
        genotype_df=genotype_df,
        variant_df=variant_df,
        phenotype_df=phenotypes,
        phenotype_pos_df=phenotypes_pos,
        covariates_df=covariates,
        group_s=group_s,
        window=args.cis_window,
        maf_threshold=args.maf_threshold,
        verbose=True
    )

    # ---- Write results ----
    # map_cis returns one row per phenotype (ungrouped) or per group (grouped
    # via group_s); the phenotype/group IDs are carried in the DataFrame index.
    if not isinstance(result, pd.DataFrame):
        # Defensive: older tensorQTL returned a dict of DataFrames per phenotype
        all_results = []
        for pheno, df in result.items():
            df['phenotype_id'] = pheno
            all_results.append(df)
        result = pd.concat(all_results, ignore_index=True)

    # Q-values on the permutation-calibrated p-value (pval_beta). Added
    # before writing so the map_cis parquet carries them; required as the
    # fdr_col input to cis.map_independent in --independent mode.
    #   storey (default, GTEx convention): tensorQTL's calculate_qvalues,
    #     which calls R's qvalue package via rpy2 (lambda estimated from data).
    #   bh: Benjamini-Hochberg escape hatch (no rpy2 dependency).
    if 'pval_beta' not in result.columns:
        if args.independent:
            sys.exit("ERROR: --independent requires a 'pval_beta' column in the "
                     "map_cis output (permutation-based correction); not found.")
    elif args.qvalue_method == 'storey':
        # GTEx convention: Storey q-values from the R qvalue package, via a
        # file-based Rscript bridge (compute_qvalues.R, alongside this
        # script). Statistically identical to tensorQTL's calculate_qvalues
        # (itself an rpy2 wrapper around qvalue::qvalue), but avoids the
        # rpy2/conda-R stack, which hits libstdc++/GLIBCXX conflicts on
        # seadragon. QVALUE_RSCRIPT sets the Rscript command;
        # 27_run_tensorqtl.sh points it at the singularity R container.
        import subprocess
        import tempfile
        rscript_cmd = os.environ.get("QVALUE_RSCRIPT", "Rscript")
        bridge = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "compute_qvalues.R")
        if not os.path.exists(bridge):
            sys.exit(f"ERROR: compute_qvalues.R not found next to {__file__} — "
                     "deploy it from gtex_conventions_xqtl_scripts.zip")
        with tempfile.TemporaryDirectory() as tmpd:
            in_tsv = os.path.join(tmpd, "pval_beta.tsv")
            out_tsv = os.path.join(tmpd, "qval.tsv")
            result[["pval_beta"]].to_csv(in_tsv, sep="\t", index=False)
            proc = subprocess.run(
                f'{rscript_cmd} "{bridge}" "{in_tsv}" "{out_tsv}"',
                shell=True, capture_output=True, text=True)
            if proc.returncode != 0 or not os.path.exists(out_tsv):
                sys.exit(
                    "ERROR: Storey q-value computation failed (Rscript bridge).\n"
                    f"  command: {rscript_cmd} {bridge} <in> <out>\n"
                    f"  stderr tail: {proc.stderr[-1500:]}\n"
                    "  Fix: install qvalue into the R library used by "
                    "QVALUE_RSCRIPT (rerun 21_install_tensorqtl.sh), or use "
                    "--qvalue-method bh (not the GTEx convention).")
            for line in proc.stdout.splitlines():
                if line.strip():
                    print(f"  {line.strip()}")
            qdf = pd.read_csv(out_tsv, sep="\t")
            if len(qdf) != len(result) or 'qval' not in qdf.columns:
                sys.exit("ERROR: compute_qvalues.R output malformed "
                         f"({len(qdf)} rows vs {len(result)} expected)")
            result["qval"] = qdf["qval"].values
        print(f"\n  Storey q-values (on pval_beta, GTEx convention via R qvalue): "
              f"{(result['qval'] <= args.independent_fdr).sum()} "
              f"phenotypes/groups at FDR <= {args.independent_fdr}")
    else:  # bh
        result['qval'] = bh_qvalues(result['pval_beta'].values)
        print(f"\n  BH q-values (on pval_beta): "
              f"{(result['qval'] <= args.independent_fdr).sum()} "
              f"phenotypes/groups at FDR <= {args.independent_fdr}")

    summary_path = os.path.join(output_dir, f"{anc}_{out_label}_cisqtl.parquet")
    result.to_parquet(summary_path)
    print(f"\n  Written: {summary_path} ({len(result)} associations)")

    # Top table: same associations sorted by nominal p-value, with the
    # phenotype/group IDs kept as a column (not lost in the index).
    top = result.copy()
    if not isinstance(top.index, pd.RangeIndex):
        top = top.reset_index()
    pcol = next((c for c in ['pval_nominal', 'pval'] if c in top.columns), None)
    if pcol is not None:
        top = top.sort_values(pcol)
    top_path = os.path.join(output_dir, f"{anc}_{out_label}_cisqtl_top.tsv")
    top.to_csv(top_path, sep='\t', index=False)
    print(f"  Written: {top_path} ({len(top)} associations, sorted by {pcol})")

    # Summary stats
    if pcol is not None:
        print(f"\n  Summary ({pcol}):")
        print(f"    p < 5e-8: {(top[pcol] < 5e-8).sum()}")
        print(f"    p < 1e-5: {(top[pcol] < 1e-5).sum()}")

    # ---- Stepwise regression for conditionally independent xQTLs ----
    if args.independent:
        print(f"\n  Running cis.map_independent (forward-backward stepwise; "
              f"FDR <= {args.independent_fdr} entry threshold)...")
        ind_path = os.path.join(
            output_dir, f"{anc}_{out_label}_cisqtl_independent.parquet")
        ind_top_path = os.path.join(
            output_dir, f"{anc}_{out_label}_cisqtl_independent_top.tsv")
        try:
            ind_result = cis.map_independent(
                genotype_df=genotype_df,
                variant_df=variant_df,
                cis_df=result,
                phenotype_df=phenotypes,
                phenotype_pos_df=phenotypes_pos,
                covariates_df=covariates,
                group_s=group_s,
                maf_threshold=args.maf_threshold,
                fdr=args.independent_fdr,
                fdr_col='qval',
                window=args.cis_window,
                verbose=True,
            )
        except ValueError as e:
            if "No significant phenotypes" in str(e):
                print(f"  WARNING: {e}")
                print(f"  No phenotypes/groups pass FDR <= "
                      f"{args.independent_fdr} — writing empty independent "
                      f"outputs.")
                ind_result = pd.DataFrame()
            else:
                raise

        ind_result.to_parquet(ind_path)
        print(f"  Written: {ind_path} ({len(ind_result)} independent "
              f"associations)")

        ind_top = ind_result.copy()
        if len(ind_top) and not isinstance(ind_top.index, pd.RangeIndex):
            ind_top = ind_top.reset_index()
        if len(ind_top) and 'pval_nominal' in ind_top.columns:
            ind_top = ind_top.sort_values('pval_nominal')
        ind_top.to_csv(ind_top_path, sep='\t', index=False)
        print(f"  Written: {ind_top_path} ({len(ind_top)} associations)")
        if len(ind_result) and 'rank' in ind_result.columns:
            if 'phenotype_id' in ind_result.columns:
                gid = ind_result['phenotype_id']
            elif not isinstance(ind_result.index, pd.RangeIndex):
                gid = pd.Series(ind_result.index)
            else:
                gid = None
            if gid is not None:
                max_rank = ind_result['rank'].groupby(gid.values).max()
                print(f"  Genes/groups with >1 conditionally independent "
                      f"xQTL: {(max_rank > 1).sum()}")

    print(f"\n  Done: {anc} / {mod}")


if __name__ == "__main__":
    main()
