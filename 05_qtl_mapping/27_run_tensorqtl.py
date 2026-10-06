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

  Chunked independent mode (--independent-only [--chunk-index K]):
     Skips map_cis and loads the existing {ANC}_{label}_cisqtl.parquet
     (which already carries q-values), then runs only the stepwise scan.
     With --chunk-index K (1-based, e.g. $LSB_JOBINDEX) and --chunk-size N
     (default 100), runs a single ~N-gene chunk of the significant
     phenotypes/groups and writes
     {output-dir}/independent_chunks/{ANC}_{label}/chunk_KKKK.parquet
     (an empty parquet if the chunk contains no significant phenotypes).
     The FULL cis_df is passed to cis.map_independent so that its internal
     significance threshold (max pval_beta among FDR-passing rows) is
     identical to the monolithic run — chunked results are statistically
     identical to a monolithic run (bit-identical with a fixed --seed).
     Without --chunk-index, --independent-only runs the whole stepwise
     scan from the saved parquet (useful for small modalities).
     Merge chunk outputs with 27c_merge_independent.py; submit arrays with
     28b_submit_independent.sh.

Usage:
  python3 27_run_tensorqtl.py \
      --qtl-dir <qtl_inputs dir> \
      --output-dir <qtl_results dir> \
      --ancestry EAS \
      --modality splicing \
      --cis-window 1000000 \
      --maf-threshold 0.05 \
      [--no-groups] [--independent]

  python3 27_run_tensorqtl.py \
      --qtl-dir <qtl_inputs dir> \
      --output-dir <qtl_results dir> \
      --ancestry EAS \
      --modality splicing \
      --independent-only --chunk-index 7 --chunk-size 100 [--seed 12345]
"""

import argparse
import json
import os
import sys
import numpy as np
import pandas as pd


def bh_qvalues(pvals):
    """Benjamini-Hochberg adjusted p-values (q-values), preserving input
    order. Non-finite p-values are treated as 1."""
    p = np.asarray(pvals, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
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


def run_independent_chunk(cis_mod, genotype_df, variant_df, result, phenotypes,
                          phenotypes_pos, covariates, group_s, args, anc,
                          out_label, output_dir):
    """Run cis.map_independent on one chunk of significant phenotypes/groups.

    Chunking contract (verified against tensorqtl 1.0.x): map_independent
    derives its per-phenotype significance threshold as max(pval_beta) over
    the FDR-passing rows of the cis_df it is given. The FULL cis_df is
    therefore passed unchanged, and only phenotype_df/phenotype_pos_df are
    subset to the chunk (grouped mode: all member phenotypes of the chunk's
    significant groups, via group_s). Per-phenotype/group computations are
    independent, so chunked results are statistically identical to the
    monolithic run (bit-identical with a fixed --seed).
    """
    fdr_col = 'qval'
    signif = result[result[fdr_col] <= args.independent_fdr]
    n_sig = len(signif)
    n_chunks = int(np.ceil(n_sig / args.chunk_size)) if n_sig else 0
    start = (args.chunk_index - 1) * args.chunk_size
    chunk = signif.iloc[start:start + args.chunk_size]

    if group_s is not None:
        # Grouped mode: chunk unit = gene/group (one cis_df row per group).
        if 'group_id' in chunk.columns:
            chunk_units = chunk['group_id'].tolist()
        else:  # defensive: older outputs may carry the group id as the index
            chunk_units = chunk.index.tolist()
        member = set(chunk_units)
        pheno_ids = [p for p in phenotypes.index
                     if group_s.get(p, None) in member]
        unit_label = 'groups'
    else:
        # Ungrouped mode: chunk unit = phenotype (cis_df index).
        chunk_units = chunk.index.tolist()
        pheno_ids = [p for p in phenotypes.index if p in set(chunk_units)]
        unit_label = 'phenotypes'

    chunk_dir = os.path.join(output_dir, 'independent_chunks',
                             f'{anc}_{out_label}')
    os.makedirs(chunk_dir, exist_ok=True)
    chunk_path = os.path.join(
        chunk_dir, f'chunk_{args.chunk_index:04d}.parquet')

    print(f"  Chunk {args.chunk_index} of ~{n_chunks}: "
          f"{len(chunk_units)}/{n_sig} significant {unit_label} "
          f"({len(pheno_ids)} member phenotypes), "
          f"chunk size {args.chunk_size}")

    # Sidecar metadata: the merge (27c) validates chunk_size/fdr of every
    # chunk against its own parameters, so stale chunks from a run with
    # different chunking parameters are caught instead of silently mixed.
    meta = {"chunk_index": args.chunk_index, "chunk_size": args.chunk_size,
            "independent_fdr": args.independent_fdr, "seed": args.seed,
            "n_significant_units": n_sig,
            "n_units_in_chunk": len(chunk_units)}
    with open(chunk_path.replace('.parquet', '.json'), 'w') as _f:
        json.dump(meta, _f, indent=2)

    if len(pheno_ids) == 0:
        # Chunk beyond the significant set (or nothing significant at all):
        # write an empty parquet so the merge can verify completeness.
        pd.DataFrame().to_parquet(chunk_path)
        print(f"  No significant {unit_label} in this chunk — wrote empty "
              f"{chunk_path}")
        return

    pheno_sub = phenotypes.loc[pheno_ids]
    pos_sub = phenotypes_pos.loc[pheno_ids]

    try:
        ind_result = cis_mod.map_independent(
            genotype_df=genotype_df,
            variant_df=variant_df,
            cis_df=result,
            phenotype_df=pheno_sub,
            phenotype_pos_df=pos_sub,
            covariates_df=covariates,
            group_s=group_s,
            maf_threshold=args.maf_threshold,
            fdr=args.independent_fdr,
            fdr_col=fdr_col,
            window=args.cis_window,
            seed=args.seed,
            verbose=True,
        )
    except ValueError as e:
        if "No significant phenotypes" in str(e):
            print(f"  WARNING: {e}")
            ind_result = pd.DataFrame()
        else:
            raise

    ind_result.to_parquet(chunk_path)
    print(f"  Written: {chunk_path} ({len(ind_result)} independent "
          f"associations)")


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
    parser.add_argument("--independent-only", action="store_true",
                        help="Skip map_cis and load the existing "
                             "{ANC}_{label}_cisqtl.parquet from --output-dir (must "
                             "already carry q-values), then run only the stepwise "
                             "independent scan. With --chunk-index, runs one chunk "
                             "of the scan (LSF array-task mode).")
    parser.add_argument("--chunk-index", type=int, default=None,
                        help="1-based chunk index for a chunked independent scan "
                             "(requires --independent-only; typically $LSB_JOBINDEX "
                             "from an LSF job array)")
    parser.add_argument("--chunk-size", type=int, default=100,
                        help="Significant phenotypes (ungrouped) or genes/groups "
                             "(grouped) per chunk (default: 100)")
    parser.add_argument("--seed", type=int, default=None,
                        help="Seed for the tensorQTL permutation RNG (passed to "
                             "map_cis and/or map_independent). Default: None "
                             "(unseeded, as before). Set a fixed seed for "
                             "bit-reproducible chunked independent scans.")
    parser.add_argument("--qvalue-method", choices=["storey", "bh"], default="storey",
                        help="Q-value method on pval_beta. 'storey' (default, GTEx "
                             "convention): R qvalue package via compute_qvalues.R "
                             "(Rscript bridge; QVALUE_RSCRIPT env var sets the R "
                             "command). 'bh': Benjamini-Hochberg escape hatch, "
                             "no R dependency.")
    parser.add_argument("--covariates-file", default=None,
                        help="Covariates TSV (tensorQTL orientation). "
                             "Default: {qtl-dir}/{ancestry}_covariates_{modality}.tsv "
                             "(per-modality optimized set, module 25b); if absent, "
                             "falls back to {qtl-dir}/{ancestry}_covariates.tsv "
                             "with a warning")
    args = parser.parse_args()

    if args.chunk_index is not None:
        if not args.independent_only:
            sys.exit("ERROR: --chunk-index requires --independent-only "
                     "(chunked independent scans run from the saved map_cis "
                     "parquet).")
        if args.chunk_index < 1:
            sys.exit("ERROR: --chunk-index is 1-based (LSF arrays pass "
                     "$LSB_JOBINDEX).")

    anc = args.ancestry
    mod = args.modality
    qtl_dir = args.qtl_dir
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # ---- File paths ----
    phenotype_bed = os.path.join(qtl_dir, f"{anc}_{mod}.bed.gz")
    groups_path = os.path.join(qtl_dir, f"{anc}_{mod}.phenotype_groups.txt")
    plink_prefix = os.path.join(qtl_dir, f"{anc}_qtl")
    if args.covariates_file:
        covariates_path = args.covariates_file
    else:
        # Per-modality optimized covariates (module 25b) are the default;
        # fall back to the shared file with a loud warning so the log
        # records exactly which covariate set was used.
        per_mod_path = os.path.join(qtl_dir, f"{anc}_covariates_{mod}.tsv")
        shared_path = os.path.join(qtl_dir, f"{anc}_covariates.tsv")
        if os.path.exists(per_mod_path):
            covariates_path = per_mod_path
        else:
            covariates_path = shared_path
            print(f"  WARN: per-modality covariates not found: {per_mod_path}")
            print(f"        falling back to shared covariates: {shared_path}")
            print(f"        (run 25b + the canonical per-modality "
                  f"25_build_covariates.py to generate the optimized set)")

    for f in [phenotype_bed, f"{plink_prefix}.pgen", covariates_path]:
        if not os.path.exists(f):
            sys.exit(f"ERROR: required file not found: {f}")

    print(f"[{anc} / {mod}] tensorQTL cis-xQTL mapping")
    print(f"  Phenotype: {phenotype_bed}")
    print(f"  Genotype:  {plink_prefix}.pgen")
    print(f"  Covariates: {covariates_path}")
    print(f"  cis window: ±{args.cis_window // 1000} kb")
    print(f"  MAF threshold: {args.maf_threshold}")
    if args.independent_only:
        print(f"  Mode: independent-only (map_cis skipped; results loaded "
              f"from parquet)" +
              (f" chunk {args.chunk_index} (size {args.chunk_size})"
               if args.chunk_index is not None else ""))

    # ---- Import tensorQTL ----
    import tensorqtl
    from tensorqtl import genotypeio, cis
    import torch
    _cuda = torch.cuda.is_available()
    print(f"  tensorqtl {tensorqtl.__version__}, torch {torch.__version__}, "
          f"CUDA available: {_cuda}" +
          (f" ({torch.cuda.get_device_name(0)})" if _cuda else " (CPU mode)"))

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

    # cis.map_cis/map_independent require covariates_df.index to exactly equal
    # phenotype_df.columns (same samples, same order), so align all three
    # inputs to the common sample set.
    common_ordered = [s for s in phenotypes.columns if s in common]
    phenotypes = phenotypes[common_ordered]
    genotype_df = genotype_df[common_ordered]
    covariates = covariates.loc[common_ordered]

    summary_path = os.path.join(output_dir, f"{anc}_{out_label}_cisqtl.parquet")
    top_path = os.path.join(output_dir, f"{anc}_{out_label}_cisqtl_top.tsv")

    if args.independent_only:
        # ---- Load existing map_cis results (skip the permutation scan) ----
        if not os.path.exists(summary_path):
            sys.exit(f"ERROR: --independent-only requires the map_cis results "
                     f"file, not found: {summary_path}\n"
                     f"  Run the map_cis job first (27_run_tensorqtl.sh without "
                     f"--independent, via 28_submit_modalities.sh), or let "
                     f"28b_submit_independent.sh chain it automatically.")
        print(f"\n  Loading map_cis results: {summary_path}")
        result = pd.read_parquet(summary_path)
        for col in ['qval', 'pval_beta']:
            if col not in result.columns:
                sys.exit(f"ERROR: {summary_path} lacks a '{col}' column — "
                         f"cannot gate the independent scan. Re-run the "
                         f"map_cis step with the current 27_run_tensorqtl.py.")
        n_sig = int((result['qval'] <= args.independent_fdr).sum())
        print(f"    {len(result)} phenotypes/groups tested; {n_sig} at "
              f"FDR <= {args.independent_fdr}")
    else:
        # ---- Run cis-xQTL mapping ----
        print(f"\n  Running cis.map_cis...")

        result = cis.map_cis(
            genotype_df=genotype_df,
            variant_df=variant_df,
            phenotype_df=phenotypes,
            phenotype_pos_df=phenotypes_pos,
            covariates_df=covariates,
            group_s=group_s,
            window=args.cis_window,
            maf_threshold=args.maf_threshold,
            seed=args.seed,
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
            rscript_cmd = os.environ.get(
                "QVALUE_RSCRIPT",
                os.path.join(os.path.dirname(os.path.dirname(
                    os.path.abspath(__file__))), "bin", "Rscript_sif"))
            bridge = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "compute_qvalues.R")
            if not os.path.exists(bridge):
                sys.exit(f"ERROR: compute_qvalues.R not found next to {__file__} — "
                         "deploy it from gtex_conventions_xqtl_scripts.zip")
            # tensorQTL's beta approximation can occasionally return NaN/Inf
            # for an otherwise completed permutation test (e.g. a numerically
            # degenerate phenotype/group at high covariate counts). Do not let a
            # single invalid p-value abort the whole mapping run. For q-value
            # estimation only, replace non-finite values with p=1 so the full
            # number of tested hypotheses remains in Storey's calculation; then
            # force those rows to q=1 below. The original pval_beta values are
            # retained in the written tensorQTL output for QC/auditability.
            pval_beta = pd.to_numeric(result["pval_beta"], errors="coerce").to_numpy(dtype=float)
            nonfinite_mask = ~np.isfinite(pval_beta)
            if nonfinite_mask.any():
                bad_ids = result.index[nonfinite_mask].astype(str).tolist()
                preview = ", ".join(bad_ids[:10])
                suffix = " ..." if len(bad_ids) > 10 else ""
                print(f"  WARNING: {nonfinite_mask.sum()} non-finite pval_beta "
                      f"value(s); using p=1 for Storey q-value estimation and "
                      f"forcing qval=1 for those rows. IDs: {preview}{suffix}")
            qvalue_input = pval_beta.copy()
            qvalue_input[nonfinite_mask] = 1.0

            with tempfile.TemporaryDirectory() as tmpd:
                in_tsv = os.path.join(tmpd, "pval_beta.tsv")
                out_tsv = os.path.join(tmpd, "qval.tsv")
                pd.DataFrame({"pval_beta": qvalue_input}).to_csv(
                    in_tsv, sep="\t", index=False)
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
                if nonfinite_mask.any():
                    result.loc[nonfinite_mask, "qval"] = 1.0
            print(f"\n  Storey q-values (on pval_beta, GTEx convention via R qvalue): "
                  f"{(result['qval'] <= args.independent_fdr).sum()} "
                  f"phenotypes/groups at FDR <= {args.independent_fdr}")
        else:  # bh
            result['qval'] = bh_qvalues(result['pval_beta'].values)
            print(f"\n  BH q-values (on pval_beta): "
                  f"{(result['qval'] <= args.independent_fdr).sum()} "
                  f"phenotypes/groups at FDR <= {args.independent_fdr}")

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
        top.to_csv(top_path, sep='\t', index=False)
        print(f"  Written: {top_path} ({len(top)} associations, sorted by {pcol})")

        # Summary stats
        if pcol is not None:
            print(f"\n  Summary ({pcol}):")
            print(f"    p < 5e-8: {(top[pcol] < 5e-8).sum()}")
            print(f"    p < 1e-5: {(top[pcol] < 1e-5).sum()}")

    # ---- Stepwise regression for conditionally independent xQTLs ----
    if args.independent or args.independent_only:
        if args.chunk_index is not None:
            # Chunked mode (LSF array task): one ~chunk-size slice of the
            # significant phenotypes/groups; writes independent_chunks/.
            print(f"\n  Running cis.map_independent chunk "
                  f"{args.chunk_index} (forward-backward stepwise; "
                  f"FDR <= {args.independent_fdr} entry threshold)...")
            run_independent_chunk(cis, genotype_df, variant_df, result,
                                  phenotypes, phenotypes_pos, covariates,
                                  group_s, args, anc, out_label, output_dir)
        else:
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
                    seed=args.seed,
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
