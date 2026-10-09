#!/usr/bin/env python3
"""39_run_nominal.py — genome-wide nominal cis-xQTL summary statistics.

Module 06 (Objective 1.6) needs full cis-window summary statistics for every
tested phenotype (coloc/colocBoost inputs); the module-05 mapping stored
lead-only map_cis results. This worker runs tensorQTL cis.map_nominal for
one ancestry x modality x chromosome, reusing the exact module-05 inputs
(harmonized BED, intersected pgen, per-modality covariates), then a merge
step concatenates chromosome shards into a bgzipped, tabix-indexed TSV for
locus slicing.

This is the nominal-mode extension of the module-05 mapping worker
(27_run_tensorqtl.py), kept as a separate script so module 05 remains the
frozen record of the project as run. Loading, chromosome-name
reconciliation, covariate fallback, and sample-alignment logic mirror that
script.

Usage
-----
  # one chromosome shard (LSF array task)
  python3 39_run_nominal.py --qtl-dir ... --output-dir ... \
      --ancestry EAS --modality isoform_expression --chrom 21

  # merge chromosome shards -> {ANC}_{MOD}.nominal.tsv.gz (+ .tbi)
  python3 39_run_nominal.py --qtl-dir ... --output-dir ... \
      --ancestry EAS --modality isoform_expression --merge

Outputs
-------
  {output_dir}/nominal/{ANC}/{ANC}_{MOD}.nominal.chr{C}.parquet   (per shard)
  {output_dir}/nominal/{ANC}/{ANC}_{MOD}.nominal.top.tsv         (per-phenotype lead)
  {output_dir}/nominal/{ANC}/{ANC}_{MOD}.nominal.tsv.gz + .tbi   (merged, sorted)
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def load_inputs(qtl_dir, anc, mod, covariates_file=None):
    """Load BED/pgen/covariates mirroring 27_run_tensorqtl.py."""
    import tensorqtl
    from tensorqtl import genotypeio

    phenotype_bed = os.path.join(qtl_dir, f"{anc}_{mod}.bed.gz")
    plink_prefix = os.path.join(qtl_dir, f"{anc}_qtl")
    if covariates_file is None:
        covariates_path = os.path.join(qtl_dir, f"{anc}_covariates_{mod}.tsv")
        if not os.path.exists(covariates_path):
            fallback = os.path.join(qtl_dir, f"{anc}_covariates.tsv")
            print(f"  WARNING: {covariates_path} not found; falling back to "
                  f"shared covariates {fallback}")
            covariates_path = fallback
    else:
        covariates_path = covariates_file.format(ANC=anc, MOD=mod)
    for f in (phenotype_bed, plink_prefix + ".pgen", covariates_path):
        if not os.path.exists(f):
            sys.exit(f"ERROR: required file not found: {f}")

    phenotypes, phenotypes_pos = tensorqtl.read_phenotype_bed(phenotype_bed)
    genotype_df, variant_df = genotypeio.load_genotypes(plink_prefix)
    covariates = pd.read_csv(covariates_path, sep="\t", index_col=0).T

    # reconcile chr prefix (same logic as 27_run_tensorqtl.py)
    pheno_has_chr = phenotypes_pos["chr"].astype(str).str.startswith("chr").any()
    var_has_chr = variant_df["chrom"].astype(str).str.startswith("chr").any()
    if pheno_has_chr and not var_has_chr:
        phenotypes_pos = phenotypes_pos.copy()
        phenotypes_pos["chr"] = phenotypes_pos["chr"].str.replace(
            "^chr", "", regex=True)
    elif var_has_chr and not pheno_has_chr:
        variant_df = variant_df.copy()
        variant_df["chrom"] = "chr" + variant_df["chrom"].astype(str)

    common = (set(phenotypes.columns) & set(genotype_df.columns)
              & set(covariates.index))
    if len(common) < 10:
        sys.exit(f"ERROR: too few common samples ({len(common)})")
    common_ordered = [s for s in phenotypes.columns if s in common]
    return (phenotypes[common_ordered], phenotypes_pos,
            genotype_df[common_ordered], variant_df,
            covariates.loc[common_ordered])


def run_chromosome(args):
    import tensorqtl
    from tensorqtl import cis
    import torch

    (phenotypes, phenotypes_pos, genotype_df, variant_df,
     covariates) = load_inputs(args.qtl_dir, args.ancestry, args.modality,
                              args.covariates_file)
    chrom = str(args.chrom)
    mask = phenotypes_pos["chr"].astype(str) == chrom
    if mask.sum() == 0:
        sys.exit(f"ERROR: no phenotypes on chr{chrom}")
    phenotypes = phenotypes.loc[mask]
    phenotypes_pos = phenotypes_pos.loc[mask]
    print(f"[{args.ancestry}/{args.modality}/chr{chrom}] "
          f"{phenotypes.shape[0]} phenotypes x {phenotypes.shape[1]} samples; "
          f"CUDA: {torch.cuda.is_available()}")

    out_dir = Path(args.output_dir) / "nominal" / args.ancestry
    out_dir.mkdir(parents=True, exist_ok=True)
    shard_path = out_dir / f"{args.ancestry}_{args.modality}.nominal.chr{chrom}.parquet"
    top_shard = out_dir / f"{args.ancestry}_{args.modality}.nominal.chr{chrom}.top.tsv"
    # both outputs must exist — a crash between the two writes (observed
    # with the top-assoc filename mismatch) must not look "done"
    if shard_path.exists() and top_shard.exists() and not args.force:
        print(f"  shard exists: {shard_path.name} (use --force to rerun)")
        return

    tmp_dir = out_dir / f"tmp_{args.modality}_chr{chrom}"
    if tmp_dir.exists():
        # clear leftovers from a crashed run so the globs below only see
        # this run's files
        for stale in tmp_dir.iterdir():
            stale.unlink()
    tmp_dir.mkdir(exist_ok=True)
    cis.map_nominal(
        genotype_df=genotype_df,
        variant_df=variant_df,
        phenotype_df=phenotypes,
        phenotype_pos_df=phenotypes_pos,
        prefix=f"{args.ancestry}_{args.modality}",
        covariates_df=covariates,
        maf_threshold=args.maf_threshold,
        window=args.cis_window,
        run_eigenmt=False,
        output_dir=str(tmp_dir),
        write_top=True,
        write_stats=True,
        verbose=True,
    )
    # tensorQTL writes {prefix}.cis_qtl_pairs.{chr}.parquet; the
    # cis_qtl_top_assoc file is only written in interaction mode, so derive
    # per-phenotype leads from the pairs table (min nominal p per
    # phenotype) — version-independent
    pair_files = list(tmp_dir.glob("*.cis_qtl_pairs.*.parquet"))
    if len(pair_files) != 1:
        sys.exit(f"ERROR: expected 1 pairs parquet in {tmp_dir}, "
                 f"found {len(pair_files)}")
    pairs = pd.read_parquet(pair_files[0])
    p = pairs.dropna(subset=["pval_nominal"])
    top = p.loc[p.groupby("phenotype_id")["pval_nominal"].idxmin()]
    pairs.to_parquet(shard_path)
    top.to_csv(top_shard, sep="\t", index=False)
    for f in tmp_dir.iterdir():
        f.unlink()
    tmp_dir.rmdir()
    print(f"  wrote {len(pairs)} pairs -> {shard_path.name}")


def merge_shards(args):
    out_dir = Path(args.output_dir) / "nominal" / args.ancestry
    shards = sorted(out_dir.glob(f"{args.ancestry}_{args.modality}.nominal.chr*.parquet"),
                    key=lambda p: int(p.stem.split(".chr")[1]))
    if not shards:
        sys.exit(f"ERROR: no shards found in {out_dir} for "
                 f"{args.ancestry}_{args.modality}")
    expect = set(str(c) for c in range(1, 23))
    have = {p.stem.split(".chr")[1] for p in shards}
    missing = expect - have
    if missing and not args.allow_missing_chroms:
        sys.exit(f"ERROR: missing chromosome shards: {sorted(missing)} "
                 f"(pass --allow-missing-chroms to merge anyway)")
    print(f"[merge] {args.ancestry}_{args.modality}: {len(shards)} shards")
    df = pd.concat([pd.read_parquet(p) for p in shards], ignore_index=True)
    # attach variant coordinates for tabix/locus slicing
    from tensorqtl import genotypeio
    _, variant_df = genotypeio.load_genotypes(
        os.path.join(args.qtl_dir, f"{args.ancestry}_qtl"))
    vmap = variant_df.reset_index()[["index", "chrom", "pos"]].rename(
        columns={"index": "variant_id"})
    df = df.merge(vmap, on="variant_id", how="left")
    df["chrom"] = df["chrom"].astype(str).str.replace("^chr", "", regex=True)
    df = df.sort_values(["chrom", "pos"],
                        key=lambda s: pd.to_numeric(s, errors="coerce"))
    # tabix needs the sequence and position as the first two columns
    lead = ["chrom", "pos"]
    df = df[lead + [c for c in df.columns if c not in lead]]
    out_tsv = out_dir / f"{args.ancestry}_{args.modality}.nominal.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    # bgzip + tabix for locus slicing
    subprocess.run(f"bgzip -f {out_tsv} && tabix -f -s 1 -b 2 -e 2 "
                   f"-S 1 {out_tsv}.gz", shell=True, check=True)
    # merged per-phenotype top table
    tops = sorted(out_dir.glob(
        f"{args.ancestry}_{args.modality}.nominal.chr*.top.tsv"))
    top_df = pd.concat([pd.read_csv(p, sep="\t") for p in tops],
                       ignore_index=True)
    top_df.to_csv(out_dir / f"{args.ancestry}_{args.modality}.nominal.top.tsv",
                  sep="\t", index=False)
    print(f"  {len(df)} pairs -> {out_tsv}.gz (+ .tbi); "
          f"{len(top_df)} phenotype leads")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--ancestry", required=True)
    p.add_argument("--modality", required=True)
    p.add_argument("--chrom", default=None, help="chromosome shard to run")
    p.add_argument("--merge", action="store_true",
                   help="merge chromosome shards instead of mapping")
    p.add_argument("--cis-window", type=int, default=1_000_000)
    p.add_argument("--maf-threshold", type=float, default=0.01)
    p.add_argument("--covariates-file", default=None,
                   help="override; supports {ANC} and {MOD} placeholders")
    p.add_argument("--allow-missing-chroms", action="store_true")
    p.add_argument("--force", action="store_true")
    args = p.parse_args()
    if args.merge:
        merge_shards(args)
    else:
        if args.chrom is None:
            sys.exit("ERROR: --chrom required unless --merge")
        run_chromosome(args)


if __name__ == "__main__":
    main()
