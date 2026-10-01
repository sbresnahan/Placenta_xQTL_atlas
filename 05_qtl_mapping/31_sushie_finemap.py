#!/usr/bin/env python3
"""31_sushie_finemap.py — cross-ancestry SuSHiE fine-mapping of xQTL loci.

Implements Objective 1.5 (cross-ancestry fine-mapping): joint multi-ancestry
SuShiE fine-mapping of every phenotype with an FDR-significant lead variant
in at least one ancestry stratum, using the exact phenotype values,
covariates, and intersected genotypes used in the tensorQTL mapping.

Subcommands
-----------
prepare-loci
    Build the per-modality locus list from phenotype-level TensorQTL
    discoveries: expression uses {ANC}_expression_cisqtl_top.tsv; all other
    modalities use {ANC}_{MOD}_ungrouped_cisqtl_top.tsv. The exact tested
    cis window is reconstructed from the phenotype BED, and per-locus L
    (max causal signals) is derived from the matching conditionally
    independent result set.

run
    Fine-map a shard of loci. For each locus, writes per-ancestry
    phenotype/covariate TSVs in SuSHiE's no-header format and calls
    `sushie finemap` on per-locus VCFs exported from the intersected pgens. Phenotypes missing from an
    ancestry's harmonized BED (e.g., detection-filtered in that stratum)
    drop that ancestry for the locus; loci with no eligible ancestry are
    recorded as skipped. Per-locus diagnostics (convergence, credible-set
    count, SNP count, wall time) are written per shard.

Conventions
-----------
- Phenotype BEDs: {qtl_dir}/{ANC}_{MOD}.bed.gz (QN+INT+ComBat values keyed
  by array_id — the exact values used in mapping).
- Genotypes: {qtl_dir}/{ANC}_qtl.pgen/.pvar/.psam (intersected, MAC>=5).
- Covariates: {qtl_dir}/{ANC}_covariates_{MOD}.tsv (tensorQTL orientation:
  rows = covariates, columns = samples); transposed here for SuSHiE.
- Outputs: {out_dir}/{MOD}/{phenotype}.sushie.{weights,cs,corr}.tsv (+ .log,
  .ancestries); shard diagnostics in {out_dir}/logs/.

Aims fallback: loci that fail convergence are flagged in the diagnostics
(status column) for the documented FINEMAP / single-signal-restriction
fallback (not implemented here).
"""

import argparse
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# prepare-loci
# ---------------------------------------------------------------------------

VAR_ID_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<pos>\d+):")


def _parse_variant_positions(variant_ids: pd.Series):
    """Parse 'chr:pos:ref:alt' variant IDs into (chrom, pos) arrays."""
    chrom = variant_ids.str.extract(VAR_ID_RE)["chrom"]
    pos = variant_ids.str.extract(VAR_ID_RE)["pos"].astype(int)
    return chrom, pos


def prepare_loci(args):
    mod = args.modality
    anc_list = args.ancestries
    # Expression already has one phenotype per gene in the primary scan.
    # For all multi-phenotype modalities, fine-map every phenotype that passes
    # phenotype-level FDR in the ungrouped TensorQTL analysis rather than only
    # the single phenotype selected to represent each grouped gene-level test.
    result_label = mod if mod == "expression" else f"{mod}_ungrouped"
    print(f"[prepare-loci] modality={mod} ancestries={anc_list} "
          f"discovery={result_label}")

    # --- union of significant phenotypes ---
    sig = {}
    for anc in anc_list:
        top_path = Path(args.results_dir) / f"{anc}_{result_label}_cisqtl_top.tsv"
        if not top_path.exists():
            print(f"  WARNING: missing {top_path.name}; skipping {anc}")
            continue
        top = pd.read_csv(top_path, sep="\t", usecols=["phenotype_id", "qval"])
        hits = set(top.loc[top["qval"] <= args.qval, "phenotype_id"])
        sig[anc] = hits
        print(f"  {anc}: {len(hits)} significant phenotypes (q <= {args.qval})")
    if not sig:
        raise SystemExit(f"No top tables found for modality {mod}")
    union = sorted(set().union(*sig.values()))
    print(f"  union: {len(union)} unique phenotypes")

    # --- independent signal counts -> per-locus L ---
    n_signals = {}
    for anc in anc_list:
        ind_path = Path(args.results_dir) / f"{anc}_{result_label}_cisqtl_independent_top.tsv"
        if not ind_path.exists():
            continue
        ind = pd.read_csv(ind_path, sep="\t", usecols=["phenotype_id", "rank"])
        per = ind.groupby("phenotype_id")["rank"].max()
        for ph, r in per.items():
            n_signals[ph] = max(n_signals.get(ph, 1), int(r))

    # --- reconstruct the cis windows actually used by tensorQTL ---
    # read_phenotype_bed() converts BED start to 1-based (start + 1), then
    # tests [phenotype_start - window, phenotype_end + window].  The cis
    # permutation parquet contains only the top variant per phenotype, so it
    # cannot be used to reconstruct the tested interval.
    qtl_dir = (Path(args.qtl_dir) if args.qtl_dir is not None
               else Path(args.results_dir).parent / "qtl_inputs")
    windows = {}
    for anc in anc_list:
        bed_path = Path(args.bed_template.format(
            qtl_dir=qtl_dir, anc=anc, mod=mod))
        if not bed_path.exists():
            print(f"  WARNING: missing {bed_path.name}; skipping {anc} windows")
            continue
        bed = pd.read_csv(
            bed_path, sep="\t", compression="gzip",
            usecols=["#chr", "start", "end", "phenotype_id"],
            dtype={"#chr": str})
        bed = bed[bed["phenotype_id"].isin(union)]
        for _, row in bed.iterrows():
            ph = row["phenotype_id"]
            chrom = str(row["#chr"])
            phenotype_start = int(row["start"]) + 1
            phenotype_end = int(row["end"])
            lo = max(1, phenotype_start - args.cis_window)
            hi = phenotype_end + args.cis_window
            if ph in windows:
                old_chrom, old_lo, old_hi = windows[ph]
                if old_chrom != chrom:
                    raise ValueError(
                        f"Chromosome mismatch for {ph}: {old_chrom} vs {chrom}")
                windows[ph] = (chrom, min(old_lo, lo), max(old_hi, hi))
            else:
                windows[ph] = (chrom, lo, hi)
        print(f"  {anc}: windows for {len(bed)} phenotypes")

    rows = []
    missing_window = 0
    for ph in union:
        if ph not in windows:
            missing_window += 1
            continue
        chrom, lo, hi = windows[ph]
        ns = n_signals.get(ph, 1)
        L = min(args.l_max, max(args.l_min, ns + args.l_buffer))
        sig_in = ",".join(a for a in anc_list if ph in sig.get(a, set()))
        rows.append({"modality": mod, "phenotype_id": ph, "chrom": chrom,
                     "start": lo, "end": hi, "n_signals": ns, "L": L,
                     "sig_in": sig_in})
    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"  wrote {len(out)} loci -> {args.out}")
    if missing_window:
        print(f"  WARNING: {missing_window} significant phenotypes had no "
              f"cis window in any phenotype BED (dropped from locus list)")


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------

def sanitize_trait(phenotype_id: str) -> str:
    """Filesystem-safe locus prefix."""
    return re.sub(r"[^A-Za-z0-9._+-]", "_", phenotype_id)


def load_bed_phenotypes(bed_path: Path, wanted: set) -> pd.DataFrame:
    """Read a harmonized BED and return rows for wanted phenotype_ids."""
    df = pd.read_csv(bed_path, sep="\t", compression="gzip")
    df = df[df["phenotype_id"].isin(wanted)]
    return df


def write_sushie_inputs(bed_row: pd.Series, cov_t: pd.DataFrame,
                        pheno_path: Path, covar_path: Path):
    """Write SuSHiE no-header phenotype (2-col) and covariate TSVs.

    bed_row: one BED row (index = #chr/start/end/phenotype_id + sample IDs).
    cov_t: covariates transposed to samples x covariates (index = sample ID).
    Samples are matched to the BED header order intersected with covariates.
    """
    samples = [c for c in bed_row.index
               if c not in ("#chr", "start", "end", "phenotype_id")]
    if cov_t is not None:
        keep = [s for s in samples if s in cov_t.index]
    else:
        keep = samples
    with open(pheno_path, "w") as f:
        for s in keep:
            f.write(f"{s}\t{bed_row[s]}\n")
    if cov_t is not None:
        with open(covar_path, "w") as f:
            for s in keep:
                vals = "\t".join(f"{v:.6g}" for v in cov_t.loc[s].values)
                f.write(f"{s}\t{vals}\n")
    return len(keep)


def parse_sushie_log(log_path: Path):
    """Extract convergence diagnostics from a SuSHiE .log file."""
    diag = {"n_snps": None, "n_individuals": None, "n_cs": None,
            "converged": None, "n_iter": None, "elbo": None}
    if not log_path.exists():
        return diag
    txt = log_path.read_text(errors="replace")
    m = re.search(r"Prepare (\d+) SNPs for (\d+) individuals", txt)
    if m:
        diag["n_snps"] = int(m.group(1))
        diag["n_individuals"] = int(m.group(2))
    m = re.search(r"(\d+) out of \d+ credible sets remain", txt)
    if m:
        diag["n_cs"] = int(m.group(1))
    m = re.search(r"concludes after (\d+) iterations\. Final ELBO score: "
                  r"(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", txt)
    if m:
        diag["n_iter"] = int(m.group(1))
        diag["elbo"] = float(m.group(2))
    if "Reach minimum tolerance threshold" in txt:
        diag["converged"] = True
    elif "maximum number of iterations" in txt.lower():
        diag["converged"] = False
    return diag


def run_shard(args):
    loci = pd.read_csv(args.loci_file, sep="\t")
    out_dir = Path(args.out_dir)
    logs_dir = out_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    shard_name = Path(args.loci_file).stem

    # sanity: sanitized trait names must be unique within the shard
    safe = loci["phenotype_id"].map(sanitize_trait)
    if safe.duplicated().any():
        dupes = loci.loc[safe.duplicated(), "phenotype_id"].tolist()
        raise SystemExit(f"Sanitized phenotype_id collision in shard: {dupes[:3]}")
    loci = loci.assign(_safe=safe)

    diagnostics = []
    bed_cache = {}        # (anc, mod) -> DataFrame indexed by phenotype_id
    cov_cache = {}        # (anc, mod) -> transposed covariate DataFrame
    cov_written = set()   # (anc, mod) with covar TSV already written

    for _, locus in loci.iterrows():
        mod, ph = locus["modality"], locus["phenotype_id"]
        trait = locus["_safe"]
        t0 = time.time()
        rec = {"modality": mod, "phenotype_id": ph, "chrom": locus["chrom"],
               "start": locus["start"], "end": locus["end"], "L": locus["L"],
               "sig_in": locus.get("sig_in", ""), "status": "ok"}
        locus_dir = out_dir / mod
        locus_dir.mkdir(parents=True, exist_ok=True)
        out_prefix = locus_dir / trait
        done_path = out_prefix.with_suffix(".done")
        weights_path = Path(str(out_prefix) + ".sushie.weights.tsv")

        if done_path.exists() and not args.force:
            if weights_path.exists() and weights_path.stat().st_size > 0:
                rec["status"] = "skipped_done"
                diagnostics.append(rec)
                continue
            print(f"  WARNING: stale .done without weights for {mod}/{trait}; rerunning")
            done_path.unlink()

        # --- gather per-ancestry inputs ---
        pheno_paths, covar_paths, vcf_paths, anc_used = [], [], [], []
        for anc in args.ancestries:
            key = (anc, mod)
            if key not in bed_cache:
                bed_path = Path(args.bed_template.format(
                    qtl_dir=args.qtl_dir, anc=anc, mod=mod))
                if not bed_path.exists():
                    bed_cache[key] = None
                else:
                    wanted = set(loci.loc[loci["modality"] == mod,
                                          "phenotype_id"])
                    df = load_bed_phenotypes(bed_path, wanted)
                    bed_cache[key] = df.set_index("phenotype_id")
                    cov_path = Path(args.covar_template.format(
                        qtl_dir=args.qtl_dir, anc=anc, mod=mod))
                    if cov_path.exists():
                        cov = pd.read_csv(cov_path, sep="\t")
                        cov_cache[key] = cov.set_index("covariate_id").T
                    else:
                        cov_cache[key] = None
            bed = bed_cache[key]
            if bed is None or ph not in bed.index:
                continue  # phenotype not measured/QC'd in this ancestry
            bed_row = bed.loc[ph]
            rec.setdefault("pheno_start", int(bed_row["start"]))
            rec.setdefault("pheno_end", int(bed_row["end"]))
            pheno_path = locus_dir / f"{trait}.{anc}.pheno.tsv"
            covar_path = locus_dir / f"{trait}.{anc}.covar.tsv"
            n = write_sushie_inputs(bed_row, cov_cache[key],
                                    pheno_path, covar_path)
            if n < args.min_samples:
                continue
            pgen_prefix = args.pgen_template.format(
                qtl_dir=args.qtl_dir, anc=anc, mod=mod)
            vcf_prefix = locus_dir / f"{trait}.{anc}.geno"
            vcf_path = Path(str(vcf_prefix) + ".vcf.gz")
            export_cmd = [
                "bash", "-lc",
                "module load plink && "
                f"plink2 --pfile {pgen_prefix} "
                f"--chr {locus['chrom']} "
                f"--from-bp {int(locus['start'])} "
                f"--to-bp {int(locus['end'])} "
                "--export vcf bgz vcf-dosage=DS-force "
                f"--out {vcf_prefix}"
            ]
            export_proc = subprocess.run(
                export_cmd, capture_output=True, text=True)
            if export_proc.returncode != 0:
                rec["status"] = "vcf_export_failed"
                rec["stderr_tail"] = (
                    export_proc.stderr or export_proc.stdout)[-500:]
                continue
            # Append ancestry-specific inputs only after the genotype export
            # succeeds, so VCF/pheno/covar lists cannot get out of sync.
            vcf_paths.append(str(vcf_path))
            pheno_paths.append(str(pheno_path))
            if cov_cache[key] is not None:
                covar_paths.append(str(covar_path))
            anc_used.append(anc)

        if not anc_used:
            rec["status"] = "no_phenotype"
            diagnostics.append(rec)
            continue

        # --- call sushie ---
        cmd = [args.sushie_bin, "finemap"]
        cmd += ["--vcf"] + vcf_paths
        cmd += ["--pheno"] + pheno_paths
        if all(c is not None for c in
               [cov_cache.get((a, mod)) for a in anc_used]):
            cmd += ["--covar"] + covar_paths
        cmd += ["--chrom", str(locus["chrom"]),
                "--start", str(int(locus["start"])),
                "--end", str(int(locus["end"])),
                "--L", str(int(locus["L"])),
                "--purity", str(args.purity),
                "--trait", ph,
                "--platform", args.platform,
                "--output", str(out_prefix)]
        env = dict(os.environ)
        env.update({"OMP_NUM_THREADS": str(args.threads),
                    "OPENBLAS_NUM_THREADS": str(args.threads),
                    "MKL_NUM_THREADS": str(args.threads),
                    "JAX_PLATFORMS": args.platform})
        proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if proc.returncode != 0:
            rec["status"] = "failed"
            rec["stderr_tail"] = (proc.stderr or proc.stdout)[-500:]
            diagnostics.append(rec)
            continue

        diag = parse_sushie_log(Path(str(out_prefix) + ".log"))
        rec.update(diag)
        if not weights_path.exists() or weights_path.stat().st_size == 0:
            rec["status"] = "failed_missing_weights"
            rec["stderr_tail"] = (
                "SuShiE exited 0 but no non-empty .sushie.weights.tsv was produced")
            diagnostics.append(rec)
            continue
        rec["ancestries"] = ",".join(anc_used)
        rec["n_ancestries"] = len(anc_used)
        rec["wall_sec"] = round(time.time() - t0, 2)
        # record ancestry order (weights columns are positional)
        Path(str(out_prefix) + ".ancestries").write_text(
            ",".join(anc_used) + "\n")
        done_path.write_text("ok\n")
        # clean per-locus input TSVs (regenerable; keeps the tree small)
        for p in pheno_paths + covar_paths + vcf_paths:
            if os.path.exists(p):
                os.remove(p)
        diagnostics.append(rec)
        print(f"  {mod}/{trait}: anc={len(anc_used)} "
              f"snps={diag['n_snps']} cs={diag['n_cs']} "
              f"conv={diag['converged']} ({rec['wall_sec']}s)")

    diag_df = pd.DataFrame(diagnostics)
    diag_path = logs_dir / f"{mod}_{shard_name}.diagnostics.tsv"
    diag_df.to_csv(diag_path, sep="\t", index=False)
    n_ok = int((diag_df["status"] == "ok").sum())
    print(f"[run] {shard_name}: {n_ok}/{len(diag_df)} loci fine-mapped; "
          f"diagnostics -> {diag_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pp = sub.add_parser("prepare-loci", help="Build per-modality locus list")
    pp.add_argument("--results-dir", required=True)
    pp.add_argument("--qtl-dir", default=None,
                    help="QTL input directory; default: sibling qtl_inputs")
    pp.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    pp.add_argument("--modality", required=True)
    pp.add_argument("--bed-template", default="{qtl_dir}/{anc}_{mod}.bed.gz")
    pp.add_argument("--cis-window", type=int, default=1_000_000,
                    help="TensorQTL cis-window in bp (default: 1,000,000)")
    pp.add_argument("--qval", type=float, default=0.05)
    pp.add_argument("--l-min", type=int, default=5)
    pp.add_argument("--l-max", type=int, default=10)
    pp.add_argument("--l-buffer", type=int, default=2)
    pp.add_argument("--out", required=True)
    pp.set_defaults(func=prepare_loci)

    pr = sub.add_parser("run", help="Fine-map a shard of loci")
    pr.add_argument("--loci-file", required=True)
    pr.add_argument("--qtl-dir", required=True)
    pr.add_argument("--out-dir", required=True)
    pr.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    pr.add_argument("--bed-template",
                    default="{qtl_dir}/{anc}_{mod}.bed.gz")
    pr.add_argument("--pgen-template", default="{qtl_dir}/{anc}_qtl")
    pr.add_argument("--covar-template",
                    default="{qtl_dir}/{anc}_covariates_{mod}.tsv")
    pr.add_argument("--sushie-bin", default="sushie")
    pr.add_argument("--platform", default="cpu")
    pr.add_argument("--purity", type=float, default=0.5)
    pr.add_argument("--threads", type=int, default=1)
    pr.add_argument("--min-samples", type=int, default=20,
                    help="Skip an ancestry for a locus below this n")
    pr.add_argument("--force", action="store_true",
                    help="Re-run loci with existing .done markers")
    pr.set_defaults(func=run_shard)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
