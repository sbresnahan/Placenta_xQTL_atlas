#!/usr/bin/env python3
"""33_aggregate_finemap.py — aggregate per-locus SuSHiE outputs.

Walks {finemap_dir}/{MOD}/ for per-locus SuSHiE outputs (from
31_sushie_finemap.py run) and the shard diagnostics in
{finemap_dir}/logs/, and writes three analysis-ready tables:

- finemap_pips.tsv.gz — one row per variant per locus: joint PIP
  (sushie_pip_all), credible-set membership, and per-ancestry posterior
  effect weights (columns named by the actual ancestry labels recorded in
  each locus's .ancestries file — SuSHiE's ancestry1/2/... output columns
  are positional).
- finemap_credible_sets.tsv.gz — one row per credible set: size, lead
  variant, span, mean/min PIP, and the cross-ancestry effect correlation
  (rho) from the .corr.tsv (multi-ancestry loci only).
- finemap_locus_summary.tsv — one row per locus: shard diagnostics
  (convergence, n_snps, wall time, status) plus aggregated signal stats
  (max PIP, lead variant, number of credible sets).

Handles SuSHiE's degenerate empty-CS output (a cs.tsv containing only a
trait header line) and single-ancestry loci (no .corr.tsv).
"""

import argparse
from pathlib import Path

import pandas as pd


def read_ancestries(locus_prefix: Path, fallback):
    f = Path(str(locus_prefix) + ".ancestries")
    if f.exists():
        return f.read_text().strip().split(",")
    return fallback


def collect(finemap_dir: Path):
    """Yield per-locus (weights, cs, corr) DataFrames with provenance."""
    weight_frames, cs_frames, corr_frames = [], [], []
    for mod_dir in sorted(p for p in finemap_dir.iterdir()
                          if p.is_dir() and p.name not in ("logs", "loci",
                                                           "shards",
                                                           "aggregated")):
        mod = mod_dir.name
        for w_path in sorted(mod_dir.glob("*.sushie.weights.tsv")):
            prefix = w_path.with_suffix("")  # strips .tsv
            prefix = Path(str(prefix).replace(".sushie.weights", ""))
            trait = prefix.name
            anc = read_ancestries(prefix, None)
            w = pd.read_csv(w_path, sep="\t")
            w["modality"] = mod
            w["phenotype_id"] = w["trait"]
            # rename positional ancestry columns with real labels
            if anc:
                ren = {}
                for i, a in enumerate(anc, start=1):
                    if f"ancestry{i}_sushie_weight" in w.columns:
                        ren[f"ancestry{i}_sushie_weight"] = f"{a}_effect_weight"
                w = w.rename(columns=ren)
            w["cs_index"] = pd.to_numeric(
                w["sushie_cs_index"].replace("No CS", None),
                errors="coerce").astype("Int64")
            w = w.drop(columns=["sushie_cs_index"])
            weight_frames.append(w)

            cs_path = Path(str(prefix) + ".sushie.cs.tsv")
            if cs_path.exists():
                raw = cs_path.read_text().strip().split("\n")
                # degenerate empty-CS file: header 'trait' + one name line
                if len(raw) <= 2 and raw[0].strip() == "trait":
                    pass
                else:
                    cs = pd.read_csv(cs_path, sep="\t")
                    cs["modality"] = mod
                    cs_frames.append(cs)

            corr_path = Path(str(prefix) + ".sushie.corr.tsv")
            if corr_path.exists():
                corr = pd.read_csv(corr_path, sep="\t")
                if len(corr):
                    corr["modality"] = mod
                    corr["phenotype_id"] = corr["trait"]
                    if anc:
                        ren = {}
                        for i, a in enumerate(anc, start=1):
                            ren[f"ancestry{i}_est_var"] = f"{a}_est_var"
                            for j, b in enumerate(anc, start=1):
                                if i < j:
                                    ren[f"ancestry{i}_ancestry{j}_est_covar"] = \
                                        f"{a}_{b}_est_covar"
                                    ren[f"ancestry{i}_ancestry{j}_est_corr"] = \
                                        f"{a}_{b}_est_corr"
                        corr = corr.rename(columns=ren)
                    corr_frames.append(corr)
    return weight_frames, cs_frames, corr_frames


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--finemap-dir", required=True,
                   help="SuSHiE output root ($RESULTS_DIR/finemap)")
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"],
                   help="Fallback ancestry labels if a .ancestries file is "
                        "missing (positional: ancestry1 -> first, ...)")
    p.add_argument("--out-dir", default=None,
                   help="Default: {finemap_dir}/aggregated")
    args = p.parse_args()

    finemap_dir = Path(args.finemap_dir)
    out_dir = Path(args.out_dir) if args.out_dir else finemap_dir / "aggregated"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Collecting per-locus outputs under {finemap_dir} ...")
    weight_frames, cs_frames, corr_frames = collect(finemap_dir)
    print(f"  {len(weight_frames)} loci with weights, "
          f"{len(cs_frames)} with credible sets, "
          f"{len(corr_frames)} with cross-ancestry correlations")

    # --- diagnostics / locus summary ---
    diag_frames = []
    for d in sorted((finemap_dir / "logs").glob("*.diagnostics.tsv")):
        diag_frames.append(pd.read_csv(d, sep="\t"))
    if not diag_frames:
        raise SystemExit("No shard diagnostics found under "
                         f"{finemap_dir}/logs — nothing to aggregate")
    diag = pd.concat(diag_frames, ignore_index=True)
    ok = diag[diag["status"] == "ok"]
    print(f"  diagnostics: {len(diag)} loci "
          f"({len(ok)} ok, {len(diag) - len(ok)} skipped/failed)")

    # --- PIPs ---
    pips = pd.concat(weight_frames, ignore_index=True)
    keep_cols = (["modality", "phenotype_id", "chrom", "pos", "snp", "a0",
                  "a1", "sushie_pip_all", "sushie_pip_cs", "cs_index"]
                 + [c for c in pips.columns if c.endswith("_effect_weight")])
    pips = pips[keep_cols].rename(columns={"sushie_pip_all": "pip_all",
                                           "sushie_pip_cs": "pip_cs"})
    pips = pips.merge(
        diag[["modality", "phenotype_id", "pheno_start", "pheno_end",
              "sig_in"]].drop_duplicates(),
        on=["modality", "phenotype_id"], how="left")
    pips_path = out_dir / "finemap_pips.tsv.gz"
    pips.to_csv(pips_path, sep="\t", index=False, compression="gzip")
    print(f"  wrote {len(pips)} variant rows -> {pips_path}")

    # --- credible sets ---
    if cs_frames:
        cs_all = pd.concat(cs_frames, ignore_index=True)
        cs_summary = (cs_all.groupby(["modality", "trait", "CSIndex"])
                      .agg(n_variants=("snp", "size"),
                           lead_snp=("snp", lambda s: s[cs_all.loc[s.index,
                                                              "pip_all"].idxmax()]),
                           lead_pos=("pos", lambda s: s[cs_all.loc[s.index,
                                                            "pip_all"].idxmax()]),
                           lead_pip=("pip_all", "max"),
                           mean_pip=("pip_all", "mean"),
                           min_pip=("pip_all", "min"),
                           chrom=("chrom", "first"),
                           cs_start=("pos", "min"),
                           cs_end=("pos", "max"))
                      .reset_index()
                      .rename(columns={"trait": "phenotype_id",
                                       "CSIndex": "cs_index"}))
        if corr_frames:
            corr_all = pd.concat(corr_frames, ignore_index=True)
            corr_cols = [c for c in corr_all.columns
                         if c.endswith("_est_corr") or c.endswith("_est_var")
                         or c.endswith("_est_covar")]
            cs_summary = cs_summary.merge(
                corr_all[["modality", "phenotype_id", "CSIndex"] + corr_cols]
                .rename(columns={"CSIndex": "cs_index"}),
                on=["modality", "phenotype_id", "cs_index"], how="left")
        cs_path = out_dir / "finemap_credible_sets.tsv.gz"
        cs_summary.to_csv(cs_path, sep="\t", index=False, compression="gzip")
        print(f"  wrote {len(cs_summary)} credible sets -> {cs_path}")
    else:
        cs_summary = pd.DataFrame()
        print("  no credible sets passed purity pruning in any locus")

    # --- locus summary ---
    if len(cs_summary):
        per_locus = (cs_summary.groupby(["modality", "phenotype_id"])
                     .agg(n_cs_out=("cs_index", "size"),
                          max_pip=("lead_pip", "max"))
                     .reset_index())
        lead = (cs_summary.loc[cs_summary.groupby(["modality", "phenotype_id"])
                               ["lead_pip"].idxmax()]
                [["modality", "phenotype_id", "lead_snp", "lead_pos"]]
                .rename(columns={"lead_snp": "top_snp",
                                 "lead_pos": "top_pos"}))
        locus_summary = diag.merge(per_locus, on=["modality", "phenotype_id"],
                                   how="left").merge(
            lead, on=["modality", "phenotype_id"], how="left")
    else:
        locus_summary = diag.copy()
        locus_summary["n_cs_out"] = pd.NA
        locus_summary["max_pip"] = pd.NA
    ls_path = out_dir / "finemap_locus_summary.tsv"
    locus_summary.to_csv(ls_path, sep="\t", index=False)
    print(f"  wrote {len(locus_summary)} loci -> {ls_path}")


if __name__ == "__main__":
    main()
