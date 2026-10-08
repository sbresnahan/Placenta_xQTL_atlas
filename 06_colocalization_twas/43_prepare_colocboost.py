#!/usr/bin/env python3
"""43_prepare_colocboost.py — build multi-trait colocBoost region manifests.

Objective 1.6: colocBoost runs jointly over many outcomes at a locus. The
analysis unit here is a *region*: the per-ancestry merged union of all
fine-mapped loci (module 05, {RESULTS_DIR}/finemap/loci/{MOD}.loci.tsv)
across the 9 per-modality layers. Each region is analyzed per ancestry with:

  - every xQTL phenotype whose fine-mapped locus overlaps the region and is
    significant in that ancestry (sig_in), using the merged genome-wide
    nominal stats (39_run_nominal.py) and in-sample LD from the intersected
    {ANC}_qtl pgen;
  - every ancestry-matched GWAS from the catalog (primary_ancestry equal to
    the stratum or "both") with a harmonized file present, using 1KG
    superpopulation LD.

Outputs (under {out_dir}, default {RESULTS_DIR}/coloc/colocboost):
  {ANC}.regions.tsv   region_id, ancestry, chrom, start, end,
                      n_xqtl_outcomes, n_gwas_traits, n_outcomes
  {ANC}.outcomes.tsv  region_id, outcome_name, type (xqtl|gwas), modality,
                      phenotype_id, trait_id, trait_type, prop_cases,
                      file, ld_side (xqtl|gwas)

Regions with no ancestry-matched GWAS or < 2 total outcomes are dropped.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def merge_intervals(df, gap):
    """Merge overlapping/nearby (chrom, start, end) intervals."""
    out = []
    for chrom, sub in df.groupby("chrom"):
        sub = sub.sort_values("start")
        cur_s = cur_e = None
        for s, e in zip(sub["start"], sub["end"]):
            if cur_s is None:
                cur_s, cur_e = s, e
            elif s <= cur_e + gap:
                cur_e = max(cur_e, e)
            else:
                out.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        if cur_s is not None:
            out.append((chrom, cur_s, cur_e))
    return pd.DataFrame(out, columns=["chrom", "start", "end"])


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-dir", required=True)
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--gwas-dir", required=True,
                   help="harmonized GWAS directory (38_harmonize_gwas.py out)")
    p.add_argument("--catalog", required=True)
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    p.add_argument("--modalities", nargs="*", default=None,
                   help="default: all 9 per-modality layers")
    p.add_argument("--merge-gap", type=int, default=100_000,
                   help="merge loci closer than this into one region")
    p.add_argument("--max-outcomes", type=int, default=500,
                   help="warn (and skip) regions with more outcomes")
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()

    results_dir = Path(args.results_dir)
    loci_dir = results_dir / "finemap" / "loci"
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    modalities = args.modalities or [
        "expression", "isoforms", "isoform_expression", "splicing",
        "intron_retention", "alt_TSS", "alt_polyA", "RNA_editing", "stability"]

    # --- catalog: ancestry-matched GWAS with harmonized files present -------
    cat = pd.read_csv(args.catalog, sep="\t")
    cat["gwas_file"] = cat["trait_id"].map(
        lambda t: str(Path(args.gwas_dir) / f"{t}.sumstats.tsv.gz"))
    cat["gwas_ok"] = cat["gwas_file"].map(lambda f: Path(f).exists())

    # --- load all loci -------------------------------------------------------
    frames = []
    for mod in modalities:
        path = loci_dir / f"{mod}.loci.tsv"
        if not path.exists():
            print(f"  WARNING: no locus list for {mod}; skipped")
            continue
        df = pd.read_csv(path, sep="\t")
        df["chrom"] = df["chrom"].astype(str).str.replace("^chr", "", regex=True)
        frames.append(df)
    if not frames:
        sys.exit("ERROR: no locus lists found under " + str(loci_dir))
    loci = pd.concat(frames, ignore_index=True)

    for anc in args.ancestries:
        # loci significant in this ancestry with nominal stats available
        loc_a = loci[loci["sig_in"].astype(str).str.contains(anc)]
        keep = []
        for _, r in loc_a.iterrows():
            xqtl = results_dir / "nominal" / anc / f"{anc}_{r['modality']}.nominal.tsv.gz"
            keep.append(xqtl.exists())
        loc_a = loc_a[keep]
        if len(loc_a) == 0:
            print(f"  {anc}: no usable loci; skipped")
            continue

        regions = merge_intervals(loc_a[["chrom", "start", "end"]],
                                  args.merge_gap)
        regions["ancestry"] = anc
        regions["region_id"] = [f"{anc}_{c}_{s}_{e}" for c, s, e in
                                zip(regions["chrom"], regions["start"],
                                    regions["end"])]

        gwas_a = cat[cat["gwas_ok"] &
                     cat["primary_ancestry"].isin([anc, "both"])]
        if len(gwas_a) == 0:
            print(f"  {anc}: no ancestry-matched GWAS files; skipped")
            continue

        outcome_rows, region_rows = [], []
        for _, reg in regions.iterrows():
            ov = loc_a[(loc_a["chrom"] == reg["chrom"]) &
                       (loc_a["start"] <= reg["end"]) &
                       (loc_a["end"] >= reg["start"])]
            rows = []
            for _, r in ov.iterrows():
                rows.append({
                    "region_id": reg["region_id"], "type": "xqtl",
                    "outcome_name": f"{r['modality']}:{r['phenotype_id']}",
                    "modality": r["modality"], "phenotype_id": r["phenotype_id"],
                    "trait_id": "", "trait_type": "quantitative",
                    "prop_cases": "",
                    "file": str(results_dir / "nominal" / anc /
                                f"{anc}_{r['modality']}.nominal.tsv.gz"),
                    "ld_side": "xqtl"})
            for _, g in gwas_a.iterrows():
                rows.append({
                    "region_id": reg["region_id"], "type": "gwas",
                    "outcome_name": f"gwas:{g['trait_id']}",
                    "modality": "", "phenotype_id": "",
                    "trait_id": g["trait_id"],
                    "trait_type": g.get("trait_type", "quantitative"),
                    "prop_cases": g.get("prop_cases", ""),
                    "file": g["gwas_file"], "ld_side": "gwas"})
            n_xqtl = sum(r["type"] == "xqtl" for r in rows)
            n_gwas = sum(r["type"] == "gwas" for r in rows)
            if n_xqtl + n_gwas < 2:
                continue
            if n_xqtl + n_gwas > args.max_outcomes:
                print(f"  WARNING: {reg['region_id']} has "
                      f"{n_xqtl + n_gwas} outcomes (> {args.max_outcomes}); "
                      f"skipped")
                continue
            outcome_rows.extend(rows)
            region_rows.append({"region_id": reg["region_id"], "ancestry": anc,
                                "chrom": reg["chrom"], "start": reg["start"],
                                "end": reg["end"], "n_xqtl_outcomes": n_xqtl,
                                "n_gwas_traits": n_gwas,
                                "n_outcomes": n_xqtl + n_gwas})

        if not region_rows:
            print(f"  {anc}: no regions with >= 2 outcomes")
            continue
        pd.DataFrame(region_rows).to_csv(out_dir / f"{anc}.regions.tsv",
                                         sep="\t", index=False)
        pd.DataFrame(outcome_rows).to_csv(out_dir / f"{anc}.outcomes.tsv",
                                          sep="\t", index=False)
        print(f"  {anc}: {len(region_rows)} regions, "
              f"{len(outcome_rows)} outcome assignments -> "
              f"{anc}.regions.tsv / {anc}.outcomes.tsv")


if __name__ == "__main__":
    main()
