#!/usr/bin/env python3
"""40_prepare_coloc_loci.py — build the colocalization task list.

Objective 1.6 colocalizes every FDR-significant xQTL locus (the module-05
fine-mapping locus lists: union of q <= 0.05 phenotypes per modality) with
each ancestry-matched GWAS from gwas_catalog.tsv. This script joins those
two inputs into a per-modality task table consumed by 41_susie_coloc.R
(one row per phenotype x ancestry x GWAS trait) and writes the 1000-Genomes
superpopulation keep-files used for GWAS-side LD.

Ancestry matching (catalog `primary_ancestry` column):
  EUR  -> colocalized in the EUR stratum
  EAS  -> colocalized in the EAS stratum
  both -> colocalized in both strata (e.g. multi-ancestry ProDiGY)
A stratum is only tasked if its locus list exists and its harmonized GWAS
file is present (missing inputs are reported and skipped, so the task list
can be built before all manual GWAS files arrive).

Outputs
-------
  {coloc_dir}/loci/{MOD}.tasks.tsv   one row per task:
      modality, phenotype_id, chrom, start, end, L, ancestry, trait_id,
      gwas_file, xqtl_file (merged nominal TSV), ld_xqtl (pgen prefix),
      ld_gwas (1KG pgen prefix + keep file)
  {coloc_dir}/loci/{ANC}.1kg.keep    sample IDs for the 1KG superpopulation
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def write_kg_keep_files(sample_map_path, ancestries, out_dir, kg_pgen=None):
    """Write {ANC}.1kg.keep files from a 1KG sample->superpopulation map.

    The map must have columns sample_id and superpop (EAS/EUR/AFR/AMR/SAS),
    e.g. the 1kGP pedigree/integrated_call_samples table.

    plink2 --keep matches on FID+IID, and a single-column file is only
    interpreted as IID when every FID in the dataset is missing/'0'. To be
    robust to the 1KG pgen's actual ID scheme, the .psam is read (when
    kg_pgen is given and readable) and matching rows are emitted with their
    true FID/IID; otherwise a two-column file with the sample ID duplicated
    is written (correct for the common FID==IID convention).
    """
    sm = pd.read_csv(sample_map_path, sep=None, engine="python")
    cols = {c.lower(): c for c in sm.columns}
    sid = next((cols[c] for c in ("sample_id", "sample", "iid", "sampleid")
                if c in cols), None)
    ssp = next((cols[c] for c in ("superpop", "superpopulation",
                                  "superpopulation_code", "population")
                if c in cols), None)
    if sid is None or ssp is None:
        sys.exit(f"ERROR: cannot identify sample/superpop columns in "
                 f"{sample_map_path} (columns: {list(sm.columns)})")

    # map sample ID -> (FID, IID) from the pgen's psam when available
    id_map = {}
    psam = Path(f"{kg_pgen}.psam") if kg_pgen else None
    if psam is not None and psam.exists():
        with open(psam) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    parts = line.split()
                if len(parts) >= 2:
                    fid, iid = parts[0], parts[1]
                    id_map.setdefault(iid, (fid, iid))
                    id_map.setdefault(fid, (fid, iid))
        print(f"  read {len(id_map)} sample IDs from {psam.name}")
    else:
        print(f"  WARNING: {psam} not readable; writing duplicated-ID "
              f"keep files (assumes FID==IID in the 1KG pgen)")

    written = {}
    for anc in ancestries:
        ids = sm.loc[sm[ssp] == anc, sid].astype(str)
        if len(ids) == 0:
            print(f"  WARNING: no 1KG samples for superpopulation {anc}")
            continue
        pairs, missing = [], 0
        for s in ids:
            if s in id_map:
                pairs.append(id_map[s])
            elif id_map:
                missing += 1
            else:
                pairs.append((s, s))
        if missing:
            print(f"  WARNING: {missing}/{len(ids)} {anc} samples absent "
                  f"from {psam.name}")
        path = Path(out_dir) / f"{anc}.1kg.keep"
        pd.DataFrame(pairs).to_csv(path, sep="\t", index=False, header=False)
        written[anc] = str(path)
        print(f"  {anc}: {len(pairs)} 1KG samples -> {path.name}")
    return written


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-dir", required=True,
                   help="qtl_results dir (holds finemap/loci/ and nominal/)")
    p.add_argument("--catalog", required=True, help="gwas_catalog.tsv")
    p.add_argument("--gwas-dir", required=True,
                   help="harmonized GWAS directory (38_harmonize_gwas.py out)")
    p.add_argument("--qtl-dir", required=True,
                   help="qtl_inputs dir (for {ANC}_qtl pgen paths)")
    p.add_argument("--kg-pgen", required=True,
                   help="1KG pgen prefix (module 01: 1kGP_hg38)")
    p.add_argument("--kg-sample-map", default=None,
                   help="1KG sample->superpopulation table; required once to "
                        "write keep files (skipped if they already exist)")
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    p.add_argument("--modalities", nargs="*", default=None,
                   help="default: all modalities with a locus list")
    p.add_argument("--out-dir", required=True, help="coloc dir")
    args = p.parse_args()

    results_dir = Path(args.results_dir)
    loci_dir = results_dir / "finemap" / "loci"
    out_dir = Path(args.out_dir)
    (out_dir / "loci").mkdir(parents=True, exist_ok=True)

    # --- 1KG keep files ---
    keep = {}
    missing_keep = [a for a in args.ancestries
                    if not (out_dir / "loci" / f"{a}.1kg.keep").exists()]
    if missing_keep:
        if args.kg_sample_map is None:
            sys.exit("ERROR: 1KG keep files missing for "
                     f"{missing_keep} and --kg-sample-map not given")
        keep = write_kg_keep_files(args.kg_sample_map, args.ancestries,
                                   out_dir / "loci", kg_pgen=args.kg_pgen)
    for a in args.ancestries:
        path = out_dir / "loci" / f"{a}.1kg.keep"
        if path.exists():
            keep[a] = str(path)

    # --- catalog ---
    cat = pd.read_csv(args.catalog, sep="\t")
    cat["gwas_file"] = cat["trait_id"].map(
        lambda t: str(Path(args.gwas_dir) / f"{t}.sumstats.tsv.gz"))
    cat["gwas_ok"] = cat["gwas_file"].map(lambda f: Path(f).exists())
    for _, r in cat[~cat["gwas_ok"]].iterrows():
        print(f"  WARNING: harmonized GWAS missing for {r['trait_id']} "
              f"({r['gwas_file']}); its tasks are skipped")

    # --- per-modality task tables ---
    modalities = args.modalities
    if not modalities:
        modalities = sorted(p.stem.replace(".loci", "")
                            for p in loci_dir.glob("*.loci.tsv"))
    summary = []
    for mod in modalities:
        loci_path = loci_dir / f"{mod}.loci.tsv"
        if not loci_path.exists():
            print(f"  WARNING: no locus list for {mod}; skipped")
            continue
        loci = pd.read_csv(loci_path, sep="\t")
        rows = []
        for _, locus in loci.iterrows():
            for anc in args.ancestries:
                xqtl_file = (results_dir / "nominal" / anc /
                             f"{anc}_{mod}.nominal.tsv.gz")
                if not xqtl_file.exists():
                    continue
                matched = cat[cat["gwas_ok"] &
                              cat["primary_ancestry"].isin([anc, "both"])]
                for _, g in matched.iterrows():
                    rows.append({
                        "modality": mod,
                        "phenotype_id": locus["phenotype_id"],
                        "chrom": str(locus["chrom"]).replace("chr", ""),
                        "start": int(locus["start"]), "end": int(locus["end"]),
                        "L": int(locus["L"]),
                        "ancestry": anc,
                        "trait_id": g["trait_id"],
                        "trait_type": g.get("trait_type", "quantitative"),
                        "prop_cases": g.get("prop_cases", ""),
                        "gwas_file": g["gwas_file"],
                        "xqtl_file": str(xqtl_file),
                        "ld_xqtl_pgen": str(Path(args.qtl_dir) / f"{anc}_qtl"),
                        "ld_gwas_pgen": args.kg_pgen,
                        "ld_gwas_keep": keep.get(anc, ""),
                    })
        tasks = pd.DataFrame(rows)
        out_path = out_dir / "loci" / f"{mod}.tasks.tsv"
        tasks.to_csv(out_path, sep="\t", index=False)
        summary.append({"modality": mod, "loci": len(loci),
                        "tasks": len(tasks)})
        print(f"  {mod}: {len(loci)} loci -> {len(tasks)} tasks "
              f"-> {out_path.name}")
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
