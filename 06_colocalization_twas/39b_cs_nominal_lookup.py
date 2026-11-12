#!/usr/bin/env python3
"""39b_cs_nominal_lookup.py — nominal slope/SE of every fine-mapped
credible-set variant against every PANTRY feature of the same gene.

For each variant in ANY SuSHiE 95% credible set (any modality, any
phenotype), look up the stage-39 nominal cis-xQTL statistics (tensorQTL
slope, slope_se, pval_nominal, af, ...) for every phenotype of the
fine-mapped phenotype's gene, across all 9 modalities and per ancestry, and
package the result as a self-describing archive.

Inputs
------
  --finemap-dir   $RESULTS_DIR/finemap   (SuSHiE per-locus outputs, 31/32)
                  CS membership is read from the raw {MOD}/*.sushie.cs.tsv
                  files (one row per CS x variant, so a variant in two CSs of
                  one locus is kept twice). `--cs-source pips` instead reads
                  aggregated/finemap_pips.tsv.gz (33), which stores a single
                  cs_index per variant.
  --results-dir   $RESULTS_DIR           (holds nominal/{ANC}/{ANC}_{MOD}.nominal.tsv.gz
                                          + .tbi from 39_run_nominal.py)
  --qtl-dir       $QTL_DIR (optional)    ({ANC}_{MOD}.bed.gz, used to list
                                          every feature of a gene so untested
                                          variant x feature pairs get explicit
                                          NA rows with a reason)

Gene assignment follows the PANTRY / 30_combine_modalities.py convention:
gene_id = phenotype_id up to the first '__' (bare ENSG IDs map to
themselves), with the Ensembl version suffix stripped for matching.

Subcommands
-----------
  index    CS membership table + per-variant query keys (run once)
  extract  one ancestry x modality lookup (LSF array task)
  package  README + MANIFEST + tar.gz of the output directory
  all      index, every ancestry x modality extract, package (single job)

Outputs (in --out-dir)
----------------------
  cs_membership.tsv.gz                 one row per credible set x variant
  {ANC}/{ANC}_{MOD}.cs_nominal.tsv.gz  one row per CS variant x feature of
                                       the CS gene (bgzip + .tbi when
                                       bgzip/tabix are on PATH)
  {ANC}/{ANC}_{MOD}.cs_nominal.done.json  per-task counts + parameters
  README.txt, MANIFEST.tsv             written by `package`
  ../{out_dir_name}_{YYYYMMDD}.tar.gz  written by `package`

Allele convention: tensorQTL `slope` is the effect per copy of the ALT
allele of variant_id (chr:pos:REF:ALT; pgen ALT dosage) on the
normalized phenotype, and `af` is the ALT allele frequency in that
ancestry. Effects from different ancestries are on the same allele.
"""

import argparse
import datetime as dt
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import pandas as pd

MODALITIES = ["expression", "isoforms", "isoform_expression", "splicing",
              "intron_retention", "alt_TSS", "alt_polyA", "RNA_editing",
              "stability"]
ANCESTRIES = ["EAS", "EUR"]
GENE_VER_RE = r"^(ENS[A-Z]*G\d+)\.\d+$"
STAT_COLS = ["slope", "slope_se", "pval_nominal", "af", "ma_count",
             "ma_samples", "start_distance", "end_distance"]
CHUNK = 2_000_000


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def gene_of(pheno_ids: pd.Series) -> pd.Series:
    """PANTRY gene for a phenotype id: prefix before '__', version stripped."""
    g = pheno_ids.astype(str).str.split("__", n=1).str[0]
    return g.str.replace(GENE_VER_RE, r"\1", regex=True)


def norm_chrom(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace(r"^chr", "", regex=True)


def chrom_sort_key(c):
    c = str(c)
    return (0, int(c)) if c.isdigit() else (1, c)


def work_dir(out_dir: Path) -> Path:
    return out_dir / "_work"


def task_paths(out_dir: Path, anc: str, mod: str):
    d = out_dir / anc
    return (d / f"{anc}_{mod}.cs_nominal.tsv",
            d / f"{anc}_{mod}.cs_nominal.done.json")


def read_gz_header(path: Path):
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")


# ---------------------------------------------------------------------------
# index
# ---------------------------------------------------------------------------

def sanitize_trait(phenotype_id: str) -> str:
    """Locus file prefix — must match 31_sushie_finemap.sanitize_trait."""
    return re.sub(r"[^A-Za-z0-9._+-]", "_", phenotype_id)


def _cs_paths(finemap_dir: Path, mod: str):
    """cs.tsv paths for one modality, built from the 32 locus list
    (finemap/loci/{MOD}.loci.tsv) so the large per-locus output directory is
    never listed; falls back to a directory scan when the list is missing."""
    mdir = finemap_dir / mod
    loci = finemap_dir / "loci" / f"{mod}.loci.tsv"
    if loci.exists():
        ph = pd.read_csv(loci, sep="\t", usecols=["phenotype_id"],
                         dtype=str)["phenotype_id"]
        return [mdir / f"{sanitize_trait(x)}.sushie.cs.tsv"
                for x in ph.drop_duplicates()], "locus list"
    print(f"  WARNING: {loci} not found; listing {mdir} (slow on GPFS)",
          flush=True)
    if not mdir.is_dir():
        return [], "missing"
    with os.scandir(mdir) as it:
        return [Path(e.path) for e in it
                if e.name.endswith(".sushie.cs.tsv")], "directory scan"


def _read_one_cs(path: Path):
    """Parse one SuSHiE cs.tsv. Returns (status, rows, header)."""
    try:
        with open(path) as fh:
            text = fh.read()
    except FileNotFoundError:
        return "missing", None, None
    lines = [l for l in text.split("\n") if l.strip()]
    # degenerate empty-CS file: header 'trait' + one name line (33)
    if len(lines) <= 2 and lines and lines[0].strip() == "trait":
        return "empty", None, None
    if len(lines) < 2:
        return "empty", None, None
    header = lines[0].split("\t")
    return "ok", [l.split("\t") for l in lines[1:]], header


def _read_cs_files(finemap_dir: Path, modalities, threads=16):
    from concurrent.futures import ThreadPoolExecutor
    need = ["trait", "CSIndex", "snp", "chrom", "pos", "pip_all"]
    frames, totals = [], {"ok": 0, "empty": 0, "missing": 0}
    for mod in modalities:
        paths, how = _cs_paths(finemap_dir, mod)
        n = len(paths)
        print(f"  {mod}: {n} loci ({how})", flush=True)
        counts = {"ok": 0, "empty": 0, "missing": 0}
        recs, header0 = [], None
        step = max(1000, n // 10)
        with ThreadPoolExecutor(max_workers=threads) as ex:
            for i, (st, rows, header) in enumerate(
                    ex.map(_read_one_cs, paths, chunksize=64), start=1):
                counts[st] += 1
                if st == "ok":
                    if header0 is None:
                        header0 = header
                        miss = [c for c in need if c not in header0]
                        if miss:
                            sys.exit(f"ERROR: cs.tsv lacks columns {miss} "
                                     f"({paths[i - 1]})")
                    if header != header0:  # column order differs: realign
                        idx = [header.index(c) for c in header0
                               if c in header]
                        rows = [[r[j] for j in idx] for r in rows]
                    recs.extend(rows)
                if i % step == 0:
                    print(f"    {mod}: {i}/{n} read", flush=True)
        print(f"    {mod}: {counts['ok']} with CS, {counts['empty']} without, "
              f"{counts['missing']} missing (failed/not run)", flush=True)
        for k in totals:
            totals[k] += counts[k]
        if recs:
            df = pd.DataFrame(recs, columns=header0)
            df["modality"] = mod
            frames.append(df)
    print(f"  total: {totals}", flush=True)
    if not frames:
        sys.exit("ERROR: no credible-set rows found")
    cs = pd.concat(frames, ignore_index=True)
    out = pd.DataFrame({
        "modality": cs["modality"],
        "phenotype_id": cs["trait"].astype(str),
        "cs_index": pd.to_numeric(cs["CSIndex"], errors="coerce").astype("Int64"),
        "variant_id": cs["snp"].astype(str),
        "chrom": norm_chrom(cs["chrom"]),
        "pos": pd.to_numeric(cs["pos"], errors="coerce").astype("Int64"),
        "pip_all": pd.to_numeric(cs["pip_all"], errors="coerce"),
        "pip_cs": (pd.to_numeric(cs["pip_cs"], errors="coerce")
                   if "pip_cs" in cs.columns else pd.NA),
    })
    return out


def _read_pips(pips_path: Path, modalities):
    pips = pd.read_csv(pips_path, sep="\t",
                       dtype={"chrom": str, "snp": str, "phenotype_id": str})
    pips = pips[pips["cs_index"].notna() & pips["modality"].isin(modalities)]
    return pd.DataFrame({
        "modality": pips["modality"],
        "phenotype_id": pips["phenotype_id"],
        "cs_index": pd.to_numeric(pips["cs_index"]).astype("Int64"),
        "variant_id": pips["snp"],
        "chrom": norm_chrom(pips["chrom"]),
        "pos": pd.to_numeric(pips["pos"]).astype("Int64"),
        "pip_all": pips["pip_all"],
        "pip_cs": pips.get("pip_cs", pd.NA),
    })


def cmd_index(args):
    out_dir = Path(args.out_dir)
    wd = work_dir(out_dir)
    wd.mkdir(parents=True, exist_ok=True)
    finemap_dir = Path(args.finemap_dir)
    print(f"[index] CS source: {args.cs_source}")
    if args.cs_source == "cs":
        cs = _read_cs_files(finemap_dir, args.modalities, args.threads)
    else:
        pips_path = Path(args.pips or finemap_dir / "aggregated" /
                         "finemap_pips.tsv.gz")
        cs = _read_pips(pips_path, args.modalities)

    bad = cs["pos"].isna() | cs["cs_index"].isna()
    if bad.any():
        print(f"  WARNING: dropping {int(bad.sum())} rows with missing "
              f"pos/CSIndex")
        cs = cs[~bad]
    cs = cs.drop_duplicates(["modality", "phenotype_id", "cs_index",
                             "variant_id"])

    # variant-ID sanity: chr:pos:ref:alt with pos matching the pos column
    parts = cs["variant_id"].str.split(":", expand=True)
    if parts.shape[1] >= 2:
        id_pos = pd.to_numeric(parts[1], errors="coerce")
        mism = (id_pos != cs["pos"].astype(float)).sum()
        if mism:
            print(f"  WARNING: {mism} variant IDs whose position field "
                  f"disagrees with the pos column")

    grp = cs.groupby(["modality", "phenotype_id", "cs_index"], sort=False)
    cs["cs_size"] = grp["variant_id"].transform("size")
    lead_idx = grp["pip_all"].idxmax()
    cs["is_cs_lead"] = False
    cs.loc[lead_idx.dropna().astype(int), "is_cs_lead"] = True
    cs["gene_id"] = gene_of(cs["phenotype_id"])

    n_all = len(cs)
    if args.min_pip > 0:
        cs = cs[cs["pip_all"] >= args.min_pip]
    if args.lead_only:
        cs = cs[cs["is_cs_lead"]]
    print(f"  {n_all} CS x variant rows; {len(cs)} after filters "
          f"(min_pip={args.min_pip}, lead_only={args.lead_only})")

    cs = cs[["chrom", "pos", "variant_id", "gene_id", "modality",
             "phenotype_id", "cs_index", "cs_size", "pip_all", "pip_cs",
             "is_cs_lead"]]
    cs["_c"] = cs["chrom"].map(lambda c: chrom_sort_key(c)[1]
                               if str(c).isdigit() else 10**6)
    cs = cs.sort_values(["_c", "chrom", "pos", "variant_id", "modality",
                         "phenotype_id", "cs_index"]).drop(columns="_c")
    cs.to_csv(out_dir / "cs_membership.tsv.gz", sep="\t", index=False,
              compression="gzip")

    keys = cs[["chrom", "pos", "variant_id", "gene_id"]].drop_duplicates()
    keys.to_csv(wd / "query_keys.tsv.gz", sep="\t", index=False,
                compression="gzip")
    own = cs[["variant_id", "modality", "phenotype_id"]].drop_duplicates()
    own.to_csv(wd / "own_cs_pairs.tsv.gz", sep="\t", index=False,
               compression="gzip")
    params = {"cs_source": args.cs_source, "min_pip": args.min_pip,
              "lead_only": args.lead_only, "modalities": args.modalities,
              "finemap_dir": str(finemap_dir),
              "n_cs_rows": int(len(cs)),
              "n_credible_sets": int(cs.groupby(["modality", "phenotype_id",
                                                 "cs_index"]).ngroups),
              "n_unique_variants": int(cs["variant_id"].nunique()),
              "n_variant_gene_keys": int(len(keys)),
              "n_genes": int(keys["gene_id"].nunique()),
              "created": dt.datetime.now().isoformat(timespec="seconds")}
    (wd / "index_params.json").write_text(json.dumps(params, indent=2))
    print(f"  {params['n_credible_sets']} credible sets, "
          f"{params['n_unique_variants']} unique variants, "
          f"{params['n_genes']} genes, {len(keys)} variant x gene keys")


# ---------------------------------------------------------------------------
# extract
# ---------------------------------------------------------------------------

def _iter_tabix(nominal: Path, header, chrom, positions, tmp_dir):
    reg = Path(tmp_dir) / f"regions.chr{chrom}.txt"
    with open(reg, "w") as fh:
        for p in positions:
            fh.write(f"{chrom}\t{p}\t{p}\n")
    proc = subprocess.Popen(["tabix", "-R", str(reg), str(nominal)],
                            stdout=subprocess.PIPE, text=True)
    try:
        yield from pd.read_csv(proc.stdout, sep="\t", header=None,
                               names=header, chunksize=CHUNK,
                               dtype={"chrom": str, "phenotype_id": str,
                                      "variant_id": str})
    except pd.errors.EmptyDataError:
        return
    finally:
        proc.stdout.close()
        if proc.wait() != 0:
            sys.exit(f"ERROR: tabix failed on {nominal} chr{chrom}")


def _iter_scan(nominal: Path, chroms):
    """Full-file stream; yields (chrom, chunk) restricted to `chroms`."""
    for chunk in pd.read_csv(nominal, sep="\t", chunksize=CHUNK,
                             dtype={"chrom": str, "phenotype_id": str,
                                    "variant_id": str}):
        chunk["chrom"] = norm_chrom(chunk["chrom"])
        chunk = chunk[chunk["chrom"].isin(chroms)]
        for c, sub in chunk.groupby("chrom", sort=False):
            yield c, sub


def _load_features(qtl_dir: Path, anc: str, mod: str):
    bed = qtl_dir / f"{anc}_{mod}.bed.gz"
    if not bed.exists():
        print(f"  WARNING: {bed} not found; untested features will not be "
              f"listed")
        return None
    ids = []
    with gzip.open(bed, "rt") as fh:
        fh.readline()
        for line in fh:
            ids.append(line.split("\t", 4)[3])
    f = pd.DataFrame({"phenotype_id": ids}).drop_duplicates()
    f["gene_id"] = gene_of(f["phenotype_id"])
    return f


def _process_chrom(rows, keys_c, feats, own, scope, anc, mod, out_cols):
    """rows: nominal rows at the query positions of one chromosome."""
    rows = rows[rows["variant_id"].isin(set(keys_c["variant_id"]))]
    present_variants = set(rows["variant_id"])
    rows = rows.assign(gene_id=gene_of(rows["phenotype_id"]))
    kg = keys_c[["variant_id", "gene_id"]].drop_duplicates()
    if scope == "gene":
        res = rows.merge(kg, on=["variant_id", "gene_id"], how="inner")
    else:
        kg2 = kg.assign(cs_gene=True)
        res = rows.merge(kg2, on=["variant_id", "gene_id"], how="left")
        res["cs_gene"] = res["cs_gene"].fillna(False).astype(bool)
    res["status"] = "tested"

    if feats is not None and scope == "gene":
        exp = kg.merge(feats, on="gene_id", how="inner")
        seen = res[["variant_id", "phenotype_id"]].assign(_s=1)
        miss = exp.merge(seen, on=["variant_id", "phenotype_id"],
                         how="left")
        miss = miss[miss["_s"].isna()].drop(columns="_s")
        if len(miss):
            miss = miss.merge(keys_c[["variant_id", "chrom", "pos"]]
                              .drop_duplicates("variant_id"),
                              on="variant_id", how="left")
            miss["status"] = miss["variant_id"].isin(present_variants).map(
                {True: "not_tested_cis_window",
                 False: "variant_not_in_store"})
            res = pd.concat([res, miss], ignore_index=True)

    res["ancestry"] = anc
    res["modality"] = mod
    o = own[own["modality"] == mod][["variant_id", "phenotype_id"]]
    res = res.merge(o.assign(variant_in_feature_cs=True),
                    on=["variant_id", "phenotype_id"], how="left")
    res["variant_in_feature_cs"] = (res["variant_in_feature_cs"]
                                    .fillna(False).astype(bool))
    vp = res["variant_id"].str.split(":", expand=True)
    res["ref"] = vp[2] if vp.shape[1] > 3 else pd.NA
    res["alt"] = vp[3] if vp.shape[1] > 3 else pd.NA
    for c in out_cols:
        if c not in res.columns:
            res[c] = pd.NA
    for c in ["pos", "ma_count", "ma_samples", "start_distance",
              "end_distance"]:
        if c in res.columns:
            res[c] = pd.to_numeric(res[c]).round().astype("Int64")
    return res[out_cols].sort_values(["pos", "variant_id", "gene_id",
                                      "phenotype_id"])


def run_extract(out_dir: Path, results_dir: Path, anc: str, mod: str,
                qtl_dir, scope: str, backend: str, force: bool):
    tsv, done = task_paths(out_dir, anc, mod)
    if done.exists() and not force:
        print(f"[extract] {anc}/{mod}: done (use --force to rerun)")
        return json.loads(done.read_text())
    wd = work_dir(out_dir)
    keys_path = wd / "query_keys.tsv.gz"
    if not keys_path.exists():
        sys.exit(f"ERROR: {keys_path} missing — run `index` first")
    nominal = results_dir / "nominal" / anc / f"{anc}_{mod}.nominal.tsv.gz"
    if not nominal.exists():
        sys.exit(f"ERROR: nominal store not found: {nominal}")
    if backend == "auto":
        backend = ("tabix" if shutil.which("tabix")
                   and Path(str(nominal) + ".tbi").exists() else "scan")
    print(f"[extract] {anc}/{mod}: backend={backend} scope={scope}")
    if backend == "scan":
        print("  NOTE: scan backend streams the whole nominal store "
              "(slow); load samtools/htslib for tabix")

    keys = pd.read_csv(keys_path, sep="\t",
                       dtype={"chrom": str, "variant_id": str,
                              "gene_id": str})
    own = pd.read_csv(wd / "own_cs_pairs.tsv.gz", sep="\t", dtype=str)
    feats = (_load_features(Path(qtl_dir), anc, mod)
             if qtl_dir and scope == "gene" else None)
    if feats is not None:
        print(f"  {len(feats)} {mod} features in {anc} BED")

    header = read_gz_header(nominal)
    missing_cols = {"chrom", "pos", "phenotype_id", "variant_id",
                    "slope", "slope_se"} - set(header)
    if missing_cols:
        sys.exit(f"ERROR: {nominal} lacks columns {sorted(missing_cols)}")
    stat_cols = [c for c in STAT_COLS if c in header]
    out_cols = (["chrom", "pos", "variant_id", "ref", "alt", "ancestry",
                 "gene_id", "modality", "phenotype_id"] + stat_cols
                + ["status", "variant_in_feature_cs"]
                + (["cs_gene"] if scope == "cis" else []))

    tsv.parent.mkdir(parents=True, exist_ok=True)
    tmp_tsv = tsv.with_suffix(".tsv.partial")
    counts = {"tested": 0, "not_tested_cis_window": 0,
              "variant_not_in_store": 0}
    chroms = sorted(keys["chrom"].unique(), key=chrom_sort_key)
    first = True

    def emit(res):
        nonlocal first
        res.to_csv(tmp_tsv, sep="\t", index=False, header=first,
                   mode="w" if first else "a", na_rep="NA")
        first = False
        for k, v in res["status"].value_counts().items():
            counts[k] = counts.get(k, 0) + int(v)

    with tempfile.TemporaryDirectory(dir=wd) as tmp:
        if backend == "tabix":
            for c in chroms:
                kc = keys[keys["chrom"] == c]
                positions = sorted(kc["pos"].unique())
                parts = [ch for ch in _iter_tabix(nominal, header, c,
                                                  positions, tmp)
                         if len(ch)]
                rows = (pd.concat(parts, ignore_index=True) if parts
                        else pd.DataFrame(columns=header))
                rows["chrom"] = norm_chrom(rows["chrom"])
                emit(_process_chrom(rows, kc, feats, own, scope, anc, mod,
                                    out_cols))
                print(f"  chr{c}: {len(positions)} positions")
        else:
            var_set = set(keys["variant_id"])
            buf = {c: [] for c in chroms}
            for c, sub in _iter_scan(nominal, set(chroms)):
                sub = sub[sub["variant_id"].isin(var_set)]
                if len(sub):
                    buf[c].append(sub)
            for c in chroms:
                rows = (pd.concat(buf[c], ignore_index=True) if buf[c]
                        else pd.DataFrame(columns=header))
                emit(_process_chrom(rows, keys[keys["chrom"] == c], feats,
                                    own, scope, anc, mod, out_cols))

    if first:  # no chromosomes at all: still write a header
        pd.DataFrame(columns=out_cols).to_csv(tmp_tsv, sep="\t",
                                              index=False)
    tmp_tsv.replace(tsv)
    final = _compress(tsv)
    info = {"ancestry": anc, "modality": mod, "scope": scope,
            "backend": backend, "nominal_store": str(nominal),
            "qtl_dir": str(qtl_dir) if qtl_dir else None,
            "output": final.name, "rows": counts,
            "n_variant_gene_keys": int(len(keys)),
            "finished": dt.datetime.now().isoformat(timespec="seconds")}
    done.write_text(json.dumps(info, indent=2))
    print(f"  wrote {final} {counts}")
    return info


def _compress(tsv: Path) -> Path:
    gz = Path(str(tsv) + ".gz")
    if shutil.which("bgzip"):
        subprocess.run(["bgzip", "-f", str(tsv)], check=True)
        if shutil.which("tabix"):
            subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2",
                            "-S", "1", str(gz)], check=True)
    else:
        with open(tsv, "rb") as fi, gzip.open(gz, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        tsv.unlink()
    return gz


def cmd_extract(args):
    run_extract(Path(args.out_dir), Path(args.results_dir), args.ancestry,
                args.modality, args.qtl_dir, args.scope, args.backend,
                args.force)


# ---------------------------------------------------------------------------
# package
# ---------------------------------------------------------------------------

README = """\
Fine-mapped credible-set variants x PANTRY gene features: nominal xQTL stats
===========================================================================
Generated by 06_colocalization_twas/39b_cs_nominal_lookup.py on {date}.

For every variant in ANY SuSHiE 95% credible set (any modality/phenotype),
this archive holds the stage-39 tensorQTL cis.map_nominal statistics of that
variant against every PANTRY feature (all 9 modalities) of the fine-mapped
phenotype's gene, separately per ancestry.

Index parameters
----------------
{params}

Files
-----
cs_membership.tsv.gz
    One row per credible set x member variant.
    chrom, pos, variant_id      variant (chr:pos:REF:ALT, chrom without 'chr')
    gene_id                     gene of the fine-mapped phenotype (prefix
                                before '__', Ensembl version stripped)
    modality, phenotype_id      fine-mapped phenotype
    cs_index, cs_size           SuSHiE credible set index / size
    pip_all, pip_cs             SuSHiE joint PIP / within-CS PIP
    is_cs_lead                  max-pip_all variant of the CS (33 convention)

{{ANC}}/{{ANC}}_{{MOD}}.cs_nominal.tsv.gz  (+ .tbi when bgzipped)
    One row per (CS variant, feature of the CS gene) in modality MOD,
    ancestry ANC. Sorted by chrom/pos; tabix: -s1 -b2 -e2 -S1.
    ref, alt                    parsed from variant_id
    gene_id, modality, phenotype_id   target feature
    slope, slope_se             tensorQTL nominal effect per ALT allele copy
                                on the normalized phenotype, and its SE
    pval_nominal, af (ALT AF in ANC), ma_count, ma_samples,
    start_distance, end_distance      tensorQTL nominal columns
    status
      tested                    statistics present
      not_tested_cis_window     variant is in this ANC x MOD store but this
                                feature's cis window (1 Mb) excludes it
      variant_not_in_store      variant has no row in this ANC x MOD store
                                (below MAF {maf_note} / not in the ANC
                                genotypes, or no MOD feature within 1 Mb)
                                -> statistics NA
    variant_in_feature_cs       TRUE if this variant is itself in a credible
                                set of this exact phenotype
    cs_gene (scope=cis only)    TRUE if the feature's gene is a CS gene of
                                the variant
    A variant in credible sets of several features of one gene appears once
    per target feature; join to cs_membership.tsv.gz on variant_id (+ gene_id)
    for which credible sets it belongs to.

MANIFEST.tsv
    Per ancestry x modality task: row counts by status, backend, inputs.
"""


def cmd_package(args):
    out_dir = Path(args.out_dir)
    params_path = work_dir(out_dir) / "index_params.json"
    params = (json.loads(params_path.read_text())
              if params_path.exists() else {})
    rows, missing = [], []
    for anc in args.ancestries:
        for mod in args.modalities:
            _, done = task_paths(out_dir, anc, mod)
            if not done.exists():
                missing.append(f"{anc}/{mod}")
                continue
            info = json.loads(done.read_text())
            r = {"ancestry": anc, "modality": mod,
                 "file": f"{anc}/{info['output']}",
                 "scope": info["scope"], "backend": info["backend"]}
            r.update({f"n_{k}": v for k, v in info["rows"].items()})
            r["nominal_store"] = info["nominal_store"]
            rows.append(r)
    if missing:
        print(f"  WARNING: {len(missing)} ancestry x modality tasks not done: "
              f"{', '.join(missing)}")
        if not args.allow_missing:
            sys.exit("ERROR: refusing to package an incomplete set "
                     "(pass --allow-missing to override)")
    man = pd.DataFrame(rows)
    man.to_csv(out_dir / "MANIFEST.tsv", sep="\t", index=False)
    date = dt.date.today().strftime("%Y%m%d")
    pjson = json.dumps({k: v for k, v in params.items()}, indent=2)
    if missing:
        pjson += f"\n\nMISSING TASKS: {', '.join(missing)}"
    (out_dir / "README.txt").write_text(README.format(
        date=date, params=pjson, maf_note=args.maf_note))
    archive = Path(args.archive) if args.archive else (
        out_dir.parent / f"{out_dir.name}_{date}.tar.gz")
    with tarfile.open(archive, "w:gz") as tar:
        for p in sorted(out_dir.rglob("*")):
            if p.is_dir() or work_dir(out_dir) in p.parents:
                continue
            tar.add(p, arcname=str(Path(out_dir.name) / p.relative_to(out_dir)))
    print(f"[package] {len(man)} tasks -> {archive} "
          f"({archive.stat().st_size / 1e6:.1f} MB)")
    if len(man):
        print(man.drop(columns=["nominal_store"]).to_string(index=False))


def cmd_all(args):
    cmd_index(args)
    for anc in args.ancestries:
        for mod in args.modalities:
            run_extract(Path(args.out_dir), Path(args.results_dir), anc, mod,
                        args.qtl_dir, args.scope, args.backend, args.force)
    cmd_package(args)


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--out-dir", required=True)
        sp.add_argument("--modalities", nargs="+", default=MODALITIES)

    def index_args(sp):
        sp.add_argument("--finemap-dir", required=True)
        sp.add_argument("--cs-source", choices=["cs", "pips"], default="cs")
        sp.add_argument("--pips", default=None,
                        help="aggregated finemap_pips.tsv.gz "
                             "(default {finemap_dir}/aggregated/...)")
        sp.add_argument("--min-pip", type=float, default=0.0,
                        help="keep CS variants with pip_all >= this")
        sp.add_argument("--lead-only", action="store_true",
                        help="keep only the lead (max PIP) variant per CS")
        sp.add_argument("--threads", type=int, default=16,
                        help="parallel cs.tsv readers (I/O-bound)")

    def extract_args(sp):
        sp.add_argument("--results-dir", required=True)
        sp.add_argument("--qtl-dir", default=None,
                        help="list untested features from {ANC}_{MOD}.bed.gz")
        sp.add_argument("--scope", choices=["gene", "cis"], default="gene",
                        help="gene: features of the CS gene only; cis: every "
                             "feature tested against the variant")
        sp.add_argument("--backend", choices=["auto", "tabix", "scan"],
                        default="auto")
        sp.add_argument("--force", action="store_true")

    def package_args(sp):
        sp.add_argument("--ancestries", nargs="+", default=ANCESTRIES)
        sp.add_argument("--archive", default=None)
        sp.add_argument("--allow-missing", action="store_true")
        sp.add_argument("--maf-note", default="0.01")

    sp = sub.add_parser("index"); common(sp); index_args(sp)
    sp.set_defaults(func=cmd_index)
    sp = sub.add_parser("extract"); common(sp); extract_args(sp)
    sp.add_argument("--ancestry", required=True)
    sp.add_argument("--modality", required=True)
    sp.set_defaults(func=cmd_extract)
    sp = sub.add_parser("package"); common(sp); package_args(sp)
    sp.set_defaults(func=cmd_package)
    sp = sub.add_parser("all"); common(sp); index_args(sp)
    extract_args(sp); package_args(sp)
    sp.set_defaults(func=cmd_all)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    main()
