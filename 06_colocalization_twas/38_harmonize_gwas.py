#!/usr/bin/env python3
"""38_harmonize_gwas.py — harmonize raw GWAS summary statistics to GRCh38
and to the pooled-pgen variant ID/allele conventions.

Objective 1.6 colocalization requires every GWAS to share the xQTL variant
key space (`chr:pos:ref:alt`, bare numeric chromosomes, GRCh38) with alleles
oriented consistently. For each catalog trait this script:

  1. Reads the raw file with a tolerant column-alias map (per-trait overrides
     via --column-map trait:field=name). If {trait_id}.txt.gz is absent but
     {trait_id}.vcf.gz exists, the file is parsed as GWAS-VCF (OpenGWAS/IEU
     spec: FORMAT keys ES/SE/LP/AF/ID; effect allele = ALT; LP = -log10 p).
  2. Lifts hg19 -> GRCh38 positions with a pure-Python chain-file mapper
     (chain from module 01; negative-strand chains complement alleles).
     Skipped for build-38 traits.
  3. Standardizes to: chr, pos, rsid, effect_allele, other_allele, eaf,
     beta, se, pval, n. Odds ratios are log-transformed (flagged in QC).
  4. Aligns alleles against the union of variants in the pooled pgens
     ({ANC}_pooled.pvar): matches on chr:pos + allele set, flips beta/eaf
     when the GWAS effect allele equals the pgen REF, complements alleles
     for negative-strand chain mappings, and drops unresolvable palindromic
     variants (A/T or C/G with 0.4 <= eaf <= 0.6).
  5. Writes {out_dir}/{trait_id}.sumstats.tsv.gz (bgzip+tabix if available)
     and a per-trait QC row in {out_dir}/harmonization_qc.tsv.

Only variants matching a pgen position are kept — coloc/TWAS only ever use
variants present in the xQTL genotype data, so unmatched GWAS variants are
dead weight (counted in QC).
"""

import argparse
import bisect
import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# column aliases (lowercase)
# ---------------------------------------------------------------------------
ALIASES = {
    "chr": ["chr", "chromosome", "#chr", "chrom", "#chrom", "hg19chrc",
            "chromosome_name", "hm_chrom"],
    "pos": ["pos", "position", "bp", "base_pair_location", "hm_pos",
            "pos_grch38", "position_grch38"],
    "rsid": ["rsid", "snp", "rs_id", "markername", "variant_id",
             "hm_rsid", "rsids", "snpid"],
    "effect_allele": ["effect_allele", "ea", "a1", "allele1", "tested_allele",
                      "hm_effect_allele", "effectallele", "alt"],
    "other_allele": ["other_allele", "nea", "a2", "allele2",
                     "non_effect_allele", "hm_other_allele", "otherallele",
                     "ref"],
    "eaf": ["eaf", "effect_allele_frequency", "freq", "eaf1", "freq1",
            "hm_effect_allele_frequency", "frq", "af"],
    "beta": ["beta", "b", "effect", "estimate", "hm_beta", "effect_size"],
    "or": ["or", "odds_ratio", "hm_odds_ratio"],
    "se": ["se", "standard_error", "stderr", "hm_se", "sebeta"],
    "pval": ["p", "pval", "p_value", "pvalue", "pval_nominal", "hm_p_value",
             "p-value", "p.value"],
    "n": ["n", "sample_size", "n_total", "samples", "hm_n", "n_samples",
          "total_sample_size"],
}

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


# ---------------------------------------------------------------------------
# GWAS-VCF reader (OpenGWAS/IEU spec)
# ---------------------------------------------------------------------------
def read_gwas_vcf(path):
    """Parse a GWAS-VCF into a DataFrame with canonical column names that
    detect_columns()/standardize() already understand.

    FORMAT keys: ES (effect size of ALT), SE, LP (-log10 p), AF (ALT allele
    frequency), ID (rsid; falls back to the VCF ID column). Multiallelic
    rows (',' in ALT) are skipped — pgen alignment is biallelic-only.
    """
    rows = []
    n_multi = n_bad = 0
    idx = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 9)
            if len(f) < 10:
                n_bad += 1
                continue
            chrom, pos, vid, ref, alt = f[0], f[1], f[2], f[3], f[4]
            if "," in alt:
                n_multi += 1
                continue
            if idx is None:
                idx = {k: i for i, k in enumerate(f[8].split(":"))}
            vals = f[9].split("\t")[0].split(":")

            def get(key):
                i = idx.get(key)
                return vals[i] if i is not None and i < len(vals) else ""

            try:
                beta = float(get("ES"))
                se = float(get("SE"))
            except ValueError:
                n_bad += 1
                continue
            lp = get("LP")
            try:
                pval = 10.0 ** (-float(lp))
            except ValueError:
                pval = float("nan")
            rows.append((chrom, pos, get("ID") or vid, alt, ref, get("AF"),
                         beta, se, pval))
    if n_multi:
        print(f"    [vcf] skipped {n_multi} multiallelic row(s)")
    if n_bad:
        print(f"    [vcf] skipped {n_bad} malformed row(s)")
    return pd.DataFrame(rows, columns=["chr", "pos", "rsid",
                                       "effect_allele", "other_allele",
                                       "eaf", "beta", "se", "pval"])


# ---------------------------------------------------------------------------
# chain-file liftover (single-position, pure python)
# ---------------------------------------------------------------------------
class ChainMapper:
    """Map single-base positions through a UCSC chain file.

    Chains map target (e.g. hg19) -> query (e.g. hg38). For each alignment
    block, target [t, t+size) maps to query [q, q+size) on qStrand; with
    qStrand == '-', q coordinates run backwards within the reverse-complement
    block, and alleles must be complemented.
    """

    def __init__(self, chain_path):
        self.blocks = {}  # chrom -> (t_starts array, records list)
        cur = None
        opener = gzip.open if str(chain_path).endswith(".gz") else open
        with opener(chain_path, "rt") as f:
            for line in f:
                if not line.strip():
                    continue
                if line.startswith("chain"):
                    # chain score tName tSize tStrand tStart tEnd
                    #       qName qSize qStrand qStart qEnd id
                    parts = line.split()
                    cur = {"t_name": parts[2], "t": int(parts[5]),
                           "q_name": parts[7], "q_size": int(parts[8]),
                           "q_strand": parts[9], "q": int(parts[10])}
                    continue
                if cur is None:
                    continue
                parts = line.split()
                size = int(parts[0])
                rec = (cur["t"], cur["q"], size, cur["q_name"],
                       cur["q_strand"], cur["q_size"])
                self.blocks.setdefault(cur["t_name"], []).append(rec)
                cur["t"] += size
                cur["q"] += size
                if len(parts) == 3:
                    # gap: dt unaligned target bases, dq unaligned query
                    # bases (q accumulates identically on both strands —
                    # minus-strand qStart is already in reverse coordinates)
                    cur["t"] += int(parts[1])
                    cur["q"] += int(parts[2])
                # len==1 terminates the chain
        for chrom, recs in self.blocks.items():
            recs.sort(key=lambda r: r[0])
            self.blocks[chrom] = (np.array([r[0] for r in recs]), recs)

    def map_pos(self, chrom, pos1):
        """Map a 1-based position. Returns (q_chrom, q_pos1, q_strand) or
        None. Tries the chromosome name as given and with/without 'chr'."""
        for name in (chrom, "chr" + chrom if not chrom.startswith("chr")
                     else chrom[3:]):
            if name in self.blocks:
                chrom_key = name
                break
        else:
            return None
        starts, recs = self.blocks[chrom_key]
        pos0 = pos1 - 1  # chain coords are 0-based half-open
        i = bisect.bisect_right(starts, pos0) - 1
        if i < 0:
            return None
        t, q, size, q_name, q_strand, q_size = recs[i]
        if not (t <= pos0 < t + size):
            return None
        off = pos0 - t
        if q_strand == "+":
            return q_name, q + off + 1, "+"
        return q_name, q_size - (q + off), "-"


def complement(allele):
    return "".join(COMPLEMENT.get(b, "N") for b in allele.upper())


# ---------------------------------------------------------------------------
# schema standardization
# ---------------------------------------------------------------------------
def detect_columns(df, overrides):
    """Map canonical fields to actual columns via aliases + overrides."""
    lower = {c.lower(): c for c in df.columns}
    mapping = {}
    for field, aliases in ALIASES.items():
        if field in overrides:
            if overrides[field] not in df.columns:
                raise ValueError(f"override column '{overrides[field]}' "
                                 f"for {field} not in file")
            mapping[field] = overrides[field]
            continue
        for a in aliases:
            if a in lower:
                mapping[field] = lower[a]
                break
    for required in ("pos", "effect_allele", "other_allele"):
        if required not in mapping:
            raise ValueError(f"cannot identify a '{required}' column; "
                             f"pass --column-map trait:{required}=NAME. "
                             f"Columns: {list(df.columns)[:15]}")
    if "beta" not in mapping and "or" not in mapping:
        raise ValueError("need a beta or odds-ratio column")
    if "chr" not in mapping and "rsid" not in mapping:
        raise ValueError("need a chromosome column (rsid-only files are "
                         "not supported)")
    return mapping


def standardize(df, mapping):
    out = pd.DataFrame()
    out["chr"] = (df[mapping["chr"]].astype(str)
                  .str.replace("^chr", "", regex=True)
                  .replace({"X": "23", "Y": "24", "M": "25", "MT": "25"}))
    out = out[out["chr"].str.fullmatch(r"\d+")]
    out["pos"] = pd.to_numeric(df.loc[out.index, mapping["pos"]],
                               errors="coerce")
    out["rsid"] = (df.loc[out.index, mapping["rsid"]].astype(str)
                   if "rsid" in mapping else "")
    out["effect_allele"] = (df.loc[out.index, mapping["effect_allele"]]
                            .astype(str).str.upper())
    out["other_allele"] = (df.loc[out.index, mapping["other_allele"]]
                           .astype(str).str.upper())
    for field in ("eaf", "se", "pval", "n"):
        out[field] = (pd.to_numeric(df.loc[out.index, mapping[field]],
                                    errors="coerce")
                      if field in mapping else np.nan)
    used_or = "beta" not in mapping
    if used_or:
        out["beta"] = np.log(pd.to_numeric(df.loc[out.index, mapping["or"]],
                                           errors="coerce"))
    else:
        out["beta"] = pd.to_numeric(df.loc[out.index, mapping["beta"]],
                                    errors="coerce")
    out = out.dropna(subset=["pos", "beta", "se"])
    out["pos"] = out["pos"].astype(int)
    # SNVs only (pgen variant space is SNV-dominated; indels rarely align)
    out = out[out["effect_allele"].str.fullmatch(r"[ACGT]") &
              out["other_allele"].str.fullmatch(r"[ACGT]")]
    return out, used_or


# ---------------------------------------------------------------------------
# pgen variant index
# ---------------------------------------------------------------------------
def load_pgen_index(pgen_dir, ancestries):
    """chr -> pos -> list of (ref, alt) across the pooled pgens."""
    index = {}
    for anc in ancestries:
        pvar = Path(pgen_dir) / f"{anc}_pooled.pvar"
        if not pvar.exists():
            print(f"  WARNING: {pvar} not found; skipping {anc}")
            continue
        # pvar files may carry extra columns beyond CHROM/POS/ID/REF/ALT
        # (QUAL/FILTER/INFO/CM); positional names= then silently misaligns
        # (pandas shifts surplus leading columns into the index). Read the
        # #CHROM header line and select columns by name instead.
        hdr = None
        with open(pvar) as fh:
            for i, line in enumerate(fh):
                if line.startswith("#CHROM"):
                    hdr = i
                    break
                if not line.startswith("##"):
                    break
        if hdr is not None:
            pv = pd.read_csv(pvar, sep="\t", skiprows=hdr, dtype=str)
            pv.columns = ["chr" if c == "#CHROM" else c.lower()
                          for c in pv.columns]
        else:  # headerless pvar: fall back to the 5-column layout
            pv = pd.read_csv(pvar, sep="\t", comment="#",
                             names=["chr", "pos", "id", "ref", "alt"],
                             dtype={"chr": str})
        pv["chr"] = pv["chr"].str.replace("^chr", "", regex=True)
        for chrom, sub in pv.groupby("chr"):
            d = index.setdefault(chrom, {})
            for pos, ref, alt in zip(sub["pos"].values, sub["ref"].values,
                                     sub["alt"].values):
                d.setdefault(int(pos), set()).add((ref, alt))
        print(f"  {anc}: {len(pv)} variants indexed")
    return index


def align_to_pgen(df, index):
    """Align GWAS alleles to pgen REF/ALT at matched positions.

    Returns (aligned_df, qc_counts). beta/eaf are flipped when the GWAS
    effect allele matches the pgen REF. Palindromic SNPs with mid-range EAF
    are dropped as strand-unresolvable.
    """
    qc = {"n_in": len(df), "n_pos_matched": 0, "n_allele_matched": 0,
          "n_flipped": 0, "n_palindromic_dropped": 0, "n_unmatched": 0}
    keep_rows = []
    for chrom, sub in df.groupby("chr", sort=False):
        pos_index = index.get(chrom)
        if pos_index is None:
            qc["n_unmatched"] += len(sub)
            continue
        for row in sub.itertuples():
            alts = pos_index.get(row.pos)
            if alts is None:
                qc["n_unmatched"] += 1
                continue
            qc["n_pos_matched"] += 1
            ea, oa = row.effect_allele, row.other_allele
            matched = False
            for ref, alt in alts:
                if ea == alt and oa == ref:
                    keep_rows.append((row, 1.0))
                    matched = True
                    break
                if ea == ref and oa == alt:
                    keep_rows.append((row, -1.0))
                    matched = True
                    break
            if not matched:
                # try strand complement (residual strand flips)
                ea_c, oa_c = complement(ea), complement(oa)
                for ref, alt in alts:
                    if ea_c == alt and oa_c == ref:
                        keep_rows.append((row, 1.0))
                        matched = True
                        break
                    if ea_c == ref and oa_c == alt:
                        keep_rows.append((row, -1.0))
                        matched = True
                        break
            if not matched:
                qc["n_unmatched"] += 1
    qc["n_allele_matched"] = len(keep_rows)
    if not keep_rows:
        return df.iloc[0:0], qc
    rows, signs = zip(*keep_rows)
    out = pd.DataFrame(rows)
    signs = np.array(signs)
    # palindromic filter
    pal = ((out["effect_allele"].isin(["A", "T"]) &
            out["other_allele"].isin(["A", "T"])) |
           (out["effect_allele"].isin(["C", "G"]) &
            out["other_allele"].isin(["C", "G"])))
    mid = out["eaf"].between(0.4, 0.6).fillna(False)
    drop = pal & mid
    qc["n_palindromic_dropped"] = int(drop.sum())
    out = out[~drop]
    signs = signs[~drop.values]
    out["beta"] = out["beta"] * signs
    out.loc[signs < 0, "eaf"] = 1 - out.loc[signs < 0, "eaf"]
    qc["n_flipped"] = int((signs < 0).sum())
    out["var_id"] = (out["chr"] + ":" + out["pos"].astype(str) + ":" +
                     out["other_allele"] + ":" + out["effect_allele"])
    return out, qc


# ---------------------------------------------------------------------------
def harmonize_trait(row, raw_dir, out_dir, mapper, index, overrides):
    tid = row["trait_id"]
    raw = Path(raw_dir) / f"{tid}.txt.gz"
    raw_vcf = Path(raw_dir) / f"{tid}.vcf.gz"
    if not raw.exists() and raw_vcf.exists():
        print(f"    [vcf] parsing GWAS-VCF: {raw_vcf.name}")
        df = read_gwas_vcf(raw_vcf)
        input_kind = "gwas-vcf"
    elif not raw.exists():
        return {"trait_id": tid, "status": "missing_raw"}
    else:
        # VCF-style tabular releases (e.g. JECS) carry a '#CHROM' header
        # line; comment='#' would silently drop it and promote the first
        # data row to column names.
        with gzip.open(raw, "rt") as fh:
            first_line = fh.readline()
        comment = None if first_line.upper().startswith("#CHROM") else "#"
        df = pd.read_csv(raw, sep=None, engine="python", nrows=None,
                         comment=comment)
        input_kind = "tabular"
    qc = {"trait_id": tid, "status": "ok", "n_raw": len(df),
          "input": input_kind}
    mapping = detect_columns(df, overrides.get(tid, {}))
    qc["column_map"] = ";".join(f"{k}={v}" for k, v in sorted(mapping.items()))
    df, used_or = standardize(df, mapping)
    qc["or_to_beta"] = used_or
    qc["n_schema"] = len(df)
    # GWAS-VCF carries no per-variant N; fill from the catalog sample_size
    # (41_susie_coloc.R takes median(n) for coloc's sample-size argument).
    if df["n"].isna().all():
        ss = str(row.get("sample_size", "")).strip()
        if ss and ss.lower() not in ("", "nan"):
            df["n"] = float(ss)
            qc["n_from_catalog"] = ss

    if str(row["build"]) == "37":
        if mapper is None:
            raise RuntimeError(f"{tid} is build 37 but no --chain given")
        mapped = [mapper.map_pos(c, p) for c, p in zip(df["chr"], df["pos"])]
        ok = [m is not None for m in mapped]
        qc["n_lifted"] = int(sum(ok))
        df = df[ok].copy()
        mapped_ok = [m for m in mapped if m is not None]
        df["pos"] = [m[1] for m in mapped_ok]
        neg = np.array([m[2] == "-" for m in mapped_ok])
        if neg.any():
            df.loc[neg, "effect_allele"] = df.loc[neg, "effect_allele"].map(complement)
            df.loc[neg, "other_allele"] = df.loc[neg, "other_allele"].map(complement)
        qc["n_neg_strand"] = int(neg.sum()) if len(mapped_ok) else 0
    else:
        qc["n_lifted"] = len(df)

    df, aqc = align_to_pgen(df, index)
    qc.update(aqc)
    df["trait_id"] = tid
    df = df.sort_values(["chr", "pos"],
                        key=lambda s: pd.to_numeric(s, errors="coerce"))
    cols = ["chr", "pos", "rsid", "effect_allele", "other_allele",
            "eaf", "beta", "se", "pval", "n", "var_id", "trait_id"]
    out_plain = Path(out_dir) / f"{tid}.sumstats.tsv"
    df.to_csv(out_plain, sep="\t", index=False, columns=cols)
    out_path = Path(str(out_plain) + ".gz")
    try:
        # bgzip + tabix so 41_susie_coloc.R can slice loci without reading
        # the whole genome (chr = col 1, pos = col 2)
        subprocess.run(["bgzip", "-f", str(out_plain)], check=True)
        subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2",
                        "-S", "1", str(out_path)], check=True)
        qc["tabix"] = True
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        print(f"  WARNING: bgzip/tabix failed ({e}); writing plain gzip "
              f"(locus slicing will fall back to full-file reads)")
        df.to_csv(out_path, sep="\t", index=False, columns=cols)
        out_plain.unlink(missing_ok=True)
        qc["tabix"] = False
    qc["file"] = str(out_path)
    return qc


def parse_column_map(entries):
    """Parse --column-map trait:field=name triples."""
    out = {}
    for e in entries or []:
        trait_field, _, name = e.partition("=")
        trait, _, field = trait_field.partition(":")
        if not (trait and field and name):
            sys.exit(f"ERROR: malformed --column-map entry: {e}")
        out.setdefault(trait, {})[field] = name
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--catalog", required=True)
    p.add_argument("--raw-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--chain", default=None,
                   help="hg19->hg38 chain file (required if any trait is "
                        "build 37). Module-01 chain: "
                        "/home/stbresnahan/bhattacharya_lab/data/"
                        "GenomicReferences/liftover/hg19ToHg38.over.chain")
    p.add_argument("--pgen-dir", required=True,
                   help="directory with {ANC}_pooled.pvar")
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    p.add_argument("--traits", nargs="*", default=None)
    p.add_argument("--column-map", nargs="*", default=None,
                   help="per-trait column overrides: trait:field=name")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cat = pd.read_csv(args.catalog, sep="\t")
    if args.traits:
        cat = cat[cat["trait_id"].isin(args.traits)]

    mapper = None
    if any(str(b) == "37" for b in cat["build"]):
        if not args.chain or not Path(args.chain).exists():
            sys.exit("ERROR: build-37 traits present; --chain required")
        print(f"[harmonize] loading chain: {args.chain}")
        mapper = ChainMapper(args.chain)

    print(f"[harmonize] indexing pooled pgen variants from {args.pgen_dir}")
    index = load_pgen_index(args.pgen_dir, args.ancestries)
    overrides = parse_column_map(args.column_map)

    qcs = []
    for _, row in cat.iterrows():
        print(f"[harmonize] {row['trait_id']} (build {row['build']})")
        qc = harmonize_trait(row, args.raw_dir, out_dir, mapper, index,
                             overrides)
        qcs.append(qc)
        print("  " + ", ".join(f"{k}={v}" for k, v in qc.items()
                               if k not in ("column_map", "file")))
    qc_df = pd.DataFrame(qcs)
    qc_path = out_dir / "harmonization_qc.tsv"
    qc_df.to_csv(qc_path, sep="\t", index=False)
    print(f"[harmonize] QC -> {qc_path}")


if __name__ == "__main__":
    main()
