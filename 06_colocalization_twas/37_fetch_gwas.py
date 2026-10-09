#!/usr/bin/env python3
"""37_fetch_gwas.py — download Table-2 GWAS summary statistics (login node).

Objective 1.6 colocalization/TWAS consumes the GWAS catalogued in
gwas_catalog.tsv. Compute nodes on seadragon have no internet access, so
this runs on a login node via 37_fetch_gwas.sh.

Access modes (catalog `access` column)
--------------------------------------
direct      : `source_url` is the file itself.
page_scrape : `source_url` is an HTML page; the first href whose target or
              link text contains `file_pattern` (case-insensitive) is
              downloaded. Used for the EGG consortium pages, which are
              Wix-hosted and do not expose stable file URLs.
gcst        : `source_url` is a GWAS Catalog accession (GCST...). The
              harmonised (GRCh38) file under the FTP tree is preferred; if
              absent, the raw submission file is taken and flagged for
              liftover in 38_harmonize_gwas.py.
manual      : not downloadable programmatically (request form / registered
              portal). Instructions are printed and the script verifies the
              expected manually-placed file in the raw directory.

Outputs
-------
{out_dir}/{trait_id}.txt.gz (+ .source_url sidecar) for fetched traits;
{out_dir}/fetch_manifest.tsv with per-trait status, bytes, and resolved URL.

Manual-placement contract: place the file at {out_dir}/{trait_id}.txt.gz
(gzipped text) or {out_dir}/{trait_id}.vcf.gz (GWAS-VCF, e.g. OpenGWAS
downloads; 38_harmonize_gwas.py parses ES/SE/LP/AF FORMAT fields directly).
For JECS age-stratified releases, add one catalog row per
age bin (trait_id jecs_childhood_bmi_2025_<agebin>) and place each file
accordingly.
"""

import argparse
import gzip
import re
import sys
import urllib.request
import urllib.error
from html.parser import HTMLParser
from pathlib import Path

import pandas as pd

GCST_FTP = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            "{bucket}/{gcst}/")
UA = {"User-Agent": "placenta-xqtl-atlas/1.0 (academic research)"}


class HrefParser(HTMLParser):
    """Collect (href, link_text) pairs from an HTML page."""

    def __init__(self):
        super().__init__()
        self.links = []
        self._cur = None

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self._cur = dict(attrs).get("href")

    def handle_data(self, data):
        if self._cur is not None:
            self.links.append((self._cur, data.strip()))

    def handle_endtag(self, tag):
        if tag == "a":
            self._cur = None


def http_get(url, dest=None, timeout=300):
    """GET url -> dest file (streaming) or return text."""
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        if dest is None:
            return r.read().decode("utf-8", errors="replace")
        n = 0
        with open(dest, "wb") as f:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
                n += len(chunk)
        return n


def resolve_page_scrape(row):
    """Find the download href on an HTML page matching file_pattern."""
    html = http_get(row["source_url"])
    parser = HrefParser()
    parser.feed(html)
    pat = row["file_pattern"].lower()
    cands = [(h, t) for h, t in parser.links
             if h and (pat in h.lower() or pat in t.lower())]
    # EGG pages link the phrase "downloaded here"; the href target carries
    # the file name. Prefer hrefs that look like files.
    file_like = [c for c in cands if re.search(r"\.(gz|txt|tsv|zip)($|\?)",
                                               c[0], re.I)]
    for href, text in file_like or cands:
        if href.startswith("http"):
            return href
        # protocol-relative or root-relative
        if href.startswith("//"):
            return "https:" + href
        if href.startswith("/"):
            m = re.match(r"(https?://[^/]+)", row["source_url"])
            return m.group(1) + href
    raise RuntimeError(
        f"no href matching '{row['file_pattern']}' on {row['source_url']}")


def resolve_gcst(row):
    """Resolve a GWAS Catalog accession to a summary-statistics file URL."""
    gcst = row["source_url"].strip()
    m = re.fullmatch(r"GCST(\d+)", gcst)
    if not m:
        raise RuntimeError(f"gcst access expects a GCST accession, got {gcst}")
    num = int(m.group(1))
    bucket = f"GCST{num // 1000 * 1000:06d}-GCST{num // 1000 * 1000 + 999:06d}"
    base = GCST_FTP.format(bucket=bucket, gcst=gcst)
    # Prefer the harmonised (GRCh38) file.
    for sub, pat in (("harmonised/", r"harmonised.*\.h\.tsv\.gz|\.h\.tsv\.gz"),
                     ("harmonised/", r"\.tsv\.gz"),
                     ("", r"\.tsv\.gz|\.txt\.gz|\.tsv\.bgz")):
        try:
            html = http_get(base + sub)
        except urllib.error.URLError:
            continue
        hrefs = re.findall(r'href="([^"]+)"', html)
        files = [h for h in hrefs if re.search(pat, h, re.I)
                 and not h.endswith(".tbi")]
        files = [h for h in files if not h.startswith("?")]
        if files:
            return base + sub + sorted(files)[0]
    raise RuntimeError(f"no summary-statistics file found under {base}")


def fetch_trait(row, out_dir, force=False):
    """Fetch one catalog row. Returns a manifest record dict."""
    tid = row["trait_id"]
    rec = {"trait_id": tid, "access": row["access"], "status": "pending",
           "resolved_url": "", "bytes": 0, "file": ""}
    dest = out_dir / f"{tid}.txt.gz"
    rec["file"] = str(dest)

    if row["access"] == "manual":
        # accept gzipped tabular ({tid}.txt.gz) or GWAS-VCF ({tid}.vcf.gz,
        # e.g. OpenGWAS downloads); 38_harmonize_gwas.py parses both
        for cand in (dest, out_dir / f"{tid}.vcf.gz"):
            if cand.exists() and cand.stat().st_size > 0:
                rec.update(status="present_manual",
                           bytes=cand.stat().st_size, file=str(cand))
                return rec
        rec["status"] = "awaiting_manual"
        return rec

    if dest.exists() and dest.stat().st_size > 0 and not force:
        rec.update(status="skipped_done", bytes=dest.stat().st_size)
        return rec

    try:
        if row["access"] == "direct":
            url = row["source_url"]
        elif row["access"] == "page_scrape":
            url = resolve_page_scrape(row)
        elif row["access"] == "gcst":
            url = resolve_gcst(row)
        else:
            raise RuntimeError(f"unknown access mode: {row['access']}")
        rec["resolved_url"] = url
        tmp = dest.with_suffix(".part")
        nbytes = http_get(url, dest=tmp)
        # sanity: must be gzip and non-trivial
        with open(tmp, "rb") as f:
            magic = f.read(2)
        if magic != b"\x1f\x8b":
            # some servers return plain text; gzip it ourselves
            raw = tmp.read_bytes()
            with gzip.open(dest, "wt") as f:
                f.write(raw.decode("utf-8", errors="replace"))
            tmp.unlink()
        else:
            tmp.rename(dest)
        rec.update(status="fetched", bytes=dest.stat().st_size)
        Path(str(dest) + ".source_url").write_text(url + "\n")
    except Exception as e:  # report, don't abort the batch
        rec["status"] = f"failed: {e}"
    return rec


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--catalog", required=True, help="gwas_catalog.tsv")
    p.add_argument("--out-dir", required=True, help="raw download directory")
    p.add_argument("--traits", nargs="*", default=None,
                   help="restrict to these trait_ids (default: all)")
    p.add_argument("--force", action="store_true",
                   help="re-download traits with existing files")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cat = pd.read_csv(args.catalog, sep="\t")
    if args.traits:
        cat = cat[cat["trait_id"].isin(args.traits)]
        if cat.empty:
            sys.exit("ERROR: no catalog rows match --traits")

    print(f"[fetch] {len(cat)} catalog traits -> {out_dir}")
    records = []
    for _, row in cat.iterrows():
        rec = fetch_trait(row, out_dir, force=args.force)
        records.append(rec)
        print(f"  {rec['trait_id']}: {rec['status']}"
              + (f" ({rec['bytes'] / 1e6:.1f} MB)" if rec["bytes"] else ""))
        if rec["status"] == "awaiting_manual":
            print(f"    MANUAL PLACEMENT NEEDED: see catalog notes;")
            print(f"    place the file at {rec['file']}")
            print(f"    source: {row['source_url']}")

    man = pd.DataFrame(records)
    man_path = out_dir / "fetch_manifest.tsv"
    man.to_csv(man_path, sep="\t", index=False)
    n_ok = man["status"].isin(["fetched", "skipped_done",
                               "present_manual"]).sum()
    print(f"[fetch] {n_ok}/{len(man)} traits available; manifest -> {man_path}")
    failed = man[man["status"].str.startswith("failed")]
    if len(failed):
        print("  FAILED traits:", ", ".join(failed["trait_id"]), file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
