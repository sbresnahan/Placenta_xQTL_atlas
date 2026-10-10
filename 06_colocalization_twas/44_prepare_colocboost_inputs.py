#!/usr/bin/env python3
"""Prepare one colocBoost shard outside R.

All external command execution lives here on the host compute node. For each
region assigned to the shard, this script uses tabix to slice summary
statistics and plink2 to export the two reference dosage panels. It writes a
small manifest plus ordinary TSV/.raw files for the pure-R colocBoost worker.

The R worker never invokes tabix, plink2, zcat, head, or any shell command.
"""

import argparse
import csv
import gzip
import io
import math
import re
import subprocess
from pathlib import Path


MANIFEST_FIELDS = [
    "region_id", "ancestry", "chrom", "start", "end",
    "prep_status", "prep_message", "sumstats_file",
    "xqtl_raw", "gwas_raw",
]


def open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_header(path):
    with open_text(path) as fh:
        line = fh.readline()
    if not line:
        raise RuntimeError(f"empty summary-statistics file: {path}")
    return line.rstrip("\r\n").split("\t")


def tabix_rows(tabix, path, chrom, start, end):
    header = read_header(path)
    region = f"{chrom}:{start}-{end}"
    proc = subprocess.run(
        [tabix, str(path), region], text=True, capture_output=True, check=False
    )
    if proc.returncode != 0:
        msg = proc.stderr.strip() or proc.stdout.strip()
        raise RuntimeError(msg or f"tabix exited {proc.returncode}")
    if not proc.stdout.strip():
        return []
    return list(csv.DictReader(io.StringIO(proc.stdout), fieldnames=header,
                               delimiter="\t"))


def fnum(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def psam_n(pgen):
    path = Path(f"{pgen}.psam")
    if not path.exists():
        return None
    n = 0
    with path.open() as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n


def sanitize_id(value):
    return re.sub(r"[^A-Za-z0-9._+-]", "_", str(value))


def slice_outcome(tabix, oc, chrom, start, end, n_xqtl):
    rows = tabix_rows(tabix, oc["file"], chrom, start, end)
    if not rows:
        return []

    out = []
    seen = set()
    if oc["type"] == "xqtl":
        for r in rows:
            if r.get("phenotype_id") != oc.get("phenotype_id"):
                continue
            variant = r.get("variant_id", "")
            beta = fnum(r.get("slope"))
            se = fnum(r.get("slope_se"))
            if (not variant or variant in seen or not math.isfinite(beta)
                    or not math.isfinite(se) or se <= 0):
                continue
            seen.add(variant)
            out.append({"outcome_name": oc["outcome_name"], "type": "xqtl",
                        "variant": variant, "beta": beta, "sebeta": se,
                        "n": n_xqtl})
    else:
        for r in rows:
            variant = r.get("var_id", "")
            beta = fnum(r.get("beta"))
            se = fnum(r.get("se"))
            n = fnum(r.get("n"))
            if (not variant or variant in seen or not math.isfinite(beta)
                    or not math.isfinite(se) or se <= 0):
                continue
            seen.add(variant)
            out.append({"outcome_name": oc["outcome_name"], "type": "gwas",
                        "variant": variant, "beta": beta, "sebeta": se,
                        "n": n})
    return out


def run_export(plink2, pgen, keep, variants, prefix):
    extract = Path(str(prefix) + ".extract.txt")
    extract.write_text("".join(f"{v}\n" for v in variants))
    cmd = [plink2, "--pfile", pgen, "--extract", str(extract),
           "--export", "A", "--threads", "1", "--silent",
           "--out", str(prefix)]
    if keep and Path(keep).exists():
        cmd[3:3] = ["--keep", keep]
    proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    raw = Path(str(prefix) + ".raw")
    if proc.returncode != 0 or not raw.exists():
        msg = proc.stderr.strip() or proc.stdout.strip()
        return None, (msg or f"plink2 exited {proc.returncode}")[-1000:]
    return raw, ""


def unique_preserving(values):
    seen = set()
    out = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def write_sumstats(rows, path):
    fields = ["outcome_name", "type", "variant", "beta", "sebeta", "n"]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def prepare_region(reg, reg_outcomes, work_dir, out_dir, tabix, plink2,
                   ld_xqtl_pgen, ld_gwas_pgen, ld_gwas_keep, min_outcomes):
    rec = {k: "" for k in MANIFEST_FIELDS}
    rec.update({k: reg.get(k, "") for k in
                ["region_id", "ancestry", "chrom", "start", "end"]})

    def finish(status, message=""):
        rec["prep_status"] = status
        rec["prep_message"] = str(message).replace("\t", " ").replace("\n", " ")[:500]
        return rec

    done = out_dir / "results" / reg["ancestry"] / f"{reg['region_id']}.done"
    if done.exists():
        return finish("already_done")

    region_dir = work_dir / sanitize_id(reg["region_id"])
    region_dir.mkdir(parents=True, exist_ok=True)
    n_xqtl = psam_n(ld_xqtl_pgen)
    if n_xqtl is None or n_xqtl < 10:
        return finish("error", "could not read xQTL N from psam")

    try:
        by_outcome = []
        for oc in reg_outcomes:
            rows = slice_outcome(tabix, oc, reg["chrom"], reg["start"],
                                 reg["end"], n_xqtl)
            if rows:
                by_outcome.append((oc, rows))
    except Exception as exc:
        return finish("error", exc)

    if len(by_outcome) < min_outcomes:
        return finish("too_few_outcomes")
    if not any(oc["type"] == "gwas" for oc, _ in by_outcome):
        return finish("too_few_outcomes", "no GWAS outcome with variants")

    rows_all = [r for _, rows in by_outcome for r in rows]
    vars_x = unique_preserving(
        r["variant"] for oc, rows in by_outcome if oc["type"] == "xqtl" for r in rows
    )
    vars_g = unique_preserving(
        r["variant"] for oc, rows in by_outcome if oc["type"] == "gwas" for r in rows
    )
    if not vars_x:
        return finish("too_few_outcomes", "no xQTL variants")
    if not vars_g:
        return finish("too_few_outcomes", "no GWAS variants")

    xraw, msg = run_export(plink2, ld_xqtl_pgen, "", vars_x,
                           region_dir / "xref_xqtl")
    if xraw is None:
        return finish("dosage_failed_xqtl", msg)
    graw, msg = run_export(plink2, ld_gwas_pgen, ld_gwas_keep, vars_g,
                           region_dir / "xref_gwas")
    if graw is None:
        return finish("dosage_failed_gwas", msg)

    sumstats = region_dir / "sumstats.tsv"
    write_sumstats(rows_all, sumstats)
    rec.update({"prep_status": "ok", "sumstats_file": str(sumstats),
                "xqtl_raw": str(xraw), "gwas_raw": str(graw)})
    return rec


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--regions", required=True)
    p.add_argument("--outcomes", required=True)
    p.add_argument("--shard-index", required=True, type=int)
    p.add_argument("--n-shards", required=True, type=int)
    p.add_argument("--outdir", required=True)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--ld-xqtl-pgen", required=True)
    p.add_argument("--ld-gwas-pgen", required=True)
    p.add_argument("--ld-gwas-keep", default="")
    p.add_argument("--tabix", default="tabix")
    p.add_argument("--plink2", default="plink2")
    p.add_argument("--min-outcomes", type=int, default=2)
    args = p.parse_args()

    regions = read_tsv(args.regions)
    outcomes = read_tsv(args.outcomes)
    selected = [r for i, r in enumerate(regions)
                if (i % args.n_shards) + 1 == args.shard_index]
    if not regions:
        raise SystemExit(f"ERROR: empty regions file: {args.regions}")

    outcomes_by_region = {}
    for oc in outcomes:
        outcomes_by_region.setdefault(oc["region_id"], []).append(oc)

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    out_dir = Path(args.outdir)
    prepared = []
    for i, reg in enumerate(selected, 1):
        print(f"[prep {i}/{len(selected)}] {reg['region_id']}", flush=True)
        prepared.append(prepare_region(
            reg, outcomes_by_region.get(reg["region_id"], []), work_dir, out_dir,
            args.tabix, args.plink2, args.ld_xqtl_pgen, args.ld_gwas_pgen,
            args.ld_gwas_keep, args.min_outcomes))

    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(prepared)
    print(f"prepared {len(prepared)} regions -> {manifest}")


if __name__ == "__main__":
    main()
