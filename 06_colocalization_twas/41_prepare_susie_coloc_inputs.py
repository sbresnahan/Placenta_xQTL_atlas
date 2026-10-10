#!/usr/bin/env python3
"""Prepare one SuSiE-coloc shard for the pure-R model worker.

All external command-line I/O belongs here, outside R.  For each task assigned
to this shard the script:
  * extracts the xQTL and GWAS locus with tabix;
  * merges/filters shared variants;
  * exports ALT-coded dosages from each PGEN with plink2;
  * computes signed LD from those dosages with NumPy;
  * writes a prepared-task manifest consumed by 41_susie_coloc.R.

The R worker never invokes tabix, plink2, zcat, head, or any other command-line
tool.  This is important because R runs inside Singularity whereas this script
runs in the host LSF environment where the plink/samtools modules are loaded.

PLINK 2.00a3.6LM (the cluster module) predates bulk --r/--r-unphased, but it
does support dosage-aware ``--export Av``.  We therefore export variant-major
ALT dosages with plink2 and calculate the signed Pearson correlation matrix in
Python/NumPy.  This keeps the workflow compatible with the installed PLINK2
without requiring PLINK 1.9.
"""

import argparse
import csv
import gzip
import io
import math
import re
import subprocess
from pathlib import Path

import numpy as np


PREP_FIELDS = [
    "prep_status", "prep_message", "merged_file",
    "ld_x_matrix", "ld_x_vars", "ld_g_matrix", "ld_g_vars", "n_xqtl",
]


def sanitize_id(value):
    return re.sub(r"[^A-Za-z0-9._+-]", "_", str(value))


def open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


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
        [tabix, str(path), region], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True, check=False
    )
    if proc.returncode != 0:
        msg = proc.stderr.strip() or f"tabix exited {proc.returncode}"
        raise RuntimeError(msg)
    if not proc.stdout.strip():
        return []
    return list(csv.DictReader(io.StringIO(proc.stdout), fieldnames=header,
                               delimiter="\t"))


def fnum(value):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return math.nan
    return x


def finite(value):
    return math.isfinite(fnum(value))


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


def merge_stats(x_rows, g_rows, phenotype_id):
    # Match 41_susie_coloc.R semantics: phenotype-specific xQTL rows, merge on
    # canonical var_id, keep the smallest xQTL nominal p for duplicate IDs,
    # and require finite effects/SE/allele frequencies.
    x_best = {}
    for r in x_rows:
        if r.get("phenotype_id") != phenotype_id:
            continue
        vid = r.get("variant_id", "")
        if not vid:
            continue
        p = fnum(r.get("pval_nominal"))
        old = x_best.get(vid)
        if old is None or p < fnum(old.get("pval_nominal")):
            x_best[vid] = r

    g_best = {}
    for r in g_rows:
        vid = r.get("var_id", "")
        if vid and vid not in g_best:
            g_best[vid] = r

    merged = []
    for vid in x_best.keys() & g_best.keys():
        x = x_best[vid]
        g = g_best[vid]
        slope = fnum(x.get("slope"))
        slope_se = fnum(x.get("slope_se"))
        af = fnum(x.get("af"))
        beta = fnum(g.get("beta"))
        se = fnum(g.get("se"))
        eaf = fnum(g.get("eaf"))
        if not (math.isfinite(slope) and math.isfinite(slope_se) and slope_se > 0
                and math.isfinite(beta) and math.isfinite(se) and se > 0
                and math.isfinite(af) and 0 < af < 1
                and math.isfinite(eaf) and 0 < eaf < 1):
            continue
        merged.append({
            "var_id": vid,
            "pos_x": int(float(x["pos"])),
            "slope": slope,
            "slope_se": slope_se,
            "af": af,
            "pval_nominal": fnum(x.get("pval_nominal")),
            "pos_g": int(float(g["pos"])),
            "beta": beta,
            "se": se,
            "eaf": eaf,
            "pval": fnum(g.get("pval")),
            "n": fnum(g.get("n")),
        })
    merged.sort(key=lambda r: r["pos_x"])
    return merged


def write_merged(rows, path):
    fields = ["var_id", "pos_x", "slope", "slope_se", "af", "pval_nominal",
              "pos_g", "beta", "se", "eaf", "pval", "n"]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def _alt_from_var_id(var_id):
    # Canonical harmonized IDs are CHROM:POS:REF:ALT.  LD signs must refer to
    # the same ALT allele used by tensorQTL/GWAS effect estimates.
    parts = str(var_id).split(":", 3)
    if len(parts) != 4 or not parts[3]:
        raise ValueError(f"cannot recover ALT allele from canonical var_id: {var_id}")
    return parts[3]


def run_ld(plink2, pgen, keep, var_ids, prefix, plink_memory_mb):
    """Export ALT-coded dosages with PLINK2 and compute signed LD in NumPy.

    The installed PLINK v2.00a3.6LM predates bulk --r/--r-unphased, while
    --export Av (variant-major additive dosage) is supported.  --export-allele
    explicitly makes the canonical ALT allele the counted allele, matching the
    effect direction in the harmonized xQTL/GWAS summary statistics.
    """
    extract = Path(str(prefix) + ".extract.txt")
    extract.write_text("".join(f"{v}\n" for v in var_ids))

    alt_path = Path(str(prefix) + ".alt.txt")
    alt_by_id = {}
    with alt_path.open("w") as fh:
        for vid in var_ids:
            alt = _alt_from_var_id(vid)
            alt_by_id[vid] = alt
            fh.write(f"{vid} {alt}\n")

    cmd = [plink2, "--pfile", pgen, "--extract", str(extract),
           "--max-alleles", "2", "--export", "Av",
           "--export-allele", str(alt_path), "--threads", "1",
           "--memory", str(plink_memory_mb),
           "--silent", "--out", str(prefix)]
    if keep and Path(keep).exists():
        cmd[3:3] = ["--keep", keep]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, check=False)
    traw = Path(str(prefix) + ".traw")
    if proc.returncode != 0 or not traw.exists():
        msg = proc.stderr.strip() or proc.stdout.strip()
        if not msg:
            msg = f"plink2 exited {proc.returncode}; expected {traw.name}"
        return None, msg[-1000:]

    ids = []
    dosage_rows = []
    try:
        with traw.open(newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader, None)
            if header is None or len(header) < 7:
                return None, f"malformed PLINK2 dosage export: {traw.name}"
            nsamples = len(header) - 6
            for fields in reader:
                if not fields or len(fields) != len(header):
                    continue
                vid = fields[1]
                counted = fields[4]
                expected_alt = alt_by_id.get(vid)
                if expected_alt is None:
                    continue
                if counted != expected_alt:
                    return None, (f"PLINK2 counted allele mismatch for {vid}: "
                                  f"expected ALT={expected_alt}, got {counted}")
                vals = []
                for value in fields[6:]:
                    if value in ("NA", "nan", ".", ""):
                        vals.append(np.nan)
                    else:
                        vals.append(float(value))
                if len(vals) != nsamples:
                    return None, f"wrong dosage column count for {vid}"
                ids.append(vid)
                dosage_rows.append(vals)
    except Exception as exc:
        return None, f"could not parse {traw.name}: {exc}"

    if len(ids) == 0:
        return None, f"no variants exported to {traw.name}"
    if len(ids) != len(set(ids)):
        return None, f"duplicate variant IDs in {traw.name}"
    if nsamples < 2:
        return None, f"fewer than 2 samples in {traw.name}"

    G = np.asarray(dosage_rows, dtype=float)  # variants x samples
    with np.errstate(invalid="ignore"):
        means = np.nanmean(G, axis=1)
    centered = G - means[:, None]
    # Mean-impute missing dosages after centering (missing -> centered value 0).
    centered[~np.isfinite(centered)] = 0.0
    ss = np.sqrt(np.sum(centered * centered, axis=1))
    denom = np.outer(ss, ss)
    with np.errstate(divide="ignore", invalid="ignore"):
        ld = centered.dot(centered.T) / denom
    valid = np.isfinite(ss) & (ss > 0)
    diag_idx = np.arange(len(ids))
    ld[diag_idx[valid], diag_idx[valid]] = 1.0
    ld[~np.isfinite(ld)] = np.nan
    finite = np.isfinite(ld)
    ld[finite] = np.clip(ld[finite], -1.0, 1.0)

    matrix = Path(str(prefix) + ".ld")
    vars_path = Path(str(prefix) + ".ld.vars")
    np.savetxt(str(matrix), ld, delimiter="\t", fmt="%.10g")
    vars_path.write_text("".join(f"{v}\n" for v in ids))
    return (matrix, vars_path), ""


def prepare_task(task, index, workdir, tabix, plink2, min_variants,
                 plink_memory_mb):
    out = dict(task)
    out.update({k: "" for k in PREP_FIELDS})
    task_dir = workdir / f"task-{index:04d}"
    task_dir.mkdir(parents=True, exist_ok=True)

    def fail(status, message=""):
        out["prep_status"] = status
        out["prep_message"] = str(message).replace("\t", " ").replace("\n", " ")[:500]
        return out

    try:
        x_rows = tabix_rows(tabix, task["xqtl_file"], task["chrom"],
                            task["start"], task["end"])
        if not x_rows:
            return fail("no_xqtl_variants")
        if not any(r.get("phenotype_id") == task["phenotype_id"] for r in x_rows):
            return fail("no_xqtl_variants", "phenotype absent from slice")

        g_rows = tabix_rows(tabix, task["gwas_file"], task["chrom"],
                            task["start"], task["end"])
        if not g_rows:
            return fail("no_gwas_variants")

        merged = merge_stats(x_rows, g_rows, task["phenotype_id"])
        if not merged:
            return fail("no_shared_variants")
        if len(merged) < min_variants:
            return fail("too_few_variants", "after merge")

        merged_path = task_dir / "merged.tsv"
        write_merged(merged, merged_path)
        var_ids = [r["var_id"] for r in merged]

        n_xqtl = psam_n(task["ld_xqtl_pgen"])
        if n_xqtl is None or n_xqtl < 10:
            return fail("error", "could not read xQTL N from psam")

        ldx, msg = run_ld(plink2, task["ld_xqtl_pgen"], "", var_ids,
                          task_dir / "ld_xqtl", plink_memory_mb)
        if ldx is None:
            return fail("ld_failed_xqtl", msg)
        ldg, msg = run_ld(plink2, task["ld_gwas_pgen"],
                          task.get("ld_gwas_keep", ""), var_ids,
                          task_dir / "ld_gwas", plink_memory_mb)
        if ldg is None:
            return fail("ld_failed_gwas", msg)

        out.update({
            "prep_status": "ok", "prep_message": "",
            "merged_file": str(merged_path),
            "ld_x_matrix": str(ldx[0]), "ld_x_vars": str(ldx[1]),
            "ld_g_matrix": str(ldg[0]), "ld_g_vars": str(ldg[1]),
            "n_xqtl": str(n_xqtl),
        })
        return out
    except Exception as exc:
        return fail("error", exc)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tasks", required=True)
    p.add_argument("--shard-index", required=True, type=int)
    p.add_argument("--n-shards", required=True, type=int)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--tabix", default="tabix")
    p.add_argument("--plink2", default="plink2")
    p.add_argument("--plink-memory-mb", type=int, default=2048,
                   help="Memory cap passed to plink2 --memory (MiB; default 2048)")
    p.add_argument("--min-variants", type=int, default=50)
    args = p.parse_args()
    if args.plink_memory_mb < 256:
        p.error("--plink-memory-mb must be at least 256 MiB")

    with open(args.tasks, newline="") as fh:
        tasks = list(csv.DictReader(fh, delimiter="\t"))
    selected = [t for i, t in enumerate(tasks)
                if (i % args.n_shards) + 1 == args.shard_index]
    if not tasks:
        raise SystemExit(f"ERROR: empty task list: {args.tasks}")

    workdir = Path(args.work_dir)
    workdir.mkdir(parents=True, exist_ok=True)
    prepared = [prepare_task(t, i + 1, workdir, args.tabix, args.plink2,
                             args.min_variants, args.plink_memory_mb)
                for i, t in enumerate(selected)]
    fields = list(tasks[0].keys()) + PREP_FIELDS
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(prepared)
    counts = {}
    for row in prepared:
        counts[row["prep_status"]] = counts.get(row["prep_status"], 0) + 1
    print(f"prepared {len(prepared)} tasks for shard {args.shard_index}/{args.n_shards}: {counts}")


if __name__ == "__main__":
    main()
