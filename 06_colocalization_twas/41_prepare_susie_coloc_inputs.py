#!/usr/bin/env python3
"""Prepare one SuSiE-coloc shard for the pure-R model worker.

All external command-line I/O belongs here, outside R.  For each task assigned
to this shard the script:
  * extracts the xQTL and GWAS locus with tabix;
  * merges/filters shared variants;
  * subsets each PGEN with plink2, then computes signed LD with PLINK 1.9;
  * writes a prepared-task manifest consumed by 41_susie_coloc.R.

The R worker never invokes tabix, plink2, zcat, head, or any other command-line
tool.  This is important because R runs inside Singularity whereas this script
runs in the host LSF environment where the plink/samtools modules are loaded.

PLINK 2.00a3.6LM (the cluster module) predates bulk --r/--r-unphased.
For compatibility, PGEN filtering/conversion is done with plink2 and the actual
signed correlation matrix is computed by PLINK 1.9 ``--r square``.
"""

import argparse
import csv
import gzip
import io
import math
import re
import subprocess
from pathlib import Path


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


def run_ld(plink2, plink1, pgen, keep, var_ids, prefix):
    extract = Path(str(prefix) + ".extract.txt")
    extract.write_text("".join(f"{v}\n" for v in var_ids))

    # The seadragon `module load plink` currently provides PLINK2
    # v2.00a3.6LM (14 Aug 2022).  That build does not implement bulk --r;
    # convert the filtered PGEN to a small PLINK1 binary fileset first.
    bed_prefix = Path(str(prefix) + ".bedtmp")
    cmd = [plink2, "--pfile", pgen, "--extract", str(extract),
           "--max-alleles", "2", "--make-bed", "--threads", "1",
           "--silent", "--out", str(bed_prefix)]
    if keep and Path(keep).exists():
        cmd[3:3] = ["--keep", keep]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, check=False)
    bim = Path(str(bed_prefix) + ".bim")
    bed = Path(str(bed_prefix) + ".bed")
    fam = Path(str(bed_prefix) + ".fam")
    if proc.returncode != 0 or not (bed.exists() and bim.exists() and fam.exists()):
        msg = proc.stderr.strip() or proc.stdout.strip()
        if not msg:
            msg = (f"plink2 exited {proc.returncode}; expected "
                   f"{bed.name}/{bim.name}/{fam.name}")
        return None, msg[-1000:]

    # Explicitly force A1 to the canonical ALT allele.  PLINK1's signed
    # --r matrix is then on the same allele scale as the summary-stat effects.
    a1_path = Path(str(prefix) + ".a1.txt")
    with a1_path.open("w") as fh:
        for vid in var_ids:
            fh.write(f"{vid} {_alt_from_var_id(vid)}\n")

    cmd = [plink1, "--bfile", str(bed_prefix), "--extract", str(extract),
           "--a1-allele", str(a1_path), "2", "1",
           "--keep-allele-order", "--r", "square", "--threads", "1",
           "--silent", "--out", str(prefix)]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, check=False)
    matrix = Path(str(prefix) + ".ld")
    vars_path = Path(str(prefix) + ".ld.vars")
    if proc.returncode != 0 or not matrix.exists():
        msg = proc.stderr.strip() or proc.stdout.strip()
        if not msg:
            msg = f"plink 1.9 exited {proc.returncode}; expected {matrix.name}"
        return None, msg[-1000:]

    # PLINK1 square matrices follow .bim variant order but do not emit an ID
    # sidecar.  Capture that exact order for the pure-R reader.
    bim_ids = []
    with bim.open() as fh:
        for line in fh:
            fields = line.rstrip("\n").split()
            if len(fields) >= 2:
                bim_ids.append(fields[1])
    if not bim_ids:
        return None, f"no variants found in {bim.name}"
    vars_path.write_text("".join(f"{v}\n" for v in bim_ids))
    return (matrix, vars_path), ""


def prepare_task(task, index, workdir, tabix, plink2, plink1, min_variants):
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

        ldx, msg = run_ld(plink2, plink1, task["ld_xqtl_pgen"], "", var_ids,
                          task_dir / "ld_xqtl")
        if ldx is None:
            return fail("ld_failed_xqtl", msg)
        ldg, msg = run_ld(plink2, plink1, task["ld_gwas_pgen"],
                          task.get("ld_gwas_keep", ""), var_ids,
                          task_dir / "ld_gwas")
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
    p.add_argument("--plink1", default="plink",
                   help="PLINK 1.9 executable used for signed --r square LD")
    p.add_argument("--min-variants", type=int, default=50)
    args = p.parse_args()

    with open(args.tasks, newline="") as fh:
        tasks = list(csv.DictReader(fh, delimiter="\t"))
    selected = [t for i, t in enumerate(tasks)
                if (i % args.n_shards) + 1 == args.shard_index]
    if not tasks:
        raise SystemExit(f"ERROR: empty task list: {args.tasks}")

    workdir = Path(args.work_dir)
    workdir.mkdir(parents=True, exist_ok=True)
    prepared = [prepare_task(t, i + 1, workdir, args.tabix, args.plink2,
                             args.plink1, args.min_variants)
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
