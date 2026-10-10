#!/usr/bin/env python3
"""Prepare one isoTWAS/TWAS training shard outside R.

The host-side preparer selects the shard's genes, extracts their expression and
isoform-expression rows from bgzipped BEDs using Python's gzip module, and
runs plink2 to export cis dosages for each ancestry panel. It writes only
ordinary TSV/PLINK .raw files plus a manifest.

46_isotwas_train.R consumes those files and performs only statistical
work in R; it never invokes command-line programs.
"""

import argparse
import csv
import gzip
import re
import subprocess
from pathlib import Path


MANIFEST_FIELDS = [
    "gene", "weight_set", "chrom", "win_start", "win_end",
    "prep_status", "prep_message", "data_dir",
]


def open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def sanitize_id(value):
    return re.sub(r"[^A-Za-z0-9._+-]", "_", str(value))


def bed_meta(path):
    meta = {}
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 4:
            raise RuntimeError(f"invalid/empty BED: {path}")
        for row in reader:
            if len(row) < 4:
                continue
            meta[row[3]] = (row[0], int(float(row[1])), int(float(row[2])))
    return header, meta


def selected_rows(path, genes, isoforms=False):
    genes = set(genes)
    rows = {g: [] for g in genes}
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise RuntimeError(f"empty BED: {path}")
        for row in reader:
            if len(row) < 4:
                continue
            pid = row[3]
            gene = pid.split("__", 1)[0] if isoforms else pid
            if gene in genes and (not isoforms or pid.startswith(gene + "__")):
                rows[gene].append(row)
    return header, rows


def write_rows(path, header, rows):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def run_export_region(plink2, pgen, chrom, start, end, prefix):
    cmd = [plink2, "--pfile", pgen,
           "--chr", str(chrom), "--from-bp", str(max(1, int(start))),
           "--to-bp", str(int(end)), "--export", "A",
           "--threads", "1", "--silent", "--out", str(prefix)]
    proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    raw = Path(str(prefix) + ".raw")
    if proc.returncode != 0 or not raw.exists():
        msg = proc.stderr.strip() or proc.stdout.strip()
        return None, (msg or f"plink2 exited {proc.returncode}")[-1000:]
    return raw, ""


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--weight-set", required=True)
    p.add_argument("--shard-index", required=True, type=int)
    p.add_argument("--n-shards", required=True, type=int)
    p.add_argument("--qtl-dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--cis-window", type=int, default=1_000_000)
    p.add_argument("--gene-list", default="")
    p.add_argument("--plink2", default="plink2")
    p.add_argument("--force", action="store_true",
                   help="rebuild selected genes even when weight files exist")
    args = p.parse_args()

    ws = args.weight_set
    ancestries = ["EAS", "EUR"] if ws == "pooled" else [ws]
    allowed = {"EAS", "EUR", "AFR", "AMR", "SAS"}
    if not set(ancestries) <= allowed:
        raise SystemExit(f"ERROR: unknown weight set: {ws}")

    qtl_dir = Path(args.qtl_dir)
    outdir = Path(args.outdir)
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    expr_headers = {}
    expr_meta = {}
    for anc in ancestries:
        path = qtl_dir / f"{anc}_expression.bed.gz"
        expr_headers[anc], expr_meta[anc] = bed_meta(path)

    gene_sets = [set(expr_meta[a]) for a in ancestries]
    genes = sorted(set.intersection(*gene_sets)) if gene_sets else []
    if args.gene_list:
        requested = {x.strip() for x in Path(args.gene_list).read_text().splitlines()
                     if x.strip()}
        genes = [g for g in genes if g in requested]
    genes = [g for i, g in enumerate(genes)
             if (i % args.n_shards) + 1 == args.shard_index]
    print(f"[prep shard {args.shard_index}/{args.n_shards}] {len(genes)} genes, "
          f"weight set {ws}", flush=True)

    # Capture only this shard's phenotype rows, keeping decompression in Python.
    expr_rows = {}
    iso_rows = {}
    iso_headers = {}
    for anc in ancestries:
        _, expr_rows[anc] = selected_rows(
            qtl_dir / f"{anc}_expression.bed.gz", genes, isoforms=False)
        iso_headers[anc], iso_rows[anc] = selected_rows(
            qtl_dir / f"{anc}_isoform_expression.bed.gz", genes, isoforms=True)

    manifest = []
    for i, gene in enumerate(genes, 1):
        chrom0, start0, end0 = expr_meta[ancestries[0]][gene]
        chrom = re.sub(r"^chr", "", str(chrom0), flags=re.IGNORECASE)
        win_start = max(0, start0 - args.cis_window)
        win_end = end0 + args.cis_window
        rec = {"gene": gene, "weight_set": ws, "chrom": chrom,
               "win_start": win_start, "win_end": win_end,
               "prep_status": "error", "prep_message": "", "data_dir": ""}

        def fail(status, message=""):
            rec["prep_status"] = status
            rec["prep_message"] = str(message).replace("\t", " ").replace("\n", " ")[:500]
            return rec

        # Match the existing worker's skip-if-done rule (gene-expression model).
        genes_dir = outdir / "weights" / ws / "genes"
        rdat = genes_dir / f"{gene}.wgt.RDat"
        if args.force:
            for old in [rdat, *genes_dir.glob(f"{gene}__*.wgt.RDat")]:
                if old.exists():
                    old.unlink()
        elif rdat.exists():
            rec["prep_status"] = "already_done"
            manifest.append(rec)
            continue

        gene_dir = work_dir / sanitize_id(gene)
        gene_dir.mkdir(parents=True, exist_ok=True)
        rec["data_dir"] = str(gene_dir)
        print(f"[prep {i}/{len(genes)}] {gene} {chrom}:{win_start}-{win_end}",
              flush=True)

        bad = None
        for anc in ancestries:
            erows = expr_rows[anc].get(gene, [])
            if not erows:
                bad = fail("error", f"{anc}: gene missing from expression BED")
                break
            write_rows(gene_dir / f"{anc}.expression.tsv", expr_headers[anc], erows[:1])
            write_rows(gene_dir / f"{anc}.isoform_expression.tsv",
                       iso_headers[anc], iso_rows[anc].get(gene, []))
            raw, msg = run_export_region(
                args.plink2, str(qtl_dir / f"{anc}_qtl"), chrom,
                win_start, win_end, gene_dir / anc)
            if raw is None:
                bad = fail("dosage_failed", f"{anc}: {msg}")
                break
        if bad is not None:
            manifest.append(bad)
            continue

        rec["prep_status"] = "ok"
        manifest.append(rec)

    manifest_path = Path(args.manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(manifest)
    print(f"prepared {len(manifest)} genes -> {manifest_path}")


if __name__ == "__main__":
    main()
