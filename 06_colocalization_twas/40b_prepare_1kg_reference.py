#!/usr/bin/env python3
"""Prepare ancestry keep files for the one-time 1KG LD-reference build.

This helper is intentionally stdlib-only and independent of step 40 task-list
construction. It reads the local sample->superpopulation map plus the actual
PLINK2 .psam schema and writes keep files under the LD-reference directory.
The keep files are consumed only by 40c_run_1kg_ld_reference.sh.
"""

import argparse
import csv
import re
from pathlib import Path


def norm(name):
    return re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")


def data_lines(path):
    kept = []
    with Path(path).open(errors="replace") as fh:
        for raw in fh:
            if not raw.strip():
                continue
            stripped = raw.lstrip()
            if stripped.startswith("#"):
                candidate = stripped.lstrip("#").strip()
                n = norm(candidate)
                if ("sample" in n or "iid" in n) and ("superpop" in n or "superpopulation" in n):
                    kept.append(candidate)
                continue
            kept.append(raw.rstrip("\r\n"))
    if not kept:
        raise RuntimeError(f"empty sample map: {path}")
    return kept


def split_table(lines):
    header = lines[0]
    if "\t" in header:
        delim = "\t"
        rows = [line.split(delim) for line in lines]
    elif "," in header:
        rows = list(csv.reader(lines))
    else:
        rows = [re.split(r"\s+", line.strip()) for line in lines]
    return rows


def read_sample_map(path):
    rows = split_table(data_lines(path))
    header = [norm(x) for x in rows[0]]
    sid_aliases = ("sample_id", "sample", "iid", "sampleid", "sample_name")
    sp_aliases = ("superpop", "superpopulation", "superpopulation_code",
                  "super_pop", "super_pop_code")
    sid_i = next((header.index(x) for x in sid_aliases if x in header), None)
    sp_i = next((header.index(x) for x in sp_aliases if x in header), None)
    if sid_i is None or sp_i is None:
        raise RuntimeError(
            f"cannot identify sample/superpopulation columns in {path}; "
            f"normalized header={header}")
    out = {}
    need = max(sid_i, sp_i)
    for fields in rows[1:]:
        if len(fields) <= need:
            continue
        sid = fields[sid_i].strip()
        sp = fields[sp_i].strip().upper()
        if not sid or not sp:
            continue
        old = out.get(sid)
        if old is not None and old != sp:
            raise RuntimeError(f"conflicting superpopulation for {sid}: {old} vs {sp}")
        out[sid] = sp
    if not out:
        raise RuntimeError(f"no usable rows in sample map: {path}")
    return out


def read_psam(prefix):
    path = Path(str(prefix) + ".psam")
    if not path.is_file():
        raise RuntimeError(f"missing PSAM: {path}")
    header = None
    rows = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line.split()
                first = parts[0].lstrip("#").upper() if parts else ""
                if first in ("FID", "IID"):
                    header = [x.lstrip("#") for x in parts]
                continue
            rows.append(line.split())
    if header is None:
        fid_i, iid_i, has_fid = 0, 1, True
    else:
        upper = [x.upper() for x in header]
        if "IID" not in upper:
            raise RuntimeError(f"PSAM header lacks IID: {header}")
        iid_i = upper.index("IID")
        has_fid = "FID" in upper
        fid_i = upper.index("FID") if has_fid else None
    need = max(iid_i, fid_i if fid_i is not None else iid_i)
    ids = {}
    for fields in rows:
        if len(fields) <= need:
            continue
        iid = fields[iid_i]
        ids[iid] = (fields[fid_i], iid) if has_fid else (iid,)
    if not ids:
        raise RuntimeError(f"no samples parsed from {path}")
    return ids, has_fid


def write_keeps(sample_map, kg_pgen, ancestries, out_dir):
    mapping = read_sample_map(sample_map)
    psam, has_fid = read_psam(kg_pgen)
    out_dir = Path(out_dir)
    keep_dir = out_dir / "keep"
    keep_dir.mkdir(parents=True, exist_ok=True)
    print(f"parsed sample map: {len(mapping)} samples")
    print(f"parsed PSAM: {len(psam)} samples ({'FID+IID' if has_fid else 'IID-only'})")
    written = {}
    for anc in ancestries:
        anc = str(anc).upper()
        source_ids = [sid for sid, sp in mapping.items() if sp == anc]
        rows = [psam[sid] for sid in source_ids if sid in psam]
        missing = len(source_ids) - len(rows)
        if not source_ids:
            raise RuntimeError(f"sample map contains no {anc} samples")
        if not rows:
            raise RuntimeError(f"none of {len(source_ids)} {anc} sample-map IDs match the PSAM")
        path = keep_dir / f"{anc}.keep"
        with path.open("w") as fh:
            for row in rows:
                fh.write("\t".join(row) + "\n")
        print(f"{anc}: {len(rows)} matched; {missing} absent from PSAM -> {path}")
        written[anc] = path
    return written


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample-map", required=True)
    p.add_argument("--kg-pgen", required=True)
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()
    write_keeps(args.sample_map, args.kg_pgen, args.ancestries, args.out_dir)


if __name__ == "__main__":
    main()
