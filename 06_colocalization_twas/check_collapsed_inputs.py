#!/usr/bin/env python3
"""check_collapsed_inputs.py — verify module 06 operates on the collapsed
replicate sample set, matching the contract module 07 enforces.

Module 05's optional replicate-collapse workflow (collapse_replicates.py)
writes replicate_collapsed/reports/ancestry_map_collapsed.tsv: one retained
RNA run per individual, with cross-protocol individuals removed entirely.
Module 07 treats that map as the retained-run authority and subsets every
input to it. Module 06 reads the module-05 qtl_inputs ({ANC}_*.bed.gz,
{ANC}_covariates_*.tsv, {ANC}_qtl.psam) directly and previously applied no
such check, so it could silently analyze samples module 07 excludes (e.g.
cross-protocol individuals, or unaveraged duplicate-run columns if the
collapse staging was not used upstream).

Per ancestry this script verifies:

  1. the collapsed ancestry map exists, has the required columns, and has no
     duplicate retained runs;
  2. exactly one retained run per array_id (joining the map to module-05
     {ANC}_metadata.tsv on rnaseq_id, as module 07 does);
  3. every sample column in the qtl_inputs BEDs / covariate files and every
     IID in {ANC}_qtl.psam is a retained array_id (extras are ERRORS);
  4. no duplicated sample columns (a signature of uncollapsed replicates);
  5. retained runs absent from module-05 metadata are WARNINGS only — module
     05's final RNA/DNA intersection legitimately drops individuals upstream,
     and module 07 inner-joins the map to the same metadata, so both modules
     resolve to the same sample set. A run found under a DIFFERENT ancestry's
     metadata is flagged separately (map/metadata ancestry disagreement);
  6. retained array_ids absent from qtl_inputs are reported as WARNINGS only
     (modality-specific coverage dropouts).

Stdlib only: runs under any python3, no conda environment required.

Usage:
  python3 check_collapsed_inputs.py --qtl-dir $QTL_DIR \
      --output-base $OUTPUT_BASE --ancestries EAS EUR
"""

import argparse
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

META_COLS = {"#chr", "chr", "start", "end", "phenotype_id", "id"}


def log(msg):
    print(msg, flush=True)


def read_header_samples(path):
    """Sample columns of a BED/covariates file (header line only)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    # drop the leading ID column, then any recognized metadata columns
    return [c for c in header[1:] if c not in META_COLS]


def read_psam_iids(path):
    """IID column of a plink2 .psam (handles #FID/IID and IID-only)."""
    iids = []
    with open(path) as fh:
        cols = None
        for line in fh:
            if line.startswith("#"):
                cols = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if cols is None:
                cols = ["IID"]  # no header: single IID column
            fields = line.rstrip("\n").split("\t")
            iids.append(fields[cols.index("IID")] if "IID" in cols else fields[0])
    return iids


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--qtl-dir", required=True,
                   help="module-05 qtl_inputs directory")
    p.add_argument("--output-base", required=True,
                   help="PANTRY output root (dirname of config.yml)")
    p.add_argument("--collapsed-ancestry-map", default=None,
                   help="default: {output-base}/replicate_collapsed/reports/"
                        "ancestry_map_collapsed.tsv")
    p.add_argument("--ancestries", nargs="+", default=["EAS", "EUR"])
    args = p.parse_args()

    qtl_dir = Path(args.qtl_dir)
    map_path = Path(args.collapsed_ancestry_map) if args.collapsed_ancestry_map else \
        Path(args.output_base) / "replicate_collapsed" / "reports" / \
        "ancestry_map_collapsed.tsv"

    n_err = 0

    def err(msg):
        nonlocal n_err
        n_err += 1
        log(f"  ERROR: {msg}")

    # --- collapsed ancestry map ----------------------------------------------
    if not map_path.exists():
        err(f"collapsed ancestry map not found: {map_path}\n"
            "      Module 06 must use the same collapsed-replicate sample set "
            "as module 07.\n      Run module-05 collapse_replicates.py first, "
            "or pass --collapsed-ancestry-map.\n      If this cohort genuinely "
            "has no technical replicates and the collapse step was "
            "intentionally skipped, rerun with SKIP_COLLAPSE_CHECK=1.")
        sys.exit(1)

    with open(map_path, newline="") as fh:
        amap = list(csv.DictReader(fh, delimiter="\t"))
    required = {"sample_id", "assigned_ancestry", "cohort"}
    missing = required - set(amap[0].keys() if amap else [])
    if missing:
        log(f"ERROR: collapsed ancestry map missing columns "
            f"{sorted(missing)}: {map_path}")
        sys.exit(1)
    seen = set()
    dups = set()
    for row in amap:
        sid = row["sample_id"]
        if sid in seen:
            dups.add(sid)
        seen.add(sid)
    if dups:
        err(f"collapsed ancestry map has duplicate sample_id rows: "
            f"{sorted(dups)[:5]}")
    log(f"collapsed ancestry map: {len(amap)} retained runs ({map_path})")

    # --- module-05 metadata for ALL ancestries (loaded up front so a retained
    # run absent from its own ancestry's metadata can be distinguished as a
    # cross-ancestry disagreement vs. a genuine intersection dropout) ---------
    meta_by_anc = {}
    rna2array_all = {}
    for anc in args.ancestries:
        meta_path = qtl_dir / f"{anc}_metadata.tsv"
        if not meta_path.exists():
            meta_by_anc[anc] = (meta_path, None)
            continue
        with open(meta_path, newline="") as fh:
            meta = list(csv.DictReader(fh, delimiter="\t"))
        if not {"rnaseq_id", "array_id"} <= set(meta[0].keys() if meta else []):
            meta_by_anc[anc] = (meta_path, None)
            continue
        meta_by_anc[anc] = (meta_path, meta)
        for m in meta:
            rna2array_all[m["rnaseq_id"]] = (anc, m["array_id"])

    # --- per-ancestry checks ---------------------------------------------------
    for anc in args.ancestries:
        log(f"[{anc}]")
        meta_path, meta = meta_by_anc[anc]
        if meta is None:
            err(f"module-05 metadata not found or missing "
                f"rnaseq_id/array_id columns: {meta_path} "
                "(needed to map retained rnaseq_id -> array_id)")
            continue
        rna2array = {m["rnaseq_id"]: m["array_id"] for m in meta}

        retained_runs = [r for r in amap if r["assigned_ancestry"] == anc]
        retained_rna = [r["sample_id"] for r in retained_runs]
        not_in_meta = [r for r in retained_rna if r not in rna2array]
        if not_in_meta:
            dropped = [r for r in not_in_meta if r not in rna2array_all]
            elsewhere = sorted(
                f"{r} (in {rna2array_all[r][0]} metadata)"
                for r in not_in_meta if r in rna2array_all)
            if dropped:
                log(f"  WARNING: {len(dropped)} retained run(s) absent from "
                    f"all module-05 metadata — dropped by the module-05 "
                    f"RNA/DNA intersection or upstream QC; module 07 "
                    f"inner-joins these away: {dropped[:5]}"
                    f"{' ...' if len(dropped) > 5 else ''}")
            if elsewhere:
                log(f"  WARNING: {len(elsewhere)} retained run(s) appear in "
                    f"a DIFFERENT ancestry's module-05 metadata (collapsed "
                    f"map vs. module-05 ancestry disagreement; excluded by "
                    f"module 07's inner join): {elsewhere[:5]}"
                    f"{' ...' if len(elsewhere) > 5 else ''}")
        retained_arrays = [rna2array[r] for r in retained_rna if r in rna2array]
        if len(set(retained_arrays)) != len(retained_arrays):
            dup_ids = [a for a, c in Counter(retained_arrays).items() if c > 1]
            err(f">1 retained run per array_id (violates the module-05 "
                f"replicate-collapse contract): {dup_ids[:5]}")
        retained = set(retained_arrays)
        log(f"  retained: {len(retained_runs)} runs in map -> "
            f"{len(retained)} array_ids with module-05 metadata "
            f"({len(not_in_meta)} dropped upstream)")

        # sample-bearing qtl_inputs for this ancestry
        targets = []
        targets += sorted(qtl_dir.glob(f"{anc}_*.bed.gz"))
        targets += sorted(qtl_dir.glob(f"{anc}_covariates_*.tsv"))
        psam = qtl_dir / f"{anc}_qtl.psam"
        if psam.exists():
            targets.append(psam)
        if not targets:
            err(f"no qtl_inputs files matched for {anc} in {qtl_dir}")
            continue

        for path in targets:
            if path.suffix == ".psam":
                samples = read_psam_iids(path)
            else:
                samples = read_header_samples(path)
            dup_cols = sorted({s for s in samples if samples.count(s) > 1})
            if dup_cols:
                err(f"{path.name}: {len(dup_cols)} duplicated sample "
                    f"column(s) (uncollapsed technical replicates?): "
                    f"{dup_cols[:5]}")
            extra = sorted(set(samples) - retained)
            if extra:
                err(f"{path.name}: {len(extra)} sample(s) NOT in the "
                    f"collapsed ancestry map (module 07 excludes these): "
                    f"{extra[:5]}"
                    f"{' ...' if len(extra) > 5 else ''}")
            missing_retained = sorted(retained - set(samples))
            if missing_retained:
                log(f"  WARNING: {path.name}: {len(missing_retained)} "
                    f"retained array_id(s) absent (dropped by the module-05 "
                    f"intersection?): {missing_retained[:5]}"
                    f"{' ...' if len(missing_retained) > 5 else ''}")
            if not dup_cols and not extra:
                log(f"  OK: {path.name}: {len(samples)} samples, all retained")

    if n_err:
        log(f"\nFAILED: {n_err} error(s). Module 06 qtl_inputs do not match "
            "the collapsed-replicate sample set module 07 uses.")
        log("Rerun module-05 steps 17-23 against the collapse staging tree "
            "(see module-05 README 'Technical replicates'), or override with "
            "SKIP_COLLAPSE_CHECK=1 if this is intentional.")
        sys.exit(1)
    log("\nOK: module 06 inputs match the collapsed-replicate contract "
        "(same retained-run authority as module 07).")


if __name__ == "__main__":
    main()
