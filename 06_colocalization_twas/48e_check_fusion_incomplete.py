#!/usr/bin/env python3
"""
48e_check_fusion_incomplete.py

Inspect adaptive FUSION Step-48 outputs and write a machine-readable recovery
manifest for a follow-up submitter.

For each weight-set x GWAS-trait pair that does NOT yet have a non-empty merged
TWAS output, this script:

  1. finds the current content-addressed feature-shard manifest for the weight set,
  2. checks every expected per-shard .twas.tsv and .done sentinel,
  3. records exactly which original LSF array indices are incomplete,
  4. writes one TSV row per incomplete pair.

Rows have action:
  resubmit_shards  -> one or more array indices are incomplete
  merge_only       -> all shards are complete, but merged output is absent/empty

The recovery TSV is designed to be consumed by a separate submission script.
It does NOT query LSF and does NOT submit jobs.

Expected current Step-48 layout:
  $ISOTWAS_DIR/weights/{WS}/fusion_shards/latest_manifest.txt
  $ISOTWAS_DIR/weights/{WS}/fusion_shards/features_{N}_{hash}/manifest.tsv
  $ISOTWAS_DIR/fusion/{WS}/per_shard/{trait}/{plan_id}/chr*.shard*.twas.tsv
  $ISOTWAS_DIR/fusion/{WS}/{WS}_{trait}.twas.tsv

Only traits with an existing harmonized
  $GWAS_DIR/{trait_id}.sumstats.tsv.gz
are considered, matching 48_fusion_twas.sh.

Python standard library only.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Shard:
    array_index: int
    chromosome: int
    shard_index: int
    n_features: int
    weights_file: Path


def nonempty(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def compress_indices(indices: Iterable[int]) -> str:
    """Compress [1,2,3,7,9,10] -> '1-3,7,9-10' (valid LSF array-index syntax)."""
    vals = sorted(set(int(x) for x in indices))
    if not vals:
        return ""
    pieces: list[str] = []
    start = prev = vals[0]
    for x in vals[1:]:
        if x == prev + 1:
            prev = x
            continue
        pieces.append(str(start) if start == prev else f"{start}-{prev}")
        start = prev = x
    pieces.append(str(start) if start == prev else f"{start}-{prev}")
    return ",".join(pieces)


def parse_weight_sets(text: str) -> list[str]:
    vals = [x for x in text.replace(",", " ").split() if x]
    if not vals:
        raise ValueError("no weight sets supplied")
    return vals


def choose_manifest(
    isotwas_dir: Path,
    weight_set: str,
    features_per_shard: int | None,
) -> Path:
    root = isotwas_dir / "weights" / weight_set / "fusion_shards"
    latest = root / "latest_manifest.txt"

    if not latest.is_file():
        raise FileNotFoundError(f"{weight_set}: missing {latest}")

    target = latest.read_text(encoding="utf-8").strip()
    if not target:
        raise ValueError(f"{weight_set}: empty {latest}")

    manifest = Path(target)
    if not manifest.is_absolute():
        manifest = (latest.parent / manifest).resolve()

    if not nonempty(manifest):
        raise FileNotFoundError(f"{weight_set}: manifest missing/empty: {manifest}")

    if features_per_shard is not None:
        expected_prefix = f"features_{features_per_shard}_"
        plan_id = manifest.parent.name
        if not plan_id.startswith(expected_prefix):
            # Prefer an unambiguous matching plan if latest points to another shard size.
            candidates = sorted(
                (p / "manifest.tsv")
                for p in root.glob(f"{expected_prefix}*")
                if (p / "manifest.tsv").is_file()
            )
            if len(candidates) == 1:
                manifest = candidates[0].resolve()
            elif len(candidates) == 0:
                raise FileNotFoundError(
                    f"{weight_set}: latest plan is {plan_id!r}, and no "
                    f"{expected_prefix}*/manifest.tsv exists under {root}"
                )
            else:
                raise RuntimeError(
                    f"{weight_set}: latest plan is {plan_id!r}, but multiple "
                    f"{expected_prefix}* plans exist; cannot safely infer which one "
                    f"was used. Point latest_manifest.txt at the intended plan or "
                    f"omit --features-per-shard."
                )
    return manifest.resolve()


def read_shards(manifest: Path) -> list[Shard]:
    shards: list[Shard] = []
    with manifest.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {
            "array_index",
            "chromosome",
            "shard_index",
            "n_features",
            "weights_file",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{manifest}: missing columns {sorted(missing)}")

        for row in reader:
            shards.append(
                Shard(
                    array_index=int(row["array_index"]),
                    chromosome=int(row["chromosome"]),
                    shard_index=int(row["shard_index"]),
                    n_features=int(row["n_features"]),
                    weights_file=Path(row["weights_file"]),
                )
            )

    if not shards:
        raise ValueError(f"{manifest}: contains no shard rows")

    indices = [s.array_index for s in shards]
    if indices != list(range(1, len(shards) + 1)):
        raise ValueError(
            f"{manifest}: array_index must be contiguous 1..N; "
            f"observed first/last={indices[:3]}...{indices[-3:]}"
        )
    return shards


def read_catalog(catalog: Path, gwas_dir: Path) -> list[tuple[str, str]]:
    traits: list[tuple[str, str]] = []
    with catalog.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"trait_id", "primary_ancestry"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{catalog}: missing columns {sorted(missing)}")

        for row in reader:
            trait_id = (row.get("trait_id") or "").strip()
            primary = (row.get("primary_ancestry") or "").strip()
            if not trait_id:
                continue
            # Match 48_fusion_twas.sh: only harmonized traits that exist are eligible.
            if (gwas_dir / f"{trait_id}.sumstats.tsv.gz").is_file():
                traits.append((trait_id, primary))
    return traits


def ld_ancestry(primary: str, weight_set: str) -> str:
    # Match 48_fusion_twas.sh exactly.
    if primary in {"EAS", "EUR"}:
        return primary
    if weight_set != "pooled":
        return weight_set
    return "EUR"


def shard_stem(chromosome: int, shard_index: int) -> str:
    return f"chr{chromosome:02d}.shard{shard_index:04d}"


def main() -> int:
    p = argparse.ArgumentParser(
        description="Write a recovery TSV for missing/incomplete adaptive FUSION Step-48 pairs."
    )
    p.add_argument(
        "--isotwas-dir",
        type=Path,
        default=Path(os.environ.get("ISOTWAS_DIR", "")) if os.environ.get("ISOTWAS_DIR") else None,
        help="Step-48 isotwas root (default: $ISOTWAS_DIR)",
    )
    p.add_argument(
        "--gwas-dir",
        type=Path,
        default=Path(os.environ.get("GWAS_DIR", "")) if os.environ.get("GWAS_DIR") else None,
        help="GWAS root containing *.sumstats.tsv.gz (default: $GWAS_DIR)",
    )
    p.add_argument(
        "--catalog",
        type=Path,
        default=None,
        help="gwas_catalog.tsv (default: $SCRIPTS_DIR/gwas_catalog.tsv or ./gwas_catalog.tsv)",
    )
    p.add_argument(
        "--weight-sets",
        default=os.environ.get("WEIGHT_SETS", "EAS EUR pooled"),
        help='Space/comma-separated weight sets (default: "EAS EUR pooled")',
    )
    p.add_argument(
        "--features-per-shard",
        type=int,
        default=None,
        help=(
            "Optional guard selecting features_N_* plans. Recommended for recovery "
            "after a known production run, e.g. --features-per-shard 100."
        ),
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("fusion_recovery.tsv"),
        help="Output recovery TSV (default: fusion_recovery.tsv)",
    )
    p.add_argument(
        "--include-complete",
        action="store_true",
        help="Also write complete pairs with action=complete (normally omitted).",
    )
    args = p.parse_args()

    if args.isotwas_dir is None:
        p.error("--isotwas-dir is required unless ISOTWAS_DIR is set")
    if args.gwas_dir is None:
        p.error("--gwas-dir is required unless GWAS_DIR is set")
    if args.features_per_shard is not None and args.features_per_shard < 1:
        p.error("--features-per-shard must be >= 1")

    isotwas_dir = args.isotwas_dir.resolve()
    gwas_dir = args.gwas_dir.resolve()

    if args.catalog is not None:
        catalog = args.catalog.resolve()
    else:
        scripts_dir = os.environ.get("SCRIPTS_DIR")
        catalog = (
            Path(scripts_dir).resolve() / "gwas_catalog.tsv"
            if scripts_dir
            else Path(__file__).resolve().parent / "gwas_catalog.tsv"
        )
    if not catalog.is_file():
        p.error(f"catalog not found: {catalog}")

    try:
        weight_sets = parse_weight_sets(args.weight_sets)
        traits = read_catalog(catalog, gwas_dir)
        if not traits:
            raise RuntimeError(
                f"no harmonized catalog traits found under {gwas_dir} "
                f"(expected <trait>.sumstats.tsv.gz)"
            )

        plans: dict[str, tuple[Path, list[Shard]]] = {}
        for ws in weight_sets:
            manifest = choose_manifest(
                isotwas_dir, ws, args.features_per_shard
            )
            plans[ws] = (manifest, read_shards(manifest))

        fields = [
            "action",
            "weight_set",
            "trait_id",
            "primary_ancestry",
            "ld_ancestry",
            "plan_id",
            "manifest",
            "job_name",
            "merge_job_name",
            "n_total_shards",
            "n_complete_shards",
            "n_missing_shards",
            "missing_indices",
            "missing_indices_count",
            "missing_feature_count",
            "final_output",
            "per_shard_dir",
        ]

        rows: list[dict[str, str]] = []
        total_pairs = 0
        merged_complete = 0
        need_shards = 0
        need_merge_only = 0

        for ws in weight_sets:
            manifest, shards = plans[ws]
            plan_id = manifest.parent.name
            for trait_id, primary in traits:
                total_pairs += 1
                final = isotwas_dir / "fusion" / ws / f"{ws}_{trait_id}.twas.tsv"
                per_shard_dir = (
                    isotwas_dir / "fusion" / ws / "per_shard" / trait_id / plan_id
                )

                missing_indices: list[int] = []
                missing_feature_count = 0
                complete_count = 0

                for shard in shards:
                    stem = shard_stem(shard.chromosome, shard.shard_index)
                    out = per_shard_dir / f"{stem}.twas.tsv"
                    done = Path(str(out) + ".done")
                    if nonempty(out) and done.is_file():
                        complete_count += 1
                    else:
                        missing_indices.append(shard.array_index)
                        missing_feature_count += shard.n_features

                if nonempty(final):
                    merged_complete += 1
                    action = "complete"
                elif missing_indices:
                    need_shards += 1
                    action = "resubmit_shards"
                else:
                    need_merge_only += 1
                    action = "merge_only"

                if action == "complete" and not args.include_complete:
                    continue

                rows.append(
                    {
                        "action": action,
                        "weight_set": ws,
                        "trait_id": trait_id,
                        "primary_ancestry": primary,
                        "ld_ancestry": ld_ancestry(primary, ws),
                        "plan_id": plan_id,
                        "manifest": str(manifest),
                        "job_name": f"fusion_{ws}_{trait_id}",
                        "merge_job_name": f"fusionmerge_{ws}_{trait_id}",
                        "n_total_shards": str(len(shards)),
                        "n_complete_shards": str(complete_count),
                        "n_missing_shards": str(len(missing_indices)),
                        "missing_indices": compress_indices(missing_indices),
                        "missing_indices_count": str(len(missing_indices)),
                        "missing_feature_count": str(missing_feature_count),
                        "final_output": str(final),
                        "per_shard_dir": str(per_shard_dir),
                    }
                )

        args.out.parent.mkdir(parents=True, exist_ok=True)
        tmp = args.out.with_name(args.out.name + f".tmp.{os.getpid()}")
        with tmp.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        tmp.replace(args.out)

        print(f"Wrote recovery manifest: {args.out.resolve()}")
        print(f"Eligible weight-set x trait pairs: {total_pairs}")
        print(f"Merged outputs complete:           {merged_complete}")
        print(f"Need shard resubmission:           {need_shards}")
        print(f"Need merge only:                   {need_merge_only}")
        print(f"Recovery rows written:             {len(rows)}")
        return 0

    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
