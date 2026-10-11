#!/usr/bin/env python3
"""Prepare deterministic chromosome x feature shards for FUSION association.

The merged FUSION .pos file contains all trained models for one weight set.
FUSION itself is single-core and chromosome-specific, so this helper splits
models within each chromosome into approximately equal-sized shards, each
containing no more than --features-per-shard models.  The resulting global
manifest is suitable for a single LSF array spanning all chromosomes.

Only Python's standard library is required.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path
import sys
import tempfile

REQUIRED = ("WGT", "ID", "CHR")
AUTOSOMES = tuple(range(1, 23))


def normalize_chr(value: str) -> int:
    text = str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    try:
        chrom = int(text)
    except ValueError as exc:
        raise ValueError(f"non-autosomal CHR value in .pos: {value!r}") from exc
    if chrom not in AUTOSOMES:
        raise ValueError(f"CHR must be 1-22 for Step 48; got {value!r}")
    return chrom


def balanced_sizes(n: int, target: int) -> list[int]:
    """Return near-equal shard sizes, each <= target, summing to n."""
    if n <= 0:
        return []
    nshards = math.ceil(n / target)
    q, r = divmod(n, nshards)
    return [q + (1 if i < r else 0) for i in range(nshards)]


def plan_id(pos_path: Path, features_per_shard: int) -> str:
    h = hashlib.sha256()
    h.update(pos_path.read_bytes())
    h.update(b"\0features-per-shard=")
    h.update(str(features_per_shard).encode("ascii"))
    return h.hexdigest()[:12]


def write_tsv_atomic(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", newline="", encoding="utf-8", dir=path.parent, delete=False
    ) as tmp:
        writer = csv.DictWriter(tmp, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        tmp_path = Path(tmp.name)
    tmp_path.replace(path)


def prepare(pos_path: Path, out_root: Path, features_per_shard: int) -> Path:
    if features_per_shard < 1:
        raise ValueError("--features-per-shard must be >= 1")
    if not pos_path.is_file():
        raise FileNotFoundError(f"merged FUSION .pos not found: {pos_path}")

    with pos_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"empty .pos file: {pos_path}")
        missing = [c for c in REQUIRED if c not in reader.fieldnames]
        if missing:
            raise ValueError(f".pos missing required columns {missing}: {pos_path}")
        header = list(reader.fieldnames)
        by_chr: dict[int, list[dict[str, str]]] = {c: [] for c in AUTOSOMES}
        for row in reader:
            chrom = normalize_chr(row["CHR"])
            row["CHR"] = str(chrom)
            by_chr[chrom].append(row)

    total = sum(map(len, by_chr.values()))
    if total == 0:
        raise ValueError(f"no models in merged FUSION .pos: {pos_path}")

    pid = plan_id(pos_path, features_per_shard)
    plan_dir = out_root / f"features_{features_per_shard}_{pid}"
    plan_dir.mkdir(parents=True, exist_ok=True)

    # A plan directory is content-addressed by the .pos bytes + target shard
    # size. Rebuilding the same plan is safe; remove only its generated shard
    # .pos files/manifest so interrupted preparation cannot leave stale rows.
    for old in plan_dir.glob("chr*.shard*.pos"):
        old.unlink()
    manifest = plan_dir / "manifest.tsv"
    if manifest.exists():
        manifest.unlink()

    manifest_rows: list[dict[str, str]] = []
    array_index = 0
    for chrom in AUTOSOMES:
        rows = by_chr[chrom]
        sizes = balanced_sizes(len(rows), features_per_shard)
        offset = 0
        for shard_index, size in enumerate(sizes, start=1):
            array_index += 1
            shard_rows = rows[offset : offset + size]
            offset += size
            shard_path = plan_dir / f"chr{chrom:02d}.shard{shard_index:04d}.pos"
            write_tsv_atomic(shard_path, header, shard_rows)
            manifest_rows.append(
                {
                    "array_index": str(array_index),
                    "chromosome": str(chrom),
                    "shard_index": str(shard_index),
                    "n_features": str(size),
                    "weights_file": str(shard_path.resolve()),
                }
            )

    manifest_header = [
        "array_index",
        "chromosome",
        "shard_index",
        "n_features",
        "weights_file",
    ]
    write_tsv_atomic(manifest, manifest_header, manifest_rows)

    # Stable convenience pointer for inspection only. Workers receive the
    # content-addressed manifest path explicitly, so concurrent/old plans do
    # not affect running jobs.
    latest = out_root / "latest_manifest.txt"
    out_root.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=out_root, delete=False
    ) as tmp:
        tmp.write(str(manifest.resolve()) + "\n")
        tmp_latest = Path(tmp.name)
    tmp_latest.replace(latest)

    counts = ", ".join(
        f"chr{c}={len(by_chr[c])}" for c in AUTOSOMES if by_chr[c]
    )
    print(
        f"prepared {len(manifest_rows)} FUSION shards from {total} models "
        f"(target <= {features_per_shard}/shard; {counts})",
        file=sys.stderr,
    )
    # stdout is intentionally machine-readable for command substitution.
    print(manifest.resolve())
    return manifest


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--pos", required=True, type=Path)
    p.add_argument("--out-root", required=True, type=Path)
    p.add_argument("--features-per-shard", type=int, default=200)
    args = p.parse_args()
    try:
        prepare(args.pos, args.out_root, args.features_per_shard)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
