import csv
import importlib.util
from pathlib import Path


def _load_module():
    here = Path(__file__).resolve().parent
    path = here / "48d_prepare_fusion_shards.py"
    spec = importlib.util.spec_from_file_location("fusion_shards", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def _write_pos(path: Path, chr_counts):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["WGT", "ID", "GENE", "CHR", "P0", "P1"])
        k = 0
        for chrom, n in chr_counts.items():
            for _ in range(n):
                k += 1
                w.writerow([f"m{k}.wgt.RDat", f"m{k}", f"g{k}", chrom, 1, 2])


def test_balanced_sizes():
    m = _load_module()
    assert m.balanced_sizes(0, 200) == []
    assert m.balanced_sizes(199, 200) == [199]
    assert m.balanced_sizes(200, 200) == [200]
    assert m.balanced_sizes(201, 200) == [101, 100]
    sizes = m.balanced_sizes(1001, 200)
    assert sum(sizes) == 1001
    assert max(sizes) <= 200
    assert max(sizes) - min(sizes) <= 1


def test_prepare_adaptive_manifest(tmp_path):
    m = _load_module()
    pos = tmp_path / "EAS.pos"
    _write_pos(pos, {1: 401, 2: 50, 22: 200})

    manifest = m.prepare(pos, tmp_path / "shards", 200)
    with manifest.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    # chr1 -> 3 balanced shards, chr2 -> 1, chr22 -> 1.
    assert len(rows) == 5
    assert [int(r["array_index"]) for r in rows] == [1, 2, 3, 4, 5]
    chr1 = [r for r in rows if r["chromosome"] == "1"]
    assert [int(r["n_features"]) for r in chr1] == [134, 134, 133]
    assert [int(r["n_features"]) for r in rows if r["chromosome"] == "2"] == [50]
    assert [int(r["n_features"]) for r in rows if r["chromosome"] == "22"] == [200]

    for row in rows:
        shard = Path(row["weights_file"])
        assert shard.is_file()
        with shard.open() as fh:
            assert sum(1 for _ in fh) - 1 == int(row["n_features"])

    # Same input + target size -> same content-addressed plan directory.
    manifest2 = m.prepare(pos, tmp_path / "shards", 200)
    assert manifest2 == manifest
