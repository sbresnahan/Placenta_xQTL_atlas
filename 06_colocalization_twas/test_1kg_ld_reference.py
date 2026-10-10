#!/usr/bin/env python3
"""Unit tests for the prebuilt 1KG ancestry/chromosome LD-reference contract."""

import importlib.util
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SUSIE = load("susie_prep", "41_prepare_susie_coloc_inputs.py")
CB = load("colocboost_prep", "44_prepare_colocboost_inputs.py")
REF = load("kg_ref_prep", "40b_prepare_1kg_reference.py")


def make_ref(tmp_path):
    prefix = tmp_path / "ref" / "EAS" / "chr1"
    prefix.parent.mkdir(parents=True)
    Path(str(prefix) + ".pgen").write_bytes(b"pgen")
    Path(str(prefix) + ".pvar").write_text("#CHROM\tPOS\tID\tREF\tALT\n")
    Path(str(prefix) + ".psam").write_text("#IID\nHG00096\n")
    Path(str(prefix) + ".done").touch()
    return prefix


def test_resolve_compact_reference(tmp_path):
    prefix = make_ref(tmp_path)
    assert SUSIE.resolve_kg_ref(tmp_path / "ref", "EAS", "chr1") == prefix
    assert CB.resolve_kg_ref(tmp_path / "ref", "eas", "1") == prefix


def test_cache_key_is_variant_set_specific_and_reference_sensitive(tmp_path):
    prefix = make_ref(tmp_path)
    a = ["1:10:A:C", "1:20:G:T"]
    b = list(reversed(a))
    m1, v1 = SUSIE.ld_cache_paths(tmp_path / "cache", "EAS", "1", prefix, a)
    m2, v2 = SUSIE.ld_cache_paths(tmp_path / "cache", "EAS", "1", prefix, b)
    assert (m1, v1) == (m2, v2), "same variant set should share one cache key"

    # Changing the reference identity must invalidate the cache key.
    pgen = Path(str(prefix) + ".pgen")
    pgen.write_bytes(b"pgen-new")
    os.utime(pgen, None)
    m3, _ = SUSIE.ld_cache_paths(tmp_path / "cache", "EAS", "1", prefix, a)
    assert m3 != m1


def test_missing_reference_fails_fast(tmp_path):
    try:
        SUSIE.resolve_kg_ref(tmp_path / "ref", "EUR", "1")
    except FileNotFoundError as exc:
        assert "40b_submit_1kg_ld_reference.sh" in str(exc)
    else:
        raise AssertionError("missing compact reference should fail")


def test_private_keep_builder_handles_iid_only_psam(tmp_path):
    sample_map = tmp_path / "map.tsv"
    sample_map.write_text(
        "sample_id\tsuperpop\n"
        "HG00096\tEUR\n"
        "HG00403\tEAS\n")
    prefix = tmp_path / "kg"
    Path(str(prefix) + ".psam").write_text(
        "#IID\tSEX\nHG00096\t0\nHG00403\t0\n")
    got = REF.write_keeps(sample_map, prefix, ["EAS", "EUR"], tmp_path / "ref")
    assert got["EAS"].read_text().strip() == "HG00403"
    assert got["EUR"].read_text().strip() == "HG00096"


def test_private_keep_builder_handles_fid_iid_psam(tmp_path):
    sample_map = tmp_path / "map.tsv"
    sample_map.write_text(
        "sample_id\tsuperpop\n"
        "HG00096\tEUR\n"
        "HG00403\tEAS\n")
    prefix = tmp_path / "kg"
    Path(str(prefix) + ".psam").write_text(
        "#FID\tIID\tSEX\nF1\tHG00096\t0\nF2\tHG00403\t0\n")
    got = REF.write_keeps(sample_map, prefix, ["EAS", "EUR"], tmp_path / "ref")
    assert got["EAS"].read_text().strip() == "F2\tHG00403"
    assert got["EUR"].read_text().strip() == "F1\tHG00096"


def test_private_keep_builder_accepts_3202_population_header(tmp_path):
    sample_map = tmp_path / "population.txt"
    sample_map.write_text(
        "FamilyID SampleID FatherID MotherID Sex Population Superpopulation\n"
        "0 HG00096 0 0 1 GBR EUR\n"
        "0 HG00403 0 0 2 CHS EAS\n")
    prefix = tmp_path / "kg"
    Path(str(prefix) + ".psam").write_text(
        "#IID\nHG00096\nHG00403\n")
    got = REF.write_keeps(sample_map, prefix, ["EAS", "EUR"], tmp_path / "ref")
    assert got["EAS"].read_text().strip() == "HG00403"
    assert got["EUR"].read_text().strip() == "HG00096"


def test_runtime_workers_do_not_subset_full_1kg():
    for name in [
        "41_prepare_susie_coloc_inputs.py",
        "42a_run_coloc_shard.sh",
        "44_prepare_colocboost_inputs.py",
        "45a_run_colocboost_shard.sh",
        "48_fusion_twas.sh",
    ]:
        text = (HERE / name).read_text()
        # Documentation may mention --keep; executable argument literals must not.
        assert '"--keep"' not in text, name
        assert "'--keep'" not in text, name
        assert not any(line.strip().startswith("--keep ")
                       for line in text.splitlines()), name
