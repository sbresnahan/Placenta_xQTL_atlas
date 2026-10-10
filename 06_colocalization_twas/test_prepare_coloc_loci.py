#!/usr/bin/env python3
"""Unit tests for robust 1KG sample-map parsing in step 40."""

import importlib.util
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "prepare_coloc_loci", HERE / "40_prepare_coloc_loci.py")
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def test_read_sample_map_skips_provenance_and_accepts_igsr_headers(tmp_path):
    path = tmp_path / "1kg_sample_superpop.tsv"
    path.write_text(
        "# downloaded/generated for the 1KG GRCh38 reference\n"
        "# this free-text line has no delimiter\n"
        "Sample name\tPopulation code\tSuperpopulation code\n"
        "HG00096\tGBR\tEUR\n"
        "HG00403\tCHS\tEAS\n"
    )

    got = MOD.read_kg_sample_map(path)
    assert got.to_dict("records") == [
        {"sample_id": "HG00096", "superpop": "EUR"},
        {"sample_id": "HG00403", "superpop": "EAS"},
    ]


def test_read_sample_map_accepts_whitespace_and_normalizes_case(tmp_path):
    path = tmp_path / "map.txt"
    path.write_text(
        "sample_id superpop\n"
        "HG00096 eur\n"
        "HG00403 EAS\n"
    )

    got = MOD.read_kg_sample_map(path)
    assert list(got["sample_id"]) == ["HG00096", "HG00403"]
    assert list(got["superpop"]) == ["EUR", "EAS"]


def test_write_keep_files_uses_psam_fid_iid(tmp_path):
    sample_map = tmp_path / "map.tsv"
    sample_map.write_text(
        "sample_id\tsuperpop\n"
        "HG00096\tEUR\n"
        "HG00403\tEAS\n"
    )
    pgen_prefix = tmp_path / "kg"
    Path(f"{pgen_prefix}.psam").write_text(
        "#FID\tIID\tSEX\n"
        "FAM1\tHG00096\t1\n"
        "FAM2\tHG00403\t2\n"
    )

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    written = MOD.write_kg_keep_files(
        sample_map, ["EAS", "EUR"], out_dir, kg_pgen=pgen_prefix)

    eas = pd.read_csv(written["EAS"], sep="\t", header=None, dtype=str)
    eur = pd.read_csv(written["EUR"], sep="\t", header=None, dtype=str)
    assert eas.values.tolist() == [["FAM2", "HG00403"]]
    assert eur.values.tolist() == [["FAM1", "HG00096"]]


def test_write_keep_files_rejects_missing_requested_superpopulation(tmp_path):
    sample_map = tmp_path / "map.tsv"
    sample_map.write_text("sample_id\tsuperpop\nHG00096\tEUR\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    try:
        MOD.write_kg_keep_files(sample_map, ["EAS", "EUR"], out_dir)
    except SystemExit as exc:
        assert "does not contain requested superpopulation" in str(exc)
        assert "EAS" in str(exc)
    else:
        raise AssertionError("expected missing EAS to terminate step 40")


def test_task_generation_respects_sig_in(tmp_path, monkeypatch):
    results = tmp_path / "qtl_results"
    loci_dir = results / "finemap" / "loci"
    loci_dir.mkdir(parents=True)
    nominal = results / "nominal"
    for anc in ("EAS", "EUR"):
        anc_dir = nominal / anc
        anc_dir.mkdir(parents=True)
        (anc_dir / f"{anc}_expression.nominal.tsv.gz").touch()

    pd.DataFrame([
        {"modality": "expression", "phenotype_id": "gene_eas", "chrom": "1",
         "start": 100, "end": 200, "n_signals": 1, "L": 5, "sig_in": "EAS"},
        {"modality": "expression", "phenotype_id": "gene_eur", "chrom": "2",
         "start": 300, "end": 400, "n_signals": 1, "L": 5, "sig_in": "EUR"},
        {"modality": "expression", "phenotype_id": "gene_both", "chrom": "3",
         "start": 500, "end": 600, "n_signals": 1, "L": 5,
         "sig_in": "EAS,EUR"},
    ]).to_csv(loci_dir / "expression.loci.tsv", sep="\t", index=False)

    gwas_dir = tmp_path / "gwas"
    gwas_dir.mkdir()
    (gwas_dir / "trait_both.sumstats.tsv.gz").touch()
    catalog = tmp_path / "gwas_catalog.tsv"
    pd.DataFrame([{
        "trait_id": "trait_both",
        "primary_ancestry": "both",
        "trait_type": "quantitative",
    }]).to_csv(catalog, sep="\t", index=False)

    out_dir = tmp_path / "coloc"
    qtl_dir = tmp_path / "qtl_inputs"
    qtl_dir.mkdir()
    argv = [
        "40_prepare_coloc_loci.py",
        "--results-dir", str(results),
        "--catalog", str(catalog),
        "--gwas-dir", str(gwas_dir),
        "--qtl-dir", str(qtl_dir),
        "--kg-pgen", str(tmp_path / "unused_1kg"),
        "--skip-kg-keep",
        "--ancestries", "EAS", "EUR",
        "--modalities", "expression",
        "--out-dir", str(out_dir),
    ]
    monkeypatch.setattr("sys.argv", argv)
    MOD.main()

    tasks = pd.read_csv(out_dir / "loci" / "expression.tasks.tsv", sep="\t")
    got = set(zip(tasks["phenotype_id"], tasks["ancestry"]))
    assert got == {
        ("gene_eas", "EAS"),
        ("gene_eur", "EUR"),
        ("gene_both", "EAS"),
        ("gene_both", "EUR"),
    }
    assert all(
        row.ancestry in row.sig_in.split(",")
        for row in tasks[["ancestry", "sig_in"]].itertuples(index=False)
    )


def test_task_generation_requires_sig_in(tmp_path, monkeypatch):
    results = tmp_path / "qtl_results"
    loci_dir = results / "finemap" / "loci"
    loci_dir.mkdir(parents=True)
    pd.DataFrame([{
        "phenotype_id": "gene1", "chrom": "1", "start": 100,
        "end": 200, "L": 5,
    }]).to_csv(loci_dir / "expression.loci.tsv", sep="\t", index=False)

    catalog = tmp_path / "gwas_catalog.tsv"
    pd.DataFrame(columns=["trait_id", "primary_ancestry"]).to_csv(
        catalog, sep="\t", index=False)
    gwas_dir = tmp_path / "gwas"
    gwas_dir.mkdir()
    qtl_dir = tmp_path / "qtl_inputs"
    qtl_dir.mkdir()
    out_dir = tmp_path / "coloc"

    argv = [
        "40_prepare_coloc_loci.py",
        "--results-dir", str(results),
        "--catalog", str(catalog),
        "--gwas-dir", str(gwas_dir),
        "--qtl-dir", str(qtl_dir),
        "--kg-pgen", str(tmp_path / "unused_1kg"),
        "--skip-kg-keep",
        "--modalities", "expression",
        "--out-dir", str(out_dir),
    ]
    monkeypatch.setattr("sys.argv", argv)
    try:
        MOD.main()
    except SystemExit as exc:
        assert "has no sig_in column" in str(exc)
    else:
        raise AssertionError("expected missing sig_in to terminate step 40")
