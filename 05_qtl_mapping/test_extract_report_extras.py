"""Fixture tests for 35_extract_report_extras.py.

Exercises the OLS core and the TargetScanner lookup engine on a synthetic
fixture with planted effects — no tensorqtl / pgen / seadragon resources
required (genotypes are supplied as an in-memory dosage DataFrame, BEDs and
covariates as tiny TSVs).

Planted model: y = beta * g + 0.3 * cov + noise, so the lookup must recover
beta for the causal pair and ~0 for a null pair, and return NaN for IDs
absent from the target.
"""

import gzip
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "mod35", HERE / "35_extract_report_extras.py")
mod35 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod35)

RNG = np.random.default_rng(42)
N = 120  # samples


def make_fixture(tmp_path, anc="EAS", mod="expression", beta=1.5):
    """Write a tiny phenotype BED + covariates; return (qtl_dir, geno_df,
    truth) where geno_df is variants x samples with variant_id index."""
    qtl_dir = tmp_path / "qtl_inputs"
    qtl_dir.mkdir(exist_ok=True)
    samples = [f"S{i:03d}" for i in range(N)]
    cov = RNG.normal(size=N)
    g_causal = RNG.binomial(2, 0.3, size=N).astype(float)
    g_null = RNG.binomial(2, 0.4, size=N).astype(float)
    y1 = beta * g_causal + 0.3 * cov + RNG.normal(scale=0.2, size=N)
    y2 = 0.0 * g_null + 0.3 * cov + RNG.normal(scale=0.2, size=N)
    # BED (bgzipped): #chr start end phenotype_id samples...
    bed = pd.DataFrame({
        "#chr": ["chr5", "chr5"], "start": [1000, 2000], "end": [2000, 3000],
        "phenotype_id": ["PHENO_A", "PHENO_B"]})
    vals = pd.DataFrame([y1, y2], columns=samples)
    bed = pd.concat([bed, vals], axis=1)
    bed_path = qtl_dir / f"{anc}_{mod}.bed.gz"
    with gzip.open(bed_path, "wt") as f:
        bed.to_csv(f, sep="\t", index=False)
    # covariates TSV: rows = covariates, cols = samples
    covdf = pd.DataFrame([cov], index=["COV1"], columns=samples)
    covdf.to_csv(qtl_dir / f"{anc}_covariates_{mod}.tsv", sep="\t")
    # genotypes as DataFrame (variant_id index like genotypeio output)
    geno = pd.DataFrame(
        [g_causal, g_null],
        index=["5:10000:A:G", "5:20000:C:T"], columns=samples)
    truth = {"causal": ("5:10000:A:G", "PHENO_A", beta),
             "null": ("5:20000:C:T", "PHENO_B", 0.0)}
    return qtl_dir, geno, truth


def test_residualize_removes_covariate():
    C = RNG.normal(size=(N, 2))
    M = np.stack([C[:, 0] * 2 + C[:, 1], RNG.normal(size=N)])
    R = mod35.residualize(M, C)
    assert np.abs(R[0]).max() < 1e-8       # fully explained -> ~0
    assert np.abs(R[1].mean()) < 0.2       # intercept removed


def test_pair_lookup_recovers_planted_slope(tmp_path):
    qtl_dir, geno, truth = make_fixture(tmp_path)
    sc = mod35.TargetScanner(str(qtl_dir), "EAS", "expression", geno)
    pairs = pd.DataFrame({
        "variant_id": [truth["causal"][0], truth["null"][0], "5:99999:X:Y"],
        "phenotype_id": [truth["causal"][1], truth["null"][1], "PHENO_A"]})
    res = sc.lookup(pairs)
    b_causal = res.loc[0, "repl_slope"]
    assert b_causal == pytest.approx(truth["causal"][2], rel=0.05)
    assert res.loc[0, "repl_pval"] < 1e-8
    assert abs(res.loc[1, "repl_slope"]) < 0.2
    assert res.loc[1, "repl_pval"] > 0.05
    assert np.isnan(res.loc[2, "repl_pval"])  # variant absent -> NaN


def test_pair_lookup_missing_phenotype(tmp_path):
    qtl_dir, geno, truth = make_fixture(tmp_path)
    sc = mod35.TargetScanner(str(qtl_dir), "EAS", "expression", geno)
    pairs = pd.DataFrame({"variant_id": [truth["causal"][0]],
                          "phenotype_id": ["PHENO_MISSING"]})
    res = sc.lookup(pairs)
    assert np.isnan(res.loc[0, "repl_pval"])


def test_load_leads_layers(tmp_path):
    res_dir = tmp_path / "qtl_results"
    res_dir.mkdir()
    base = pd.DataFrame({
        "phenotype_id": ["P1", "P2"], "group_id": ["G1", "G2"],
        "variant_id": ["1:1:A:T", "1:2:A:T"],
        "slope": [1.0, 0.5], "slope_se": [0.1, 0.1],
        "pval_nominal": [1e-9, 1e-3], "qval": [0.001, 0.4]})
    base.to_csv(res_dir / "EAS_expression_cisqtl_top.tsv", sep="\t", index=False)
    base.to_csv(res_dir / "EAS_splicing_ungrouped_cisqtl_top.tsv", sep="\t",
                index=False)
    leads = mod35.load_leads(str(res_dir), "EAS", 0.05)
    assert set(leads["layer"]) == {"grouped", "ungrouped"}
    assert leads["is_lead"].sum() == 2  # one per layer
    assert set(leads["modality"]) == {"expression", "splicing"}


def test_pick_showcase_loci(tmp_path):
    res_dir = tmp_path / "qtl_results"
    res_dir.mkdir()
    comb = pd.DataFrame({
        "group_id": ["ENSG00000164308", "ENSG00000164307", "GENE_X", "GENE_Y"],
        "phenotype_id": ["p1", "p2", "px", "py"],
        "variant_id": ["5:1:A:T"] * 4,
        "modality": ["expression", "alt_polyA", "splicing", "expression"],
        "qval": [1e-20, 1e-30, 1e-8, 1e-6]})
    for anc in ["EAS", "EUR"]:
        comb.to_csv(res_dir / f"{anc}_combined_cisqtl_top.tsv",
                    sep="\t", index=False)
    loci = mod35.pick_showcase_loci(str(res_dir), ["EAS", "EUR"], 0.05)
    assert loci["ERAP2"] == "ENSG00000164308"
    # auto pick: shared, non-expression driver, not ERAP1/2 -> GENE_X
    assert "GENE_X" in loci.values()
    assert "ENSG00000164307" not in loci.values()


def test_pick_showcase_loci_versioned_ids(tmp_path):
    """Real combined tables carry Ensembl version suffixes: ERAP1/ERAP2
    exclusion and cross-ancestry intersection must be version-insensitive,
    and the returned id must be unversioned."""
    res_dir = tmp_path / "qtl_results"
    res_dir.mkdir()
    comb = pd.DataFrame({
        "group_id": ["ENSG00000164308.10", "ENSG00000164307.13",
                     "ENSG0000099999.5"],
        "phenotype_id": ["p1", "p2", "px"],
        "variant_id": ["5:1:A:T"] * 3,
        "modality": ["expression", "splicing", "splicing"],
        "qval": [1e-20, 1e-30, 1e-8]})
    for anc in ["EAS", "EUR"]:
        comb.to_csv(res_dir / f"{anc}_combined_cisqtl_top.tsv",
                    sep="\t", index=False)
    loci = mod35.pick_showcase_loci(str(res_dir), ["EAS", "EUR"], 0.05)
    assert loci["ERAP2"] == "ENSG00000164308"
    # versioned ERAP1 must be excluded; auto pick is the unversioned GENE
    assert "ENSG0000099999" in loci.values()
    assert all("." not in v for v in loci.values())


def test_strip_ver():
    s = pd.Series(["ENSG00000164307.13", "ENSG00000164308", "HPLRG_0001"])
    out = mod35._strip_ver(s)
    assert list(out) == ["ENSG00000164307", "ENSG00000164308", "HPLRG_0001"]
    assert mod35._unver("ENSG00000164307.13") == "ENSG00000164307"
