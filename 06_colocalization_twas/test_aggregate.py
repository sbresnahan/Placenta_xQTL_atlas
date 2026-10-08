"""Smoke test for 49_aggregate_coloc_twas.py.

Builds a fake results tree (coloc results + diagnostics, colocBoost clusters,
FUSION TWAS outputs, .pos files, expression BED + covariates + deconvolution)
and checks the aggregated tables.
Run:  pytest test_aggregate.py -v
"""

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "49_aggregate_coloc_twas.py"
RNG = np.random.default_rng(5)


@pytest.fixture(scope="module")
def results_tree(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("agg")
    results = tmp / "qtl_results"
    qtl = tmp / "qtl_inputs"
    coloc = results / "coloc"

    # --- coloc results (2 tasks, one with a PP.H4 call) ----------------------
    rdir = coloc / "results" / "expression"
    rdir.mkdir(parents=True)
    base_cols = dict(modality="expression", ancestry="EAS", chrom="1",
                     start=100, end=200, nsnps=300, n_xqtl=280, n_gwas=1000,
                     converged_xqtl=True, converged_gwas=True, n_cs_xqtl=2,
                     n_cs_gwas=1)
    pd.DataFrame([
        {**base_cols, "phenotype_id": "GENE1", "trait_id": "traitA",
         "hit1": "v1", "hit2": "v2", "PP.H0.abf": 0.01, "PP.H1.abf": 0.02,
         "PP.H2.abf": 0.02, "PP.H3.abf": 0.10, "PP.H4.abf": 0.85,
         "idx1": 1, "idx2": 1, "coloc_call": True},
        {**base_cols, "phenotype_id": "GENE2", "trait_id": "traitA",
         "hit1": "v3", "hit2": "v4", "PP.H0.abf": 0.01, "PP.H1.abf": 0.02,
         "PP.H2.abf": 0.02, "PP.H3.abf": 0.80, "PP.H4.abf": 0.15,
         "idx1": 1, "idx2": 1, "coloc_call": False},
    ]).to_csv(rdir / "x.coloc.tsv", sep="\t", index=False)
    ddir = coloc / "diagnostics"
    ddir.mkdir(parents=True)
    pd.DataFrame({"modality": ["expression"] * 2,
                  "phenotype_id": ["GENE1", "GENE2"], "trait_id": ["traitA"] * 2,
                  "ancestry": ["EAS"] * 2, "status": ["ok", "ok"]}
                 ).to_csv(ddir / "expression.shard-0001.diagnostics.tsv",
                          sep="\t", index=False)

    # --- colocBoost clusters --------------------------------------------------
    cdir = coloc / "colocboost" / "results" / "EAS"
    cdir.mkdir(parents=True)
    pd.DataFrame({"region_id": ["EAS_1_1_2"], "ancestry": ["EAS"],
                  "cos_id": ["cos1"], "colocalized_outcomes": ["a; b"],
                  "top_variable": ["1:100:A:G"], "top_variable_vcp": [0.9]}
                 ).to_csv(cdir / "EAS_1_1_2.clusters.tsv", sep="\t", index=False)

    # --- FUSION TWAS outputs + pos ---------------------------------------------
    fdir = results / "isotwas" / "fusion" / "EAS"
    fdir.mkdir(parents=True)
    pd.DataFrame({
        "PANEL": ["NA"] * 3, "FILE": ["f"] * 3,
        "ID": ["GENE1", "GENE1__T1", "GENE3"], "CHR": [1] * 3,
        "P0": [0] * 3, "P1": [1] * 3, "HSQ": ["NA"] * 3,
        "BEST.GWAS.ID": ["v"] * 3, "BEST.GWAS.Z": [3.0] * 3,
        "EQTL.ID": ["v"] * 3, "EQTL.R2": [0.1] * 3, "EQTL.Z": [3.0] * 3,
        "EQTL.GWAS.Z": [3.0] * 3, "NSNP": [100] * 3, "NWGT": [5] * 3,
        "MODEL": ["m"] * 3, "MODELCV.R2": [0.1] * 3, "MODELCV.PV": [1e-4] * 3,
        "TWAS.Z": [4.0, 3.5, 0.1], "TWAS.P": [6e-5, 5e-4, 0.92],
    }).to_csv(fdir / "EAS_traitA.twas.tsv", sep=" ", index=False)
    wdir = results / "isotwas" / "weights"
    wdir.mkdir(parents=True)
    pd.DataFrame({"WGT": ["GENE1.wgt.RDat", "GENE1__T1.wgt.RDat",
                          "GENE3.wgt.RDat"],
                  "ID": ["GENE1", "GENE1__T1", "GENE3"],
                  "GENE": ["GENE1", "GENE1", "GENE3"],
                  "CHR": [1] * 3, "P0": [0] * 3, "P1": [1] * 3}
                 ).to_csv(wdir / "EAS.pos", sep="\t", index=False)

    # --- expression BED + covariates + deconvolution ---------------------------
    qtl.mkdir(parents=True)
    samples = [f"S{i}" for i in range(60)]
    prop_a = RNG.uniform(0, 1, 60)
    expr_g1 = 2.0 * prop_a + RNG.normal(0, 0.5, 60)  # correlated with CT_A
    bed = pd.DataFrame({"#chr": ["1"] * 3, "start": [1] * 3, "end": [2] * 3,
                        "phenotype_id": ["GENE1", "GENE2", "GENE3"]})
    for i, s in enumerate(samples):
        bed[s] = [expr_g1[i], RNG.normal(0, 1), RNG.normal(0, 1)]
    bed.to_csv(qtl / "EAS_expression.bed.gz", sep="\t", index=False,
               compression="gzip")
    cov = pd.DataFrame({"cov": ["PC1"]})
    for s in samples:
        cov[s] = RNG.normal(0, 1, 1)
    cov.to_csv(qtl / "EAS_covariates_expression.tsv", sep="\t", index=False)
    pd.DataFrame({"sample_id": samples, "CT_A": prop_a,
                  "CT_B": 1 - prop_a}).to_csv(
        qtl / "EAS_deconvolution_harmonized.tsv", sep="\t", index=False)
    return results, qtl


def test_aggregation(results_tree):
    results, qtl = results_tree
    proc = subprocess.run(
        ["python3", str(SCRIPT), "--results-dir", str(results),
         "--qtl-dir", str(qtl), "--ancestries", "EAS"],
        capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, f"aggregator failed:\n{proc.stdout}\n{proc.stderr}"
    out = results / "coloc" / "aggregated"

    coloc = pd.read_csv(out / "coloc_results.tsv.gz", sep="\t")
    assert len(coloc) == 2
    best = pd.read_csv(out / "coloc_best.tsv.gz", sep="\t")
    assert best.loc[best["phenotype_id"] == "GENE1", "PP.H4.abf"].iloc[0] == 0.85

    clu = pd.read_csv(out / "colocboost_clusters.tsv.gz", sep="\t")
    assert len(clu) == 1

    twas = pd.read_csv(out / "twas_results.tsv.gz", sep="\t")
    assert {"GENE1", "GENE1__T1", "GENE3"} == set(twas["ID"])
    assert "TWAS.Q" in twas.columns
    genes = pd.read_csv(out / "twas_gene_results.tsv.gz", sep="\t")
    g1 = genes[genes["GENE"] == "GENE1"].iloc[0]
    assert g1["n_models"] == 2
    assert g1["acat_p"] < 1e-3
    g3 = genes[genes["GENE"] == "GENE3"].iloc[0]
    assert g3["acat_p"] > 0.5

    ann = pd.read_csv(out / "gene_celltype_annotation.tsv", sep="\t")
    row = ann[ann["gene"] == "GENE1"].iloc[0]
    assert row["primary_cell_type_EAS"] == "CT_A"
    assert row["max_abs_spearman_EAS"] > 0.2
