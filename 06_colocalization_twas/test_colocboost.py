"""Fixture test for 44_prepare_colocboost_inputs.py -> 44_colocboost.R.

Reuses the simulated panels from test_susie_coloc.py. Builds one region with
three outcomes sharing the same causal variant (two xQTL phenotypes + one
GWAS) and checks that colocBoost reports a cluster containing all three with
the causal variant on top. A second region whose GWAS has a distinct causal
variant should yield no cross-outcome cluster containing the GWAS.

Requires plink2, bgzip, tabix on PATH and R with colocboost installed.
Run:  pytest test_colocboost.py -v
"""

import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from test_susie_coloc import (CHROM, POSITIONS, VAR_IDS, make_pgen, regress,
                              simulate_genotypes, write_tabix_tsv, write_vcf)

HERE = Path(__file__).resolve().parent
PREP_SCRIPT = HERE / "44_prepare_colocboost_inputs.py"
R_SCRIPT = HERE / "44_colocboost.R"
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")
RNG = np.random.default_rng(7)


@pytest.fixture(scope="module")
def cb_fixture(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("colocboost")

    g_xqtl = simulate_genotypes(300, seed=11)
    g_1kg = simulate_genotypes(500, seed=12)
    xqtl_samples = [f"XQ{i}" for i in range(300)]
    kg_samples = [f"KG{i}" for i in range(500)]
    write_vcf(g_xqtl, VAR_IDS, xqtl_samples, tmp / "xqtl.vcf")
    write_vcf(g_1kg, VAR_IDS, kg_samples, tmp / "kg.vcf")
    make_pgen(tmp / "xqtl.vcf", tmp / "EAS_qtl")
    make_pgen(tmp / "kg.vcf", tmp / "1kg")
    pd.DataFrame({"FID": kg_samples, "IID": kg_samples}).to_csv(
        tmp / "EAS.1kg.keep", sep="\t", index=False, header=False)

    causal = 150
    causal_other = 320

    def xqtl_rows(pheno, mod, y):
        slope, se, p = regress(g_xqtl, y)
        return pd.DataFrame({
            "chrom": CHROM, "pos": POSITIONS, "phenotype_id": pheno,
            "variant_id": VAR_IDS, "tss_distance": POSITIONS - POSITIONS[0],
            "af": g_xqtl.mean(0) / 2, "ma_samples": 300, "ma_count": 100,
            "pval_nominal": p, "slope": slope, "slope_se": se})

    def gwas_rows(trait, y):
        slope, se, p = regress(g_1kg, y)
        ref = [v.split(":")[2] for v in VAR_IDS]
        alt = [v.split(":")[3] for v in VAR_IDS]
        return pd.DataFrame({
            "chr": CHROM, "pos": POSITIONS, "rsid": ".",
            "effect_allele": alt, "other_allele": ref,
            "eaf": g_1kg.mean(0) / 2, "beta": slope, "se": se, "pval": p,
            "n": 500, "var_id": VAR_IDS, "trait_id": trait})

    # region 1: all three outcomes share causal variant 150
    y_e1 = 0.6 * g_xqtl[:, causal] + RNG.normal(0, 1, 300)
    y_s1 = 0.5 * g_xqtl[:, causal] + RNG.normal(0, 1, 300)
    y_g1 = 0.5 * g_1kg[:, causal] + RNG.normal(0, 1, 500)
    # region 2 (same coordinates in the fixture; separate region_id): GWAS
    # driven by an independent variant
    y_g2 = 0.5 * g_1kg[:, causal_other] + RNG.normal(0, 1, 500)

    write_tabix_tsv(pd.concat([xqtl_rows("GENEA", "expression", y_e1)]),
                    tmp / "EAS_expression.nominal.tsv", "chrom", "pos")
    write_tabix_tsv(pd.concat([xqtl_rows("GENEA_clu", "splicing", y_s1)]),
                    tmp / "EAS_splicing.nominal.tsv", "chrom", "pos")
    write_tabix_tsv(gwas_rows("trait1", y_g1), tmp / "trait1.sumstats.tsv",
                    "chr", "pos")
    write_tabix_tsv(gwas_rows("trait2", y_g2), tmp / "trait2.sumstats.tsv",
                    "chr", "pos")

    start, end = int(POSITIONS[0] - 1000), int(POSITIONS[-1] + 1000)
    regions = pd.DataFrame([
        dict(region_id="EAS_1_r1", ancestry="EAS", chrom=CHROM, start=start,
             end=end, n_xqtl_outcomes=2, n_gwas_traits=1, n_outcomes=3),
        dict(region_id="EAS_1_r2", ancestry="EAS", chrom=CHROM, start=start,
             end=end, n_xqtl_outcomes=2, n_gwas_traits=1, n_outcomes=3),
    ])
    regions.to_csv(tmp / "EAS.regions.tsv", sep="\t", index=False)

    def oc(region, name, typ, mod, pheno, trait, file, side):
        return dict(region_id=region, outcome_name=name, type=typ,
                    modality=mod, phenotype_id=pheno, trait_id=trait,
                    trait_type="quantitative", prop_cases="", file=file,
                    ld_side=side)

    outcomes = pd.DataFrame([
        oc("EAS_1_r1", "expression:GENEA", "xqtl", "expression", "GENEA", "",
           str(tmp / "EAS_expression.nominal.tsv.gz"), "xqtl"),
        oc("EAS_1_r1", "splicing:GENEA_clu", "xqtl", "splicing", "GENEA_clu",
           "", str(tmp / "EAS_splicing.nominal.tsv.gz"), "xqtl"),
        oc("EAS_1_r1", "gwas:trait1", "gwas", "", "", "trait1",
           str(tmp / "trait1.sumstats.tsv.gz"), "gwas"),
        oc("EAS_1_r2", "expression:GENEA", "xqtl", "expression", "GENEA", "",
           str(tmp / "EAS_expression.nominal.tsv.gz"), "xqtl"),
        oc("EAS_1_r2", "splicing:GENEA_clu", "xqtl", "splicing", "GENEA_clu",
           "", str(tmp / "EAS_splicing.nominal.tsv.gz"), "xqtl"),
        oc("EAS_1_r2", "gwas:trait2", "gwas", "", "", "trait2",
           str(tmp / "trait2.sumstats.tsv.gz"), "gwas"),
    ])
    outcomes.to_csv(tmp / "EAS.outcomes.tsv", sep="\t", index=False)
    return tmp


def run_worker(tmp):
    prep_dir = tmp / "prepared_colocboost"
    prep_dir.mkdir(exist_ok=True)
    manifest = prep_dir / "prepared.regions.tsv"
    env = {**os.environ, "SCRIPTS_DIR": str(HERE)}
    prep = subprocess.run(
        ["python3", str(PREP_SCRIPT),
         "--regions", str(tmp / "EAS.regions.tsv"),
         "--outcomes", str(tmp / "EAS.outcomes.tsv"),
         "--shard-index", "1", "--n-shards", "1",
         "--outdir", str(tmp / "cb"),
         "--work-dir", str(prep_dir), "--manifest", str(manifest),
         "--ld-xqtl-pgen", str(tmp / "EAS_qtl"),
         "--ld-gwas-pgen", str(tmp / "1kg"),
         "--ld-gwas-keep", str(tmp / "EAS.1kg.keep")],
        capture_output=True, text=True, timeout=3600, env=env)
    assert prep.returncode == 0, f"preparer failed:\n{prep.stdout}\n{prep.stderr}"
    proc = subprocess.run(
        [RSCRIPT, str(R_SCRIPT),
         "--prepared-manifest", str(manifest),
         "--outdir", str(tmp / "cb"), "--min-variants", "50"],
        capture_output=True, text=True, timeout=3600, env=env)
    assert proc.returncode == 0, f"worker failed:\n{proc.stdout}\n{proc.stderr}"
    return proc


def test_shared_causal_cluster(cb_fixture):
    run_worker(cb_fixture)
    out = cb_fixture / "cb" / "results" / "EAS" / "EAS_1_r1.clusters.tsv"
    assert out.exists(), "clusters.tsv missing for region 1"
    df = pd.read_csv(out, sep="\t")
    assert len(df) >= 1 and df["cos_id"].notna().any(), "no cluster reported"
    best = df.loc[df["top_variable_vcp"].idxmax()]
    outcomes = set(str(best["colocalized_outcomes"]).split("; "))
    assert {"expression:GENEA", "splicing:GENEA_clu", "gwas:trait1"} <= outcomes
    # top variable should be the causal variant or a tight LD proxy
    top_pos = int(str(best["top_variable"]).split(":")[1])
    assert abs(top_pos - POSITIONS[150]) <= 5000


def test_distinct_causal_no_joint_cluster(cb_fixture):
    out = cb_fixture / "cb" / "results" / "EAS" / "EAS_1_r2.clusters.tsv"
    diag = pd.read_csv(cb_fixture / "cb" / "diagnostics" /
                       "EAS.shard-0001.diagnostics.tsv", sep="\t")
    row = diag[diag["region_id"] == "EAS_1_r2"].iloc[0]
    assert row["status"] == "ok"
    if not out.exists():
        return  # no clusters at all is acceptable
    df = pd.read_csv(out, sep="\t")
    df = df[df["cos_id"].notna()]
    for _, r in df.iterrows():
        outcomes = set(str(r["colocalized_outcomes"]).split("; "))
        # the GWAS with an independent causal variant must not cluster with
        # the xQTL outcomes
        assert not ("gwas:trait2" in outcomes and len(outcomes) > 1), \
            f"unexpected joint cluster: {outcomes}"


def test_done_sentinels(cb_fixture):
    res_dir = cb_fixture / "cb" / "results" / "EAS"
    assert (res_dir / "EAS_1_r1.done").exists()
    assert (res_dir / "EAS_1_r2.done").exists()
