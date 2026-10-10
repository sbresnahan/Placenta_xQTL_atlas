"""Fixture test for 46_isotwas_train.R + FUSION round-trip.

Simulates two ancestry panels (EAS 200, EUR 150) over one locus. GENE1 has
two isoforms with distinct causal variants (idx 100 and 200) and its
gene-level expression is driven by variant 100; GENE2 is null. Checks:
  1. GENE1 .wgt.RDat loads, contains gene + isoform model columns with
     CV R^2 > 0.01, and the causal variants carry non-zero weights.
  2. GENE2 is gated out (no_heritable_model).
  3. Pooled weight set trains and stacks both ancestries (N.tot = 350).
  4. FUSION round-trip: FUSION.assoc_test.R runs on the trained weights with
     a 1KG-like LD reference and reports a TWAS association for GENE1.

Requires plink2, bgzip, tabix, R with glmnet, and the FUSION repo clone
(path via FUSION_DIR env or /workspace/fusion_twas).
Run:  pytest test_isotwas_train.py -v
"""

import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from test_susie_coloc import (CHROM, POSITIONS, VAR_IDS, make_pgen, regress,
                              simulate_genotypes, write_vcf)

HERE = Path(__file__).resolve().parent
PREP_SCRIPT = HERE / "46_prepare_isotwas_inputs.py"
R_SCRIPT = HERE / "46_isotwas_train.R"
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")
FUSION_DIR = Path(os.environ.get("FUSION_DIR", "/workspace/fusion_twas"))
RNG = np.random.default_rng(13)

N_VARIANTS = len(VAR_IDS)


def write_bed(path, rows, samples):
    """rows: list of (phenotype_id, values array)."""
    df = pd.DataFrame({"#chr": CHROM, "start": 100_000, "end": 110_000,
                       "phenotype_id": [r[0] for r in rows]})
    for i, s in enumerate(samples):
        df[s] = [r[1][i] for r in rows]
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def write_covs(path, samples, seed=0):
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({"covariate": ["PC1", "PC2"]})
    for s in samples:
        df[s] = rng.normal(0, 1, 2)
    df.to_csv(path, sep="\t", index=False)


@pytest.fixture(scope="module")
def twas_fixture(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("isotwas")
    qtl = tmp / "qtl"
    qtl.mkdir()

    panels = {}
    for anc, n, seed in [("EAS", 200, 21), ("EUR", 150, 22)]:
        g = simulate_genotypes(n, seed=seed)
        samples = [f"{anc}{i}" for i in range(n)]
        write_vcf(g, VAR_IDS, samples, tmp / f"{anc}.vcf")
        make_pgen(tmp / f"{anc}.vcf", qtl / f"{anc}_qtl")
        panels[anc] = (g, samples)

    # phenotypes: GENE1 isoform1 ~ variant 100, isoform2 ~ variant 200,
    # expression ~ variant 100; GENE2 null
    for anc, (g, samples) in panels.items():
        n = len(samples)
        iso1 = 0.7 * g[:, 100] + RNG.normal(0, 1, n)
        iso2 = 0.7 * g[:, 200] + RNG.normal(0, 1, n)
        expr = 0.6 * g[:, 100] + RNG.normal(0, 1, n)
        null = RNG.normal(0, 1, n)
        write_bed(qtl / f"{anc}_expression.bed.gz",
                  [("GENE1", expr), ("GENE2", null)], samples)
        write_bed(qtl / f"{anc}_isoform_expression.bed.gz",
                  [("GENE1__T1", iso1), ("GENE1__T2", iso2),
                   ("GENE2__T1", null)], samples)
        write_covs(qtl / f"{anc}_covariates_expression.tsv", samples, seed=1)
        write_covs(qtl / f"{anc}_covariates_isoform_expression.tsv", samples,
                   seed=2)

    return tmp


def run_trainer(tmp, weight_set):
    outdir = tmp / "isotwas"
    prep_dir = outdir / "prepared_test" / weight_set
    prep_dir.mkdir(parents=True, exist_ok=True)
    manifest = prep_dir / "shard-0001.tsv"
    env = {**os.environ, "SCRIPTS_DIR": str(HERE)}
    prep = subprocess.run(
        ["python3", str(PREP_SCRIPT), "--weight-set", weight_set,
         "--shard-index", "1", "--n-shards", "1",
         "--qtl-dir", str(tmp / "qtl"), "--outdir", str(outdir),
         "--work-dir", str(prep_dir), "--manifest", str(manifest),
         "--force"],
        capture_output=True, text=True, timeout=3600, env=env)
    assert prep.returncode == 0, f"preparer failed:\n{prep.stdout}\n{prep.stderr}"
    cmd = [RSCRIPT, str(R_SCRIPT), "--weight-set", weight_set,
           "--shard-index", "1", "--prepared-manifest", str(manifest),
           "--qtl-dir", str(tmp / "qtl"), "--outdir", str(outdir),
           "--min-variants", "20", "--r2-min", "0.05"]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, env=env)
    assert proc.returncode == 0, f"trainer failed:\n{proc.stdout}\n{proc.stderr}"
    return proc


def load_rdat(path):
    """Load a FUSION .wgt.RDat via R and return wgt.matrix + cv.performance."""
    r = HERE / "_dump_rdat.R"
    r.write_text(
        'load(commandArgs(TRUE)[1])\n'
        'out <- commandArgs(TRUE)[2]\n'
        'write.csv(as.data.frame(wgt.matrix), paste0(out, ".wgt.csv"))\n'
        'write.csv(as.data.frame(cv.performance), paste0(out, ".cvp.csv"))\n'
        'write.csv(data.frame(N = N.tot), paste0(out, ".n.csv"))\n')
    subprocess.run([RSCRIPT, str(r), str(path), str(path)], check=True)
    wgt = pd.read_csv(str(path) + ".wgt.csv", index_col=0)
    cvp = pd.read_csv(str(path) + ".cvp.csv", index_col=0)
    n = pd.read_csv(str(path) + ".n.csv")["N"].iloc[0]
    return wgt, cvp, n


def test_gene1_weights_trained(twas_fixture):
    run_trainer(twas_fixture, "EAS")
    genes_dir = twas_fixture / "isotwas" / "weights" / "EAS" / "genes"
    # per-model files: expression model + at least one isoform model
    expr_rdat = genes_dir / "GENE1.wgt.RDat"
    assert expr_rdat.exists(), "GENE1 expression weights missing"
    iso_rdats = sorted(genes_dir.glob("GENE1__*.wgt.RDat"))
    assert iso_rdats, "no isoform weight files for GENE1"
    wgt, cvp, n = load_rdat(expr_rdat)
    assert n == 200
    assert list(wgt.columns) == ["GENE1"]
    assert cvp.loc["rsq", "GENE1"] > 0.05
    # causal variants carry non-zero weights in the relevant models
    v100, v200 = VAR_IDS[100], VAR_IDS[200]
    assert wgt.loc[v100, "GENE1"] != 0
    # isoform models: T1 driven by v100, T2 by v200
    by_model = {}
    for rd in iso_rdats:
        w, c, _ = load_rdat(rd)
        by_model[rd.name.replace(".wgt.RDat", "")] = (w, c)
        assert (c.loc["rsq"] > 0.05).all()
    if "GENE1__T1" in by_model:
        assert by_model["GENE1__T1"][0].loc[v100].iloc[0] != 0
    if "GENE1__T2" in by_model:
        assert by_model["GENE1__T2"][0].loc[v200].iloc[0] != 0


def test_gene2_gated_out(twas_fixture):
    diag = pd.read_csv(twas_fixture / "isotwas" / "diagnostics" /
                       "isotwas_EAS.shard-0001.diagnostics.tsv", sep="\t")
    row = diag[diag["gene"] == "GENE2"].iloc[0]
    assert row["status"] == "no_heritable_model"
    assert not (twas_fixture / "isotwas" / "weights" / "EAS" / "genes" /
                "GENE2.wgt.RDat").exists()


def test_pooled_weight_set(twas_fixture):
    run_trainer(twas_fixture, "pooled")
    rdat = (twas_fixture / "isotwas" / "weights" / "pooled" / "genes" /
            "GENE1.wgt.RDat")
    assert rdat.exists()
    wgt, cvp, n = load_rdat(rdat)
    assert n == 350  # 200 EAS + 150 EUR stacked


def test_fusion_round_trip(twas_fixture):
    """FUSION.assoc_test.R must load our weights and report GENE1."""
    if not (FUSION_DIR / "FUSION.assoc_test.R").exists():
        pytest.skip("FUSION clone not available")
    tmp = twas_fixture
    # LD reference: plink1 BED of the EAS panel, chr-prefixed as FUSION wants
    subprocess.run(["plink2", "--pfile", str(tmp / "qtl" / "EAS_qtl"),
                    "--chr", "1", "--make-bed", "--out",
                    str(tmp / "ldref.chr1"), "--silent"], check=True)
    # FUSION-format GWAS sumstats: trait driven by variant 100
    g = simulate_genotypes(200, seed=21)  # same panel as EAS training
    y = 0.5 * g[:, 100] + RNG.normal(0, 1, 200)
    slope, se, p = regress(g, y)
    ref = [v.split(":")[2] for v in VAR_IDS]
    alt = [v.split(":")[3] for v in VAR_IDS]
    pd.DataFrame({"SNP": VAR_IDS, "A1": alt, "A2": ref,
                  "Z": slope / se}).to_csv(tmp / "fusion.sumstats",
                                           sep=" ", index=False)
    # .pos file from the shard fragment
    pos_src = (tmp / "isotwas" / "weights" / "EAS" / "shard-0001.pos")
    pos = pd.read_csv(pos_src, sep="\t")
    pos = pos[pos["GENE"] == "GENE1"]
    pos_path = tmp / "fusion.pos"
    pos.to_csv(pos_path, sep="\t", index=False)
    out_prefix = tmp / "fusion_out"
    # FUSION.assoc_test.R sources utils/ via here(): run from FUSION_DIR
    proc = subprocess.run(
        [RSCRIPT, str(FUSION_DIR / "FUSION.assoc_test.R"),
         "--sumstats", str(tmp / "fusion.sumstats"),
         "--weights", str(pos_path),
         "--weights_dir", str(tmp / "isotwas" / "weights" / "EAS" / "genes"),
         "--ref_ld_chr", str(tmp / "ldref.chr"), "--chr", "1",
         "--min_r2pred", "0", "--out", str(out_prefix)],
        capture_output=True, text=True, timeout=1800, cwd=str(FUSION_DIR))
    assert proc.returncode == 0, \
        f"FUSION.assoc_test failed:\n{proc.stdout}\n{proc.stderr}"
    res = pd.read_csv(out_prefix, sep=r"\s+")  # FUSION writes to --out verbatim
    assert "GENE1" in set(res["ID"])
    row = res[res["ID"] == "GENE1"].iloc[0]
    assert row["TWAS.P"] < 0.01, f"TWAS.P = {row['TWAS.P']}"
    # isoform 1 shares the GWAS causal variant; isoform 2 does not
    if "GENE1__T1" in set(res["ID"]):
        assert res.loc[res["ID"] == "GENE1__T1", "TWAS.P"].iloc[0] < 0.01
    if "GENE1__T2" in set(res["ID"]):
        assert res.loc[res["ID"] == "GENE1__T2", "TWAS.P"].iloc[0] > 0.01
