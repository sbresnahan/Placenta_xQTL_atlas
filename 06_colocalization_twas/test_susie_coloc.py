"""Fixture test for 41_prepare_susie_coloc_inputs.py -> 41_susie_coloc.R.

Simulates two small genotype panels (an in-sample "xQTL" pgen and a
reference "1KG" pgen) over one locus with blocky LD, then generates:
  - task A: xQTL phenotype and GWAS trait driven by the SAME causal variant
            -> expect PP.H4 >= 0.7 (colocalized)
  - task B: xQTL phenotype and GWAS trait driven by two DIFFERENT,
            low-LD causal variants -> expect PP.H4 < 0.7

Requires plink2, bgzip, tabix on PATH, Python/NumPy, and R with coloc installed.
Run:  pytest test_susie_coloc.py -v
"""

import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
PREP_SCRIPT = HERE / "41_prepare_susie_coloc_inputs.py"
R_SCRIPT = HERE / "41_susie_coloc.R"
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")

RNG = np.random.default_rng(42)
N_VARIANTS = 400
CHROM = "1"
POSITIONS = 100_000 + 100 * np.arange(N_VARIANTS)
ALLELES = RNG.choice(["A", "C", "G", "T"], size=(N_VARIANTS, 2), p=None)
# ensure ref != alt
for i in range(N_VARIANTS):
    while ALLELES[i, 0] == ALLELES[i, 1]:
        ALLELES[i, 1] = RNG.choice(["A", "C", "G", "T"])
VAR_IDS = np.array([f"{CHROM}:{p}:{a[0]}:{a[1]}"
                    for p, a in zip(POSITIONS, ALLELES)])
MAFS = RNG.uniform(0.1, 0.4, N_VARIANTS)


def simulate_genotypes(n_samples, rho=0.85, seed=0):
    """Markov-chain haplotypes -> diploid dosages with decaying LD."""
    rng = np.random.default_rng(seed)
    haps = np.zeros((n_samples * 2, N_VARIANTS), dtype=np.int8)
    haps[:, 0] = rng.binomial(1, MAFS[0], n_samples * 2)
    for j in range(1, N_VARIANTS):
        flip = rng.binomial(1, 1 - rho, n_samples * 2)
        draw = rng.binomial(1, MAFS[j], n_samples * 2)
        haps[:, j] = np.where(flip, draw, haps[:, j - 1])
    g = haps[0::2] + haps[1::2]
    return g.astype(np.float64)  # samples x variants


def write_vcf(geno, var_ids, sample_ids, path):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(sample_ids) + "\n")
        for j, vid in enumerate(var_ids):
            c, p, ref, alt = vid.split(":")
            gts = "\t".join({0: "0|0", 1: "0|1", 2: "1|1"}[int(x)]
                            for x in geno[:, j])
            f.write(f"{c}\t{p}\t{vid}\t{ref}\t{alt}\t.\t.\t.\tGT\t{gts}\n")


def make_pgen(vcf, out_prefix):
    subprocess.run(["plink2", "--vcf", str(vcf), "--double-id",
                    "--make-pgen", "--out", str(out_prefix), "--silent"],
                   check=True)


def regress(geno, y):
    """Per-variant simple linear regression -> slope, se, p (Wald)."""
    from scipy import stats
    n = len(y)
    y0 = y - y.mean()
    g0 = geno - geno.mean(0)
    denom = (g0 ** 2).sum(0)
    denom[denom == 0] = np.nan
    slope = (g0 * y0[:, None]).sum(0) / denom
    resid = y0[:, None] - slope[None, :] * g0
    s2 = (resid ** 2).sum(0) / (n - 2)
    se = np.sqrt(s2 / denom)
    t = slope / se
    p = 2 * stats.t.sf(np.abs(t), n - 2)
    return slope, se, p


def write_tabix_tsv(df, path, chrom_col, pos_col):
    df = df.sort_values([chrom_col, pos_col])
    df.to_csv(path, sep="\t", index=False)
    subprocess.run(f"bgzip -f {path} && tabix -f -s 1 -b 2 -e 2 -S 1 {path}.gz",
                   shell=True, check=True)


@pytest.fixture(scope="module")
def coloc_fixture(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("coloc")

    # --- genotype panels -----------------------------------------------------
    g_xqtl = simulate_genotypes(300, seed=1)
    g_1kg = simulate_genotypes(500, seed=2)
    xqtl_samples = [f"XQ{i}" for i in range(300)]
    kg_samples = [f"KG{i}" for i in range(500)]
    write_vcf(g_xqtl, VAR_IDS, xqtl_samples, tmp / "xqtl.vcf")
    write_vcf(g_1kg, VAR_IDS, kg_samples, tmp / "kg.vcf")
    make_pgen(tmp / "xqtl.vcf", tmp / "EAS_qtl")
    make_pgen(tmp / "kg.vcf", tmp / "1kg")
    # keep file: first 400 1KG samples (tests the --keep path); plink2
    # matches FID+IID, and --double-id sets FID==IID, so write two columns
    pd.DataFrame({"FID": kg_samples[:400], "IID": kg_samples[:400]}).to_csv(
        tmp / "EAS.1kg.keep", sep="\t", index=False, header=False)

    # --- phenotypes / summary stats ------------------------------------------
    causal_shared = 150
    causal_distinct_gwas = 320  # low LD with 150 (Markov rho^170 ~ 0)
    beta_x, beta_g = 0.6, 0.5

    y_xqtl_a = beta_x * g_xqtl[:, causal_shared] + RNG.normal(0, 1, 300)
    y_gwas_a = beta_g * g_1kg[:, causal_shared] + RNG.normal(0, 1, 500)
    y_xqtl_b = beta_x * g_xqtl[:, causal_shared] + RNG.normal(0, 1, 300)
    y_gwas_b = beta_g * g_1kg[:, causal_distinct_gwas] + RNG.normal(0, 1, 500)

    def xqtl_rows(pheno, y):
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

    xqtl = pd.concat([xqtl_rows("GENEA", y_xqtl_a), xqtl_rows("GENEB", y_xqtl_b)])
    write_tabix_tsv(xqtl, tmp / "EAS_expression.nominal.tsv", "chrom", "pos")
    for trait, y in [("trait_shared", y_gwas_a), ("trait_distinct", y_gwas_b)]:
        write_tabix_tsv(gwas_rows(trait, y), tmp / f"{trait}.sumstats.tsv",
                        "chr", "pos")

    # --- task list -------------------------------------------------------------
    start, end = POSITIONS[0] - 1000, POSITIONS[-1] + 1000
    tasks = pd.DataFrame([
        dict(modality="expression", phenotype_id="GENEA", chrom=CHROM,
             start=start, end=end, L=5, ancestry="EAS", trait_id="trait_shared",
             trait_type="quantitative", prop_cases="",
             gwas_file=str(tmp / "trait_shared.sumstats.tsv.gz"),
             xqtl_file=str(tmp / "EAS_expression.nominal.tsv.gz"),
             ld_xqtl_pgen=str(tmp / "EAS_qtl"), ld_gwas_pgen=str(tmp / "1kg"),
             ld_gwas_keep=str(tmp / "EAS.1kg.keep")),
        dict(modality="expression", phenotype_id="GENEB", chrom=CHROM,
             start=start, end=end, L=5, ancestry="EAS",
             trait_id="trait_distinct", trait_type="quantitative", prop_cases="",
             gwas_file=str(tmp / "trait_distinct.sumstats.tsv.gz"),
             xqtl_file=str(tmp / "EAS_expression.nominal.tsv.gz"),
             ld_xqtl_pgen=str(tmp / "EAS_qtl"), ld_gwas_pgen=str(tmp / "1kg"),
             ld_gwas_keep=str(tmp / "EAS.1kg.keep")),
    ])
    tasks.to_csv(tmp / "expression.tasks.tsv", sep="\t", index=False)
    return tmp


def run_worker(tmp):
    env = dict(os.environ)
    prep_dir = tmp / "prepared_susie"
    prep_dir.mkdir(exist_ok=True)
    manifest = prep_dir / "prepared.tasks.tsv"
    prep = subprocess.run(
        ["python3", str(PREP_SCRIPT),
         "--tasks", str(tmp / "expression.tasks.tsv"),
         "--shard-index", "1", "--n-shards", "1",
         "--work-dir", str(prep_dir / "work"),
         "--out", str(manifest), "--min-variants", "50"],
        capture_output=True, text=True, env=env, timeout=1800)
    assert prep.returncode == 0, f"preparer failed:\n{prep.stdout}\n{prep.stderr}"
    proc = subprocess.run(
        [RSCRIPT, str(R_SCRIPT), "--prepared-tasks", str(manifest),
         "--shard-index", "1", "--n-shards", "1",
         "--outdir", str(tmp / "coloc"), "--min-variants", "50"],
        capture_output=True, text=True, env=env, timeout=1800)
    assert proc.returncode == 0, f"worker failed:\n{proc.stdout}\n{proc.stderr}"
    return proc


def test_shared_causal_colocalizes(coloc_fixture):
    run_worker(coloc_fixture)
    out = (coloc_fixture / "coloc" / "results" / "expression" /
           "EAS_GENEA_trait_shared.coloc.tsv")
    assert out.exists(), "coloc.tsv missing for shared-causal task"
    df = pd.read_csv(out, sep="\t")
    assert df["PP.H4.abf"].max() >= 0.7, \
        f"expected PP.H4 >= 0.7, got {df['PP.H4.abf'].max():.3f}"
    assert bool(df.loc[df["PP.H4.abf"].idxmax(), "coloc_call"])
    # per-variant table written with SNP.PP.H4
    vpath = out.parent / out.name.replace(".coloc.tsv", ".variants.tsv")
    assert vpath.exists()
    v = pd.read_csv(vpath, sep="\t")
    assert "SNP.PP.H4_best" in v.columns and len(v) >= 50
    # the causal variant should carry high posterior for the shared task
    causal_id = VAR_IDS[150]
    assert v.loc[v["var_id"] == causal_id, "SNP.PP.H4_best"].iloc[0] > 0.01


def test_distinct_causal_does_not(coloc_fixture):
    out = (coloc_fixture / "coloc" / "results" / "expression" /
           "EAS_GENEB_trait_distinct.coloc.tsv")
    diag = pd.read_csv(coloc_fixture / "coloc" / "diagnostics" /
                       "expression.shard-0001.diagnostics.tsv", sep="\t")
    row = diag[diag["phenotype_id"] == "GENEB"].iloc[0]
    if not out.exists():
        # acceptable: one side found no credible set at all
        assert row["status"] in {"no_cs_xqtl", "no_cs_gwas", "no_cs_pair"}
        return
    df = pd.read_csv(out, sep="\t")
    assert row["status"] == "ok"
    assert df["PP.H4.abf"].max() < 0.7, \
        f"distinct causal variants should not colocalize (PP.H4 " \
        f"{df['PP.H4.abf'].max():.3f})"


def test_done_sentinels_and_rerun_skip(coloc_fixture):
    res_dir = coloc_fixture / "coloc" / "results" / "expression"
    assert (res_dir / "EAS_GENEA_trait_shared.done").exists()
    # rerun: worker must succeed and mark tasks 'already done' in diagnostics
    proc = run_worker(coloc_fixture)
    assert proc.returncode == 0
    diag = pd.read_csv(coloc_fixture / "coloc" / "diagnostics" /
                       "expression.shard-0001.diagnostics.tsv", sep="\t")
    assert set(diag["status"]) == {"ok"}
    assert (diag["message"] == "already done").all()
