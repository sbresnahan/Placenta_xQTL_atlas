"""End-to-end tests for the Objective 1.5 fine-mapping modules
(31_sushie_finemap.py, 33_aggregate_finemap.py,
34_pip_annotation_enrichment.py) on a synthetic two-ancestry fixture with
planted causal variants.

Fixture: ANC1 (n=100) + ANC2 (n=80), 300 variants on chr1 in LD blocks of
5 (founder dosage + 10% allele flips). Variant idx 150 (pos 1,300,000) is
causal in both ancestries (beta 1.2); idx 225 (pos 1,450,000) is causal in
ANC1 only (beta 1.0); noise sd 0.5; one covariate with effect 0.3.

Requires plink2, bgzip, and the sushie CLI (pip install sushie); the
report-render test additionally requires Rscript with rmarkdown/ggprism.
All are skipped gracefully when absent.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pyarrow")

HERE = Path(__file__).resolve().parent
MOD31 = HERE / "31_sushie_finemap.py"
MOD33 = HERE / "33_aggregate_finemap.py"
MOD34 = HERE / "34_pip_annotation_enrichment.py"
REPORT = HERE.parent / "reports" / "report_finemap.Rmd"

N_VAR = 300
POS0 = 1_000_000
STEP = 2000
BLOCK = 5
FLIP = 0.10
CAUSAL_SHARED = 150          # pos 1,300,000
CAUSAL_ANC1 = 225            # pos 1,450,000
POS_SHARED = POS0 + STEP * CAUSAL_SHARED
POS_ANC1 = POS0 + STEP * CAUSAL_ANC1
BETA_SHARED, BETA_ANC1 = 1.2, 1.0
ANC = {"ANC1": 100, "ANC2": 80}

SUSHIE_BIN = shutil.which("sushie") or str(
    Path(sys.executable).parent / "sushie")
RSCRIPT = shutil.which("Rscript") or "/opt/conda/bin/Rscript"

needs_tools = pytest.mark.skipif(
    not (shutil.which("plink2") and shutil.which("bgzip")
         and Path(SUSHIE_BIN).exists()),
    reason="requires plink2, bgzip, and sushie CLI")


# ---------------------------------------------------------------------------
# fixture builder
# ---------------------------------------------------------------------------

def _write_plink1(prefix: Path, dosage, positions, samples):
    n_var, n_samp = dosage.shape
    with open(prefix.with_suffix(".bim"), "w") as f:
        for p in positions:
            f.write(f"1\t1:{p}:A:G\t0\t{p}\tA\tG\n")
    with open(prefix.with_suffix(".fam"), "w") as f:
        for s in samples:
            f.write(f"{s}\t{s}\t0\t0\t0\t-9\n")
    code = {0: 0b00, 1: 0b10, 2: 0b11}
    with open(prefix.with_suffix(".bed"), "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))
        for i in range(n_var):
            packed = bytearray()
            for j in range(0, n_samp, 4):
                byte = 0
                for k, s in enumerate(range(j, min(j + 4, n_samp))):
                    byte |= code[int(dosage[i, s])] << (2 * k)
                packed.append(byte)
            f.write(bytes(packed))


def _sim_dosage(n_samp, seed):
    r = np.random.default_rng(seed)
    dosage = np.zeros((N_VAR, n_samp), dtype=np.int8)
    for b0 in range(0, N_VAR, BLOCK):
        founder = r.binomial(2, r.uniform(0.05, 0.5), size=n_samp)
        for i in range(b0, min(b0 + BLOCK, N_VAR)):
            flip = r.random(n_samp) < FLIP
            d = founder.copy()
            d[flip & (founder == 0)] = 2
            d[flip & (founder == 2)] = 0
            dosage[i] = d
    return dosage


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    root = tmp_path_factory.mktemp("finemap_fx")
    (root / "geno").mkdir()
    (root / "pheno").mkdir()
    (root / "results").mkdir()
    positions = np.arange(POS0, POS0 + STEP * N_VAR, STEP)
    var_ids = np.array([f"1:{p}:A:G" for p in positions])

    for ai, (anc, n) in enumerate(ANC.items()):
        samples = [f"{anc}_S{j:03d}" for j in range(n)]
        dosage = _sim_dosage(n, 42 + ai)
        prefix = root / "geno" / f"{anc}_qtl"
        _write_plink1(prefix, dosage, positions, samples)
        subprocess.run(["plink2", "--bfile", str(prefix), "--make-pgen",
                        "--out", str(prefix), "--silent"], check=True)
        cov1 = np.random.default_rng(7 + ai).normal(size=n)
        y = (BETA_SHARED * dosage[CAUSAL_SHARED]
             + (BETA_ANC1 * dosage[CAUSAL_ANC1] if anc == "ANC1" else 0)
             + 0.3 * cov1
             + np.random.default_rng(11 + ai).normal(scale=0.5, size=n))
        cov = pd.DataFrame({"covariate_id": ["COV1"],
                            **{s: [c] for s, c in zip(samples, cov1)}})
        cov.to_csv(root / "pheno" / f"{anc}_covariates_expression.tsv",
                   sep="\t", index=False)
        bed = pd.DataFrame(
            [["1", POS_SHARED, POS_SHARED + 100, "GENEX"]
             + list(np.round(y, 6))],
            columns=["#chr", "start", "end", "phenotype_id"] + samples)
        bed.to_csv(root / "pheno" / f"{anc}_expression.bed",
                   sep="\t", index=False)
        subprocess.run(["bgzip", "-f",
                        str(root / "pheno" / f"{anc}_expression.bed")],
                       check=True)
        pd.DataFrame({"phenotype_id": "GENEX", "variant_id": var_ids,
                      "slope": 0.0, "slope_se": 1.0, "tstat": 0.0,
                      "pval_nominal": 0.5}).to_parquet(
            root / "results" / f"{anc}_expression_cisqtl.parquet")
        pd.DataFrame({"phenotype_id": ["GENEX"],
                      "variant_id": [var_ids[CAUSAL_SHARED]],
                      "pval_nominal": 1e-12, "pval_beta": 1e-10,
                      "qval": 1e-6}).to_csv(
            root / "results" / f"{anc}_expression_cisqtl_top.tsv",
            sep="\t", index=False)
    pd.DataFrame({"phenotype_id": ["GENEX", "GENEX"],
                  "variant_id": [var_ids[CAUSAL_SHARED], var_ids[CAUSAL_ANC1]],
                  "rank": [1, 2]}).to_csv(
        root / "results" / "ANC1_expression_cisqtl_independent_top.tsv",
        sep="\t", index=False)
    pd.DataFrame({"phenotype_id": ["GENEX"],
                  "variant_id": [var_ids[CAUSAL_SHARED]],
                  "rank": [1]}).to_csv(
        root / "results" / "ANC2_expression_cisqtl_independent_top.tsv",
        sep="\t", index=False)
    return root


def _run_module(mod, *cli):
    res = subprocess.run([sys.executable, str(mod), *cli],
                         capture_output=True, text=True)
    assert res.returncode == 0, f"{mod.name} failed:\n{res.stderr[-2000:]}"
    return res


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

@needs_tools
def test_prepare_loci(fx):
    out = fx / "loci.expression.tsv"
    _run_module(MOD31, "prepare-loci", "--results-dir", str(fx / "results"),
                "--qtl-dir", str(fx / "pheno"),
                "--ancestries", "ANC1", "ANC2", "--modality", "expression",
                "--out", str(out))
    loci = pd.read_csv(out, sep="\t")
    assert len(loci) == 1
    row = loci.iloc[0]
    assert row["chrom"] == 1
    # windows are BED start/end +/- cis-window (1-based start), not clipped
    # to the tested-variant range
    assert row["start"] == POS_SHARED + 1 - 1_000_000
    assert row["end"] == POS_SHARED + 100 + 1_000_000
    assert row["n_signals"] == 2          # max rank across ancestries
    assert row["L"] == 5                  # max(5, 2 + 2)
    assert row["sig_in"] == "ANC1,ANC2"


@needs_tools
def test_run_recovers_causal(fx):
    loci = fx / "loci.expression.tsv"
    if not loci.exists():
        test_prepare_loci(fx)
    _run_module(MOD31, "run", "--loci-file", str(loci),
                "--qtl-dir", str(fx / "pheno"),
                "--pgen-template", str(fx / "geno" / "{anc}_qtl"),
                "--out-dir", str(fx / "finemap"),
                "--ancestries", "ANC1", "ANC2", "--sushie-bin", SUSHIE_BIN)
    out = fx / "finemap" / "expression"
    w = pd.read_csv(out / "GENEX.sushie.weights.tsv", sep="\t")
    assert w.loc[w["pos"] == POS_SHARED, "sushie_pip_all"].item() >= 0.9
    assert w.loc[w["pos"] == POS_ANC1, "sushie_pip_all"].item() >= 0.9
    assert (out / "GENEX.ancestries").read_text().strip() == "ANC1,ANC2"
    assert (out / "GENEX.done").exists()
    diag = pd.read_csv(fx / "finemap" / "logs"
                       / "expression_loci.expression.diagnostics.tsv", sep="\t")
    rec = diag.iloc[0]
    assert rec["status"] == "ok" and rec["converged"] == True  # noqa: E712
    assert rec["n_snps"] == N_VAR and rec["n_cs"] == 2
    assert rec["n_ancestries"] == 2


@needs_tools
def test_aggregate(fx):
    finemap = fx / "finemap"
    if not (finemap / "expression" / "GENEX.done").exists():
        test_run_recovers_causal(fx)
    # degenerate empty-CS locus + failed-locus diagnostics row
    (finemap / "expression" / "GENEY.sushie.cs.tsv").write_text(
        "trait\nGENEY\n")
    w = pd.read_csv(finemap / "expression" / "GENEX.sushie.weights.tsv",
                    sep="\t")
    w["trait"] = "GENEY"
    w["sushie_pip_all"] = 0.003
    w["sushie_cs_index"] = "No CS"
    w.to_csv(finemap / "expression" / "GENEY.sushie.weights.tsv",
             sep="\t", index=False)
    (finemap / "expression" / "GENEY.ancestries").write_text("ANC1,ANC2\n")
    dpath = finemap / "logs" / "expression_loci.expression.diagnostics.tsv"
    d = pd.read_csv(dpath, sep="\t")
    fail = d.iloc[0].copy()
    fail["phenotype_id"], fail["status"] = "GENEFAIL", "failed"
    degen = d.iloc[0].copy()
    degen["phenotype_id"], degen["n_cs"] = "GENEY", 0
    pd.concat([d, pd.DataFrame([fail, degen])]).to_csv(
        dpath, sep="\t", index=False)

    _run_module(MOD33, "--finemap-dir", str(finemap),
                "--ancestries", "ANC1", "ANC2")
    agg = finemap / "aggregated"
    pips = pd.read_csv(agg / "finemap_pips.tsv.gz", sep="\t")
    assert {"ANC1_effect_weight", "ANC2_effect_weight"} <= set(pips.columns)
    top = pips[pips["phenotype_id"] == "GENEX"].nlargest(1, "pip_all")
    assert top["pos"].item() in (POS_SHARED, POS_ANC1)
    cs = pd.read_csv(agg / "finemap_credible_sets.tsv.gz", sep="\t")
    assert len(cs) == 2
    assert set(cs["lead_pos"]) == {POS_SHARED, POS_ANC1}
    assert "ANC1_ANC2_est_corr" in cs.columns
    ls = pd.read_csv(agg / "finemap_locus_summary.tsv", sep="\t")
    assert len(ls) == 3
    assert ls.loc[ls["phenotype_id"] == "GENEFAIL", "max_pip"].isna().all()
    assert ls.loc[ls["phenotype_id"] == "GENEY", "n_cs_out"].isna().all()


@needs_tools
def test_enrichment(fx):
    agg = fx / "finemap" / "aggregated"
    if not (agg / "finemap_pips.tsv.gz").exists():
        test_aggregate(fx)
    pips = pd.read_csv(agg / "finemap_pips.tsv.gz", sep="\t")
    snps = sorted(pips["snp"].unique())
    # fastVEP --output-format tab: 17 columns, VEP-default-compatible layout
    rows = ["## fastVEP output",
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature"
            "\tFeature_type\tConsequence\tcDNA_position\tCDS_position"
            "\tProtein_position\tAmino_acids\tCodons\tExisting_variation"
            "\tIMPACT\tDISTANCE\tSTRAND\tFLAGS"]
    for i, s in enumerate(snps):
        chrom, pos, ref, alt = s.split(":")
        if pos == str(POS_SHARED):
            cons = "missense_variant"
        elif pos == str(POS_ANC1):
            cons = "splice_donor_variant"
        elif i % 10 == 0:
            cons = "intron_variant"
        else:
            cons = "upstream_gene_variant"
        rows.append(f"{s}\t{chrom}:{pos}\t{alt}\tG\tT\tTranscript\t{cons}"
                    "\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t1\t-")
    fv = fx / "fastvep_fixture.txt"
    fv.write_text("\n".join(rows) + "\n")
    ccre = fx / "ccre_fixture.bed"
    ccre.write_text(f"chr1\t{POS_SHARED - 1000}\t{POS_SHARED + 1000}\n"
                    f"chr1\t{POS0}\t{POS0 + 5000}\n")
    ocr = fx / "ocr_fixture.bed"
    ocr.write_text(f"chr1\t{POS_ANC1 - 1000}\t{POS_ANC1 + 1000}\n")

    out = fx / "enrichment.tsv"
    _run_module(MOD34, "--pips", str(agg / "finemap_pips.tsv.gz"),
                "--fastvep", str(fv), "--ccre", f"{ccre}:PLS",
                "--placenta-ocr", str(ocr), "--out", str(out))
    res = pd.read_csv(out, sep="\t")
    by_anno = res.set_index("annotation")
    # planted: causal variants carry these annotations -> enriched
    assert by_anno.loc["consequence_coding_nonsynonymous", "n_high_pip"] == 1
    assert by_anno.loc["consequence_splice_site", "n_high_pip"] == 1
    assert by_anno.loc["placenta_ocr", "n_high_pip"] == 1
    assert by_anno.loc["ccre_PLS", "n_high_pip"] == 1
    assert by_anno.loc["placenta_ocr", "odds_ratio"] > 1
    assert "fisher_fdr" in res.columns
    # fastVEP input generator: minimal VCF, ID = original snp string
    fv_in = fx / "fastvep_input.vcf"
    _run_module(MOD34, "--make-fastvep-input",
                "--pips", str(agg / "finemap_pips.tsv.gz"),
                "--out", str(fv_in))
    lines = fv_in.read_text().splitlines()
    assert lines[0] == "##fileformat=VCFv4.2"
    assert any(l.startswith("##contig=<ID=") for l in lines)
    hdr = lines.index("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    first = lines[hdr + 1].split("\t")
    assert len(first) == 8 and first[2].count(":") == 3
    assert not first[0].startswith("chr")


@needs_tools
def test_report_renders(fx):
    agg = fx / "finemap" / "aggregated"
    enr = fx / "enrichment.tsv"
    if not enr.exists():
        test_enrichment(fx)
    if not Path(RSCRIPT).exists():
        pytest.skip("Rscript not available")
    fig_dir = fx / "figs"
    res = subprocess.run(
        [RSCRIPT, "-e",
         f"rmarkdown::render('{REPORT}', "
         f"params=list(agg_dir='{agg}', enrichment_path='{enr}', "
         f"fig_dir='{fig_dir}', table_dir='{fx}/tabs'))"],
        capture_output=True, text=True, cwd=REPORT.parent)
    assert res.returncode == 0, res.stderr[-2000:]
    html = REPORT.parent / "report_finemap.html"
    assert html.exists() and html.stat().st_size > 100_000
    for fig in ("fig_fm_maxpip.png", "fig_fm_cs_size.png", "fig_fm_rho.png",
                "fig_fm_locus_tracks.png", "fig_fm_enrichment.png"):
        assert (fig_dir / fig).exists(), f"missing {fig}"
