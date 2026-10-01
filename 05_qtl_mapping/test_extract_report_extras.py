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


def test_bed_pheno_id():
    """Combined-layer compound ids must resolve to per-modality BED ids."""
    f = mod35._bed_pheno_id
    # ungrouped: '{modality}__{gene_id}'
    assert f("expression__ENSG00000164308.17", "expression",
             "ENSG00000164308.17") == "ENSG00000164308.17"
    # grouped: '{modality}__{group_id}__{phenotype_id}'
    assert f("isoform_expression__ENSG00000164308.17__ENST00000437043.8",
             "isoform_expression",
             "ENSG00000164308.17") == "ENST00000437043.8"
    # bare BED ids (per-modality tables) pass through unchanged
    assert f("ENST00000437043.8", "isoform_expression") == "ENST00000437043.8"
    assert f("chr5:100:200:clu_1", "splicing") == "chr5:100:200:clu_1"


# ---------------------------------------------------------------------------
# QC section fixtures
# ---------------------------------------------------------------------------

def _write_bed(df, path):
    """df: features x samples -> BED with dummy coordinates."""
    bed = pd.DataFrame({"#chr": "chr1", "start": 1000, "end": 2000,
                        "phenotype_id": df.index})
    bed = pd.concat([bed, df.reset_index(drop=True)], axis=1)
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as f:
            bed.to_csv(f, sep="\t", index=False)
    else:
        bed.to_csv(path, sep="\t", index=False)


def _write_qc_fixture(tmp_path, n_feat=200, n_per_cohort=40, shift=3.0):
    """Pooled (unnorm) + final BEDs for EAS/expression with a planted cohort
    shift that the final BED removes. Returns (qtl_dir, pooled_dir)."""
    qtl_dir = tmp_path / "qtl_inputs"
    qtl_dir.mkdir(exist_ok=True)
    pooled_dir = tmp_path / "pooled"
    pooled_dir.mkdir(exist_ok=True)
    feats = [f"PHENO_{i}" for i in range(n_feat)]
    s1 = [f"COH1_S{i:03d}" for i in range(n_per_cohort)]
    s2 = [f"COH2_S{i:03d}" for i in range(n_per_cohort)]
    X1 = RNG.normal(0, 1, size=(n_feat, n_per_cohort))
    X2_shift = RNG.normal(shift, 1, size=(n_feat, n_per_cohort))
    X2_clean = RNG.normal(0, 1, size=(n_feat, n_per_cohort))
    _write_bed(pd.DataFrame(np.hstack([X1, X2_shift]),
                            index=feats, columns=s1 + s2),
               pooled_dir / "EAS_expression_pooled.bed")
    _write_bed(pd.DataFrame(np.hstack([X1, X2_clean]),
                            index=feats, columns=s1 + s2),
               qtl_dir / "EAS_expression.bed.gz")
    return qtl_dir, pooled_dir


def test_qc_phenotype_summaries(tmp_path):
    qtl_dir, pooled_dir = _write_qc_fixture(tmp_path)
    out_dir = tmp_path / "extras"
    mod35.qc_phenotype_summaries(str(qtl_dir), ["EAS"], str(out_dir),
                                 pooled_bed_dir=str(pooled_dir),
                                 n_values=500, n_features=100)
    vals = pd.read_csv(out_dir / "qc_pheno_values.tsv.gz", sep="\t")
    assert set(vals["stage"]) == {"before", "after"}
    # cohort labels parsed from namespaced sample ids
    assert set(vals["cohort"]) == {"COH1", "COH2"}
    assert (vals.groupby(["stage", "cohort"]).size() <= 500).all()
    qs = pd.read_csv(out_dir / "qc_pheno_quantiles.tsv", sep="\t")
    assert (qs.groupby(["stage", "cohort"]).size() == 99).all()
    # planted shift visible in before-stage quantiles, gone after
    med = qs[qs["quantile"] == 50].set_index(["stage", "cohort"])["value"]
    assert abs(med[("before", "COH2")] - med[("before", "COH1")]) > 2
    assert abs(med[("after", "COH2")] - med[("after", "COH1")]) < 0.5
    pca = pd.read_csv(out_dir / "qc_pheno_pca.tsv.gz", sep="\t")
    assert set(pca["stage"]) == {"before", "after"}
    assert len(pca[pca["stage"] == "before"]) == 80
    b = pca[pca["stage"] == "before"].groupby("cohort")["PC1"].mean()
    a = pca[pca["stage"] == "after"].groupby("cohort")["PC1"].mean()
    assert abs(b["COH2"] - b["COH1"]) > 2 * abs(a["COH2"] - a["COH1"])
    cnt = pd.read_csv(out_dir / "qc_pheno_feature_counts.tsv", sep="\t")
    assert set(cnt["n_features"]) == {200}
    assert set(cnt["n_samples"]) == {80}


def test_qc_phenotype_summaries_missing_before(tmp_path):
    """Without --pooled-bed-dir only the after stage is produced."""
    qtl_dir, _ = _write_qc_fixture(tmp_path)
    out_dir = tmp_path / "extras"
    mod35.qc_phenotype_summaries(str(qtl_dir), ["EAS"], str(out_dir),
                                 n_values=100, n_features=50)
    cnt = pd.read_csv(out_dir / "qc_pheno_feature_counts.tsv", sep="\t")
    assert set(cnt["stage"]) == {"after"}


def test_qc_cohort_labels_geno_space(tmp_path):
    """After-stage array ids resolve via the genotype-space map; namespaced
    before-stage ids ({label}_{sample}) resolve via prefix stripping + the
    RNAseq-space map."""
    qtl_dir = tmp_path / "qtl_inputs"
    qtl_dir.mkdir()
    pooled_dir = tmp_path / "pooled"
    pooled_dir.mkdir()
    feats = [f"PHENO_{i}" for i in range(50)]
    rna_ids = [f"J{1000 + i}" for i in range(10)]
    arr_ids = [f"010-000{i:02d}_B010-000{i:02d}" for i in range(10)]
    X = pd.DataFrame(RNG.normal(0, 1, size=(50, 10)), index=feats)
    _write_bed(X.set_axis([f"cohort2_{s}" for s in rna_ids], axis=1),
               pooled_dir / "EAS_expression_pooled.bed")
    _write_bed(X.set_axis(arr_ids, axis=1), qtl_dir / "EAS_expression.bed.gz")
    amap_path = tmp_path / "rna_map.tsv"
    pd.DataFrame({"sample_id": rna_ids, "assigned_ancestry": "EAS",
                  "cohort": "GUSTO"}).to_csv(amap_path, sep="\t", index=False)
    gmap_path = tmp_path / "geno_map.tsv"
    pd.DataFrame({"sample_id": arr_ids, "assigned_ancestry": "EAS",
                  "cohort": "GUSTO"}).to_csv(gmap_path, sep="\t", index=False)
    out_dir = tmp_path / "extras"
    mod35.qc_phenotype_summaries(str(qtl_dir), ["EAS"], str(out_dir),
                                 pooled_bed_dir=str(pooled_dir),
                                 ancestry_map_path=str(amap_path),
                                 geno_map_path=str(gmap_path),
                                 n_values=50, n_features=50)
    pca = pd.read_csv(out_dir / "qc_pheno_pca.tsv.gz", sep="\t")
    assert set(pca["stage"]) == {"before", "after"}
    assert set(pca["cohort"]) == {"GUSTO"}  # both id spaces resolved
    vals = pd.read_csv(out_dir / "qc_pheno_values.tsv.gz", sep="\t")
    assert set(vals["cohort"]) == {"GUSTO"}


def test_qc_genotype_stage_counts(tmp_path):
    d = tmp_path / "cohortqc"
    d.mkdir()
    (d / "COHX.rsq_pass.pvar").write_text("#CHROM\tPOS\n1\t100\n1\t200\n")
    (d / "COHX.rsq_pass.psam").write_text("#IID\nS1\nS2\nS3\n")
    pooled = tmp_path / "pooled_geno"
    pooled.mkdir()
    (pooled / "EAS_pooled.pvar").write_text("#CHROM\tPOS\n1\t1\n")
    (pooled / "EAS_pooled.psam").write_text("#IID\nS1\n")
    out_dir = tmp_path / "extras"
    mod35.qc_genotype_stage_counts(str(out_dir), cohort_qc_glob=str(d),
                                   geno_pooled_dir=str(pooled),
                                   ancestries=["EAS"])
    t = pd.read_csv(out_dir / "qc_genotype_stage_counts.tsv", sep="\t")
    assert set(t["stage"]) == {"imputed_rsq_pass", "pooled_qc_pass"}
    assert t.loc[t["entity"] == "COHX", "n_variants"].iloc[0] == 2
    assert t.loc[t["entity"] == "COHX", "n_samples"].iloc[0] == 3
    assert t.loc[t["entity"] == "EAS", "n_variants"].iloc[0] == 1


def test_qc_picard_metrics(tmp_path):
    d = tmp_path / "picard"
    d.mkdir()
    for c in ["COH1", "COH2"]:
        pd.DataFrame({"sample": [f"{c}_S1"], "PCT_PF_READS_ALIGNED": [0.9],
                      "SubjectBias.GC": [0.05]}).to_csv(
            d / f"{c}_qc_metrics.tsv", sep="\t", index=False)
    out_dir = tmp_path / "extras"
    mod35.qc_picard_metrics(str(d / "*_qc_metrics.tsv"), str(out_dir))
    t = pd.read_csv(out_dir / "qc_picard_metrics.tsv", sep="\t")
    assert set(t["cohort"]) == {"COH1", "COH2"}
    assert len(t) == 2


def test_qc_picard_metrics_relabel(tmp_path):
    """cohortN file labels are relabeled to the modal ancestry-map cohort."""
    d = tmp_path / "picard"
    d.mkdir()
    for c in ["cohort1", "cohort2"]:
        pd.DataFrame({"sample": [f"S{i:03d}" for i in range(10)],
                      "SubjectBias.GC": [0.05] * 10}).to_csv(
            d / f"{c}_qc_metrics.tsv", sep="\t", index=False)
    amap = tmp_path / "amap.tsv"
    pd.DataFrame({"sample_id": [f"S{i:03d}" for i in range(10)],
                  "assigned_ancestry": ["EAS"] * 10,
                  "cohort": ["GUSTO"] * 10}).to_csv(amap, sep="\t",
                                                    index=False)
    out_dir = tmp_path / "extras"
    mod35.qc_picard_metrics(str(d / "*_qc_metrics.tsv"), str(out_dir),
                            ancestry_map_path=str(amap))
    t = pd.read_csv(out_dir / "qc_picard_metrics.tsv", sep="\t")
    assert set(t["cohort"]) == {"GUSTO"}  # both files' samples are GUSTO
    assert len(t) == 20
