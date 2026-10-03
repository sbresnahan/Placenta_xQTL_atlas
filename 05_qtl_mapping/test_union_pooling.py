"""Tests for union pooling (devBrain-style) across cohorts.

Covers:
  - pool_expression_within_ancestry.py / pool_modalities_within_ancestry.py:
    union keeps cohort-specific features with NaN fill; metadata come from a
    cohort that has the feature; intersection mode reproduces legacy behavior;
    pooling statistics are reported.
  - assemble_bed.py --skip-feature-filter: per-cohort prefilters (mean-count
    floor, usage [min_frac, max_frac], >50%-missing) are relaxed to
    drop-only-all-zero/all-NaN; per-value floors are untouched.
"""
import os
import sys

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "03_phenotyping"))

import pool_expression_within_ancestry as pe
import pool_modalities_within_ancestry as pm
import assemble_bed as ab


# ---------------------------------------------------------------------------
# BED / ancestry-map fixtures
# ---------------------------------------------------------------------------

def make_bed(path, features, samples, start=100):
    rows = []
    for i, f in enumerate(features):
        row = {"#chr": "chr1", "start": start + i * 100, "end": start + i * 100 + 1,
               "phenotype_id": f}
        for j, s in enumerate(samples):
            row[s] = float(i * 10 + j + 1)
        rows.append(row)
    df = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return df


@pytest.fixture
def two_cohort_tree(tmp_path):
    """cohortA has features g1,g2,g3 (samples s1,s2); cohortB has g2,g3,g4
    (samples s3,s4). g1 is cohortA-specific, g4 cohortB-specific."""
    dirs = {}
    for cohort, feats, samples in [("cohortA", ["g1", "g2", "g3"], ["s1", "s2"]),
                                   ("cohortB", ["g2", "g3", "g4"], ["s3", "s4"])]:
        cdir = tmp_path / cohort
        make_bed(cdir / "output" / "unnorm" / "expression.bed", feats, samples)
        make_bed(cdir / "output" / "unnorm" / "alt_TSS.bed", feats, samples)
        dirs[cohort] = str(cdir)
    amap = pd.DataFrame({
        "sample_id": ["s1", "s2", "s3", "s4"],
        "assigned_ancestry": ["EUR"] * 4,
        "cohort": ["cohortA", "cohortA", "cohortB", "cohortB"],
    })
    return dirs, amap


# ---------------------------------------------------------------------------
# Expression pooling (19 stage 3 path)
# ---------------------------------------------------------------------------

def test_expression_union_keeps_cohort_specific(two_cohort_tree, tmp_path):
    dirs, amap = two_cohort_tree
    (tmp_path / "out").mkdir()
    stats = pe.pool_ancestry_stratum(
        "EUR", amap, dirs, "expression", str(tmp_path / "out"), pool_mode="union")
    pooled = pd.read_csv(stats["output"], sep="\t")
    assert set(pooled["phenotype_id"]) == {"g1", "g2", "g3", "g4"}
    # NaN fill where a cohort lacks the feature
    g1 = pooled[pooled["phenotype_id"] == "g1"].iloc[0]
    assert g1[["s1", "s2"]].notna().all()
    assert g1[["s3", "s4"]].isna().all()
    g4 = pooled[pooled["phenotype_id"] == "g4"].iloc[0]
    assert g4[["s1", "s2"]].isna().all()
    assert g4[["s3", "s4"]].notna().all()
    # metadata for the cohortB-specific feature comes from cohortB
    assert g4["start"] == 100 + 2 * 100  # g4 is index 2 in cohortB's BED
    # stats
    assert stats["n_genes"] == 4
    assert stats["n_genes_intersection"] == 2
    assert stats["n_genes_union"] == 4
    assert stats["n_gained_vs_intersection"] == 2
    assert stats["n_na_cells_union"] == 4  # 2 features x 2 samples
    assert stats["pool_mode"] == "union"


def test_expression_intersection_legacy(two_cohort_tree, tmp_path):
    dirs, amap = two_cohort_tree
    (tmp_path / "out").mkdir()
    stats = pe.pool_ancestry_stratum(
        "EUR", amap, dirs, "expression", str(tmp_path / "out"),
        pool_mode="intersection")
    pooled = pd.read_csv(stats["output"], sep="\t")
    assert set(pooled["phenotype_id"]) == {"g2", "g3"}
    assert not pooled[["s1", "s2", "s3", "s4"]].isna().any().any()


def test_modality_union_and_report(two_cohort_tree, tmp_path):
    dirs, amap = two_cohort_tree
    (tmp_path / "outm").mkdir()
    stats = pm.pool_intersection_stratum(
        "EUR", amap, dirs, "alt_TSS", str(tmp_path / "outm"), pool_mode="union")
    pooled = pd.read_csv(stats["output"], sep="\t")
    assert set(pooled["phenotype_id"]) == {"g1", "g2", "g3", "g4"}
    assert stats["pooling"] == "union"
    assert stats["n_gained_vs_intersection"] == 2
    assert stats["n_na_cells_union"] == 4


# ---------------------------------------------------------------------------
# assemble_bed.py --skip-feature-filter
# ---------------------------------------------------------------------------

GTF = """chr1\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id "gA";
chr1\ttest\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "gA"; transcript_id "tA1";
chr1\ttest\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "gA"; transcript_id "tA2";
chr1\ttest\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "gA"; transcript_id "tA3";
chr1\ttest\tgene\t5000\t6000\t.\t-\t.\tgene_id "gB";
chr1\ttest\ttranscript\t5000\t6000\t.\t-\t.\tgene_id "gB"; transcript_id "tB1";
"""


def _write_quant(d, counts):
    """counts: dict transcript -> NumReads. TPM proportional to counts."""
    d.mkdir(parents=True, exist_ok=True)
    total = float(sum(counts.values()))
    rows = []
    for name, n in counts.items():
        rows.append({"Name": name, "Length": 1000, "EffectiveLength": 900,
                     "TPM": 1e6 * n / total if total else 0.0, "NumReads": n})
    pd.DataFrame(rows).to_csv(d / "quant.sf", sep="\t", index=False)


@pytest.fixture
def salmon_tree(tmp_path):
    """Two samples. tA2 is low-count (mean 5 < 10 floor) but non-zero;
    tB1 is all-zero in both samples."""
    gtf = tmp_path / "ref.gtf"
    gtf.write_text(GTF)
    sdir = tmp_path / "salmon"
    _write_quant(sdir / "s1", {"tA1": 600.0, "tA2": 5.0, "tA3": 395.0, "tB1": 0.0})
    _write_quant(sdir / "s2", {"tA1": 600.0, "tA2": 5.0, "tA3": 395.0, "tB1": 0.0})
    return str(sdir), str(gtf)


def test_expression_prefilter_default_drops_low_count(salmon_tree, tmp_path):
    sdir, gtf = salmon_tree
    bed_iso = tmp_path / "iso.bed"
    ab.assemble_expression(["s1", "s2"], Path(sdir), "tpm", Path(gtf),
                           str(bed_iso), None, skip_feature_filter=False)
    df = pd.read_csv(bed_iso, sep="\t")
    assert all("tA2" not in p for p in df["phenotype_id"])  # count floor + usage
    assert any("tA1" in p for p in df["phenotype_id"])
    assert any("tA3" in p for p in df["phenotype_id"])


def test_expression_skip_filter_keeps_low_count_drops_allzero(salmon_tree, tmp_path):
    sdir, gtf = salmon_tree
    bed_iso = tmp_path / "iso.bed"
    ab.assemble_expression(["s1", "s2"], Path(sdir), "tpm", Path(gtf),
                           str(bed_iso), None, skip_feature_filter=True)
    df = pd.read_csv(bed_iso, sep="\t")
    # low-count (but non-zero) isoform survives under union policy
    assert any("tA2" in p for p in df["phenotype_id"])
    assert any("tA1" in p for p in df["phenotype_id"])
    assert any("tA3" in p for p in df["phenotype_id"])
    # all-zero isoform is still dropped (minimal floor)
    assert all("tB1" not in p for p in df["phenotype_id"])


def _write_featurecounts(path, counts):
    """counts: dict gene -> count."""
    lines = ["# Program:featureCounts v2.0",
             "Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample.bam"]
    for g, c in counts.items():
        lines.append(f"{g}\tchr1\t1000\t2000\t+\t1000\t{c}")
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def stab_tree(tmp_path):
    """4 samples. gA: quantified everywhere. gB: exonic < 10 in 3/4 samples
    (75% NaN ratio — legacy >50%-missing filter drops it). gC: below the
    per-value floor everywhere (all-NaN ratio — dropped in both modes)."""
    gtf = tmp_path / "ref.gtf"
    gtf.write_text(
        'chr1\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id "gA";\n'
        'chr1\ttest\tgene\t3000\t4000\t.\t+\t.\tgene_id "gB";\n'
        'chr1\ttest\tgene\t5000\t6000\t.\t+\t.\tgene_id "gC";\n')
    sdir = tmp_path / "stab"
    sdir.mkdir()
    ex = {"s1": (100, 5, 5), "s2": (100, 5, 5), "s3": (100, 5, 5), "s4": (100, 50, 5)}
    inr = {"s1": (50, 100, 5), "s2": (50, 100, 5), "s3": (50, 100, 5), "s4": (50, 100, 5)}
    for s in ("s1", "s2", "s3", "s4"):
        e = ex[s]
        _write_featurecounts(sdir / f"{s}.exonic.counts.txt",
                             {"gA": e[0], "gB": e[1], "gC": e[2]})
        i = inr[s]
        _write_featurecounts(sdir / f"{s}.intronic.counts.txt",
                             {"gA": i[0], "gB": i[1], "gC": i[2]})
    return str(sdir), str(gtf)


def test_stability_legacy_drops_mostly_missing(stab_tree, tmp_path):
    sdir, gtf = stab_tree
    bed = tmp_path / "stab.bed"
    ab.assemble_stability(["s1", "s2", "s3", "s4"], Path(sdir), Path(gtf), str(bed),
                          skip_feature_filter=False)
    df = pd.read_csv(bed, sep="\t")
    assert "gA" in set(df["phenotype_id"])
    assert "gB" not in set(df["phenotype_id"])  # 75% NaN > 50% legacy cutoff
    assert "gC" not in set(df["phenotype_id"])


def test_stability_skip_filter_keeps_mostly_missing(stab_tree, tmp_path):
    sdir, gtf = stab_tree
    bed = tmp_path / "stab.bed"
    ab.assemble_stability(["s1", "s2", "s3", "s4"], Path(sdir), Path(gtf), str(bed),
                          skip_feature_filter=True)
    df = pd.read_csv(bed, sep="\t")
    assert "gA" in set(df["phenotype_id"])
    assert "gB" in set(df["phenotype_id"])   # kept: has 1 non-NaN value
    assert "gC" not in set(df["phenotype_id"])  # dropped: all-NaN
    # per-VALUE floor is untouched: gB's low-count cells are NaN
    gb = df[df["phenotype_id"] == "gB"].iloc[0]
    assert gb[["s1", "s2", "s3"]].isna().all()
    assert gb["s4"] == pytest.approx(0.5)  # exon 50 / intron 100


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
