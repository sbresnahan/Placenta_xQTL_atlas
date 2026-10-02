#!/usr/bin/env python3
"""
test_assemble_bed_tpm.py — Tests for the tpm_from_counts units mode of
assemble_bed.py's expression subcommand.

Background: qu_correct_salmon.R writes count-proportional CPM values into the
'TPM' column of the adjusted quant.sf files (drop-in quant.sf compatibility),
so reading them with units='tpm' yields CPM-scale (read-fraction) values.
units='tpm_from_counts' instead recomputes length-normalized TPM per sample as
    (NumReads / EffectiveLength) / sum(NumReads / EffectiveLength) * 1e6
which changes within-sample ranks (per-feature length factor) and therefore
changes both the isoform-expression BED and the isoform usage ratios
(molecule fractions instead of read fractions).

Fixtures use DISTINCT effective lengths per transcript — with uniform lengths
tpm_from_counts would be proportional to raw counts and could not discriminate
the two modes.

Verified here:
  1. load_salmon(..., 'tpm_from_counts') matches hand-computed TPM;
  2. the quant.sf 'TPM' column is ignored in this mode (fixtures carry
     CPM-scale values there, mimicking qu_correct_salmon.R output);
  3. usage ratios are molecule fractions (and differ from read fractions);
  4. the isoform-expression BED carries the recomputed TPM;
  5. the min_count filter still uses raw NumReads (a short transcript with
     high recomputed TPM but low counts is excluded);
  6. edge cases: zero-count transcript, non-positive / missing effective
     length, and an all-zero sample (no NaNs).

Run:  python3 -m pytest test_assemble_bed_tpm.py -v
  or: python3 test_assemble_bed_tpm.py
"""

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assemble_bed import assemble_expression, load_salmon  # noqa: E402

# ---- Synthetic annotation (same layout as test_assemble_bed_isoform_expr) --
# g1 (+ strand, chr1, TSS at 1000): tx1, tx2
# g2 (- strand, chr2, TSS at 5000): tx3, tx4  (tx4 fails the min_count floor)
GTF = """\
chr1\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id "g1";
chr1\ttest\ttranscript\t1000\t1500\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";
chr1\ttest\ttranscript\t1200\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "tx2";
chr2\ttest\tgene\t4000\t5000\t.\t-\t.\tgene_id "g2";
chr2\ttest\ttranscript\t4500\t5000\t.\t-\t.\tgene_id "g2"; transcript_id "tx3";
chr2\ttest\ttranscript\t4000\t4200\t.\t-\t.\tgene_id "g2"; transcript_id "tx4";
"""

SAMPLES = ["s1", "s2", "s3"]

# Distinct effective lengths so length normalization is actually exercised
EFF_LENGTH = {"tx1": 900.0, "tx2": 300.0, "tx3": 900.0, "tx4": 450.0}

NUMREADS = {
    "tx1": [90.0, 180.0, 270.0],
    "tx2": [60.0, 120.0, 30.0],
    "tx3": [60.0, 60.0, 60.0],
    "tx4": [2.0, 2.0, 2.0],  # mean 2 < min_count=10 -> excluded
}


def _expected_x():
    """Length-normalized abundances x = NumReads / EffectiveLength."""
    return {t: [c / EFF_LENGTH[t] for c in NUMREADS[t]] for t in NUMREADS}


def _expected_tpm():
    x = _expected_x()
    totals = [sum(x[t][i] for t in x) for i in range(len(SAMPLES))]
    return {t: [x[t][i] / totals[i] * 1e6 for i in range(len(SAMPLES))]
            for t in x}


EXPECTED_X = _expected_x()
EXPECTED_TPM = _expected_tpm()

# Molecule fractions: x_tx / sum_{tx in gene} x_tx (per-sample scale cancels)
EXPECTED_RATIO = {
    "g1__tx1": [1.0 / 3.0, 1.0 / 3.0, 0.75],
    "g1__tx2": [2.0 / 3.0, 2.0 / 3.0, 0.25],
    "g2__tx3": [15.0 / 16.0] * 3,
    "g2__tx4": [1.0 / 16.0] * 3,  # passes frac filter but fails min_count
}
# Read fractions (count proportions) — what the CPM-scale 'tpm' mode gives;
# used as a contrast that must NOT match the tpm_from_counts ratios.
READ_FRACTION_G1_TX1 = [0.6, 0.6, 0.9]

# BED values are written with float_format='%g' (6 significant digits)
BED_TOL = dict(rel=1e-5)


def _write_tree(root: Path):
    """Write GTF + per-sample quant.sf. The 'TPM' column deliberately holds
    CPM-scale (count-proportional) values, mimicking qu_correct_salmon.R's
    adjusted quant.sf, to prove tpm_from_counts ignores that column."""
    gtf = root / "anno.gtf"
    gtf.write_text(GTF)
    salmon = root / "salmon"
    for i, s in enumerate(SAMPLES):
        d = salmon / s
        d.mkdir(parents=True)
        counts = pd.Series({t: NUMREADS[t][i] for t in NUMREADS})
        cpm = counts / counts.sum() * 1e6
        df = pd.DataFrame({
            "Name": list(NUMREADS),
            "Length": [int(EFF_LENGTH[t]) for t in NUMREADS],
            "EffectiveLength": [EFF_LENGTH[t] for t in NUMREADS],
            "TPM": [cpm[t] for t in NUMREADS],
            "NumReads": [NUMREADS[t][i] for t in NUMREADS],
        })
        df.to_csv(d / "quant.sf", sep="\t", index=False)
    return gtf, salmon


def _read_bed(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def test_load_salmon_tpm_from_counts_values(tmp_path):
    _, salmon = _write_tree(tmp_path)
    mat = load_salmon(SAMPLES, salmon, "tpm_from_counts")
    assert mat.index.tolist() == list(NUMREADS)
    assert mat.columns.tolist() == SAMPLES
    for tx in NUMREADS:
        assert mat.loc[tx, SAMPLES].tolist() == pytest.approx(
            EXPECTED_TPM[tx]), tx


def test_tpm_column_is_ignored(tmp_path):
    """The fixture's 'TPM' column holds CPM-scale values; tpm_from_counts
    must recompute from NumReads x EffectiveLength, not read that column."""
    _, salmon = _write_tree(tmp_path)
    mat = load_salmon(SAMPLES, salmon, "tpm_from_counts")
    # CPM of tx1 in s1: 90 / 212 * 1e6 ~= 424528; recomputed TPM ~= 269461
    cpm_tx1_s1 = 90.0 / 212.0 * 1e6
    assert mat.loc["tx1", "s1"] != pytest.approx(cpm_tx1_s1, **BED_TOL)
    assert mat.loc["tx1", "s1"] == pytest.approx(EXPECTED_TPM["tx1"][0])


def test_usage_ratios_are_molecule_fractions(tmp_path):
    gtf, salmon = _write_tree(tmp_path)
    bed_ratio = tmp_path / "isoforms.bed"
    assemble_expression(SAMPLES, salmon, "tpm_from_counts", gtf, bed_ratio,
                        None, min_count=10)
    ratio = _read_bed(bed_ratio).set_index("phenotype_id")
    assert set(ratio.index) == {"g1__tx1", "g1__tx2", "g2__tx3"}
    for pid, exp in EXPECTED_RATIO.items():
        if pid in ratio.index:
            assert ratio.loc[pid, SAMPLES].tolist() == pytest.approx(
                exp, **BED_TOL), pid
    # Discriminating check: molecule fraction (1/3) != read fraction (0.6)
    assert ratio.loc["g1__tx1", SAMPLES].tolist() != pytest.approx(
        READ_FRACTION_G1_TX1, **BED_TOL)

    # Contrast: the same tree read with units='tpm' (CPM column) yields the
    # read fractions — confirming the fixture discriminates the two modes.
    bed_ratio_cpm = tmp_path / "isoforms_cpm.bed"
    assemble_expression(SAMPLES, salmon, "tpm", gtf, bed_ratio_cpm,
                        None, min_count=10)
    ratio_cpm = _read_bed(bed_ratio_cpm).set_index("phenotype_id")
    assert ratio_cpm.loc["g1__tx1", SAMPLES].tolist() == pytest.approx(
        READ_FRACTION_G1_TX1, **BED_TOL)


def test_isoform_expression_bed_carries_recomputed_tpm(tmp_path):
    gtf, salmon = _write_tree(tmp_path)
    bed_ratio = tmp_path / "isoforms.bed"
    bed_expr = tmp_path / "isoform_expression.bed"
    assemble_expression(SAMPLES, salmon, "tpm_from_counts", gtf, bed_ratio,
                        None, min_count=10, bed_iso_expr=bed_expr)
    expr = _read_bed(bed_expr).set_index("phenotype_id")
    assert set(expr.index) == {"g1__tx1", "g1__tx2", "g2__tx3"}
    for tx, gene in [("tx1", "g1"), ("tx2", "g1"), ("tx3", "g2")]:
        pid = f"{gene}__{tx}"
        assert expr.loc[pid, SAMPLES].tolist() == pytest.approx(
            EXPECTED_TPM[tx], **BED_TOL), pid


def test_min_count_filter_uses_raw_counts(tmp_path):
    """tx4's recomputed TPM is high (~1.2e4, short transcript) but its mean
    NumReads is 2 < min_count — exclusion proves the floor is on raw counts."""
    assert min(EXPECTED_TPM["tx4"]) > 1000.0  # fixture sanity check
    gtf, salmon = _write_tree(tmp_path)
    bed_ratio = tmp_path / "isoforms.bed"
    bed_expr = tmp_path / "isoform_expression.bed"
    assemble_expression(SAMPLES, salmon, "tpm_from_counts", gtf, bed_ratio,
                        None, min_count=10, bed_iso_expr=bed_expr)
    assert "g2__tx4" not in _read_bed(bed_ratio)["phenotype_id"].tolist()
    assert "g2__tx4" not in _read_bed(bed_expr)["phenotype_id"].tolist()


def test_edge_cases(tmp_path):
    """Zero counts, non-positive / missing EffectiveLength, and an all-zero
    sample must yield zeros, never NaN/inf."""
    salmon = tmp_path / "salmon"
    txs = ["ta", "tb", "tc", "td"]
    count = {"ta": 100.0, "tb": 50.0, "tc": 25.0, "td": 0.0}
    eff = {"ta": 1000.0, "tb": 0.0, "tc": np.nan, "td": 500.0}
    for s in ["e1", "e2"]:
        d = salmon / s
        d.mkdir(parents=True)
        scale = 1.0 if s == "e1" else 0.0  # e2: all-zero sample
        df = pd.DataFrame({
            "Name": txs,
            "Length": [1000] * 4,
            "EffectiveLength": [eff[t] for t in txs],
            "TPM": [0.0] * 4,
            "NumReads": [count[t] * scale for t in txs],
        })
        df.to_csv(d / "quant.sf", sep="\t", index=False)

    mat = load_salmon(["e1", "e2"], salmon, "tpm_from_counts")
    # e1: only ta has defined nonzero abundance -> carries the full 1e6
    assert mat.loc["ta", "e1"] == pytest.approx(1e6)
    assert mat.loc[["tb", "tc", "td"], "e1"].tolist() == [0.0, 0.0, 0.0]
    # e2: all-zero sample -> all-zero TPM (no division by zero)
    assert mat["e2"].tolist() == [0.0, 0.0, 0.0, 0.0]
    assert not mat.isna().any().any()
    assert np.isfinite(mat.to_numpy()).all()


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as td:
        tp = Path(td)
        test_load_salmon_tpm_from_counts_values(tp)
        test_tpm_column_is_ignored(tp)
        test_usage_ratios_are_molecule_fractions(tp)
        test_isoform_expression_bed_carries_recomputed_tpm(tp)
        test_min_count_filter_uses_raw_counts(tp)
        test_edge_cases(tp)
    print("all tests passed")
