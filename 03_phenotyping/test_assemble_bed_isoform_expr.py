#!/usr/bin/env python3
"""
test_assemble_bed_isoform_expr.py — Synthetic test for the isoform-expression
(--output-isoform-expr) path added to assemble_bed.py's expression subcommand.

Builds a tiny Salmon tree (3 samples) + tiny GTF (2 genes, 4 transcripts) and
verifies:
  1. isoform_expression BED holds PRE-ratio TPM values with the min_count
     floor applied (and no min/max-frac filter);
  2. phenotype_id scheme is gene__transcript (same as the ratio BED);
  3. regression: the ratio (isoforms) BED is byte-identical whether or not
     --output-isoform-expr is requested;
  4. gene-level BED is unaffected.

Run:  python3 -m pytest test_assemble_bed_isoform_expr.py -v
  or: python3 test_assemble_bed_isoform_expr.py
"""

import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assemble_bed import assemble_expression  # noqa: E402

# ---- Synthetic annotation -------------------------------------------------
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

# TPM values (per transcript, per sample)
TPM = {
    "tx1": [10.0, 20.0, 30.0],
    "tx2": [30.0, 20.0, 10.0],
    "tx3": [5.0, 5.0, 5.0],
    "tx4": [1.0, 1.0, 1.0],
}
# NumReads: tx4 mean = 2 < min_count=10 -> excluded everywhere
NUMREADS = {
    "tx1": [100.0, 200.0, 300.0],
    "tx2": [300.0, 200.0, 100.0],
    "tx3": [20.0, 20.0, 20.0],
    "tx4": [2.0, 2.0, 2.0],
}

EXPECTED_RATIO = {  # TPM / gene-total TPM, per sample
    "g1__tx1": [0.25, 0.5, 0.75],
    "g1__tx2": [0.75, 0.5, 0.25],
    "g2__tx3": [5.0 / 6.0] * 3,
    "g2__tx4": [1.0 / 6.0] * 3,  # passes ratio filter but fails min_count
}


def _write_tree(root: Path):
    gtf = root / "anno.gtf"
    gtf.write_text(GTF)
    salmon = root / "salmon"
    for s in SAMPLES:
        d = salmon / s
        d.mkdir(parents=True)
        df = pd.DataFrame({
            "Name": list(TPM),
            "Length": [1000] * 4,
            "EffectiveLength": [1000.0] * 4,
            "TPM": [TPM[t][SAMPLES.index(s)] for t in TPM],
            "NumReads": [NUMREADS[t][SAMPLES.index(s)] for t in TPM],
        })
        df.to_csv(d / "quant.sf", sep="\t", index=False)
    return gtf, salmon


def _read_bed(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def test_isoform_expression_bed(tmp_path):
    gtf, salmon = _write_tree(tmp_path)
    bed_ratio = tmp_path / "isoforms.bed"
    bed_expr = tmp_path / "isoform_expression.bed"
    bed_gene = tmp_path / "expression.bed"

    assemble_expression(SAMPLES, salmon, "tpm", gtf, bed_ratio, bed_gene,
                        min_count=10, bed_iso_expr=bed_expr)

    expr = _read_bed(bed_expr).set_index("phenotype_id")

    # 1. tx4 excluded (min_count floor); survivors carry raw TPM values
    assert set(expr.index) == {"g1__tx1", "g1__tx2", "g2__tx3"}
    for tx, gene in [("tx1", "g1"), ("tx2", "g1"), ("tx3", "g2")]:
        pid = f"{gene}__{tx}"
        assert expr.loc[pid, SAMPLES].tolist() == pytest.approx(TPM[tx]), pid

    # 2. BED coordinates come from the gene TSS (0-based half-open)
    assert expr.loc["g1__tx1", "#chr"] == "chr1"
    assert expr.loc["g1__tx1", "start"] == 999 and expr.loc["g1__tx1", "end"] == 1000
    assert expr.loc["g2__tx3", "#chr"] == "chr2"
    assert expr.loc["g2__tx3", "start"] == 4999 and expr.loc["g2__tx3", "end"] == 5000

    # 3. ratio BED unchanged in content and correct
    ratio = _read_bed(bed_ratio).set_index("phenotype_id")
    assert set(ratio.index) == {"g1__tx1", "g1__tx2", "g2__tx3"}
    for pid, exp in EXPECTED_RATIO.items():
        if pid in ratio.index:
            assert ratio.loc[pid, SAMPLES].tolist() == pytest.approx(exp), pid

    # 4. gene-level BED = per-gene TPM sums
    gene = _read_bed(bed_gene).set_index("phenotype_id")
    assert gene.loc["g1", SAMPLES].tolist() == pytest.approx([40.0, 40.0, 40.0])
    assert gene.loc["g2", SAMPLES].tolist() == pytest.approx([6.0, 6.0, 6.0])


def test_ratio_bed_byte_identical_regression(tmp_path):
    """Requesting --output-isoform-expr must not change the ratio BED."""
    gtf, salmon = _write_tree(tmp_path)

    bed_a = tmp_path / "ratio_without.bed"
    assemble_expression(SAMPLES, salmon, "tpm", gtf, bed_a, None,
                        min_count=10)

    bed_b = tmp_path / "ratio_with.bed"
    bed_expr = tmp_path / "expr.bed"
    assemble_expression(SAMPLES, salmon, "tpm", gtf, bed_b, None,
                        min_count=10, bed_iso_expr=bed_expr)

    assert bed_a.read_bytes() == bed_b.read_bytes()


def test_validation_requires_some_output(tmp_path):
    gtf, salmon = _write_tree(tmp_path)
    with pytest.raises(ValueError, match="at least one of"):
        assemble_expression(SAMPLES, salmon, "tpm", gtf, None, None)


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as td:
        tp = Path(td)
        test_isoform_expression_bed(tp)
        test_ratio_bed_byte_identical_regression(tp)
        test_validation_requires_some_output(tp)
    print("all tests passed")
