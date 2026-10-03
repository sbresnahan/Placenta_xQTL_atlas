"""Tests for the devBrain-style full-panel Picard extraction in picard_qc.py.

Covers:
  - parse_picard_metrics_all: multi-row tables, histogram-section termination
  - extract_full_panel: PAIR-row primary (unsuffixed), category suffixes,
    ID/categorical exclusions, FIELD_RENAME backward-compatible names
  - the original 24 whitelist column names survive with library-level values
  - parse-only collector mode (no Picard execution; reads persisted raw files)
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import picard_qc as pq


# ---------------------------------------------------------------------------
# Synthetic Picard metrics fixtures (real file layout: comment lines,
# ## METRICS CLASS line, header row, data rows, then blank/## sections)
# ---------------------------------------------------------------------------

ALIGNMENT = """## picard.analysis.AlignmentSummaryMetrics
## METRICS CLASS\tpicard.analysis.AlignmentSummaryMetrics
CATEGORY\tTOTAL_READS\tPF_READS\tPF_READS_ALIGNED\tPCT_PF_READS_ALIGNED\tMEAN_READ_LENGTH\tSTRAND_BALANCE\tSAMPLE\tLIBRARY\tREAD_GROUP
FIRST_OF_PAIR\t1000\t1000\t900\t0.9\t150\t0.51\t\t\t
SECOND_OF_PAIR\t1000\t1000\t880\t0.88\t150\t0.49\t\t\t
PAIR\t2000\t2000\t1780\t0.89\t150\t0.5\t\t\t
"""

INSERT_SIZE = """## METRICS CLASS\tpicard.analysis.InsertSizeMetrics
MEDIAN_INSERT_SIZE\tMEAN_INSERT_SIZE\tSTANDARD_DEVIATION\tMEDIAN_ABSOLUTE_DEVIATION\tMIN_INSERT_SIZE\tMAX_INSERT_SIZE\tPAIR_ORIENTATION\tREADS_USED\tSAMPLE
299\t301.5\t55.2\t42\t50\t800\tFR\tALL\t
## HISTOGRAM_CLASS\tjava.lang.Integer
insert_size\tAll_Reads.fr_count
100\t5
101\t7
"""

RNA_METRICS = """## METRICS CLASS\tpicard.analysis.CollectRnaSeqMetrics$RnaSeqMetrics
PF_BASES\tPF_ALIGNED_BASES\tRIBOSOMAL_BASES\tCODING_BASES\tUTR_BASES\tINTRONIC_BASES\tINTERGENIC_BASES\tPCT_RIBOSOMAL_BASES\tPCT_CODING_BASES\tPCT_UTR_BASES\tPCT_INTRONIC_BASES\tPCT_INTERGENIC_BASES\tMEDIAN_5PRIME_BIAS\tMEDIAN_3PRIME_BIAS\tMEDIAN_5PRIME_TO_3PRIME_BIAS\tPCT_MRNA_BASES\tPCT_USABLE_BASES\tSAMPLE\tLIBRARY\tREAD_GROUP
100000\t90000\t1000\t40000\t20000\t15000\t14000\t0.011\t0.444\t0.222\t0.167\t0.156\t0.31\t0.28\t1.1\t0.666\t0.8\t\t\t
"""

GC_SUMMARY = """## METRICS CLASS\tpicard.analysis.CollectGcBiasMetrics$GcBiasSummaryMetrics
ACCUMULATION_LEVEL\tREADS_USED\tWINDOW_SIZE\tTOTAL_CLUSTERS\tALIGNED_READS\tAT_DROPOUT\tGC_DROPOUT\tMEAN_COVERAGE\tSAMPLE\tLIBRARY\tREAD_GROUP
ALL\tALL\t100\t10\t1780\t1.5\t2.5\t35.5\t\t\t
"""

DUP_METRICS = """## METRICS CLASS\tpicard.sam.DuplicationMetrics
LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE
lib1\t50000\t900000\t5000\t90000\t0.11\t5000000

## HISTOGRAM_CLASS\tjava.lang.Double
BIN\tvalue
1.0\t100.0
"""


def _write(tmp_path, name, content):
    p = tmp_path / name
    p.write_text(content)
    return str(p)


@pytest.fixture
def raw_dir(tmp_path):
    """A persisted per-sample raw-metrics directory with all five tools."""
    d = tmp_path / "S1"
    d.mkdir()
    (d / "alignment_metrics.txt").write_text(ALIGNMENT)
    (d / "insert_size_metrics.txt").write_text(INSERT_SIZE)
    (d / "rna_metrics.txt").write_text(RNA_METRICS)
    (d / "gc_bias_summary.txt").write_text(GC_SUMMARY)
    (d / "dup_metrics.txt").write_text(DUP_METRICS)
    return str(d)


# ---------------------------------------------------------------------------
# parse_picard_metrics_all
# ---------------------------------------------------------------------------

def test_parse_all_rows_multi_category(tmp_path):
    f = _write(tmp_path, "alignment_metrics.txt", ALIGNMENT)
    rows = pq.parse_picard_metrics_all(f)
    assert len(rows) == 3
    assert [r["CATEGORY"] for r in rows] == ["FIRST_OF_PAIR", "SECOND_OF_PAIR", "PAIR"]


def test_parse_stops_at_histogram(tmp_path):
    f = _write(tmp_path, "insert_size_metrics.txt", INSERT_SIZE)
    rows = pq.parse_picard_metrics_all(f)
    assert len(rows) == 1  # histogram rows must not leak in
    assert rows[0]["MEAN_INSERT_SIZE"] == "301.5"


def test_parse_missing_file(tmp_path):
    assert pq.parse_picard_metrics_all(str(tmp_path / "nope.txt")) == []


# ---------------------------------------------------------------------------
# extract_full_panel
# ---------------------------------------------------------------------------

def test_primary_row_is_pair_unsuffixed(tmp_path):
    f = _write(tmp_path, "alignment_metrics.txt", ALIGNMENT)
    out = pq.extract_full_panel(pq.parse_picard_metrics_all(f), "AlignMetrics")
    # Library-level PAIR row is primary and unsuffixed
    assert out["AlignMetrics.TOTAL_READS"] == 2000.0
    assert out["AlignMetrics.MEAN_READ_LENGTH"] == 150.0
    # Other categories suffixed
    assert out["AlignMetrics.TOTAL_READS.FIRST_OF_PAIR"] == 1000.0
    assert out["AlignMetrics.TOTAL_READS.SECOND_OF_PAIR"] == 1000.0


def test_field_rename_preserves_legacy_names(tmp_path):
    f = _write(tmp_path, "alignment_metrics.txt", ALIGNMENT)
    out = pq.extract_full_panel(pq.parse_picard_metrics_all(f), "AlignMetrics")
    # Picard >=2.27 PF_* names map back to the pipeline's historical names
    assert out["AlignMetrics.READS_ALIGNED"] == 1780.0
    assert out["AlignMetrics.PCT_READS_ALIGNED"] == pytest.approx(0.89)
    assert "AlignMetrics.PF_READS_ALIGNED" not in out
    # Rename applies to the primary row only
    assert out["AlignMetrics.PF_READS_ALIGNED.FIRST_OF_PAIR"] == 900.0


def test_exclusions_and_nonnumeric_dropped(tmp_path):
    f = _write(tmp_path, "alignment_metrics.txt", ALIGNMENT)
    out = pq.extract_full_panel(pq.parse_picard_metrics_all(f), "AlignMetrics")
    for key in out:
        assert not any(ex in key for ex in
                       ("SAMPLE", "LIBRARY", "READ_GROUP", "CATEGORY"))


def test_single_row_table_unsuffixed(tmp_path):
    f = _write(tmp_path, "rna_metrics.txt", RNA_METRICS)
    out = pq.extract_full_panel(pq.parse_picard_metrics_all(f), "RnaMetrics")
    assert out["RnaMetrics.PCT_CODING_BASES"] == pytest.approx(0.444)
    # Full panel: fields beyond the old whitelist are present
    assert out["RnaMetrics.PCT_MRNA_BASES"] == pytest.approx(0.666)
    assert out["RnaMetrics.PCT_USABLE_BASES"] == pytest.approx(0.8)
    assert out["RnaMetrics.PF_ALIGNED_BASES"] == 90000.0
    assert not any(k.endswith(".PAIR") for k in out)


def test_empty_rows():
    assert pq.extract_full_panel([], "X") == {}


# ---------------------------------------------------------------------------
# Collectors in parse-only mode (no Picard execution)
# ---------------------------------------------------------------------------

def test_collectors_parse_only_full_panel(raw_dir):
    m = {}
    m.update(pq.collect_alignment_summary("picard", None, raw_dir, "S1",
                                          parse_only=True))
    m.update(pq.collect_insert_size("picard", None, raw_dir, "S1",
                                    parse_only=True))
    m.update(pq.collect_rna_seq_metrics("picard", None, raw_dir, "S1",
                                        "unused.refFlat", parse_only=True))
    m.update(pq.collect_gc_bias("picard", None, raw_dir, "S1",
                                parse_only=True))
    m.update(pq.mark_duplicates("picard", None, raw_dir, "S1",
                                parse_only=True))

    # The original 24 whitelist column names survive (library-level values)
    legacy_24 = [
        "AlignMetrics.TOTAL_READS", "AlignMetrics.READS_ALIGNED",
        "AlignMetrics.MEAN_READ_LENGTH", "AlignMetrics.STRAND_BALANCE",
        "AlignMetrics.PCT_READS_ALIGNED",
        "InsertSize.MEAN_INSERT_SIZE", "InsertSize.STANDARD_DEVIATION",
        "InsertSize.MEDIAN_INSERT_SIZE",
        "DupMetrics.PERCENT_DUPLICATION", "DupMetrics.UNPAIRED_READ_DUPLICATES",
        "DupMetrics.READ_PAIR_DUPLICATES",
        "RnaMetrics.PCT_RIBOSOMAL_BASES", "RnaMetrics.PCT_CODING_BASES",
        "RnaMetrics.PCT_UTR_BASES", "RnaMetrics.PCT_INTRONIC_BASES",
        "RnaMetrics.PCT_INTERGENIC_BASES", "RnaMetrics.MEDIAN_5PRIME_BIAS",
        "RnaMetrics.MEDIAN_3PRIME_BIAS", "RnaMetrics.MEDIAN_5PRIME_TO_3PRIME_BIAS",
        "GcBias.AT_DROPOUT", "GcBias.GC_DROPOUT", "GcBias.MEAN_COVERAGE",
    ]
    for name in legacy_24:
        assert name in m, f"legacy column missing: {name}"
        assert m[name] is not None

    # Full panel goes beyond the whitelist
    assert "InsertSize.MEDIAN_ABSOLUTE_DEVIATION" in m
    assert "DupMetrics.ESTIMATED_LIBRARY_SIZE" in m
    assert "DupMetrics.READ_PAIRS_EXAMINED" in m
    assert "GcBias.TOTAL_CLUSTERS" in m
    assert len(m) > 30

    # Spot-check values
    assert m["DupMetrics.PERCENT_DUPLICATION"] == pytest.approx(0.11)
    assert m["GcBias.MEAN_COVERAGE"] == pytest.approx(35.5)
    assert m["InsertSize.MEDIAN_INSERT_SIZE"] == 299.0


def test_collectors_parse_only_missing_raw_dir(tmp_path):
    empty = str(tmp_path / "no_such_sample")
    assert pq.collect_alignment_summary("picard", None, empty, "S2",
                                        parse_only=True) == {}
    assert pq.mark_duplicates("picard", None, empty, "S2", parse_only=True) == {}


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
