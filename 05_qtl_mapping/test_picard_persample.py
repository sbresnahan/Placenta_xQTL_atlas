"""Tests for the sample-parallel Picard QC rework (2026-10-04).

Covers:
  - picard_qc.py --no-impute: single-/few-sample runs leave failed-tool cells
    as NaN instead of filling them with a degenerate (1-sample) median
  - picard_persample.scan: done-set = legacy chunk rows UNION valid
    per-sample outputs; corrupt/partial per-sample files are not done
  - picard_persample.merge: chunk rows kept verbatim (chunk wins on
    conflict), new rows appended, NaN in new rows imputed with cohort-level
    medians — including against the real cohort4.chunk5.qcmetrics.tsv fixture
  - per-sample parse-only output equals the same sample's row from a
    multi-sample (chunk-style) parse-only run
"""
import os
import subprocess
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import picard_persample as pp

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PICARD_QC = os.path.join(SCRIPT_DIR, "picard_qc.py")
PERSAMPLE = os.path.join(SCRIPT_DIR, "picard_persample.py")
REAL_CHUNK = "/mnt/user-uploads/cohort4.chunk5.qcmetrics.tsv"

# ---------------------------------------------------------------------------
# Synthetic Picard metrics fixtures (same layout as test_picard_qc_fullpanel)
# ---------------------------------------------------------------------------

ALIGNMENT = """## METRICS CLASS\tpicard.analysis.AlignmentSummaryMetrics
CATEGORY\tTOTAL_READS\tPF_READS\tPF_READS_ALIGNED\tPCT_PF_READS_ALIGNED\tMEAN_READ_LENGTH\tSTRAND_BALANCE\tSAMPLE\tLIBRARY\tREAD_GROUP
FIRST_OF_PAIR\t1000\t1000\t900\t0.9\t150\t0.51\t\t\t
SECOND_OF_PAIR\t1000\t1000\t880\t0.88\t150\t0.49\t\t\t
PAIR\t2000\t2000\t1780\t0.89\t150\t0.5\t\t\t
"""

INSERT_SIZE_A = """## METRICS CLASS\tpicard.analysis.InsertSizeMetrics
MEDIAN_INSERT_SIZE\tMEAN_INSERT_SIZE\tSTANDARD_DEVIATION\tMEDIAN_ABSOLUTE_DEVIATION\tMIN_INSERT_SIZE\tMAX_INSERT_SIZE\tPAIR_ORIENTATION\tREADS_USED\tSAMPLE
299\t301.5\t55.2\t42\t50\t800\tFR\tALL\t
"""

INSERT_SIZE_B = """## METRICS CLASS\tpicard.analysis.InsertSizeMetrics
MEDIAN_INSERT_SIZE\tMEAN_INSERT_SIZE\tSTANDARD_DEVIATION\tMEDIAN_ABSOLUTE_DEVIATION\tMIN_INSERT_SIZE\tMAX_INSERT_SIZE\tPAIR_ORIENTATION\tREADS_USED\tSAMPLE
350\t355.0\t60.1\t45\t60\t900\tFR\tALL\t
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
"""


def _write_raw_sample(raw_root, sample, insert_size=INSERT_SIZE_A,
                      with_dup=True):
    d = os.path.join(raw_root, sample)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "alignment_metrics.txt"), "w") as f:
        f.write(ALIGNMENT)
    with open(os.path.join(d, "insert_size_metrics.txt"), "w") as f:
        f.write(insert_size)
    with open(os.path.join(d, "rna_metrics.txt"), "w") as f:
        f.write(RNA_METRICS)
    with open(os.path.join(d, "gc_bias_summary.txt"), "w") as f:
        f.write(GC_SUMMARY)
    if with_dup:
        with open(os.path.join(d, "dup_metrics.txt"), "w") as f:
            f.write(DUP_METRICS)
    return d


def _write_config(tmp_path):
    """Minimal config.yml satisfying config_get.py + picard_qc.main()."""
    config = tmp_path / "config.yml"
    config.write_text(
        f"output_base: {tmp_path}/output\n"
        f"reference_dir: {tmp_path}/ref\n"
        f"normalized_gtf: {tmp_path}/ref/annot.gtf\n"
        f"ref_genome: {tmp_path}/ref/genome.fa\n"
        "cohorts:\n"
        "  cohort1:\n"
        f"    samples_file: {tmp_path}/output/cohort1/samples.txt\n"
    )
    return str(config)


def _run_picard_qc_parse_only(tmp_path, samples, extra_args=()):
    """Run picard_qc.py in --parse-only mode; return the output DataFrame."""
    config = _write_config(tmp_path)
    samples_file = tmp_path / "samples.txt"
    samples_file.write_text("\n".join(samples) + "\n")
    out = tmp_path / "out.tsv"
    raw_dir = tmp_path / "raw"
    cmd = [sys.executable, PICARD_QC,
           "--config", config, "--cohort", "cohort1",
           "--samples", str(samples_file),
           "--output", str(out),
           "--refflat", str(tmp_path / "dummy.refFlat"),
           "--raw-dir", str(raw_dir),
           "--parse-only", *extra_args]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"picard_qc.py failed:\n{res.stderr}"
    return pd.read_csv(out, sep="\t", index_col=0)


def _write_chunk(qc_dir, cohort, idx, df):
    path = os.path.join(qc_dir, f"{cohort}.chunk{idx}.qcmetrics.tsv")
    df.to_csv(path, sep="\t")
    return path


def _write_per_sample(qc_dir, cohort, sample, df=None, raw_text=None):
    ps_dir = os.path.join(qc_dir, "per_sample", cohort)
    os.makedirs(ps_dir, exist_ok=True)
    path = os.path.join(ps_dir, f"{sample}.qcmetrics.tsv")
    if raw_text is not None:
        with open(path, "w") as f:
            f.write(raw_text)
    else:
        df.to_csv(path, sep="\t")
    return path


# ---------------------------------------------------------------------------
# --no-impute
# ---------------------------------------------------------------------------

def test_no_impute_leaves_nan_for_failed_tool(tmp_path):
    raw = tmp_path / "raw"
    _write_raw_sample(str(raw), "A")
    _write_raw_sample(str(raw), "B", with_dup=False)  # MarkDuplicates "failed"

    df = _run_picard_qc_parse_only(tmp_path, ["A", "B"], ["--no-impute"])
    assert pd.isna(df.loc["B", "DupMetrics.PERCENT_DUPLICATION"])
    assert df.loc["A", "DupMetrics.PERCENT_DUPLICATION"] == pytest.approx(0.11)


def test_default_impute_unchanged(tmp_path):
    """Without --no-impute the legacy median imputation still applies."""
    raw = tmp_path / "raw"
    _write_raw_sample(str(raw), "A")
    _write_raw_sample(str(raw), "B", with_dup=False)

    df = _run_picard_qc_parse_only(tmp_path, ["A", "B"])
    # 2-sample column median of [0.11, NaN] is 0.11
    assert df.loc["B", "DupMetrics.PERCENT_DUPLICATION"] == pytest.approx(0.11)


# ---------------------------------------------------------------------------
# scan: done-set and todo list
# ---------------------------------------------------------------------------

def test_scan_done_set_and_todo(tmp_path):
    qc_dir = str(tmp_path / "qc")
    os.makedirs(qc_dir)
    chunk_df = pd.DataFrame({"m1": [1.0, 2.0]}, index=["A", "B"])
    chunk_df.index.name = "sample"
    _write_chunk(qc_dir, "cohort1", 1, chunk_df)
    _write_chunk(qc_dir, "cohort1", 2,
                 pd.DataFrame({"m1": [3.0]}, index=pd.Index(["C"], name="sample")))
    _write_per_sample(qc_dir, "cohort1", "D",
                      pd.DataFrame({"m1": [4.0]}, index=pd.Index(["D"], name="sample")))
    # Corrupt per-sample outputs must NOT count as done
    _write_per_sample(qc_dir, "cohort1", "E", raw_text="sample\tm1\n")  # no data row
    _write_per_sample(qc_dir, "cohort1", "F", raw_text="")               # empty

    samples = tmp_path / "samples.txt"
    samples.write_text("A\nB\nC\nD\nE\nF\nG\nG\n")  # G duplicated

    assert pp.done_samples(qc_dir, "cohort1") == {"A", "B", "C", "D"}

    todo_file = tmp_path / "todo.txt"
    todo, n_done = pp.write_todo(str(samples), qc_dir, "cohort1", str(todo_file))
    assert todo == ["E", "F", "G"]
    assert n_done == 4
    assert todo_file.read_text().splitlines() == ["E", "F", "G"]


def test_scan_ignores_other_cohorts(tmp_path):
    qc_dir = str(tmp_path / "qc")
    os.makedirs(qc_dir)
    _write_chunk(qc_dir, "cohort2", 1,
                 pd.DataFrame({"m1": [1.0]}, index=pd.Index(["A"], name="sample")))
    samples = tmp_path / "samples.txt"
    samples.write_text("A\n")
    todo, n_done = pp.write_todo(str(samples), qc_dir, "cohort1",
                                 str(tmp_path / "todo.txt"))
    assert todo == ["A"] and n_done == 0


def test_scan_cli_stdout_is_count(tmp_path):
    """The submit driver captures scan stdout as the array size."""
    qc_dir = tmp_path / "qc"
    qc_dir.mkdir()
    _write_chunk(str(qc_dir), "cohort1", 1,
                 pd.DataFrame({"m1": [1.0]}, index=pd.Index(["A"], name="sample")))
    samples = tmp_path / "samples.txt"
    samples.write_text("A\nB\nC\n")
    res = subprocess.run(
        [sys.executable, PERSAMPLE, "scan", "--samples", str(samples),
         "--qc-dir", str(qc_dir), "--cohort", "cohort1",
         "--todo-out", str(tmp_path / "todo.txt")],
        capture_output=True, text=True)
    assert res.returncode == 0
    assert res.stdout.strip() == "2"


# ---------------------------------------------------------------------------
# merge: chunk rows verbatim, new rows imputed at cohort level
# ---------------------------------------------------------------------------

def test_merge_chunk_verbatim_new_imputed(tmp_path):
    qc_dir = str(tmp_path / "qc")
    os.makedirs(qc_dir)
    chunk = pd.DataFrame({"m1": [1.0, 3.0], "m2": [10.0, 30.0]},
                         index=pd.Index(["A", "B"], name="sample"))
    _write_chunk(qc_dir, "cohort1", 1, chunk)

    # New sample C with one missing cell
    new = pd.DataFrame({"m1": [5.0], "m2": [float("nan")]},
                       index=pd.Index(["C"], name="sample"))
    _write_per_sample(qc_dir, "cohort1", "C", new)
    # Per-sample file for A (already in chunk) must be ignored — chunk wins
    bogus = pd.DataFrame({"m1": [999.0], "m2": [999.0]},
                         index=pd.Index(["A"], name="sample"))
    _write_per_sample(qc_dir, "cohort1", "A", bogus)

    merged, stats = pp.merge_cohort(qc_dir, "cohort1")

    assert list(merged.index) == ["A", "B", "C"]
    # Chunk rows verbatim
    assert merged.loc["A", "m1"] == 1.0 and merged.loc["A", "m2"] == 10.0
    assert merged.loc["B", "m1"] == 3.0 and merged.loc["B", "m2"] == 30.0
    # New row: present value kept, NaN imputed with cohort median of
    # [10, 30, NaN] = 20.0
    assert merged.loc["C", "m1"] == 5.0
    assert merged.loc["C", "m2"] == pytest.approx(20.0)
    assert stats["n_chunk_rows"] == 2
    assert stats["n_new_rows"] == 1
    assert stats["n_chunk_wins_skipped"] == 1
    assert stats["n_imputed_cells"] == 1


def test_merge_cli_writes_cohort_file_and_missing_report(tmp_path):
    qc_dir = tmp_path / "qc"
    qc_dir.mkdir()
    chunk = pd.DataFrame({"m1": [1.0]}, index=pd.Index(["A"], name="sample"))
    _write_chunk(str(qc_dir), "cohort1", 1, chunk)
    _write_per_sample(str(qc_dir), "cohort1", "B",
                      pd.DataFrame({"m1": [2.0]}, index=pd.Index(["B"], name="sample")))
    # samples.txt lists a third sample that never ran -> missing report
    out_base = tmp_path / "output"
    (out_base / "cohort1").mkdir(parents=True)
    (out_base / "cohort1" / "samples.txt").write_text("A\nB\nC\n")

    res = subprocess.run(
        [sys.executable, PERSAMPLE, "merge", "--qc-dir", str(qc_dir),
         "--output-base", str(out_base)],
        capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    merged = pd.read_csv(qc_dir / "cohort1_qc_metrics.tsv", sep="\t", index_col=0)
    assert list(merged.index) == ["A", "B"]
    missing = (qc_dir / "cohort1.missing_samples.txt").read_text().splitlines()
    assert missing == ["C"]


@pytest.mark.skipif(not os.path.exists(REAL_CHUNK),
                    reason="real chunk fixture not attached")
def test_merge_against_real_chunk_fixture(tmp_path):
    """cohort4.chunk5.qcmetrics.tsv (10 real samples) + one synthetic new
    sample with blanked cells + one bogus duplicate of a completed sample."""
    qc_dir = str(tmp_path / "qc")
    os.makedirs(qc_dir)
    fixture = pd.read_csv(REAL_CHUNK, sep="\t", index_col=0)
    assert fixture.shape[0] == 10
    _write_chunk(qc_dir, "cohort4", 5, fixture)

    # New sample: first fixture row with two cells blanked
    new_row = fixture.iloc[[0]].copy()
    new_row.index = pd.Index(["TESTNEW1"], name="sample")
    blank_cols = ["DupMetrics.PERCENT_DUPLICATION", "InsertSize.MEAN_INSERT_SIZE"]
    new_row.loc["TESTNEW1", blank_cols] = float("nan")
    _write_per_sample(qc_dir, "cohort4", "TESTNEW1", new_row)

    # Bogus per-sample file for an already-completed sample (chunk must win)
    bogus = fixture.iloc[[0]].copy()
    bogus.iloc[0] = 999.0
    _write_per_sample(qc_dir, "cohort4", fixture.index[0], bogus)

    merged, stats = pp.merge_cohort(qc_dir, "cohort4")

    assert merged.shape == (11, fixture.shape[1])
    # Legacy rows byte-for-byte equal to the chunk file
    pd.testing.assert_frame_equal(
        merged.loc[fixture.index, fixture.columns], fixture,
        check_dtype=False)
    # New row imputed with cohort-level medians (median of the 10 fixture
    # values; the NaN in the new row does not enter the median)
    for col in blank_cols:
        assert merged.loc["TESTNEW1", col] == pytest.approx(
            fixture[col].median())
    # Non-blanked cells of the new row kept as computed
    assert merged.loc["TESTNEW1", "AlignMetrics.TOTAL_READS"] == \
        pytest.approx(fixture.iloc[0]["AlignMetrics.TOTAL_READS"])
    assert stats["n_chunk_rows"] == 10
    assert stats["n_new_rows"] == 1
    assert stats["n_chunk_wins_skipped"] == 1
    assert stats["n_imputed_cells"] == 2


# ---------------------------------------------------------------------------
# Per-sample run == same row of a multi-sample (chunk-style) run
# ---------------------------------------------------------------------------

def test_persample_row_matches_multisample_row(tmp_path):
    raw = tmp_path / "raw"
    _write_raw_sample(str(raw), "A", insert_size=INSERT_SIZE_A)
    _write_raw_sample(str(raw), "B", insert_size=INSERT_SIZE_B)

    df_ab = _run_picard_qc_parse_only(tmp_path, ["A", "B"], ["--no-impute"])
    df_a = _run_picard_qc_parse_only(tmp_path, ["A"], ["--no-impute"])

    assert df_a.shape[0] == 1
    pd.testing.assert_series_equal(df_ab.loc["A"], df_a.loc["A"],
                                   check_names=False)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
