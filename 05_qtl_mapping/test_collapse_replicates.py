#!/usr/bin/env python3
"""
test_collapse_replicates.py — Tests for collapse_replicates.py.

Builds a synthetic three-cohort fixture with:
  * A1: same-cohort concordant pair (S1+S2 in cohortA) -> collapsed
  * A3: same-cohort discordant pair (S7+S8 in cohortA) -> collapsed
        (discordance flagged in the report; never a veto)
  * A2: cross-protocol pair (S3 in cohortA, S5 in cohortB) -> dropped
        entirely (both runs, all cohorts + collapsed ancestry map)
  * cohortC: no pairs -> symlinked, untouched

Run: pytest test_collapse_replicates.py
"""

import gzip
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import collapse_replicates as cr

META_COLS = ["#chr", "start", "end", "phenotype_id"]


# ---------------------------------------------------------------------------
# Fixture construction
# ---------------------------------------------------------------------------

def _write_bed(path, pids, samples):
    df = pd.DataFrame({
        "#chr": ["chr1"] * len(pids),
        "start": range(100, 100 + len(pids)),
        "end": range(101, 101 + len(pids)),
        "phenotype_id": pids,
    })
    for s, vals in samples.items():
        assert len(vals) == len(pids)
        df[s] = vals
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, float_format="%g")


def _write_quant(path, counts, eff_lengths=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame({
        "Name": list(counts.keys()),
        "Length": [1000] * len(counts),
        "EffectiveLength": ([900.0] * len(counts) if eff_lengths is None
                            else [eff_lengths[k] for k in counts.keys()]),
        "TPM": [1.0] * len(counts),
        "NumReads": list(counts.values()),
    })
    df.to_csv(path, sep="\t", index=False)


def _write_featurecounts(path, counts):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("# Program:featureCounts v2.0.0\n")
        f.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample.bam\n")
        for g, c in counts.items():
            f.write(f"{g}\tchr1\t1\t100\t+\t100\t{c}\n")


def _write_numers(path, samples):
    """samples: dict sample -> list of counts; 5 fixed junctions."""
    rows = ["chr1:100:200:clu_1_+", "chr1:150:200:clu_1_+",
            "chr1:300:400:clu_2_-", "chr1:500:600:clu_3_+",
            "chr1:550:600:clu_3_+"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as f:
        f.write(" ".join(samples.keys()) + "\n")  # samples-only header
        for i, r in enumerate(rows):
            f.write(r + " " + " ".join(str(v[i]) for v in samples.values()) + "\n")


def _write_psi(path, samples):
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame({
        "gene_id": ["G1", "G2"], "seqid": ["chr1", "chr1"],
        "start": [100, 300], "end": [200, 400], "strand": ["+", "-"],
    })
    for s, vals in samples.items():
        df[s] = vals
    df.to_csv(path, sep="\t", index=False, compression="gzip")


@pytest.fixture(scope="module")
def fixture(tmp_path_factory):
    root = tmp_path_factory.mktemp("collapse_fixture")
    base = root / "base"
    for c in ("cohortA", "cohortB", "cohortC"):
        (base / c).mkdir(parents=True)

    # config.yml
    (root / "config.yml").write_text(
        f"output_base: {base}\n"
        "cohorts:\n"
        "  cohortA:\n"
        "    fastq_dir: x\n"
        "  cohortB:\n"
        "    fastq_dir: x\n"
        "  cohortC:\n"
        "    fastq_dir: x\n")

    # metadata + ancestry map
    (root / "metadata.tsv").write_text(
        "rnaseq_id\tarray_id\n"
        "S1\tA1\nS2\tA1\nS7\tA3\nS8\tA3\nS3\tA2\nS5\tA2\n"
        "S4\tA4\nS6\tA5\nS9\tA6\n")
    (root / "ancestry_map.tsv").write_text(
        "sample_id\tassigned_ancestry\tcohort\n"
        "S1\tEAS\tcohortA\nS3\tEAS\tcohortA\nS5\tEAS\tcohortB\n"
        "S7\tEAS\tcohortA\nS4\tEAS\tcohortA\n"
        "S6\tEUR\tcohortB\nS9\tEUR\tcohortC\n")

    # samples.txt
    (base / "cohortA" / "samples.txt").write_text("S1\nS2\nS3\nS4\nS7\nS8\n")
    (base / "cohortB" / "samples.txt").write_text("S5\nS6\n")
    (base / "cohortC" / "samples.txt").write_text("S9\n")

    # normalized GTF: tx1,tx2,tx3 -> G1; tx4 -> G2; G3 bare
    gtf = root / "ref.gtf"
    with open(gtf, "w") as f:
        f.write('chr1\ttest\tgene\t1\t1000\t.\t+\t.\tgene_id "G1";\n')
        for tx in ("tx1", "tx2", "tx3"):
            f.write(f'chr1\ttest\ttranscript\t1\t500\t.\t+\t.\t'
                    f'gene_id "G1"; transcript_id "{tx}";\n')
        f.write('chr1\ttest\tgene\t2000\t3000\t.\t-\t.\tgene_id "G2";\n')
        f.write('chr1\ttest\ttranscript\t2000\t2500\t.\t-\t.\t'
                'gene_id "G2"; transcript_id "tx4";\n')
        f.write('chr1\ttest\tgene\t4000\t5000\t.\t+\t.\tgene_id "G3";\n')

    unnormA = base / "cohortA" / "output" / "unnorm"

    # ---- expression (counts) ----
    _write_bed(unnormA / "expression.bed", ["G1", "G2", "G3", "G4"], {
        "S1": [100, 10, 500, 30], "S2": [50, 20, 450, 35],
        "S3": [200, 40, 100, 60], "S4": [300, 60, 700, 45],
        "S7": [100, 50, 10, 5], "S8": [90, 2000, 3, 500],
    })

    # ---- isoform_expression (counts) ----
    _write_bed(unnormA / "isoform_expression.bed",
               ["G1__tx1", "G1__tx2", "G2__tx4"], {
        "S1": [60, 40, 100], "S2": [30, 20, 90],
        "S3": [70, 30, 80], "S4": [80, 20, 70],
        "S7": [10, 50, 30], "S8": [50, 10, 30],
    })

    # ---- isoforms (ratios) + QU-corrected quants ----
    _write_bed(unnormA / "isoforms.bed", ["G1__tx1", "G1__tx2", "G2__tx4"], {
        "S1": [60 / 110, 40 / 110, 1.0], "S2": [0.3, 0.2, 1.0],
        "S3": [0.6, 0.4, 1.0], "S4": [0.7, 0.3, 1.0],
        "S7": [0.1, 0.5, 0.9], "S8": [0.9, 0.5, 0.1],
    })
    _write_quant(base / "cohortA" / "intermediate" / "expression_qu" / "S1" / "quant.sf",
                 {"tx1": 60, "tx2": 40, "tx3": 10, "tx4": 100},
                 eff_lengths={"tx1": 900, "tx2": 300, "tx3": 900, "tx4": 450})
    _write_quant(base / "cohortA" / "intermediate" / "expression_qu" / "S2" / "quant.sf",
                 {"tx1": 30, "tx2": 20, "tx3": 50, "tx4": 90},
                 eff_lengths={"tx1": 900, "tx2": 300, "tx3": 900, "tx4": 450})

    # ---- alt_TSS (ratios) + txrevise group quants ----
    _write_bed(unnormA / "alt_TSS.bed",
               ["G1__grp_1_upstream_u1", "G1__grp_1_upstream_u2",
                "G2__grp_2_upstream_v1", "G2__grp_2_upstream_v2"], {
        "S1": [0.7, 0.3, 0.55, 0.45], "S2": [0.65, 0.35, 0.52, 0.48],
        "S3": [0.6, 0.4, 0.5, 0.5], "S4": [0.5, 0.5, 0.6, 0.4],
        "S7": [0.1, 0.9, 0.2, 0.8], "S8": [0.9, 0.1, 0.8, 0.2],
    })
    alt = base / "cohortA" / "intermediate" / "alt_TSS_polyA"
    _write_quant(alt / "grp_1.upstream" / "S1" / "quant.sf",
                 {"G1.grp_1.upstream.u1": 70, "G1.grp_1.upstream.u2": 30})
    _write_quant(alt / "grp_1.upstream" / "S2" / "quant.sf",
                 {"G1.grp_1.upstream.u1": 65, "G1.grp_1.upstream.u2": 35})
    _write_quant(alt / "grp_2.upstream" / "S1" / "quant.sf",
                 {"G2.grp_2.upstream.v1": 55, "G2.grp_2.upstream.v2": 45})
    _write_quant(alt / "grp_2.upstream" / "S2" / "quant.sf",
                 {"G2.grp_2.upstream.v1": 52, "G2.grp_2.upstream.v2": 48})

    # ---- alt_polyA (2 features -> insufficient data for concordance) ----
    _write_bed(unnormA / "alt_polyA.bed",
               ["G1__grp_1_downstream_w1", "G1__grp_1_downstream_w2"], {
        "S1": [0.8, 0.2], "S2": [0.75, 0.25],
        "S3": [0.6, 0.4], "S4": [0.5, 0.5],
        "S7": [0.1, 0.9], "S8": [0.9, 0.1],
    })

    # ---- splicing BED + numers ----
    _write_bed(unnormA / "splicing.bed",
               ["G1__chr1_100_200_clu_1_+", "G1__chr1_150_200_clu_1_+",
                "G2__chr1_300_400_clu_2_-"], {
        "S1": [0.8, 0.2, 1.0], "S2": [0.8, 0.2, 1.0],
        "S3": [0.7, 0.3, 1.0], "S4": [0.6, 0.4, 1.0],
        "S7": [0.9, 0.1, 0.5], "S8": [0.1, 0.9, 0.5],
    })
    _write_numers(base / "cohortA" / "intermediate" / "splicing" /
                  "leafcutter_perind_numers.counts.gz", {
        "S1": [80, 20, 50, 30, 70], "S2": [40, 10, 25, 10, 30],
        "S3": [70, 30, 100, 20, 80], "S4": [60, 40, 80, 20, 60],
        "S7": [90, 10, 50, 10, 90], "S8": [10, 90, 50, 10, 90],
    })

    # ---- stability BED + featureCounts ----
    _write_bed(unnormA / "stability.bed", ["G1", "G2", "G3", "G4"], {
        "S1": [2.0, 2.5, 3.0, 1.5],
        "S2": [2.0, 135 / 65, 95 / 45, 0.25],
        "S3": [2.1, 2.2, 2.3, 2.4], "S4": [1.9, 2.1, 2.2, 2.3],
        "S7": [4.0, 0.5, 3.0, 0.25], "S8": [0.5, 4.0, 0.25, 3.0],
    })
    stab = base / "cohortA" / "intermediate" / "stability"
    _write_featurecounts(stab / "S1.exonic.counts.txt",
                         {"G1": 100, "G2": 250, "G3": 150, "G4": 60})
    _write_featurecounts(stab / "S1.intronic.counts.txt",
                         {"G1": 50, "G2": 100, "G3": 50, "G4": 40})
    _write_featurecounts(stab / "S2.exonic.counts.txt",
                         {"G1": 60, "G2": 135, "G3": 95, "G4": 5})
    _write_featurecounts(stab / "S2.intronic.counts.txt",
                         {"G1": 30, "G2": 65, "G3": 45, "G4": 20})

    # ---- RNA_editing BED + site matrix + site map ----
    _write_bed(unnormA / "RNA_editing.bed",
               ["G1__RNAedit_Cluster_1", "G2__e3", "G3__e4"], {
        "S1": [0.26378, 0.24678, 0.27273],
        "S2": [0.25466, 0.24678, 0.25773],
        "S3": [0.21469, 0.21311, 0.23810],
        "S4": [0.22399, 0.22581, 0.27273],
        "S7": [0.52381, 0.52381, 0.33333],
        "S8": [0.51336, 0.02439, 0.71429],
    })
    em = base / "cohortA" / "intermediate" / "RNA_editing" / "edit_site_matrix.tsv"
    em.parent.mkdir(parents=True, exist_ok=True)
    em.write_text(
        "site\tS1\tS2\tS3\tS4\tS7\tS8\n"
        "e1\t5/20\t15/60\t8/40\t2/10\t1/10\t9/10\n"
        "e2\t10/40\t30/120\t4/20\t8/40\t9/10\t2/20\n"
        "e3\t0/0\t0/0\t6/30\t3/15\t5/10\t0/20\n"
        "e4\t4/16\t12/48\t2/10\t1/5\t3/10\t7/10\n")
    (base / "cohortA" / "output" / "RNA_editing.site_to_phenotype.tsv").write_text(
        "site\tgene_id\tphenotype_id\n"
        "e1\tG1\tG1__RNAedit_Cluster_1\n"
        "e2\tG1\tG1__RNAedit_Cluster_1\n"
        "e3\tG2\tG2__e3\n"
        "e4\tG3\tG3__e4\n")

    # ---- intron_retention BED + PSI ----
    _write_bed(unnormA / "intron_retention.bed", ["IR1", "IR2", "IR3"], {
        "S1": [0.1, 0.2, 0.3], "S2": [0.15, 0.18, 0.32],
        "S3": [0.2, 0.3, 0.4], "S4": [0.3, 0.4, 0.5],
        "S7": [0.1, 0.5, 0.9], "S8": [0.9, 0.5, 0.1],
    })
    _write_psi(base / "cohortA" / "intermediate" / "intron_retention" /
               "retained_intron_psi.tsv.gz", {
        "S1": [0.1, 0.2], "S2": [0.15, 0.18], "S3": [0.2, 0.3],
        "S4": [0.3, 0.4], "S7": [0.1, 0.5], "S8": [0.9, 0.5],
    })

    # ---- cohortB: cross-pair run S5 + singleton S6 ----
    unnormB = base / "cohortB" / "output" / "unnorm"
    for mod, pids in [
            ("expression", ["G1", "G2", "G3", "G4"]),
            ("isoforms", ["G1__tx1", "G1__tx2", "G2__tx4"]),
            ("isoform_expression", ["G1__tx1", "G1__tx2", "G2__tx4"]),
            ("alt_TSS", ["G1__grp_1_upstream_u1", "G1__grp_1_upstream_u2",
                         "G2__grp_2_upstream_v1", "G2__grp_2_upstream_v2"]),
            ("alt_polyA", ["G1__grp_1_downstream_w1", "G1__grp_1_downstream_w2"]),
            ("splicing", ["G1__chr1_100_200_clu_1_+", "G1__chr1_150_200_clu_1_+",
                          "G2__chr1_300_400_clu_2_-"]),
            ("stability", ["G1", "G2", "G3", "G4"]),
            ("RNA_editing", ["G1__RNAedit_Cluster_1", "G2__e3", "G3__e4"]),
            ("intron_retention", ["IR1", "IR2", "IR3"])]:
        _write_bed(unnormB / f"{mod}.bed", pids, {
            "S5": [0.5] * len(pids), "S6": [0.25] * len(pids)})
    _write_numers(base / "cohortB" / "intermediate" / "splicing" /
                  "leafcutter_perind_numers.counts.gz", {
        "S5": [10, 20, 30, 40, 50], "S6": [5, 5, 5, 5, 5]})
    _write_psi(base / "cohortB" / "intermediate" / "intron_retention" /
               "retained_intron_psi.tsv.gz", {
        "S5": [0.5, 0.5], "S6": [0.25, 0.25]})

    return root


def _run_collapse(root, extra_args=()):
    out = root / "staging"
    cmd = [
        sys.executable,
        str(Path(__file__).parent / "collapse_replicates.py"),
        "--metadata", str(root / "metadata.tsv"),
        "--ancestry-map", str(root / "ancestry_map.tsv"),
        "--config", str(root / "config.yml"),
        "--ref-anno", str(root / "ref.gtf"),
        "--out-dir", str(out),
        "--concordance-min-features", "3",
        *extra_args,
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"collapse failed:\n{res.stderr}\n{res.stdout}"
    return out, res


@pytest.fixture(scope="module")
def collapsed(fixture):
    out, res = _run_collapse(fixture)
    return fixture, out, res


def _staged_bed(out, cohort, modality):
    return pd.read_csv(out / cohort / "output" / "unnorm" / f"{modality}.bed",
                       sep="\t")


# ---------------------------------------------------------------------------
# Unit tests: parsers
# ---------------------------------------------------------------------------

def test_parse_numers_header_conventions(tmp_path):
    rows = {"chr1:100:200:clu_1_+": [5, 10]}
    # samples-only header
    p1 = tmp_path / "a.counts.gz"
    with gzip.open(p1, "wt") as f:
        f.write("S1 S2\n")
        for r, v in rows.items():
            f.write(f"{r} {v[0]} {v[1]}\n")
    ids, _, mat = cr.parse_numers(p1)
    assert ids == ["S1", "S2"]
    assert mat.loc["chr1:100:200:clu_1_+", "S2"] == 10
    # 'chrom'-token header
    p2 = tmp_path / "b.counts.gz"
    with gzip.open(p2, "wt") as f:
        f.write("chrom S1 S2\n")
        for r, v in rows.items():
            f.write(f"{r} {v[0]} {v[1]}\n")
    ids2, _, mat2 = cr.parse_numers(p2)
    assert ids2 == ["S1", "S2"]
    assert mat2.loc["chr1:100:200:clu_1_+", "S1"] == 5


def test_parse_splicing_pid():
    assert cr._parse_splicing_pid("G1__chr1_100_200_clu_1_+") == \
        ("chr1", 100, 200, 1, "+")
    assert cr._parse_splicing_pid("G2__chr1_300_400_clu_2_-") == \
        ("chr1", 300, 400, 2, "-")
    assert cr._parse_splicing_pid("G1__notajunction") is None
    assert cr._parse_splicing_pid("noSEP") is None


# ---------------------------------------------------------------------------
# End-to-end tests
# ---------------------------------------------------------------------------

def test_expression_count_sum(collapsed):
    _, out, _ = collapsed
    bed = _staged_bed(out, "cohortA", "expression")
    assert "S2" not in bed.columns and "S8" not in bed.columns
    assert list(bed["S1"]) == [150, 30, 950, 65]          # summed
    assert list(bed["S7"]) == [190, 2050, 13, 505]        # discordant: summed
    assert list(bed["S4"]) == [300, 60, 700, 45]          # untouched
    assert "S3" not in bed.columns                        # cross-protocol: dropped


def test_isoform_expression_count_sum(collapsed):
    _, out, _ = collapsed
    bed = _staged_bed(out, "cohortA", "isoform_expression")
    assert list(bed["S1"]) == [90, 60, 190]
    assert list(bed["S7"]) == [60, 60, 60]                # discordant: summed


def test_isoforms_ratio_recomputed_full_denominator(collapsed):
    """Collapsed ratios are molecule fractions (length-normalized) and must
    use the full transcript denominator (including tx3, which is absent from
    the BED). Collapsed counts tx1..tx3 = 90/60/60 with lengths 900/300/900:
    x = (0.1, 0.2, 1/15), G1 total = 11/30, so tx1 = 3/11, tx2 = 6/11."""
    _, out, _ = collapsed
    bed = _staged_bed(out, "cohortA", "isoforms")
    assert bed["S1"].iloc[0] == pytest.approx(3 / 11, rel=1e-5)
    assert bed["S1"].iloc[1] == pytest.approx(6 / 11, rel=1e-5)
    assert bed["S1"].iloc[2] == pytest.approx(1.0, rel=1e-5)
    # discordant: still collapsed; the fixture has no S7/S8 quants, so the
    # fallback averages the BED ratio columns
    assert list(bed["S7"]) == [0.5, 0.5, 0.5]


def test_alt_tss_ratio_recomputed(collapsed):
    _, out, _ = collapsed
    bed = _staged_bed(out, "cohortA", "alt_TSS")
    assert bed["S1"].iloc[0] == pytest.approx(135 / 200, rel=1e-5)
    assert bed["S1"].iloc[1] == pytest.approx(65 / 200, rel=1e-5)
    assert bed["S1"].iloc[2] == pytest.approx(107 / 200, rel=1e-5)
    assert bed["S1"].iloc[3] == pytest.approx(93 / 200, rel=1e-5)
    # discordant pair, no per-run group quants -> BED-column mean
    assert list(bed["S7"]) == [0.5, 0.5, 0.5, 0.5]


def test_alt_polyA_insufficient_data_collapsed(collapsed):
    """Too few features for concordance -> flagged, but still collapsed.
    The fixture has no per-run downstream quants -> BED-column mean."""
    _, out, _ = collapsed
    bed = _staged_bed(out, "cohortA", "alt_polyA")
    assert "S2" not in bed.columns
    assert list(bed["S1"]) == [0.775, 0.225]              # mean of S1, S2


def test_stability_floor_after_sum(collapsed):
    """G4: S2 exonic count 5 is below the per-feature floor; summing first
    (60+5=65 >= 10) rescues it: ratio = 65/60, not 60/60."""
    _, out, _ = collapsed
    bed = _staged_bed(out, "cohortA", "stability")
    assert bed["S1"].iloc[0] == pytest.approx(2.0, rel=1e-5)
    assert bed["S1"].iloc[1] == pytest.approx(385 / 165, rel=1e-5)
    assert bed["S1"].iloc[2] == pytest.approx(245 / 95, rel=1e-5)
    assert bed["S1"].iloc[3] == pytest.approx(65 / 60, rel=1e-5)
    # discordant pair, no per-run featureCounts -> BED-column mean
    assert list(bed["S7"]) == [2.25, 2.25, 1.625, 1.625]


def test_rna_editing_recomputed(collapsed):
    _, out, _ = collapsed
    bed = _staged_bed(out, "cohortA", "RNA_editing")
    assert bed["S1"].iloc[0] == pytest.approx(
        (20.5 / 80.5 + 40.5 / 160.5) / 2, rel=1e-5)
    # e3 is 0/0 in both runs -> site row-mean fraction across other samples
    e3_mean = np.mean([6.5 / 30.5, 3.5 / 15.5, 5.5 / 10.5, 0.5 / 20.5])
    assert bed["S1"].iloc[1] == pytest.approx(e3_mean, rel=1e-5)
    assert bed["S1"].iloc[2] == pytest.approx(16.5 / 64.5, rel=1e-5)
    # discordant pair S7+S8: recomputed from the edit matrix
    assert bed["S7"].iloc[0] == pytest.approx(
        (10.5 / 20.5 + 11.5 / 30.5) / 2, rel=1e-5)
    assert bed["S7"].iloc[1] == pytest.approx(5.5 / 30.5, rel=1e-5)
    assert bed["S7"].iloc[2] == pytest.approx(10.5 / 20.5, rel=1e-5)


def test_splicing_numers_collapsed(collapsed):
    _, out, _ = collapsed
    # No staged per-cohort splicing BED (flows through harmonization)
    assert not (out / "cohortA" / "output" / "unnorm" / "splicing.bed").exists()
    num = out / "cohortA" / "intermediate" / "splicing" / \
        "leafcutter_perind_numers.counts.gz"
    assert num.exists()
    ids, row_ids, mat = cr.parse_numers(num)
    assert ids == ["S1", "S4", "S7"]  # S2/S8 collapsed away; S3 dropped (cross)
    assert list(mat["S1"]) == [120, 30, 75, 40, 100]
    assert list(mat["S7"]) == [100, 100, 100, 20, 180]  # discordant: summed
    # samples-only header
    with gzip.open(num, "rt") as f:
        header = f.readline().strip().split(" ")
        first = f.readline().strip().split(" ")
    assert len(header) == len(first) - 1


def test_intron_retention_averaged(collapsed):
    _, out, _ = collapsed
    assert not (out / "cohortA" / "output" / "unnorm" /
                "intron_retention.bed").exists()
    psi = pd.read_csv(out / "cohortA" / "intermediate" / "intron_retention" /
                      "retained_intron_psi.tsv.gz", sep="\t")
    assert "S2" not in psi.columns and "S8" not in psi.columns
    assert "S3" not in psi.columns                  # cross-protocol: dropped
    assert list(psi["S1"]) == [0.125, 0.19]         # mean of S1, S2
    assert list(psi["S7"]) == [0.5, 0.5]            # discordant: averaged


def test_cross_protocol_dropped(collapsed):
    """A2 (S3 in cohortA, S5 in cohortB): a sample spanning cohorts is a
    labeling error -- the individual is removed from every cohort file."""
    _, out, _ = collapsed
    bedA = _staged_bed(out, "cohortA", "expression")
    assert "S3" not in bedA.columns
    bedB = _staged_bed(out, "cohortB", "expression")
    assert "S5" not in bedB.columns
    assert "S6" in bedB.columns
    # cohortB numers + PSI also drop S5
    ids, _, _ = cr.parse_numers(out / "cohortB" / "intermediate" / "splicing" /
                                "leafcutter_perind_numers.counts.gz")
    assert ids == ["S6"]
    psiB = pd.read_csv(out / "cohortB" / "intermediate" / "intron_retention" /
                       "retained_intron_psi.tsv.gz", sep="\t")
    assert "S5" not in psiB.columns


def test_untouched_cohort_symlinked(collapsed):
    _, out, _ = collapsed
    c = out / "cohortC"
    assert c.is_symlink()


def test_reports(collapsed):
    root, out, _ = collapsed
    coll = pd.read_csv(out / "reports" / "replicate_collapses.tsv", sep="\t")
    conc = pd.read_csv(out / "reports" / "replicate_concordance.tsv", sep="\t")

    def action(aid, mod):
        rows = coll[(coll.array_id == aid) & (coll.modality == mod)]
        return set(rows["action"])

    assert action("A1", "expression") == {"collapsed"}
    assert action("A1", "isoforms") == {"collapsed"}
    assert action("A1", "splicing") == {"collapsed"}
    assert action("A1", "alt_polyA") == {"collapsed"}   # insufficient: flag-only
    assert action("A1", "intron_retention") == {"collapsed"}
    assert action("A3", "expression") == {"collapsed"}  # discordant: flag-only
    assert action("A3", "splicing") == {"collapsed"}
    assert action("A2", "expression") == {"dropped_cross_protocol"}
    # A2 kept column: none in either cohort (individual removed)
    a2 = coll[(coll.array_id == "A2") & (coll.modality == "expression")]
    assert set(a2["kept_column"].fillna("")) == {""}

    def flag(aid, mod):
        rows = conc[(conc.array_id == aid) & (conc.modality == mod)]
        assert len(rows) == 1
        return rows["flag"].iloc[0]

    assert flag("A1", "expression") == "ok"
    assert flag("A3", "expression") == "discordant"
    assert flag("A1", "alt_polyA") == "insufficient_data"
    assert flag("A2", "expression") == "cross_protocol"
    s = conc[(conc.array_id == "A1") & (conc.modality == "expression")]["spearman"].iloc[0]
    assert s == pytest.approx(1.0)

    # per-ancestry splits exist
    assert (out / "reports" / "EAS_replicate_collapses.tsv").exists()
    # collapsed ancestry map drops non-primary runs (S2/S8 not in map) and
    # cross-protocol individuals entirely (S3 and S5)
    amap = pd.read_csv(out / "reports" / "ancestry_map_collapsed.tsv", sep="\t")
    assert set(amap["sample_id"]) == {"S1", "S4", "S6", "S7", "S9"}


def test_dry_run_writes_no_staging(fixture):
    out = fixture / "staging_dry"
    cmd = [
        sys.executable,
        str(Path(__file__).parent / "collapse_replicates.py"),
        "--metadata", str(fixture / "metadata.tsv"),
        "--ancestry-map", str(fixture / "ancestry_map.tsv"),
        "--config", str(fixture / "config.yml"),
        "--ref-anno", str(fixture / "ref.gtf"),
        "--out-dir", str(out),
        "--concordance-min-features", "3",
        "--dry-run",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert "DRY RUN" in res.stdout
    assert not (out / "cohortA").exists()
    assert (out / "reports" / "replicate_collapses.tsv").exists()


def test_isoforms_requires_ref_anno(fixture):
    cmd = [
        sys.executable,
        str(Path(__file__).parent / "collapse_replicates.py"),
        "--metadata", str(fixture / "metadata.tsv"),
        "--ancestry-map", str(fixture / "ancestry_map.tsv"),
        "--config", str(fixture / "config.yml"),
        "--out-dir", str(fixture / "staging_x"),
        "--modalities", "isoforms",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 1
    assert "--ref-anno" in res.stderr
