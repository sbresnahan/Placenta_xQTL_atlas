#!/usr/bin/env python3
"""
test_harmonize.py — Synthetic test for harmonize_within_ancestry.py

Creates two fake cohorts with:
  - Splicing: leafCutter counts with known shared junctions (should merge
    into the same meta-cluster) and cohort-specific junctions.
  - IR: MAJIQ PSI tables with known shared events (exact coordinate match)
    and cohort-specific events.

Then runs the harmonization logic directly (importing the module functions)
and verifies the output structure and harmonization correctness.

Does NOT call assemble_bed.py (which needs a real GTF); instead it checks the
pooled intermediates that feed into assemble_bed.py.

Sample IDs in samples.txt are BARE (e.g. s1, s2) — matching real pipeline
behavior where the script namespaces them as {cohort}_{sample}.
"""

import gzip
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

# Import the module under test
sys.path.insert(0, "/mnt/results")
import harmonize_within_ancestry as hwa


def make_leafcutter_counts(path, sample_ids, rows):
    """Write a leafcutter_perind_numers.counts.gz file.

    rows: list of (chrom, start, end, clu_n, strand, [counts...])
    """
    with gzip.open(path, "wt") as f:
        f.write("chrom " + " ".join(sample_ids) + "\n")
        for chrom, start, end, clu_n, strand, counts in rows:
            row_id = f"{chrom}:{start}:{end}:clu_{clu_n}_{strand}"
            counts_str = " ".join(str(int(c)) for c in counts)
            f.write(row_id + " " + counts_str + "\n")


def make_ir_psi(path, sample_ids, events):
    """Write a retained_intron_psi.tsv.gz file.

    events: list of dict with keys:
        seqid, start, end, strand, gene_id_base, gene_id, gene_name,
        majiq_ir_id, ec_idx, + per-sample PSI values
    """
    metadata_cols = [
        "majiq_ir_id", "ec_idx", "seqid", "start", "end", "strand",
        "gene_id", "gene_id_base", "gene_name", "event_type",
        "is_denovo", "event_denovo", "ref_exon_start", "ref_exon_end",
        "other_exon_start", "other_exon_end",
    ]
    cols = metadata_cols + sample_ids
    df = pd.DataFrame(events, columns=cols)
    df.to_csv(path, sep="\t", index=False, compression="gzip", float_format="%g")


def test_splicing_harmonization():
    """Two cohorts share junctions that should merge into one meta-cluster."""
    print("\n=== TEST: Splicing harmonization ===")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        output_base = tmpdir / "PANTRY"
        # Cohort structure — bare sample IDs in samples.txt
        for cohort in ["cohortA", "cohortB"]:
            (output_base / cohort / "intermediate" / "splicing").mkdir(parents=True)
            (output_base / cohort).mkdir(parents=True, exist_ok=True)
            (output_base / cohort / "samples.txt").write_text("s1\ns2\n")

        # Cohort A: 2 samples (s1, s2), clusters 1 and 2
        #   cluster 1: junctions chr1:100:200:+, chr1:100:300:+  (shared with B cluster 1)
        #   cluster 2: junctions chr1:500:600:+                   (A-only)
        make_leafcutter_counts(
            output_base / "cohortA" / "intermediate" / "splicing" / "leafcutter_perind_numers.counts.gz",
            ["s1", "s2"],
            [
                ("chr1", 100, 200, 1, "+", [10, 20]),
                ("chr1", 100, 300, 1, "+", [5, 15]),
                ("chr1", 500, 600, 2, "+", [8, 12]),
            ],
        )
        # Cohort B: 2 samples (s1, s2), clusters 1 and 5
        #   cluster 1: junctions chr1:100:200:+, chr1:100:400:+  (shares chr1:100:200:+ with A cluster 1)
        #   cluster 5: junctions chr2:700:800:-                  (B-only)
        make_leafcutter_counts(
            output_base / "cohortB" / "intermediate" / "splicing" / "leafcutter_perind_numers.counts.gz",
            ["s1", "s2"],
            [
                ("chr1", 100, 200, 1, "+", [30, 40]),
                ("chr1", 100, 400, 1, "+", [25, 35]),
                ("chr2", 700, 800, 5, "-", [7, 9]),
            ],
        )

        # Ancestry map: all samples are EUR (bare IDs match samples.txt)
        ancestry_map = tmpdir / "ancestry.tsv"
        lines = ["sample_id\tancestry"]
        for cohort in ["cohortA", "cohortB"]:
            for s in ["s1", "s2"]:
                lines.append(f"{s}\tEUR")
        ancestry_map.write_text("\n".join(lines) + "\n")

        # Mock config
        config = {
            "OUTPUT_BASE": str(output_base),
            "PANTRY_SCRIPTS": "/fake/pantry_scripts",
            "NORMALIZED_GTF": "/fake/gtf",
        }
        cohorts = ["cohortA", "cohortB"]

        # Load ancestry samples
        ancestry_df = hwa.load_ancestry_map(str(ancestry_map))
        ancestry_samples = {}
        for cohort in cohorts:
            sf = output_base / cohort / "samples.txt"
            cs = hwa.get_cohort_samples_in_ancestry(str(sf), ancestry_df, "EUR")
            ancestry_samples[cohort] = cs

        # Verify ancestry sample loading
        assert ancestry_samples["cohortA"] == ["s1", "s2"]
        assert ancestry_samples["cohortB"] == ["s1", "s2"]
        print(f"  PASS: ancestry sample loading (bare IDs)")

        # --- Parse per-cohort counts ---
        cohort_data = {}
        junction_to_clusters = {}
        for cohort in cohorts:
            counts_path = output_base / cohort / "intermediate" / "splicing" / "leafcutter_perind_numers.counts.gz"
            sample_ids, rows = hwa.parse_leafcutter_counts(str(counts_path))
            cohort_samples = ancestry_samples[cohort]
            col_idx = np.array([sample_ids.index(s) for s in cohort_samples])
            namespaced = [f"{cohort}_{s}" for s in cohort_samples]
            subset_rows = []
            for jk, clu_n, counts in rows:
                cs = counts[col_idx]
                subset_rows.append((jk, clu_n, cs))
                junction_to_clusters.setdefault(jk, set()).add((cohort, clu_n))
            cohort_data[cohort] = {"sample_ids": namespaced, "rows": subset_rows}

        # Verify namespacing
        assert cohort_data["cohortA"]["sample_ids"] == ["cohortA_s1", "cohortA_s2"]
        assert cohort_data["cohortB"]["sample_ids"] == ["cohortB_s1", "cohortB_s2"]
        print(f"  PASS: sample namespacing (cohortA_s1, etc.)")

        # --- Build connected components ---
        all_clusters = sorted({(c, cl) for c, d in cohort_data.items() for _, cl, _ in d["rows"]})
        cluster_index = {c: i for i, c in enumerate(all_clusters)}
        n = len(all_clusters)
        edges_src, edges_dst = [], []
        for jk, cset in junction_to_clusters.items():
            clist = sorted(cset)
            for i in range(len(clist)):
                for j in range(i + 1, len(clist)):
                    edges_src.append(cluster_index[clist[i]])
                    edges_dst.append(cluster_index[clist[j]])

        n_comp, labels = hwa._union_find(n, edges_src, edges_dst)

        # --- Assertions ---
        # CohortA cluster 1 and CohortB cluster 1 share chr1:100:200:+
        # => they should be in the SAME meta-cluster
        a_c1 = cluster_index[("cohortA", 1)]
        b_c1 = cluster_index[("cohortB", 1)]
        assert labels[a_c1] == labels[b_c1], (
            f"FAIL: cohortA clu_1 and cohortB clu_1 should be in the same meta-cluster, "
            f"got {labels[a_c1]} vs {labels[b_c1]}"
        )
        print(f"  PASS: cohortA clu_1 and cohortB clu_1 merged (meta-cluster {labels[a_c1]})")

        # CohortA cluster 2 (A-only) should be in a DIFFERENT meta-cluster
        a_c2 = cluster_index[("cohortA", 2)]
        assert labels[a_c1] != labels[a_c2], (
            f"FAIL: cohortA clu_2 should be in a different meta-cluster than the shared one"
        )
        print(f"  PASS: cohortA clu_2 is separate (meta-cluster {labels[a_c2]})")

        # CohortB cluster 5 (B-only) should be in a DIFFERENT meta-cluster
        b_c5 = cluster_index[("cohortB", 5)]
        assert labels[a_c1] != labels[b_c5], (
            f"FAIL: cohortB clu_5 should be in a different meta-cluster"
        )
        print(f"  PASS: cohortB clu_5 is separate (meta-cluster {labels[b_c5]})")

        # Total meta-clusters: 3 (shared, A-only, B-only)
        assert n_comp == 3, f"FAIL: expected 3 meta-clusters, got {n_comp}"
        print(f"  PASS: {n_comp} meta-clusters (expected 3)")

        # Junction chr1:100:200:+ should appear in both cohorts
        shared_jk = ("chr1", 100, 200, "+")
        assert shared_jk in junction_to_clusters, "FAIL: shared junction not found"
        assert ("cohortA", 1) in junction_to_clusters[shared_jk]
        assert ("cohortB", 1) in junction_to_clusters[shared_jk]
        print(f"  PASS: shared junction chr1:100:200:+ present in both cohorts")

        print("=== Splicing harmonization: ALL TESTS PASSED ===")


def test_ir_harmonization():
    """Two cohorts share IR events by exact coordinate match."""
    print("\n=== TEST: IR harmonization ===")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        output_base = tmpdir / "PANTRY"
        for cohort in ["cohortA", "cohortB"]:
            (output_base / cohort / "intermediate" / "intron_retention").mkdir(parents=True)
            (output_base / cohort).mkdir(parents=True, exist_ok=True)
            (output_base / cohort / "samples.txt").write_text("s1\ns2\n")

        metadata_cols = [
            "majiq_ir_id", "ec_idx", "seqid", "start", "end", "strand",
            "gene_id", "gene_id_base", "gene_name", "event_type",
            "is_denovo", "event_denovo", "ref_exon_start", "ref_exon_end",
            "other_exon_start", "other_exon_end",
        ]

        def make_event(seqid, start, end, strand, gene_id_base, gene_id,
                       gene_name, ir_id, ec_idx, psi_s1, psi_s2):
            e = {mc: "" for mc in metadata_cols}
            e.update({
                "majiq_ir_id": ir_id, "ec_idx": ec_idx,
                "seqid": seqid, "start": start, "end": end, "strand": strand,
                "gene_id": gene_id, "gene_id_base": gene_id_base,
                "gene_name": gene_name, "event_type": "intron_retention",
                "is_denovo": "false",
            })
            e["s1"] = psi_s1
            e["s2"] = psi_s2
            return e

        # Cohort A events:
        #   Event 1: chr1:1000:2000:+ ENSG00000001 (shared with B)
        #   Event 2: chr1:5000:6000:+ ENSG00000002 (A-only)
        events_a = [
            make_event("chr1", 1000, 2000, "+", "ENSG00000001", "ENSG00000001.1",
                       "GENE1", "IR:ENSG00000001:chr1:1000:2000:+:0", 0,
                       0.8, 0.7),
            make_event("chr1", 5000, 6000, "+", "ENSG00000002", "ENSG00000002.1",
                       "GENE2", "IR:ENSG00000002:chr1:5000:6000:+:1", 1,
                       0.5, 0.6),
        ]
        make_ir_psi(
            output_base / "cohortA" / "intermediate" / "intron_retention" / "retained_intron_psi.tsv.gz",
            ["s1", "s2"],
            events_a,
        )

        # Cohort B events:
        #   Event 1: chr1:1000:2000:+ ENSG00000001 (shared with A, same coords)
        #   Event 3: chr2:3000:4000:- ENSG00000003 (B-only)
        events_b = [
            make_event("chr1", 1000, 2000, "+", "ENSG00000001", "ENSG00000001.1",
                       "GENE1", "IR:ENSG00000001:chr1:1000:2000:+:5", 5,
                       0.9, 0.85),
            make_event("chr2", 3000, 4000, "-", "ENSG00000003", "ENSG00000003.1",
                       "GENE3", "IR:ENSG00000003:chr2:3000:4000:-:6", 6,
                       0.4, 0.3),
        ]
        make_ir_psi(
            output_base / "cohortB" / "intermediate" / "intron_retention" / "retained_intron_psi.tsv.gz",
            ["s1", "s2"],
            events_b,
        )

        # Ancestry map (bare IDs)
        ancestry_map = tmpdir / "ancestry.tsv"
        lines = ["sample_id\tancestry"]
        for cohort in ["cohortA", "cohortB"]:
            for s in ["s1", "s2"]:
                lines.append(f"{s}\tEUR")
        ancestry_map.write_text("\n".join(lines) + "\n")

        config = {
            "OUTPUT_BASE": str(output_base),
            "PANTRY_SCRIPTS": "/fake/pantry_scripts",
            "NORMALIZED_GTF": "/fake/gtf",
        }
        cohorts = ["cohortA", "cohortB"]

        ancestry_df = hwa.load_ancestry_map(str(ancestry_map))
        ancestry_samples = {}
        for cohort in cohorts:
            sf = output_base / cohort / "samples.txt"
            ancestry_samples[cohort] = hwa.get_cohort_samples_in_ancestry(str(sf), ancestry_df, "EUR")

        # --- Replicate the IR loading logic (without assemble_bed.py call) ---
        event_data = {}
        all_samples = []
        for cohort in cohorts:
            psi_path = output_base / cohort / "intermediate" / "intron_retention" / "retained_intron_psi.tsv.gz"
            df = pd.read_csv(psi_path, sep="\t", dtype={"seqid": str})
            sample_cols = [c for c in df.columns if c not in metadata_cols]
            kept = [c for c in ancestry_samples[cohort] if c in sample_cols]
            namespaced = [f"{cohort}_{s}" for s in kept]
            all_samples.extend(namespaced)
            for _, row in df.iterrows():
                gidb = str(row.get("gene_id_base", "")).split(".")[0]
                ek = (str(row["seqid"]), int(row["start"]), int(row["end"]),
                      str(row["strand"]), gidb)
                if ek not in event_data:
                    meta = {mc: row[mc] for mc in metadata_cols if mc in row.index}
                    event_data[ek] = {"metadata": meta, "cohort_meta": {}, "cohort_psi": {}}
                cm = {mc: row[mc] for mc in metadata_cols if mc in row.index}
                event_data[ek]["cohort_meta"][cohort] = cm
                pv = {}
                for s, ns in zip(kept, namespaced):
                    v = row[s]
                    pv[ns] = v if pd.notna(v) else np.nan
                event_data[ek]["cohort_psi"][cohort] = pv

        # --- Assertions ---
        # 3 unique events total
        assert len(event_data) == 3, f"FAIL: expected 3 unique events, got {len(event_data)}"
        print(f"  PASS: {len(event_data)} unique events (expected 3)")

        # Shared event chr1:1000:2000:+ ENSG00000001 present in both cohorts
        shared_key = ("chr1", 1000, 2000, "+", "ENSG00000001")
        assert shared_key in event_data, "FAIL: shared event not found"
        assert "cohortA" in event_data[shared_key]["cohort_psi"], "FAIL: shared event missing cohortA"
        assert "cohortB" in event_data[shared_key]["cohort_psi"], "FAIL: shared event missing cohortB"
        print(f"  PASS: shared event chr1:1000:2000:+ ENSG00000001 present in both cohorts")

        # A-only event
        a_only = ("chr1", 5000, 6000, "+", "ENSG00000002")
        assert a_only in event_data, "FAIL: A-only event not found"
        assert "cohortA" in event_data[a_only]["cohort_psi"]
        assert "cohortB" not in event_data[a_only]["cohort_psi"]
        print(f"  PASS: A-only event chr1:5000:6000:+ present in cohortA only")

        # B-only event
        b_only = ("chr2", 3000, 4000, "-", "ENSG00000003")
        assert b_only in event_data, "FAIL: B-only event not found"
        assert "cohortB" in event_data[b_only]["cohort_psi"]
        assert "cohortA" not in event_data[b_only]["cohort_psi"]
        print(f"  PASS: B-only event chr2:3000:4000:- present in cohortB only")

        # Per-cohort original majiq_ir_id preserved in cohort_meta
        assert event_data[shared_key]["cohort_meta"]["cohortA"]["majiq_ir_id"] == "IR:ENSG00000001:chr1:1000:2000:+:0"
        assert event_data[shared_key]["cohort_meta"]["cohortB"]["majiq_ir_id"] == "IR:ENSG00000001:chr1:1000:2000:+:5"
        print(f"  PASS: per-cohort original majiq_ir_id preserved in cohort_meta")

        # PSI values correct (namespaced keys: cohortA_s1, cohortB_s1)
        assert event_data[shared_key]["cohort_psi"]["cohortA"]["cohortA_s1"] == 0.8
        assert event_data[shared_key]["cohort_psi"]["cohortB"]["cohortB_s1"] == 0.9
        print(f"  PASS: PSI values correct for shared event (namespaced keys)")

        print("=== IR harmonization: ALL TESTS PASSED ===")


def test_union_find():
    """Unit test for the union-find fallback."""
    print("\n=== TEST: Union-find fallback ===")
    # 5 nodes: 0-1-2 connected, 3-4 connected, no edges between groups
    n_comp, labels = hwa._union_find(5, [0, 1, 3], [1, 2, 4])
    assert n_comp == 2, f"FAIL: expected 2 components, got {n_comp}"
    assert labels[0] == labels[1] == labels[2], "FAIL: 0,1,2 should be connected"
    assert labels[3] == labels[4], "FAIL: 3,4 should be connected"
    assert labels[0] != labels[3], "FAIL: the two groups should be separate"
    print(f"  PASS: 2 components, correct grouping")
    print("=== Union-find: ALL TESTS PASSED ===")


def test_leafcutter_parsing():
    """Unit test for leafCutter counts parsing."""
    print("\n=== TEST: leafCutter counts parsing ===")
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "test.counts.gz"
        make_leafcutter_counts(
            path, ["s1", "s2"],
            [("chr1", 100, 200, 3, "+", [10, 20]),
             ("chrX", 500, 600, 7, "-", [5, 8])],
        )
        sample_ids, rows = hwa.parse_leafcutter_counts(str(path))
        assert sample_ids == ["s1", "s2"], f"FAIL: sample_ids={sample_ids}"
        assert len(rows) == 2, f"FAIL: expected 2 rows, got {len(rows)}"
        jk0, clu0, counts0 = rows[0]
        assert jk0 == ("chr1", 100, 200, "+"), f"FAIL: junction_key={jk0}"
        assert clu0 == 3, f"FAIL: clu_n={clu0}"
        assert list(counts0) == [10.0, 20.0], f"FAIL: counts={list(counts0)}"
        jk1, clu1, counts1 = rows[1]
        assert jk1 == ("chrX", 500, 600, "-"), f"FAIL: junction_key={jk1}"
        assert clu1 == 7
        print(f"  PASS: parsed 2 junctions correctly")
        print("=== leafCutter parsing: ALL TESTS PASSED ===")


if __name__ == "__main__":
    test_union_find()
    test_leafcutter_parsing()
    test_splicing_harmonization()
    test_ir_harmonization()
    print("\n\n>>> ALL TEST SUITES PASSED <<<")
