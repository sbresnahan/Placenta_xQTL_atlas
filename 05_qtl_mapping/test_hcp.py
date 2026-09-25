#!/usr/bin/env python3
"""
test_hcp.py — Synthetic test suite for the HCP latent factor extraction pipeline.

Tests all components with synthetic data:
  1. test_picard_qc: mock PicardTools metric files → verify parsing
  2. test_pool_expression: synthetic BEDs for 2 cohorts, 2 ancestries → verify pooling
  3. test_combat_normalize_hcp: synthetic expression + QC → verify QN, INT,
     ComBat-last batch removal, HCP
  3b. test_connectivity_outliers: corrupted sample → verify bicor
     connectivity outlier removal (z < -3)
  4. test_reflat_generation: verify GTF → refFlat conversion
  5. test_full_pipeline: small synthetic dataset through all scripts

Usage:
  python3 test_hcp.py
  python3 test_hcp.py --test test_pool_expression  # run single test
"""

import argparse
import os
import subprocess
import sys
import tempfile
import shutil

import numpy as np
import pandas as pd

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))


# =============================================================================
# Test helpers
# =============================================================================

def make_synthetic_bed(path, n_genes, samples, chrom="chr1", seed=42):
    """Create a synthetic expression BED file.

    BED format: #chr, start, end, phenotype_id, sample1, sample2, ...
    Values are TPM-like (non-negative floats).
    """
    rng = np.random.RandomState(seed)
    meta = []
    data = []
    for i in range(n_genes):
        gene_id = f"GENE{i:05d}"
        start = i * 1000
        end = start + 500
        meta.append((chrom, start, end, gene_id))
        # TPM values: log-normal distribution with gene-specific mean
        gene_mean = rng.exponential(scale=10) + 1
        values = rng.lognormal(mean=np.log(gene_mean), sigma=0.5, size=len(samples))
        data.append(values)

    df = pd.DataFrame(meta, columns=["#chr", "start", "end", "phenotype_id"])
    data_df = pd.DataFrame(np.array(data), columns=samples)
    result = pd.concat([df, data_df], axis=1)
    result.to_csv(path, sep="\t", index=False, float_format="%g")
    return result


def make_synthetic_qc(path, samples, n_metrics=10, seed=42):
    """Create a synthetic QC metrics TSV."""
    rng = np.random.RandomState(seed)
    metric_names = [f"QC_metric_{i}" for i in range(n_metrics)]
    data = rng.normal(size=(len(samples), n_metrics))
    df = pd.DataFrame(data, index=samples, columns=metric_names)
    df.index.name = "sample"
    df.to_csv(path, sep="\t")
    return df


def make_synthetic_ancestry_map(path, samples_by_cohort_ancestry):
    """Create a synthetic ancestry map TSV.

    Args:
        path: output file path
        samples_by_cohort_ancestry: dict of {(cohort, ancestry): [sample_ids]}
    """
    rows = []
    for (cohort, ancestry), samples in samples_by_cohort_ancestry.items():
        for s in samples:
            rows.append({"sample_id": s, "assigned_ancestry": ancestry, "cohort": cohort})
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)
    return df


def make_synthetic_gtf(path, n_genes=5, n_isoforms_per_gene=2):
    """Create a minimal synthetic GTF with gene, transcript, and exon features."""
    with open(path, "w") as f:
        f.write("##gff-version 2\n")
        for i in range(n_genes):
            gene_id = f"GENE{i:05d}"
            gene_name = f"GENE{i:05d}"
            chrom = "chr1"
            strand = "+" if i % 2 == 0 else "-"
            gene_start = i * 10000 + 1
            gene_end = gene_start + 5000

            f.write(f'{chrom}\tHPLR\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t'
                    f'gene_id "{gene_id}"; gene_name "{gene_name}";\n')

            for j in range(n_isoforms_per_gene):
                tx_id = f"{gene_id}.T{j+1}"
                tx_start = gene_start
                tx_end = gene_end
                f.write(f'{chrom}\tHPLR\ttranscript\t{tx_start}\t{tx_end}\t.\t{strand}\t.\t'
                        f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_name "{gene_name}";\n')

                # Two exons per transcript
                exon1_start = tx_start
                exon1_end = tx_start + 2000
                exon2_start = tx_start + 3000
                exon2_end = tx_end
                f.write(f'{chrom}\tHPLR\texon\t{exon1_start}\t{exon1_end}\t.\t{strand}\t.\t'
                        f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_name "{gene_name}";\n')
                f.write(f'{chrom}\tHPLR\texon\t{exon2_start}\t{exon2_end}\t.\t{strand}\t.\t'
                        f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_name "{gene_name}";\n')

                # CDS on first exon
                cds_start = exon1_start + 100
                cds_end = exon1_end - 100
                f.write(f'{chrom}\tHPLR\tCDS\t{cds_start}\t{cds_end}\t.\t{strand}\t0\t'
                        f'gene_id "{gene_id}"; transcript_id "{tx_id}"; gene_name "{gene_name}";\n')


def assert_close(a, b, tol=1e-6, msg=""):
    if abs(a - b) > tol:
        raise AssertionError(f"Expected {b}, got {a} (tol={tol}) {msg}")


def assert_true(cond, msg=""):
    if not cond:
        raise AssertionError(msg)


# =============================================================================
# Test 1: PicardTools metric parsing
# =============================================================================

def test_picard_qc():
    """Test PicardTools metric file parsing with mock metric files."""
    print("\n--- test_picard_qc ---")
    sys.path.insert(0, SCRIPTS_DIR)
    from picard_qc import parse_picard_metrics, coerce_numeric

    # Create a mock PicardTools metrics file.
    # PicardTools uses "## METRICS CLASS" (space, not underscore) as the header.
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write("## htsjdk.samtools.metrics.StringHeader\n")
        f.write("# Started on Mon Jan 01 00:00:00 UTC 2024\n")
        f.write("## METRICS CLASS\tpicard.analysis.AlignmentSummaryMetrics\n")
        f.write("CATEGORY\tTOTAL_READS\tREADS_ALIGNED\tMEAN_READ_LENGTH\tSTRAND_BALANCE\tPCT_READS_ALIGNED\n")
        f.write("FIRST_OF_PAIR\t1000000\t950000\t150\t0.5012\t0.95\n")
        f.write("SECOND_OF_PAIR\t1000000\t948000\t150\t0.4988\t0.948\n")
        mock_file = f.name

    metrics = parse_picard_metrics(mock_file)
    os.unlink(mock_file)

    # Should parse the first data row
    assert_true("TOTAL_READS" in metrics, "TOTAL_READS not in parsed metrics")
    assert_close(coerce_numeric(metrics["TOTAL_READS"]), 1000000, msg="TOTAL_READS")
    assert_close(coerce_numeric(metrics["PCT_READS_ALIGNED"]), 0.95, msg="PCT_READS_ALIGNED")
    assert_close(coerce_numeric(metrics["STRAND_BALANCE"]), 0.5012, msg="STRAND_BALANCE")

    # Test coerce_numeric edge cases
    assert_true(coerce_numeric("") is None, "empty string should be None")
    assert_true(coerce_numeric("?") is None, "? should be None")
    assert_true(coerce_numeric(None) is None, "None should be None")
    assert_close(coerce_numeric("3.14"), 3.14, msg="3.14")

    # Test subject-specific bias computation
    from picard_qc import compute_subject_specific_bias
    n_genes = 100
    n_samples = 10
    rng = np.random.RandomState(42)
    tpm_matrix = pd.DataFrame(
        rng.lognormal(mean=5, sigma=1, size=(n_genes, n_samples)),
        index=[f"GENE{i:05d}" for i in range(n_genes)],
        columns=[f"S{i}" for i in range(n_samples)],
    )
    gene_annot = pd.DataFrame({
        "gc_content": rng.uniform(0.3, 0.7, n_genes),
        "gene_length": rng.uniform(500, 5000, n_genes),
    }, index=tpm_matrix.index)

    bias_df = compute_subject_specific_bias(tpm_matrix, gene_annot)
    assert_true(len(bias_df) == n_samples, f"Expected {n_samples} rows, got {len(bias_df)}")
    assert_true("SubjectBias.GC" in bias_df.columns, "SubjectBias.GC column missing")
    assert_true("SubjectBias.GENE_LENGTH" in bias_df.columns, "SubjectBias.GENE_LENGTH column missing")
    # Correlations should be in [-1, 1]
    assert_true(bias_df["SubjectBias.GC"].between(-1, 1).all(), "GC bias out of range [-1,1]")
    assert_true(bias_df["SubjectBias.GENE_LENGTH"].between(-1, 1).all(),
                "Gene length bias out of range [-1,1]")

    print("  PASSED: metric parsing + subject-specific bias")


# =============================================================================
# Test 2: Expression pooling
# =============================================================================

def test_pool_expression():
    """Test pooling of expression BEDs across cohorts within ancestry strata."""
    print("\n--- test_pool_expression ---")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create 2 cohorts with 2 ancestries
        cohort1_samples = ["C1_S1", "C1_S2", "C1_S3", "C1_S4"]
        cohort2_samples = ["C2_S1", "C2_S2", "C2_S3", "C2_S4"]

        # Cohort 1: EUR samples C1_S1, C1_S2; EAS samples C1_S3, C1_S4
        # Cohort 2: EUR samples C2_S1, C2_S2; EAS samples C2_S3, C2_S4
        ancestry_map_path = os.path.join(tmpdir, "ancestry_map.tsv")
        make_synthetic_ancestry_map(ancestry_map_path, {
            ("cohort1", "EUR"): ["C1_S1", "C1_S2"],
            ("cohort1", "EAS"): ["C1_S3", "C1_S4"],
            ("cohort2", "EUR"): ["C2_S1", "C2_S2"],
            ("cohort2", "EAS"): ["C2_S3", "C2_S4"],
        })

        # Create cohort directories with unnorm expression BEDs
        for cohort, samples in [("cohort1", cohort1_samples), ("cohort2", cohort2_samples)]:
            cohort_dir = os.path.join(tmpdir, cohort)
            unnorm_dir = os.path.join(cohort_dir, "output", "unnorm")
            os.makedirs(unnorm_dir)
            make_synthetic_bed(
                os.path.join(unnorm_dir, "expression.bed"),
                n_genes=50, samples=samples, seed=hash(cohort) % 1000)

        # Run pooling
        output_dir = os.path.join(tmpdir, "pooled")
        cmd = [
            sys.executable, os.path.join(SCRIPTS_DIR, "pool_expression_within_ancestry.py"),
            "--ancestry-map", ancestry_map_path,
            "--cohort-dirs", f"cohort1={tmpdir}/cohort1", f"cohort2={tmpdir}/cohort2",
            "--output-dir", output_dir,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert_true(result.returncode == 0,
                    f"Pooling failed:\n{result.stderr}")

        # Check EUR pooled output
        eur_path = os.path.join(output_dir, "EUR_pooled_expression.bed")
        assert_true(os.path.exists(eur_path), f"EUR pooled BED not found: {eur_path}")

        eur_df = pd.read_csv(eur_path, sep="\t")
        eur_samples = [c for c in eur_df.columns if c not in ["#chr", "start", "end", "phenotype_id"]]
        assert_true(len(eur_samples) == 4, f"Expected 4 EUR samples, got {len(eur_samples)}")
        assert_true(set(eur_samples) == {"C1_S1", "C1_S2", "C2_S1", "C2_S2"},
                    f"Wrong EUR samples: {eur_samples}")

        # Check EAS pooled output
        eas_path = os.path.join(output_dir, "EAS_pooled_expression.bed")
        assert_true(os.path.exists(eas_path), f"EAS pooled BED not found: {eas_path}")

        eas_df = pd.read_csv(eas_path, sep="\t")
        eas_samples = [c for c in eas_df.columns if c not in ["#chr", "start", "end", "phenotype_id"]]
        assert_true(len(eas_samples) == 4, f"Expected 4 EAS samples, got {len(eas_samples)}")
        assert_true(set(eas_samples) == {"C1_S3", "C1_S4", "C2_S3", "C2_S4"},
                    f"Wrong EAS samples: {eas_samples}")

        # Check gene intersection (both cohorts have same genes → all 50 retained)
        assert_true(len(eur_df) == 50, f"Expected 50 genes, got {len(eur_df)}")

        # Check summary file
        summary_path = os.path.join(output_dir, "pooling_summary.tsv")
        assert_true(os.path.exists(summary_path), "Summary file not found")

        print("  PASSED: 2 cohorts, 2 ancestries, gene intersection, sample grouping")


# =============================================================================
# Test 3: ComBat + INT + HCP (requires R)
# =============================================================================

def test_combat_normalize_hcp():
    """Test ComBat + INT + HCP with synthetic data (requires R + sva + Rhcpp)."""
    print("\n--- test_combat_normalize_hcp ---")

    # Check if R is available
    r_available = shutil.which("Rscript") is not None
    if not r_available:
        print("  SKIPPED: Rscript not available")
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        n_genes = 200
        n_samples_cohort1 = 20
        n_samples_cohort2 = 20
        cohort1_samples = [f"C1_S{i}" for i in range(n_samples_cohort1)]
        cohort2_samples = [f"C2_S{i}" for i in range(n_samples_cohort2)]
        all_samples = cohort1_samples + cohort2_samples

        # Create pooled expression BED with a batch effect
        rng = np.random.RandomState(42)
        meta = []
        data = []
        for i in range(n_genes):
            gene_id = f"GENE{i:05d}"
            start = i * 1000
            end = start + 500
            meta.append(("chr1", start, end, gene_id))

            # Base expression
            base = rng.lognormal(mean=5, sigma=0.5, size=len(all_samples))
            # Inject batch effect: cohort2 gets an 8x multiplicative shift on
            # half the genes. Feature-specific (a uniform all-gene shift is
            # erased by QN itself) and strong relative to the lognormal
            # background, so ComBat has something to remove.
            if i < n_genes // 2:
                base[n_samples_cohort1:] = base[n_samples_cohort1:] * 8
            data.append(base)

        df = pd.DataFrame(meta, columns=["#chr", "start", "end", "phenotype_id"])
        data_df = pd.DataFrame(np.array(data), columns=all_samples)
        expr_df = pd.concat([df, data_df], axis=1)
        expr_path = os.path.join(tmpdir, "EUR_pooled_expression.bed")
        expr_df.to_csv(expr_path, sep="\t", index=False, float_format="%g")

        # Create QC metrics
        qc_path = os.path.join(tmpdir, "qc_metrics.tsv")
        make_synthetic_qc(qc_path, all_samples, n_metrics=8, seed=99)

        # Create ancestry map
        anc_map_path = os.path.join(tmpdir, "ancestry_map.tsv")
        make_synthetic_ancestry_map(anc_map_path, {
            ("cohort1", "EUR"): cohort1_samples,
            ("cohort2", "EUR"): cohort2_samples,
        })

        # Run QN + INT + ComBat + HCP. --outlier-z -999 disables connectivity
        # outlier removal so the sample set is deterministic in this test
        # (outlier detection itself is covered by test_connectivity_outliers).
        output_dir = os.path.join(tmpdir, "hcp_output")
        cmd = [
            "Rscript", os.path.join(SCRIPTS_DIR, "combat_normalize_hcp.R"),
            "--expression", expr_path,
            "--qc-metrics", qc_path,
            "--ancestry-map", anc_map_path,
            "--ancestry", "EUR",
            "--k", "5",
            "--outlier-z", "-999",
            "--output-dir", output_dir,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            # Check if it's a missing package issue
            if "there is no package called" in result.stderr or \
               "could not find function" in result.stderr:
                print(f"  SKIPPED: missing R package")
                print(f"    {result.stderr.strip().split(chr(10))[-1]}")
                return
            raise AssertionError(
                f"ComBat+HCP failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")

        # Check HCP factors output
        hcp_path = os.path.join(output_dir, "EUR_hcp_factors.tsv")
        assert_true(os.path.exists(hcp_path), f"HCP factors file not found: {hcp_path}")

        hcp_df = pd.read_csv(hcp_path, sep="\t")
        # tensorQTL format: first column = covariate name, rest = samples
        n_cov_cols = len(hcp_df.columns) - 1
        assert_true(n_cov_cols == len(all_samples),
                    f"Expected {len(all_samples)} sample columns, got {n_cov_cols}")
        assert_true(len(hcp_df) == 5, f"Expected 5 HCP factors, got {len(hcp_df)}")
        assert_true(hcp_df.iloc[:, 0].str.startswith("HCP_").all(),
                    "HCP factor names should start with HCP_")

        # Check expression output
        expr_out_path = os.path.join(output_dir, "EUR_combat_int_expression.bed")
        assert_true(os.path.exists(expr_out_path), "ComBat+INT expression BED not found")

        expr_out = pd.read_csv(expr_out_path, sep="\t")
        expr_samples = [c for c in expr_out.columns
                        if c not in ["#chr", "start", "end", "phenotype_id"]]
        assert_true(len(expr_samples) == len(all_samples),
                    f"Expected {len(all_samples)} samples in output, got {len(expr_samples)}")

        # Check that values are approximately standard normal. With ComBat
        # LAST (2026-09 schema) the output is approximately but not exactly
        # N(0,1) per gene — use a looser tolerance than the old ComBat->INT
        # ordering allowed.
        data_values = expr_out[expr_samples].values.flatten()
        data_values = data_values[~np.isnan(data_values)]
        mean_val = np.mean(data_values)
        std_val = np.std(data_values)
        assert_true(abs(mean_val) < 0.2, f"Final mean should be ~0, got {mean_val:.4f}")
        assert_true(abs(std_val - 1.0) < 0.2, f"Final std should be ~1, got {std_val:.4f}")

        # Outlier list: written even when empty (0 outliers with -999 threshold)
        outlier_path = os.path.join(output_dir, "EUR_expression_outliers.tsv")
        assert_true(os.path.exists(outlier_path), "Expression outliers file not found")
        outliers_df = pd.read_csv(outlier_path, sep="\t")
        assert_true(len(outliers_df) == 0,
                    f"Expected 0 outliers with --outlier-z -999, got {len(outliers_df)}")

        # ComBat-last batch removal: paired run with --skip-combat gives the
        # pre-ComBat (QN+INT) reference; the cohort mean difference on the
        # batch-affected genes must shrink after ComBat.
        nocombat_dir = os.path.join(tmpdir, "hcp_nocombat")
        cmd_nc = [
            "Rscript", os.path.join(SCRIPTS_DIR, "combat_normalize_hcp.R"),
            "--expression", expr_path,
            "--qc-metrics", qc_path,
            "--ancestry-map", anc_map_path,
            "--ancestry", "EUR",
            "--k", "5",
            "--outlier-z", "-999",
            "--skip-combat", "--skip-hcp",
            "--output-dir", nocombat_dir,
        ]
        result_nc = subprocess.run(cmd_nc, capture_output=True, text=True)
        assert_true(result_nc.returncode == 0,
                    f"--skip-combat run failed:\n{result_nc.stderr}")
        expr_nc = pd.read_csv(
            os.path.join(nocombat_dir, "EUR_combat_int_expression.bed"), sep="\t")
        c1, c2 = cohort1_samples, cohort2_samples
        affected = [f"GENE{i:05d}" for i in range(n_genes // 2)]
        def cohort_gap(df):
            sub = df[df["phenotype_id"].isin(affected)]
            return (sub[c1].mean(axis=1) - sub[c2].mean(axis=1)).abs().mean()
        gap_before = cohort_gap(expr_nc)
        gap_after = cohort_gap(expr_out)
        assert_true(gap_after < gap_before * 0.5,
                    f"ComBat should remove cohort shift: before={gap_before:.3f}, "
                    f"after={gap_after:.3f}")

        # Check diagnostics PDF
        diag_path = os.path.join(output_dir, "EUR_hcp_diagnostics.pdf")
        assert_true(os.path.exists(diag_path), "Diagnostics PDF not found")
        assert_true(os.path.getsize(diag_path) > 1000, "Diagnostics PDF too small (likely empty)")

        print("  PASSED: ComBat removes batch effect, INT produces ~N(0,1), "
              "HCP factors have correct dimensions")


# =============================================================================
# Test 3b: Connectivity outlier removal (devBrain §3.3)
# =============================================================================

def test_connectivity_outliers():
    """Test bicor connectivity outlier detection (z < -3) with one corrupted
    sample. Requires R + sva + Rhcpp."""
    print("\n--- test_connectivity_outliers ---")

    r_available = shutil.which("Rscript") is not None
    if not r_available:
        print("  SKIPPED: Rscript not available")
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        rng = np.random.RandomState(7)
        n_genes = 500
        n_good = 40
        # Co-varying gene module shared by all good samples: a same-sign
        # per-sample score shifts the first 250 genes, so good samples have
        # homogeneous, mutually correlated expression profiles while the
        # corrupted sample is uncorrelated with everyone. (Connectivity
        # outlier detection needs this shared structure — fully independent
        # samples are all equally "disconnected".)
        values = 2 + rng.rand(n_genes, n_good + 1)  # background in [2, 3]
        module = np.arange(n_genes // 2)
        grp = 1 + 0.3 * rng.normal(size=n_good + 1)
        values[module, :] += 2.0 * grp

        good_samples = [f"S{i}" for i in range(n_good)]
        bad_sample = "S_BAD"
        all_samples = good_samples + [bad_sample]

        # Corrupt the last sample: permute values across genes, destroying
        # its correlation with every other sample
        bad_col = all_samples.index(bad_sample)
        values[:, bad_col] = values[rng.permutation(n_genes), bad_col]

        meta = [("chr1", i * 1000, i * 1000 + 500, f"GENE{i:05d}")
                for i in range(n_genes)]
        expr_df = pd.concat([
            pd.DataFrame(meta, columns=["#chr", "start", "end", "phenotype_id"]),
            pd.DataFrame(values, columns=all_samples),
        ], axis=1)
        expr_path = os.path.join(tmpdir, "EUR_pooled_expression.bed")
        expr_df.to_csv(expr_path, sep="\t", index=False, float_format="%g")

        qc_path = os.path.join(tmpdir, "qc_metrics.tsv")
        make_synthetic_qc(qc_path, all_samples, n_metrics=6, seed=11)

        # Single cohort: ComBat skipped, isolating the outlier-removal step
        anc_map_path = os.path.join(tmpdir, "ancestry_map.tsv")
        make_synthetic_ancestry_map(anc_map_path, {
            ("cohort1", "EUR"): all_samples,
        })

        output_dir = os.path.join(tmpdir, "hcp_output")
        cmd = [
            "Rscript", os.path.join(SCRIPTS_DIR, "combat_normalize_hcp.R"),
            "--expression", expr_path,
            "--qc-metrics", qc_path,
            "--ancestry-map", anc_map_path,
            "--ancestry", "EUR",
            "--k", "3",
            "--output-dir", output_dir,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            if "there is no package called" in result.stderr or \
               "could not find function" in result.stderr:
                print("  SKIPPED: missing R package")
                return
            raise AssertionError(
                f"combat_normalize_hcp.R failed:\nSTDOUT:\n{result.stdout}\n"
                f"STDERR:\n{result.stderr}")

        # The corrupted sample must be flagged and removed
        outlier_path = os.path.join(output_dir, "EUR_expression_outliers.tsv")
        assert_true(os.path.exists(outlier_path), "Outliers file not found")
        outliers = pd.read_csv(outlier_path, sep="\t")["sample_id"].tolist()
        assert_true(bad_sample in outliers,
                    f"Corrupted sample {bad_sample} not flagged; outliers: {outliers}")

        expr_out = pd.read_csv(
            os.path.join(output_dir, "EUR_combat_int_expression.bed"), sep="\t")
        out_samples = [c for c in expr_out.columns
                       if c not in ["#chr", "start", "end", "phenotype_id"]]
        assert_true(bad_sample not in out_samples,
                    "Corrupted sample still in expression BED")
        assert_true(len(out_samples) == n_good,
                    f"Expected {n_good} samples after outlier removal, "
                    f"got {len(out_samples)}")

        hcp_df = pd.read_csv(os.path.join(output_dir, "EUR_hcp_factors.tsv"),
                             sep="\t")
        assert_true(len(hcp_df.columns) - 1 == n_good,
                    f"HCP file should have {n_good} sample columns, "
                    f"got {len(hcp_df.columns) - 1}")

        print(f"  PASSED: connectivity outlier removal flagged {outliers} "
              f"and excluded them from BED + HCP outputs")


# =============================================================================
# Test 4: refFlat generation
# =============================================================================

def test_reflat_generation():
    """Test GTF → refFlat conversion."""
    print("\n--- test_reflat_generation ---")
    sys.path.insert(0, SCRIPTS_DIR)
    from picard_qc import generate_refflat

    with tempfile.TemporaryDirectory() as tmpdir:
        gtf_path = os.path.join(tmpdir, "test.gtf")
        make_synthetic_gtf(gtf_path, n_genes=5, n_isoforms_per_gene=2)

        refflat_path = os.path.join(tmpdir, "test.refFlat")
        generate_refflat(gtf_path, refflat_path)

        assert_true(os.path.exists(refflat_path), "refFlat file not created")

        # Parse refFlat and verify
        with open(refflat_path) as f:
            lines = f.readlines()

        # 5 genes × 2 isoforms = 10 transcript lines
        assert_true(len(lines) == 10, f"Expected 10 refFlat lines, got {len(lines)}")

        # Verify format: 11 tab-delimited fields (UCSC refFlat)
        for line in lines:
            fields = line.strip().split("\t")
            assert_true(len(fields) == 11,
                        f"Expected 11 fields, got {len(fields)}: {line.strip()}")
            # geneName, name, chrom, strand, txStart, txEnd, cdsStart,
            # cdsEnd, exonCount, exonStarts, exonEnds
            assert_true(fields[2].startswith("chr"), f"Chrom should start with 'chr': {fields[2]}")
            assert_true(fields[3] in ["+", "-"], f"Strand should be + or -: {fields[3]}")
            assert_true(int(fields[8]) == 2, f"Expected 2 exons, got {fields[8]}")
            # exonStarts and exonEnds should end with comma
            assert_true(fields[9].endswith(","), "exonStarts should end with comma")
            assert_true(fields[10].endswith(","), "exonEnds should end with comma")

        print("  PASSED: GTF → refFlat conversion (5 genes, 10 transcripts, correct format)")


# =============================================================================
# Test 5: Full pipeline (Python components only)
# =============================================================================

def test_full_pipeline():
    """Test the full Python pipeline with synthetic data (pooling + QC parsing)."""
    print("\n--- test_full_pipeline ---")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create 2 cohorts, 1 ancestry
        cohort1_samples = [f"C1_S{i}" for i in range(10)]
        cohort2_samples = [f"C2_S{i}" for i in range(10)]
        all_samples = cohort1_samples + cohort2_samples

        # Ancestry map
        anc_map_path = os.path.join(tmpdir, "ancestry_map.tsv")
        make_synthetic_ancestry_map(anc_map_path, {
            ("cohort1", "EUR"): cohort1_samples,
            ("cohort2", "EUR"): cohort2_samples,
        })

        # Cohort dirs with unnorm BEDs
        for cohort, samples in [("cohort1", cohort1_samples), ("cohort2", cohort2_samples)]:
            cohort_dir = os.path.join(tmpdir, cohort)
            unnorm_dir = os.path.join(cohort_dir, "output", "unnorm")
            os.makedirs(unnorm_dir)
            make_synthetic_bed(
                os.path.join(unnorm_dir, "expression.bed"),
                n_genes=100, samples=samples, seed=hash(cohort) % 1000)

        # QC metrics
        qc_path = os.path.join(tmpdir, "all_qc_metrics.tsv")
        make_synthetic_qc(qc_path, all_samples, n_metrics=10, seed=42)

        # Pool expression
        pooled_dir = os.path.join(tmpdir, "pooled")
        cmd = [
            sys.executable, os.path.join(SCRIPTS_DIR, "pool_expression_within_ancestry.py"),
            "--ancestry-map", anc_map_path,
            "--cohort-dirs", f"cohort1={tmpdir}/cohort1", f"cohort2={tmpdir}/cohort2",
            "--output-dir", pooled_dir,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert_true(result.returncode == 0, f"Pooling failed:\n{result.stderr}")

        # Verify pooled output
        eur_path = os.path.join(pooled_dir, "EUR_pooled_expression.bed")
        assert_true(os.path.exists(eur_path), "Pooled EUR BED not found")
        eur_df = pd.read_csv(eur_path, sep="\t")
        eur_samples = [c for c in eur_df.columns
                       if c not in ["#chr", "start", "end", "phenotype_id"]]
        assert_true(len(eur_samples) == 20, f"Expected 20 samples, got {len(eur_samples)}")

        # Verify the output is usable as tensorQTL input format
        # (BED with #chr, start, end, phenotype_id + sample columns)
        assert_true(eur_df.columns[0] == "#chr", "First column should be #chr")
        assert_true(eur_df.columns[3] == "phenotype_id", "4th column should be phenotype_id")

        print("  PASSED: full Python pipeline (pooling → tensorQTL-ready BED format)")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Test suite for HCP pipeline")
    parser.add_argument("--test", help="Run a single test by name", default=None)
    args = parser.parse_args()

    tests = {
        "test_picard_qc": test_picard_qc,
        "test_pool_expression": test_pool_expression,
        "test_combat_normalize_hcp": test_combat_normalize_hcp,
        "test_connectivity_outliers": test_connectivity_outliers,
        "test_reflat_generation": test_reflat_generation,
        "test_full_pipeline": test_full_pipeline,
    }

    if args.test:
        if args.test not in tests:
            print(f"Unknown test: {args.test}")
            print(f"Available: {', '.join(tests.keys())}")
            sys.exit(1)
        tests[args.test]()
        print("\nAll requested tests passed.")
        return

    passed = 0
    failed = 0
    for name, test_fn in tests.items():
        try:
            test_fn()
            passed += 1
        except Exception as e:
            print(f"  FAILED: {e}")
            failed += 1

    print(f"\n{'='*60}")
    print(f"Results: {passed} passed, {failed} failed")
    if failed > 0:
        sys.exit(1)
    print("All tests passed.")


if __name__ == "__main__":
    main()
