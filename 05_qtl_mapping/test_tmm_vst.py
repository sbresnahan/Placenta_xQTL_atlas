#!/usr/bin/env python3
"""
test_tmm_vst.py — Regression tests for the TMM -> VST bridge used by the
PsychENCODE/isoTWAS expression normalization schema.

The bridge (tmm_vst() in normalize_expression_hcp.R and
normalize_modalities.R) injects edgeR TMM effective library sizes as DESeq2
size factors and runs vst() on the RAW counts:

    eff_lib <- colSums(counts) * edgeR::calcNormFactors(counts, "TMM")
    sizeFactors(dds) <- eff_lib / exp(mean(log(eff_lib)))
    vst(dds, blind = TRUE)

The failure mode this guards against: pre-dividing counts by the TMM factors
and handing the result to vst() with default size factors — a double
normalization that silently corrupts the scale.

Tests (R scripts run as subprocesses on synthetic fixtures):
  T1: bridge correctness — size factors equal TMM effective library sizes
      scaled to geometric mean 1; pipeline VST output matches an independent
      manual edgeR+DESeq2 computation; no NA/Inf
  T2: double-normalization guard — pipeline output differs MATERIALLY from
      the wrong "pre-divided + default size factors" computation
  T3: TMM recovers a planted composition bias (cohort2 inflated in 40% of
      genes -> TMM factors shift in the compensating direction)
  T4: feature-set parity — the VST output feature set is exactly the
      TPM-filtered set (counts BED emits all features; low-TPM genes are
      dropped), and the zero-library-size guard errors clearly
  T5: normalize_modalities.R --counts-input path matches the manual bridge
  T6: VST mode appends k-1 cohort dummies to the HCP known-covariate matrix
      (skipped gracefully when Rhcpp is unavailable)

Run:
  python3 test_tmm_vst.py        # or: pytest test_tmm_vst.py
"""

import os
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EXPR_SCRIPT = HERE / "normalize_expression_hcp.R"
MOD_SCRIPT = HERE / "normalize_modalities.R"
RSCRIPT_SIF = HERE.parent / "bin" / "Rscript_sif"

PASS = 0
FAIL = 0


def report(name, ok, detail=""):
    global PASS, FAIL
    status = "PASS" if ok else "FAIL"
    if ok:
        PASS += 1
    else:
        FAIL += 1
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")


def run_rcode(code, cwd):
    """Run an R snippet; return CompletedProcess."""
    rfile = Path(cwd) / "verify.R"
    rfile.write_text(code)
    return subprocess.run([str(RSCRIPT_SIF), "--vanilla", str(rfile)],
                          capture_output=True, text=True, cwd=cwd)


# ---------------------------------------------------------------------------
# Fixture: 2 cohorts x 20 samples, 300 genes, planted composition bias
# (cohort2: 40% of genes 5x). 3 extra genes are near-zero (TPM-filtered out).
# ---------------------------------------------------------------------------
N1, N2, N_GENES, N_LOW = 20, 20, 300, 3
C1 = [f"C1_S{i}" for i in range(N1)]
C2 = [f"C2_S{i}" for i in range(N2)]
SAMPLES = C1 + C2
N_ALL = N_GENES + N_LOW


def write_expression_fixture(root):
    """Write pooled TPM BED + counts BED + QC + ancestry map. Returns dict."""
    rng = np.random.RandomState(2026)
    counts = np.empty((N_ALL, len(SAMPLES)))
    for i in range(N_ALL):
        mu = rng.gamma(2.0, 300.0, size=len(SAMPLES)) + 100
        if i < int(0.4 * N_GENES):
            mu[N1:] *= 5.0  # planted composition bias in cohort2
        counts[i] = rng.negative_binomial(10, 10 / (10 + mu))
    # Near-zero genes: detected in 2/40 samples -> TPM filter drops them
    counts[N_GENES:, :] = 0.0
    counts[N_GENES:, 0] = 5.0
    counts[N_GENES:, N1] = 7.0

    gene_ids = [f"GENE{i:05d}" for i in range(N_ALL)]
    meta = pd.DataFrame({
        "#chr": "chr1",
        "start": [i * 1000 for i in range(N_ALL)],
        "end": [i * 1000 + 500 for i in range(N_ALL)],
        "phenotype_id": gene_ids,
    })
    pd.concat([meta, pd.DataFrame(counts, columns=SAMPLES)], axis=1).to_csv(
        Path(root) / "EUR_pooled_expression_counts.bed",
        sep="\t", index=False, float_format="%g")
    tpm = counts / counts.sum(axis=0, keepdims=True) * 1e6
    pd.concat([meta, pd.DataFrame(tpm, columns=SAMPLES)], axis=1).to_csv(
        Path(root) / "EUR_pooled_expression.bed",
        sep="\t", index=False, float_format="%g")

    qc = pd.DataFrame(rng.normal(size=(len(SAMPLES), 5)),
                      index=SAMPLES, columns=[f"QC{j}" for j in range(5)])
    qc.index.name = "sample_id"
    qc.to_csv(Path(root) / "qc.tsv", sep="\t")
    pd.DataFrame({"sample_id": SAMPLES, "assigned_ancestry": "EUR",
                  "cohort": ["cohort1"] * N1 + ["cohort2"] * N2}
                 ).to_csv(Path(root) / "anc_map.tsv", sep="\t", index=False)
    return {"counts": counts, "gene_ids": gene_ids}


# Independent re-implementation of the bridge for verification. Prints
# key=value lines consumed by the Python assertions.
VERIFY_R = textwrap.dedent(r"""
    suppressPackageStartupMessages({library(edgeR); library(DESeq2)})
    read_bed <- function(p) {
      hdr <- strsplit(readLines(p, n = 1), "\t", fixed = TRUE)[[1]]
      cc <- rep("numeric", length(hdr)); cc[1:4] <- c("character", "integer", "integer", "character")
      read.delim(p, sep = "\t", check.names = FALSE, colClasses = cc)
    }
    counts_bed <- read_bed("EUR_pooled_expression_counts.bed")
    tpm_bed    <- read_bed("EUR_pooled_expression.bed")
    out_bed    <- read_bed(file.path("OUTDIR", "EUR_vst_expression.bed"))
    samples <- setdiff(colnames(counts_bed), c("#chr", "start", "end", "phenotype_id"))

    # Reproduce the TPM detection filter (devBrain: TPM > 0.1 in > 25%)
    tpm_mat <- as.matrix(tpm_bed[, samples])
    detect <- colMeans(t(tpm_mat) > 0.1)          # per gene
    novar <- apply(tpm_mat, 1, function(x) length(unique(x)) <= 1)
    keep <- tpm_bed$phenotype_id[(detect > 0.25) & !novar]

    gidx <- match(out_bed$phenotype_id, counts_bed$phenotype_id)
    counts <- t(as.matrix(counts_bed[gidx, samples]))   # samples x genes

    # Manual bridge (independent of the pipeline's tmm_vst)
    cg <- t(counts)
    lib <- colSums(cg)
    tmm <- edgeR::calcNormFactors(cg, method = "TMM")
    eff <- lib * tmm
    sf <- eff / exp(mean(log(eff)))
    dds <- DESeq2::DESeqDataSetFromMatrix(round(cg),
          colData = data.frame(row.names = colnames(cg)), design = ~ 1)
    DESeq2::sizeFactors(dds) <- sf
    manual <- t(SummarizedExperiment::assay(
      DESeq2::vst(dds, blind = TRUE, nsub = min(1000L, nrow(cg)))))

    pipe <- t(as.matrix(out_bed[, samples]))   # samples x genes, like manual
    cat("SF_GEOMEAN", exp(mean(log(sf))), "\n")
    cat("VST_MAXDIFF", max(abs(pipe - manual)), "\n")
    cat("N_NA_INF", sum(!is.finite(pipe)), "\n")

    # WRONG approach: pre-divide by TMM factors, default size factors
    wrong_counts <- sweep(cg, 2, tmm, "/")
    dds_w <- DESeq2::DESeqDataSetFromMatrix(round(wrong_counts),
            colData = data.frame(row.names = colnames(cg)), design = ~ 1)
    wrong <- t(SummarizedExperiment::assay(
      DESeq2::vst(dds_w, blind = TRUE, nsub = min(1000L, nrow(cg)))))
    cat("WRONG_MAXDIFF", max(abs(pipe - wrong)), "\n")

    # TMM factor direction per cohort
    cat("TMM_C1_MEAN", mean(tmm[1:20]), "\n")
    cat("TMM_C2_MEAN", mean(tmm[21:40]), "\n")

    # Feature-set parity
    cat("FEATURE_PARITY", identical(sort(out_bed$phenotype_id), sort(keep)), "\n")
""")


def test_bridge_correctness():
    """T1+T2+T3+T4(parity): run the pipeline in VST mode, verify against an
    independent manual computation."""
    print("\n=== T1-T4: TMM->VST bridge (normalize_expression_hcp.R) ===")
    if not all(shutil_which_r()):
        report("Rscript available", False, "bin/Rscript_sif or singularity not found")
        return
    tmp = Path(tempfile.mkdtemp())
    write_expression_fixture(tmp)
    outdir = tmp / "out"
    cmd = [str(RSCRIPT_SIF), str(EXPR_SCRIPT),
           "--expression", str(tmp / "EUR_pooled_expression.bed"),
           "--counts", str(tmp / "EUR_pooled_expression_counts.bed"),
           "--mode", "vst",
           "--qc-metrics", str(tmp / "qc.tsv"),
           "--ancestry-map", str(tmp / "anc_map.tsv"),
           "--ancestry", "EUR", "--k", "3", "--outlier-z", "-999",
           "--skip-hcp", "--output-dir", str(outdir)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    report("pipeline VST run exits 0", r.returncode == 0,
           (r.stdout + r.stderr)[-600:])
    if r.returncode != 0:
        return
    report("vst BED written", (outdir / "EUR_vst_expression.bed").exists())

    rv = run_rcode(VERIFY_R.replace("OUTDIR", str(outdir)), tmp)
    kv = {}
    for line in rv.stdout.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].isupper():
            kv[parts[0]] = parts[1]
    report("verification R ran", rv.returncode == 0, rv.stderr[-600:])

    # T1a: size factors are TMM effective library sizes, geometric mean 1
    report("size factors geometric mean == 1",
           abs(float(kv.get("SF_GEOMEAN", "nan")) - 1) < 1e-8,
           f"got {kv.get('SF_GEOMEAN')}")
    # T1b: pipeline VST == manual bridge
    report("pipeline VST matches manual edgeR+DESeq2 bridge",
           float(kv.get("VST_MAXDIFF", "nan")) < 1e-6,
           f"max abs diff {kv.get('VST_MAXDIFF')}")
    # T1c: no NA/Inf
    report("no NA/Inf in VST output", kv.get("N_NA_INF") == "0",
           f"got {kv.get('N_NA_INF')}")
    # T2: double-normalization guard — wrong result must differ materially
    report("output differs from pre-divided double-normalization",
           float(kv.get("WRONG_MAXDIFF", "0")) > 0.05,
           f"max abs diff {kv.get('WRONG_MAXDIFF')}")
    # T3: TMM recovers the planted composition bias. cohort2 is inflated in
    # 40% of genes, so its library share is dominated by those genes and the
    # 60% non-inflated majority appears relatively DEPLETED vs the reference.
    # TMM normalizes to the majority, so cohort2 gets the SMALLER factors
    # (counts scaled UP relative to library size), which preserves the
    # planted inflation instead of dividing it away.
    tmm1, tmm2 = float(kv.get("TMM_C1_MEAN", "nan")), float(kv.get("TMM_C2_MEAN", "nan"))
    report("TMM factors shift against the planted bias (C2 < C1)",
           tmm2 < tmm1, f"C1 mean {tmm1:.4f}, C2 mean {tmm2:.4f}")
    # T4: feature-set parity with the TPM filter
    report("feature set identical to TPM-filter path",
           kv.get("FEATURE_PARITY") == "TRUE", f"got {kv.get('FEATURE_PARITY')}")


def shutil_which_r():
    import shutil
    return [RSCRIPT_SIF.exists() and shutil.which("singularity") is not None]


def test_zero_library_guard():
    """T4b: an all-zero sample must trigger the clear library-size guard."""
    print("\n=== T4b: zero-library-size guard ===")
    tmp = Path(tempfile.mkdtemp())
    fx = write_expression_fixture(tmp)
    # Zero out one sample in BOTH BEDs
    for name in ["EUR_pooled_expression_counts.bed", "EUR_pooled_expression.bed"]:
        df = pd.read_csv(tmp / name, sep="\t")
        df["C1_S0"] = 0.0
        df.to_csv(tmp / name, sep="\t", index=False)
    cmd = [str(RSCRIPT_SIF), str(EXPR_SCRIPT),
           "--expression", str(tmp / "EUR_pooled_expression.bed"),
           "--counts", str(tmp / "EUR_pooled_expression_counts.bed"),
           "--mode", "vst",
           "--qc-metrics", str(tmp / "qc.tsv"),
           "--ancestry-map", str(tmp / "anc_map.tsv"),
           "--ancestry", "EUR", "--k", "3", "--outlier-z", "-999",
           "--skip-hcp", "--output-dir", str(tmp / "out")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    report("zero-library sample fails loudly", r.returncode != 0)
    report("guard message is clear",
           "zero/non-finite library size" in (r.stdout + r.stderr),
           (r.stdout + r.stderr)[-400:])


def test_modalities_vst_path():
    """T5: normalize_modalities.R --counts-input matches the manual bridge."""
    print("\n=== T5: normalize_modalities.R VST path ===")
    tmp = Path(tempfile.mkdtemp())
    rng = np.random.RandomState(55)
    n_tx = 250
    tx_ids = [f"GENE{i:05d}__TX{t}" for i in range(n_tx // 2) for t in range(2)]
    counts = np.empty((n_tx, len(SAMPLES)))
    for i in range(n_tx):
        mu = rng.gamma(2.0, 200.0, size=len(SAMPLES)) + 60
        if i % 3 == 0:
            mu[N1:] *= 4.0
        counts[i] = rng.negative_binomial(8, 8 / (8 + mu))
    meta = pd.DataFrame({"#chr": "chr1",
                         "start": [i * 500 for i in range(n_tx)],
                         "end": [i * 500 + 200 for i in range(n_tx)],
                         "phenotype_id": tx_ids})
    pd.concat([meta, pd.DataFrame(counts, columns=SAMPLES)], axis=1).to_csv(
        tmp / "EUR_isoform_expression_counts_pooled.bed",
        sep="\t", index=False, float_format="%g")
    tpm = counts / counts.sum(axis=0, keepdims=True) * 1e6
    pd.concat([meta, pd.DataFrame(tpm, columns=SAMPLES)], axis=1).to_csv(
        tmp / "EUR_isoform_expression_pooled.bed",
        sep="\t", index=False, float_format="%g")
    pd.DataFrame({"sample_id": SAMPLES, "assigned_ancestry": "EUR",
                  "cohort": ["cohort1"] * N1 + ["cohort2"] * N2}
                 ).to_csv(tmp / "anc_map.tsv", sep="\t", index=False)

    outdir = tmp / "out"
    cmd = [str(RSCRIPT_SIF), str(MOD_SCRIPT),
           "--input", str(tmp / "EUR_isoform_expression_pooled.bed"),
           "--counts-input", str(tmp / "EUR_isoform_expression_counts_pooled.bed"),
           "--ancestry-map", str(tmp / "anc_map.tsv"),
           "--ancestry", "EUR", "--modality", "isoform_expression",
           "--output-dir", str(outdir)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    report("modalities VST run exits 0", r.returncode == 0,
           (r.stdout + r.stderr)[-600:])
    if r.returncode != 0:
        return
    report("isoform_expression_vst.bed written",
           (outdir / "EUR_isoform_expression_vst.bed").exists())

    code = textwrap.dedent(r"""
        suppressPackageStartupMessages({library(edgeR); library(DESeq2)})
        read_bed <- function(p) {
          hdr <- strsplit(readLines(p, n = 1), "\t", fixed = TRUE)[[1]]
          cc <- rep("numeric", length(hdr)); cc[1:4] <- c("character", "integer", "integer", "character")
          read.delim(p, sep = "\t", check.names = FALSE, colClasses = cc)
        }
        counts_bed <- read_bed("EUR_isoform_expression_counts_pooled.bed")
        out_bed <- read_bed(file.path("OUTDIR", "EUR_isoform_expression_vst.bed"))
        samples <- setdiff(colnames(counts_bed), c("#chr", "start", "end", "phenotype_id"))
        gidx <- match(out_bed$phenotype_id, counts_bed$phenotype_id)
        cg <- as.matrix(counts_bed[gidx, samples])     # features x samples
        lib <- colSums(cg)
        sf <- (lib * edgeR::calcNormFactors(cg, method = "TMM"))
        sf <- sf / exp(mean(log(sf)))
        dds <- DESeq2::DESeqDataSetFromMatrix(round(cg),
              colData = data.frame(row.names = colnames(cg)), design = ~ 1)
        DESeq2::sizeFactors(dds) <- sf
        manual <- t(SummarizedExperiment::assay(
          DESeq2::vst(dds, blind = TRUE, nsub = min(1000L, nrow(cg)))))
        pipe <- t(as.matrix(out_bed[, samples]))
        cat("VST_MAXDIFF", max(abs(pipe - manual)), "\n")
        cat("N_NA_INF", sum(!is.finite(pipe)), "\n")
    """).replace("OUTDIR", str(outdir))
    rv = run_rcode(code, tmp)
    kv = {}
    for line in rv.stdout.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].isupper():
            kv[parts[0]] = parts[1]
    report("modalities VST matches manual bridge",
           float(kv.get("VST_MAXDIFF", "nan")) < 1e-6,
           f"max abs diff {kv.get('VST_MAXDIFF')} (stderr: {rv.stderr[-300:]})")
    report("modalities VST: no NA/Inf", kv.get("N_NA_INF") == "0")


def test_hcp_cohort_dummies():
    """T6: VST mode appends k-1 cohort dummies to the HCP Z matrix.
    Skipped gracefully when Rhcpp is unavailable."""
    print("\n=== T6: cohort dummies in HCP Z (VST mode) ===")
    tmp = Path(tempfile.mkdtemp())
    write_expression_fixture(tmp)
    cmd = [str(RSCRIPT_SIF), str(EXPR_SCRIPT),
           "--expression", str(tmp / "EUR_pooled_expression.bed"),
           "--counts", str(tmp / "EUR_pooled_expression_counts.bed"),
           "--mode", "vst",
           "--qc-metrics", str(tmp / "qc.tsv"),
           "--ancestry-map", str(tmp / "anc_map.tsv"),
           "--ancestry", "EUR", "--k", "3", "--outlier-z", "-999",
           "--output-dir", str(tmp / "out")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 and ("there is no package called" in r.stderr
                              or "could not find function" in r.stderr):
        report("T6 SKIPPED: Rhcpp not available", True)
        return
    report("HCP run exits 0", r.returncode == 0, (r.stdout + r.stderr)[-500:])
    report("cohort dummies added to HCP Z (k-1 = 1, reference cohort1)",
           "HCP known covariates: + 1 cohort dummies (reference: cohort1)"
           in r.stdout)
    report("HCP factors written",
           (tmp / "out" / "EUR_hcp_factors.tsv").exists())


if __name__ == "__main__":
    test_bridge_correctness()
    test_zero_library_guard()
    test_modalities_vst_path()
    test_hcp_cohort_dummies()
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)
