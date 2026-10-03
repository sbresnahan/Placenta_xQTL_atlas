#!/usr/bin/env python3
"""
test_normalize_modalities.py — Synthetic test suite for the cross-cohort
pooling + within-cohort INT normalization helpers (roadmap step 4).

Tests:
  T1: detection filter + QN/INT correctness (R script, single cohort)
  T2: pooling phenotype_id intersection (Python pooler)
  T3: namespaced pass-through for pre-pooled modalities (Python pooler)
  T4: within-cohort INT end-to-end on synthetic cohort-shifted data (R script)
  T5: single-sample-cohort robustness (R script)

Run:
  python3 test_normalize_modalities.py
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
POOLER = HERE / "pool_modalities_within_ancestry.py"
R_SCRIPT = HERE / "normalize_modalities.R"
RSCRIPT_SIF = HERE.parent / "bin" / "Rscript_sif"

BED_META = ["#chr", "start", "end", "phenotype_id"]

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


def run(cmd, **kw):
    """Run a command, return CompletedProcess. Raise on failure."""
    return subprocess.run(cmd, capture_output=True, text=True, check=True, **kw)


# ---------------------------------------------------------------------------
# T1: detection filter + QN/INT correctness (R, single cohort)
# ---------------------------------------------------------------------------

def test_detection_filter_and_int():
    print("\n=== T1: detection filter + QN/INT (single cohort) ===")
    tmp = Path(tempfile.mkdtemp())

    # 202 features x 40 samples, ONE cohort (within-cohort INT over the
    # whole stratum -> output is exactly QN + INT). Proportion-scale values
    # including exact 0s.
    # (QN needs a realistic feature count to be smooth — with only a handful
    # of features every sample collapses onto the same few quantile means and
    # the per-feature INT sd is artificially deflated by ties.)
    rng = np.random.RandomState(1)
    n_samp = 40
    n_bg = 200
    samples = [f"s{i}" for i in range(n_samp)]
    feats = {f"BG{i}__f{i}": rng.rand(n_samp).tolist() for i in range(n_bg)}
    feats["G1__ok1"] = rng.rand(n_samp).tolist()
    # 75% exact zeros but detected in all samples -> KEPT (zeros are
    # real PSI values; a zeros fraction is not a detection filter)
    feats["G2__zeros"] = [0.0] * 30 + rng.rand(10).tolist()
    # detected (non-NA) in only 10/40 = 25% < 40% -> dropped (devBrain
    # detection filter for proportion modalities)
    feats["G3__lowdetect"] = rng.rand(10).tolist() + [np.nan] * 30
    feat_names = list(feats.keys())
    bed = pd.DataFrame({
        "#chr": ["chr1"] * len(feat_names),
        "start": [100 + 10 * i for i in range(len(feat_names))],
        "end": [101 + 10 * i for i in range(len(feat_names))],
        "phenotype_id": feat_names,
        **{s: [feats[f][j] for f in feat_names] for j, s in enumerate(samples)},
    })
    bed_path = tmp / "test_detect.bed"
    bed.to_csv(bed_path, sep="\t", index=False, float_format="%g")

    anc = pd.DataFrame({
        "sample_id": samples,
        "assigned_ancestry": ["EUR"] * n_samp,
        "cohort": ["cA"] * n_samp,  # single cohort -> INT over full stratum
    })
    anc_path = tmp / "ancestry.tsv"
    anc.to_csv(anc_path, sep="\t", index=False)

    out_dir = tmp / "out"
    out_dir.mkdir()

    r = run([
        str(RSCRIPT_SIF), str(R_SCRIPT),
        "--input", str(bed_path),
        "--ancestry-map", str(anc_path),
        "--ancestry", "EUR",
        "--modality", "splicing",
        "--output-dir", str(out_dir),
    ])
    if r.returncode != 0:
        report("T1 Rscript runs", False, r.stderr[-500:])
        return

    out_bed = out_dir / "EUR_splicing_int.bed"
    if not out_bed.exists():
        report("T1 output BED exists", False, str(out_bed))
        return
    report("T1 Rscript runs + output exists", True)

    df = pd.read_csv(out_bed, sep="\t")
    out_feats = list(df["phenotype_id"])
    report("T1 low-detection feature dropped (<40% detected)",
           "G3__lowdetect" not in out_feats, f"features: {out_feats}")
    report("T1 zero-heavy feature kept (zeros are real values)",
           "G2__zeros" in out_feats, f"features: {out_feats}")

    data = df.drop(columns=BED_META)
    report("T1 no NaN/Inf in output",
           not data.isna().any().any() and not np.isinf(data.values).any())

    # Single cohort => within-cohort INT over the full stratum => output is
    # exactly QN+INT.
    # Check the continuous features are ~N(0,1) (the zero-heavy feature has
    # tied ranks, so it is checked only for presence above).
    cont = df[df["phenotype_id"].str.startswith(("BG", "G1__"))]
    cont_data = cont.drop(columns=BED_META)
    means = cont_data.mean(axis=1).abs()
    sds = cont_data.std(axis=1)
    report("T1 INT means ~0", (means < 0.05).all(),
           f"max|mean|={means.max():.4f}")
    report("T1 INT sds ~1", ((sds > 0.9) & (sds < 1.1)).all(),
           f"sd range=[{sds.min():.3f},{sds.max():.3f}]")


# ---------------------------------------------------------------------------
# T2: pooling phenotype_id intersection (Python)
# ---------------------------------------------------------------------------

def test_pooling_intersection():
    print("\n=== T2: pooling phenotype_id intersection ===")
    tmp = Path(tempfile.mkdtemp())

    # Two cohorts with partial feature overlap
    # Cohort A: genes G1, G2, G3 ; samples a1, a2
    # Cohort B: genes G2, G3, G4 ; samples b1, b2
    # Intersection: G2, G3
    def make_bed(features, samples):
        n = len(features)
        return pd.DataFrame({
            "#chr": ["chr1"] * n,
            "start": list(range(100, 100 + 10 * n, 10)),
            "end": [s + 1 for s in range(100, 100 + 10 * n, 10)],
            "phenotype_id": features,
            **{s: np.random.RandomState(hash(s) % 2**31).rand(n).tolist()
               for s in samples},
        })

    cohortA = tmp / "cohortA"
    cohortB = tmp / "cohortB"
    for c in (cohortA, cohortB):
        (c / "output" / "unnorm").mkdir(parents=True)

    make_bed(["G1", "G2", "G3"], ["a1", "a2"]).to_csv(
        cohortA / "output" / "unnorm" / "alt_TSS.bed", sep="\t", index=False)
    make_bed(["G2", "G3", "G4"], ["b1", "b2"]).to_csv(
        cohortB / "output" / "unnorm" / "alt_TSS.bed", sep="\t", index=False)

    anc = pd.DataFrame({
        "sample_id": ["a1", "a2", "b1", "b2"],
        "assigned_ancestry": ["EUR"] * 4,
        "cohort": ["cohortA", "cohortA", "cohortB", "cohortB"],
    })
    anc_path = tmp / "ancestry.tsv"
    anc.to_csv(anc_path, sep="\t", index=False)

    out_dir = tmp / "pooled"
    out_dir.mkdir()

    r = run([
        "python3", str(POOLER),
        "--ancestry-map", str(anc_path),
        "--cohort-dirs", f"cohortA={cohortA}", f"cohortB={cohortB}",
        "--modality", "alt_TSS",
        "--output-dir", str(out_dir),
    ])

    out_bed = out_dir / "EUR_alt_TSS_pooled.bed"
    if not out_bed.exists():
        report("T2 output exists", False, str(out_bed))
        return
    report("T2 pooler runs + output exists", True)

    df = pd.read_csv(out_bed, sep="\t")
    features = list(df["phenotype_id"])
    samples = [c for c in df.columns if c not in BED_META]

    # Intersection should be G2, G3 (sorted)
    report("T2 feature intersection correct",
           sorted(features) == ["G2", "G3"],
           f"got {sorted(features)}")
    # All 4 samples present
    report("T2 all samples present",
           sorted(samples) == ["a1", "a2", "b1", "b2"],
           f"got {sorted(samples)}")


# ---------------------------------------------------------------------------
# T3: namespaced pass-through (Python)
# ---------------------------------------------------------------------------

def test_prepooled_passthrough():
    print("\n=== T3: namespaced pass-through ===")
    tmp = Path(tempfile.mkdtemp())

    # Simulate harmonize output: <root>/<ancestry>/<modality>/unnorm/<modality>.bed
    root = tmp / "harmonize"
    mod_dir = root / "EUR" / "splicing" / "unnorm"
    mod_dir.mkdir(parents=True)

    bed = pd.DataFrame({
        "#chr": ["chr1"] * 3,
        "start": [100, 200, 300],
        "end": [101, 201, 301],
        "phenotype_id": ["G1__j1", "G1__j2", "G2__j3"],
        "cohortA_a1": [0.2, 0.3, 0.5],
        "cohortA_a2": [0.1, 0.4, 0.5],
        "cohortB_b1": [0.3, 0.2, 0.5],
        "cohortB_b2": [0.25, 0.35, 0.4],
    })
    bed.to_csv(mod_dir / "splicing.bed", sep="\t", index=False, float_format="%g")

    anc = pd.DataFrame({
        "sample_id": ["a1", "a2", "b1", "b2"],
        "assigned_ancestry": ["EUR"] * 4,
        "cohort": ["cohortA", "cohortA", "cohortB", "cohortB"],
    })
    anc_path = tmp / "ancestry.tsv"
    anc.to_csv(anc_path, sep="\t", index=False)

    out_dir = tmp / "pooled"
    out_dir.mkdir()

    r = run([
        "python3", str(POOLER),
        "--ancestry-map", str(anc_path),
        "--modality", "splicing",
        "--pre-pooled-dir", str(root),
        "--output-dir", str(out_dir),
    ])

    out_bed = out_dir / "EUR_splicing_pooled.bed"
    label_file = out_dir / "EUR_splicing_cohort_labels.tsv"
    if not out_bed.exists():
        report("T3 output BED exists", False, str(out_bed))
        return
    if not label_file.exists():
        report("T3 cohort-label sidecar exists", False, str(label_file))
        return
    report("T3 pooler runs + BED + sidecar exist", True)

    # BED should be identical to input (pass-through)
    df_out = pd.read_csv(out_bed, sep="\t")
    report("T3 BED pass-through identical",
           df_out.shape == bed.shape and
           list(df_out.columns) == list(bed.columns),
           f"shape {df_out.shape} vs {bed.shape}")

    # Cohort labels parsed correctly from namespaced IDs
    labels = pd.read_csv(label_file, sep="\t")
    label_map = dict(zip(labels["sample_id"], labels["cohort"]))
    expected = {
        "cohortA_a1": "cohortA", "cohortA_a2": "cohortA",
        "cohortB_b1": "cohortB", "cohortB_b2": "cohortB",
    }
    report("T3 cohort labels parsed from namespaced IDs",
           label_map == expected, f"got {label_map}")


# ---------------------------------------------------------------------------
# T4: within-cohort INT end-to-end on synthetic cohort-shifted data (R)
# ---------------------------------------------------------------------------

def test_withincohort_int_endtoend():
    print("\n=== T4: within-cohort INT end-to-end (cohort-shift removal) ===")
    tmp = Path(tempfile.mkdtemp())

    # Synthetic data: 50 features x 40 samples, 2 cohorts (20 each).
    # Inject a cohort-specific mean shift so the within-cohort INT has
    # something to remove.
    rng = np.random.RandomState(42)
    n_feat, n_samp = 50, 40
    # Proportions in [0.2, 0.5]. Cohort B gets a +0.35 shift on the FIRST HALF
    # of features only (-> [0.55, 0.85], no boundary clipping). The shift must
    # be feature-specific: a uniform all-feature shift is invisible to
    # within-sample ranks and is erased by QN itself.
    base = 0.2 + 0.3 * rng.rand(n_feat, n_samp)
    affected = np.arange(n_feat // 2)
    base[affected, 20:] += 0.35

    samples = [f"s{i}" for i in range(n_samp)]
    cohorts = ["cA"] * 20 + ["cB"] * 20

    bed = pd.DataFrame({
        "#chr": ["chr1"] * n_feat,
        "start": list(range(100, 100 + 10 * n_feat, 10)),
        "end": list(range(101, 101 + 10 * n_feat, 10)),
        "phenotype_id": [f"G{i}__e{i}" for i in range(n_feat)],
    })
    for j, s in enumerate(samples):
        bed[s] = base[:, j]
    bed_path = tmp / "test_int.bed"
    bed.to_csv(bed_path, sep="\t", index=False, float_format="%g")

    # Paired design: run the SAME input twice.
    #   Run A: single-cohort ancestry map -> pooled INT (shift preserved)
    #   Run B: two-cohort ancestry map -> within-cohort INT (shift removed
    #          by construction)
    anc_single = pd.DataFrame({
        "sample_id": samples,
        "assigned_ancestry": ["EUR"] * n_samp,
        "cohort": ["cA"] * n_samp,
    })
    anc_single_path = tmp / "ancestry_single.tsv"
    anc_single.to_csv(anc_single_path, sep="\t", index=False)

    anc_two = pd.DataFrame({
        "sample_id": samples,
        "assigned_ancestry": ["EUR"] * n_samp,
        "cohort": cohorts,
    })
    anc_two_path = tmp / "ancestry_two.tsv"
    anc_two.to_csv(anc_two_path, sep="\t", index=False)

    datas = {}
    for tag, anc_path in [("pooled", anc_single_path),
                          ("withincohort", anc_two_path)]:
        out_dir = tmp / f"out_{tag}"
        out_dir.mkdir()
        r = run([
            str(RSCRIPT_SIF), str(R_SCRIPT),
            "--input", str(bed_path),
            "--ancestry-map", str(anc_path),
            "--ancestry", "EUR",
            "--modality", "splicing",
            "--output-dir", str(out_dir),
        ])
        if r.returncode != 0:
            report(f"T4 Rscript runs ({tag})", False, r.stderr[-800:])
            return
        datas[tag] = pd.read_csv(out_dir / "EUR_splicing_int.bed",
                                 sep="\t").drop(columns=BED_META)
    report("T4 Rscript runs (both)", True)

    from sklearn.decomposition import PCA
    ca_idx, cb_idx = np.arange(20), np.arange(20, 40)

    # Cohort separation on PC1: pooled INT vs within-cohort INT
    centroids = {}
    for tag, data in datas.items():
        pca = PCA(n_components=2).fit_transform(data.T.values)
        centroids[tag] = abs(pca[ca_idx, 0].mean() - pca[cb_idx, 0].mean())
    report("T4 cohort separation removed by within-cohort INT (PC1 gap)",
           centroids["withincohort"] < centroids["pooled"] * 0.5,
           f"pooled={centroids['pooled']:.3f}, "
           f"withincohort={centroids['withincohort']:.3f}")

    # Within-cohort INT guarantee: every feature is ~mean-0 within each
    # cohort, so the injected shift is removed BY CONSTRUCTION (no
    # batch-correction step). Not exactly 0 here: QN assigns the same
    # rank-mean value to a feature whenever it ties on within-sample rank
    # across samples (frequent with only 50 features), and tie-averaged
    # ranks break the exact symmetry of the INT quantiles. The residual is
    # ~1e-2 at this toy feature count and shrinks as features >> samples.
    # The pooled-INT arm must retain a positive gap, else the test is
    # vacuous.
    wc = datas["withincohort"]
    d_pooled = (datas["pooled"].iloc[affected][samples[:20]].mean(axis=1)
                - datas["pooled"].iloc[affected][samples[20:]].mean(axis=1)).abs().mean()
    gap = (wc.iloc[affected][samples[:20]].mean(axis=1)
           - wc.iloc[affected][samples[20:]].mean(axis=1)).abs().max()
    report("T4 per-feature cohort mean gap removed by construction",
           gap < 0.05 and d_pooled > 0.1,
           f"pooled gap={d_pooled:.3f}, within-cohort max|gap|={gap:.2e}")

    # Per-cohort sd equals the rank-INT quantile sd for that cohort's n
    # (identical across equal-n cohorts here).
    sd_ca = wc[samples[:20]].values.std()
    sd_cb = wc[samples[20:]].values.std()
    report("T4 per-cohort sd matches across equal-n cohorts",
           abs(sd_ca - sd_cb) < 0.05,
           f"sd cA={sd_ca:.3f}, sd cB={sd_cb:.3f}")

    # Pooled output is approximately but not exactly N(0,1) per feature -
    # loose sanity bounds only.
    means = wc.mean(axis=1).abs()
    sds = wc.std(axis=1)
    report("T4 final means ~0 (loose)", (means < 0.5).all(),
           f"max|mean|={means.max():.3f}")
    report("T4 final sds ~1 (loose)", ((sds > 0.5) & (sds < 1.5)).all(),
           f"sd range=[{sds.min():.3f},{sds.max():.3f}]")

    # Diagnostics PDF produced
    diag = tmp / "out_withincohort" / "EUR_splicing_diagnostics.pdf"
    report("T4 diagnostics PDF produced", diag.exists(), str(diag))


# ---------------------------------------------------------------------------
# T5: single-sample-cohort robustness (R)
# ---------------------------------------------------------------------------

def test_single_sample_cohort():
    print("\n=== T5: single-sample-cohort robustness ===")
    tmp = Path(tempfile.mkdtemp())

    # 2 cohorts, but cohort B has only 1 sample. No batch merging: every
    # cohort is INT'd on its own ranks regardless of size (n=1 ->
    # qnorm(1/2) = 0 for every feature).
    rng = np.random.RandomState(7)
    n_feat = 20
    bed = pd.DataFrame({
        "#chr": ["chr1"] * n_feat,
        "start": list(range(100, 100 + 10 * n_feat, 10)),
        "end": list(range(101, 101 + 10 * n_feat, 10)),
        "phenotype_id": [f"G{i}__e{i}" for i in range(n_feat)],
        "a1": rng.rand(n_feat),
        "a2": rng.rand(n_feat),
        "a3": rng.rand(n_feat),
        "b1": rng.rand(n_feat),  # single-sample cohort
    })
    bed_path = tmp / "test_single.bed"
    bed.to_csv(bed_path, sep="\t", index=False, float_format="%g")

    anc = pd.DataFrame({
        "sample_id": ["a1", "a2", "a3", "b1"],
        "assigned_ancestry": ["EUR"] * 4,
        "cohort": ["cA", "cA", "cA", "cB"],
    })
    anc_path = tmp / "ancestry.tsv"
    anc.to_csv(anc_path, sep="\t", index=False)

    out_dir = tmp / "out"
    out_dir.mkdir()

    r = run([
        str(RSCRIPT_SIF), str(R_SCRIPT),
        "--input", str(bed_path),
        "--ancestry-map", str(anc_path),
        "--ancestry", "EUR",
        "--modality", "RNA_editing",
        "--output-dir", str(out_dir),
    ])
    if r.returncode != 0:
        report("T5 Rscript runs (no crash)", False, r.stderr[-800:])
        return
    report("T5 Rscript runs (no crash with single-sample cohort)", True)

    out_bed = out_dir / "EUR_RNA_editing_int.bed"
    report("T5 output BED produced", out_bed.exists(), str(out_bed))

    if out_bed.exists():
        df = pd.read_csv(out_bed, sep="\t")
        data = df.drop(columns=BED_META)
        report("T5 no NaN in output", not data.isna().any().any())
        report("T5 single-sample cohort INTs to constant 0",
               (data["b1"].abs() < 1e-12).all(),
               f"b1 range=[{data['b1'].min()},{data['b1'].max()}]")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("test_normalize_modalities.py — synthetic test suite")
    print("=" * 60)

    if not POOLER.exists():
        print(f"ERROR: pooler not found at {POOLER}")
        sys.exit(1)
    if not R_SCRIPT.exists():
        print(f"ERROR: R script not found at {R_SCRIPT}")
        sys.exit(1)

    test_detection_filter_and_int()
    test_pooling_intersection()
    test_prepooled_passthrough()
    test_withincohort_int_endtoend()
    test_single_sample_cohort()

    print(f"\n{'=' * 60}")
    print(f"Results: {PASS} passed, {FAIL} failed")
    print("=" * 60)
    sys.exit(0 if FAIL == 0 else 1)


if __name__ == "__main__":
    main()
