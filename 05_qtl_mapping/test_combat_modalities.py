#!/usr/bin/env python3
"""
test_combat_modalities.py — Synthetic test suite for the cross-cohort ComBat +
pooling helpers (roadmap step 4).

Tests:
  T1: logit transform correctness (R script, via subprocess)
  T2: pooling phenotype_id intersection (Python pooler)
  T3: namespaced pass-through for pre-pooled modalities (Python pooler)
  T4: ComBat end-to-end on synthetic batched data (R script)
  T5: single-batch fallback (R script)

Run:
  python3 test_combat_modalities.py
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
R_SCRIPT = HERE / "combat_normalize_modalities.R"

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
# T1: logit transform correctness (R)
# ---------------------------------------------------------------------------

def test_logit_transform():
    print("\n=== T1: logit transform correctness ===")
    tmp = Path(tempfile.mkdtemp())

    # Build a tiny proportion BED: 3 features x 4 samples, values in [0,1]
    # including exact 0 and 1 to test clipping.
    bed = pd.DataFrame({
        "#chr": ["chr1"] * 3,
        "start": [100, 200, 300],
        "end": [101, 201, 301],
        "phenotype_id": ["G1__e1", "G1__e2", "G2__e3"],
        "s1": [0.0, 0.5, 1.0],
        "s2": [0.1, 0.9, 0.3],
        "s3": [0.2, 0.0, 0.8],
        "s4": [1.0, 0.4, 0.6],
    })
    bed_path = tmp / "test_logit.bed"
    bed.to_csv(bed_path, sep="\t", index=False, float_format="%g")

    # Ancestry map: all 4 samples in one ancestry, 2 cohorts
    anc = pd.DataFrame({
        "sample_id": ["s1", "s2", "s3", "s4"],
        "assigned_ancestry": ["EUR"] * 4,
        "cohort": ["cA", "cA", "cB", "cB"],
    })
    anc_path = tmp / "ancestry.tsv"
    anc.to_csv(anc_path, sep="\t", index=False)

    out_dir = tmp / "out"
    out_dir.mkdir()

    r = run([
        "Rscript", str(R_SCRIPT),
        "--input", str(bed_path),
        "--ancestry-map", str(anc_path),
        "--ancestry", "EUR",
        "--modality", "alt_TSS",  # logit modality
        "--output-dir", str(out_dir),
    ])
    if r.returncode != 0:
        report("T1 Rscript runs", False, r.stderr[-500:])
        return

    out_bed = out_dir / "EUR_alt_TSS_combat_int.bed"
    if not out_bed.exists():
        report("T1 output BED exists", False, str(out_bed))
        return
    report("T1 Rscript runs + output exists", True)

    # Check: no Inf / NaN in output (logit clipping worked)
    df = pd.read_csv(out_bed, sep="\t")
    data = df.drop(columns=BED_META)
    has_inf = np.isinf(data.values).any()
    has_nan = data.isna().any().any()
    report("T1 no Inf in output (clipping works)", not has_inf,
           f"has_inf={has_inf}")
    report("T1 no NaN in output", not has_nan, f"has_nan={has_nan}")

    # Check: INT output is approximately standard normal per feature
    # (mean ~0, sd ~1). With only 4 samples this is loose, so check |mean|<1
    # and 0.3<sd<2 as a sanity bound.
    means = data.mean(axis=1).abs()
    sds = data.std(axis=1)
    mean_ok = (means < 1.0).all()
    sd_ok = ((sds > 0.3) & (sds < 2.0)).all()
    report("T1 INT means ~0", mean_ok, f"max|mean|={means.max():.3f}")
    report("T1 INT sds ~1", sd_ok, f"sd range=[{sds.min():.3f},{sds.max():.3f}]")


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
# T4: ComBat end-to-end on synthetic batched data (R)
# ---------------------------------------------------------------------------

def test_combat_endtoend():
    print("\n=== T4: ComBat end-to-end (batch removal) ===")
    tmp = Path(tempfile.mkdtemp())

    # Synthetic data: 50 features x 40 samples, 2 cohorts (20 each).
    # Inject a cohort-specific mean shift so ComBat has something to remove.
    rng = np.random.RandomState(42)
    n_feat, n_samp = 50, 40
    base = rng.randn(n_feat, n_samp)
    # Cohort B (cols 20:40) gets +3 shift on all features
    base[:, 20:] += 3.0
    # Clip to [0.01, 0.99] to simulate proportions
    base = np.clip(np.abs(base) / (np.abs(base).max() + 1), 0.01, 0.99)

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
    bed_path = tmp / "test_combat.bed"
    bed.to_csv(bed_path, sep="\t", index=False, float_format="%g")

    anc = pd.DataFrame({
        "sample_id": samples,
        "assigned_ancestry": ["EUR"] * n_samp,
        "cohort": cohorts,
    })
    anc_path = tmp / "ancestry.tsv"
    anc.to_csv(anc_path, sep="\t", index=False)

    out_dir = tmp / "out"
    out_dir.mkdir()

    r = run([
        "Rscript", str(R_SCRIPT),
        "--input", str(bed_path),
        "--ancestry-map", str(anc_path),
        "--ancestry", "EUR",
        "--modality", "splicing",  # logit
        "--output-dir", str(out_dir),
    ])
    if r.returncode != 0:
        report("T4 Rscript runs", False, r.stderr[-800:])
        return
    report("T4 Rscript runs", True)

    out_bed = out_dir / "EUR_splicing_combat_int.bed"
    df = pd.read_csv(out_bed, sep="\t")
    data = df.drop(columns=BED_META)

    # Check: after ComBat, cohort separation in PCA space should be reduced.
    # We measure the centroid distance between cohorts on PC1, before (logit
    # of input proportions) and after (ComBat+INT output). This is scale-
    # invariant, unlike comparing raw vs INT-transformed means directly.
    from sklearn.decomposition import PCA
    ca, cb = samples[:20], samples[20:]

    # Before: logit-transform the input proportions (matching what ComBat sees)
    in_data = bed.set_index("phenotype_id")[samples]  # features x samples
    eps = 1e-4
    in_logit = np.log((in_data.clip(eps, 1 - eps)) / (1 - in_data.clip(eps, 1 - eps)))
    pca_before = PCA(n_components=2).fit_transform(in_logit.T)  # samples x PC
    before_centroid = abs(pca_before[:20, 0].mean() - pca_before[20:, 0].mean())

    # After: ComBat+INT output (features x samples)
    pca_after = PCA(n_components=2).fit_transform(data.T)
    after_centroid = abs(pca_after[:20, 0].mean() - pca_after[20:, 0].mean())

    report("T4 cohort separation reduced by ComBat (PC1 centroid distance)",
           after_centroid < before_centroid * 0.5,
           f"before={before_centroid:.3f}, after={after_centroid:.3f}")

    # INT output ~ standard normal (per feature = per row in BED orientation)
    means = data.mean(axis=1).abs()
    sds = data.std(axis=1)
    report("T4 INT means ~0", (means < 0.5).all(), f"max|mean|={means.max():.3f}")
    report("T4 INT sds ~1", ((sds > 0.7) & (sds < 1.3)).all(),
           f"sd range=[{sds.min():.3f},{sds.max():.3f}]")

    # Diagnostics PDF produced
    diag = out_dir / "EUR_splicing_combat_diagnostics.pdf"
    report("T4 diagnostics PDF produced", diag.exists(), str(diag))


# ---------------------------------------------------------------------------
# T5: single-batch fallback (R)
# ---------------------------------------------------------------------------

def test_single_batch_fallback():
    print("\n=== T5: single-batch fallback ===")
    tmp = Path(tempfile.mkdtemp())

    # 2 cohorts, but cohort B has only 1 sample -> merged into 'other'
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
        "Rscript", str(R_SCRIPT),
        "--input", str(bed_path),
        "--ancestry-map", str(anc_path),
        "--ancestry", "EUR",
        "--modality", "RNA_editing",  # logit
        "--output-dir", str(out_dir),
    ])
    if r.returncode != 0:
        report("T5 Rscript runs (no crash)", False, r.stderr[-800:])
        return
    report("T5 Rscript runs (no crash with single-sample batch)", True)

    out_bed = out_dir / "EUR_RNA_editing_combat_int.bed"
    report("T5 output BED produced", out_bed.exists(), str(out_bed))

    if out_bed.exists():
        df = pd.read_csv(out_bed, sep="\t")
        data = df.drop(columns=BED_META)
        report("T5 no NaN in output", not data.isna().any().any())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("test_combat_modalities.py — synthetic test suite")
    print("=" * 60)

    if not POOLER.exists():
        print(f"ERROR: pooler not found at {POOLER}")
        sys.exit(1)
    if not R_SCRIPT.exists():
        print(f"ERROR: R script not found at {R_SCRIPT}")
        sys.exit(1)

    test_logit_transform()
    test_pooling_intersection()
    test_prepooled_passthrough()
    test_combat_endtoend()
    test_single_batch_fallback()

    print(f"\n{'=' * 60}")
    print(f"Results: {PASS} passed, {FAIL} failed")
    print("=" * 60)
    sys.exit(0 if FAIL == 0 else 1)


if __name__ == "__main__":
    main()
