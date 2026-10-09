"""Fixture tests for the module-07 GxE scanner (gxe_core.py + 52_gxe_scan.py).

Layers:
  1. test_nominal_matches_tensorqtl — gxe_core's window fit reproduces
     tensorqtl.core.calculate_interaction_nominal (b, b_se, tstat) on
     identical inputs.
  2. test_residualize_permute_commutes — the linear-algebra identity the
     fixed-Xinv permutation scheme relies on.
  3. test_interaction_recovery_and_null — planted interaction recovered;
     main-effect-only phenotype stays null.
  4. test_permutation_scan_properties — adaptive schedule: strong signal
     burns all 10k permutations and hits the pval_perm floor; null stops
     early; closed-form rss matches the explicit residual rss.
  5. test_beta_approx_wrapper — finite on sane input, NaN on degenerate.
  6. test_ivw_meta — hand-computed IVW values + edge cases.
  7. test_tnt_table — exhaustive 3x3 mother/child genotype table.
  8. test_scanner_cli_end_to_end — full 52_gxe_scan.py cis-perm run on tiny
     two-ancestry pgens (exercises pgen reading, ID intersection, ancestry
     stacking, exposure/covariate alignment).

Run:  cd 07_gxe_mapping && pytest test_gxe_scan.py -v
"""

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "06_colocalization_twas"))

import gxe_core

_spec = importlib.util.spec_from_file_location("gxe_scan_52", HERE / "52_gxe_scan.py")
scanner = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(scanner)

from test_susie_coloc import (simulate_genotypes, write_vcf, make_pgen,
                              VAR_IDS, POSITIONS, CHROM)

CPU = torch.device("cpu")
RNG = np.random.default_rng(7)
N_SAMPLES = 200


@pytest.fixture(scope="module")
def inputs():
    """Genotypes (samples x variants), exposure, covariates, phenotypes."""
    G = simulate_genotypes(N_SAMPLES, seed=3)          # [ns x nv]
    e = RNG.normal(size=N_SAMPLES)                     # continuous exposure
    cov = RNG.normal(size=(N_SAMPLES, 3))              # 3 covariates
    return G, e, cov


def tq_reference(G, e, cov, P):
    """tensorQTL's own interaction scan (the reference implementation)."""
    from tensorqtl.core import Residualizer, calculate_interaction_nominal
    res = Residualizer(torch.tensor(cov, dtype=torch.float32))
    genotypes_t = torch.tensor(G.T, dtype=torch.float32)   # [variants x samples]
    phenotypes_t = torch.tensor(P, dtype=torch.float32)    # [phenotypes x samples]
    interaction_t = torch.tensor(e.reshape(-1, 1), dtype=torch.float32)
    tstat_t, b_t, b_se_t, af_t, ma_s, ma_c = calculate_interaction_nominal(
        genotypes_t, phenotypes_t, interaction_t, residualizer=res)
    return tstat_t.numpy(), b_t.numpy(), b_se_t.numpy(), res.dof


def core_scan(G, e, cov, P):
    """gxe_core equivalent over ALL variants (one big window)."""
    Q_t, dof_base = gxe_core.make_residualizer(cov, device=CPU)
    dof = dof_base - 2
    genotypes_t = torch.tensor(G.T, dtype=torch.float32)
    e_t = torch.tensor(e, dtype=torch.float32)
    g0, e0, ge0 = gxe_core.prepare_interaction_terms(genotypes_t, e_t, Q_t)
    Xt, Xinv, ok = gxe_core.build_window_design(g0, e0, ge0)
    bs, bses, ts = [], [], []
    for j in range(P.shape[0]):
        p0 = gxe_core.prepare_phenotype(
            torch.tensor(P[j], dtype=torch.float32), Q_t)
        b, b_se, t = gxe_core.window_fit(Xt, Xinv, p0, dof)
        bs.append(b.numpy())
        bses.append(b_se.numpy())
        ts.append(t.numpy())
    # -> [np x nv x 3]
    return np.stack(bs, 0), np.stack(bses, 0), np.stack(ts, 0), dof, Xt, Xinv



def test_interaction_maf_filter_matches_tensorqtl(inputs):
    """Our continuous-exposure filter must match tensorQTL exactly."""
    from tensorqtl.core import filter_maf_interaction

    G, e, _ = inputs
    Gv = G.T.astype(np.float32)
    threshold = 0.05
    ours = scanner.interaction_maf_mask(Gv, e, threshold)

    order = np.argsort(e, kind="mergesort")
    upper = np.ones(len(e), dtype=bool)
    upper[order[:len(e) // 2]] = False
    _, mask_t = filter_maf_interaction(
        torch.tensor(Gv, dtype=torch.float32),
        interaction_mask_t=torch.tensor(upper, dtype=torch.bool),
        maf_threshold_interaction=threshold)
    expected = mask_t.cpu().numpy().astype(bool)
    assert np.array_equal(ours, expected)


def test_interaction_maf_filter_rejects_exposure_imbalance():
    """A variant common overall but confined to one exposure half is removed."""
    e = np.arange(20, dtype=float)
    G = np.zeros((2, 20), dtype=np.float32)
    G[0, :10] = 1.0            # MAF 0.25 overall; absent in upper half
    G[1, [0, 10]] = 1.0       # MAF 0.05 in both halves
    mask = scanner.interaction_maf_mask(G, e, 0.05)
    assert mask.tolist() == [False, True]

def test_nominal_matches_tensorqtl(inputs):
    G, e, cov = inputs
    P = RNG.normal(size=(4, N_SAMPLES)).astype(np.float64)
    # add signal: phenotype 0 has main + interaction effects at variant 50
    P[0] += 0.6 * G[:, 50] + 0.8 * G[:, 50] * e
    P[1] += 0.6 * G[:, 100]

    t_ref, b_ref, bse_ref, dof_ref = tq_reference(G, e, cov, P)
    b, b_se, t, dof, _, _ = core_scan(G, e, cov, P)
    assert dof == dof_ref - 2  # tensorQTL: residualizer.dof - 2*ni

    # tensorQTL returns [nv x 3 x np]; core returns [np x nv x 3]
    b_ref = np.transpose(b_ref, (2, 0, 1))
    bse_ref = np.transpose(bse_ref, (2, 0, 1))
    t_ref = np.transpose(t_ref, (2, 0, 1))

    for name, mine, ref in [("b", b, b_ref), ("b_se", b_se, bse_ref),
                            ("tstat", t, t_ref)]:
        finite = np.isfinite(mine) & np.isfinite(ref)
        assert finite.all()
        maxrel = np.max(np.abs(mine - ref) / (np.abs(ref) + 1e-6))
        assert maxrel < 5e-3, f"{name}: max relative diff {maxrel:.2e}"


def test_permutation_matches_brute_force(inputs):
    """The fixed-Xinv fast path must equal a full refit of each permuted,
    re-residualized phenotype (the tensorQTL map_cis permutation scheme)."""
    G, e, cov = inputs
    y = 0.5 * G[:, 50] + RNG.normal(size=N_SAMPLES)
    Q_t, dof_base = gxe_core.make_residualizer(cov, device=CPU)
    dof = dof_base - 2
    genotypes_t = torch.tensor(G.T, dtype=torch.float32)
    e_t = torch.tensor(e, dtype=torch.float32)
    g0, e0, ge0 = gxe_core.prepare_interaction_terms(genotypes_t, e_t, Q_t)
    Xt, Xinv, ok = gxe_core.build_window_design(g0, e0, ge0)
    p_raw = torch.tensor(y, dtype=torch.float32)
    p0 = gxe_core.prepare_phenotype(p_raw, Q_t)
    _, _, tstat = gxe_core.window_fit(Xt, Xinv, p0, dof)
    r2_max = float(np.max(gxe_core.t_to_r2(tstat[:, 2].numpy(), dof)))

    rng = np.random.default_rng(99)
    r2_perm, nperm, _ = gxe_core.adaptive_permutation_scan(
        Xt, Xinv, p_raw, Q_t, dof, r2_max, rng, blocks=(25,), stop_p=2.0)
    # brute force: same permutations, full center+residualize+refit
    rng2 = np.random.default_rng(99)
    bf = []
    for _ in range(25):
        ix = torch.from_numpy(rng2.permutation(N_SAMPLES))
        p0p = gxe_core.prepare_phenotype(p_raw[ix], Q_t)
        _, _, tp = gxe_core.window_fit(Xt, Xinv, p0p, dof)
        bf.append(np.nanmax(gxe_core.t_to_r2(tp[:, 2].numpy(), dof)))
    assert np.allclose(r2_perm, bf, rtol=1e-4, atol=1e-6)


def test_interaction_recovery_and_null(inputs):
    G, e, cov = inputs
    causal = 50
    y_int = 0.5 * G[:, causal] + 0.9 * G[:, causal] * e \
        + RNG.normal(scale=1.0, size=N_SAMPLES)
    y_null = 0.8 * G[:, causal] + RNG.normal(scale=1.0, size=N_SAMPLES)
    P = np.stack([y_int, y_null])

    b, b_se, t, dof, _, _ = core_scan(G, e, cov, P)
    # interaction phenotype: top |t_gi| variant should be the causal one
    # (or a tight-LD proxy); b_gi at the causal variant recovers ~0.9
    top = int(np.argmax(np.abs(t[0, :, 2])))
    assert abs(b[0, top, 2] - 0.9) < 0.25
    assert abs(t[0, top, 2]) > 4
    # null phenotype: clearly separated from the signal phenotype
    assert np.max(np.abs(t[1, :, 2])) < 4
    assert abs(t[0, top, 2]) > np.max(np.abs(t[1, :, 2])) + 1


def test_permutation_scan_properties(inputs):
    G, e, cov = inputs
    causal = 50
    y_int = 0.5 * G[:, causal] + 1.2 * G[:, causal] * e \
        + RNG.normal(scale=1.0, size=N_SAMPLES)
    y_null = 0.8 * G[:, causal] + RNG.normal(scale=1.0, size=N_SAMPLES)

    Q_t, dof_base = gxe_core.make_residualizer(cov, device=CPU)
    dof = dof_base - 2
    genotypes_t = torch.tensor(G.T, dtype=torch.float32)
    e_t = torch.tensor(e, dtype=torch.float32)
    g0, e0, ge0 = gxe_core.prepare_interaction_terms(genotypes_t, e_t, Q_t)
    Xt, Xinv, ok = gxe_core.build_window_design(g0, e0, ge0)

    results = {}
    for name, y in [("int", y_int), ("null", y_null)]:
        p0 = gxe_core.prepare_phenotype(torch.tensor(y, dtype=torch.float32), Q_t)
        b, b_se, tstat = gxe_core.window_fit(Xt, Xinv, p0, dof)
        # closed-form rss == explicit residual rss (the permutation speedup)
        w = torch.matmul(Xt, p0)
        bb = torch.matmul(Xinv, w.unsqueeze(-1)).squeeze(-1)
        rss_closed = p0.dot(p0) - (bb * w).sum(1)
        X_t = torch.stack([g0, e0.expand(G.shape[1], -1), ge0], dim=2)
        r = torch.matmul(X_t, bb.unsqueeze(-1)).squeeze(-1) - p0
        rss_explicit = (r * r).sum(1)
        assert torch.allclose(rss_closed, rss_explicit, rtol=1e-3, atol=1e-3)

        r2 = gxe_core.t_to_r2(tstat[:, 2].numpy(), dof)
        r2_max = float(np.max(r2))
        rng = np.random.default_rng(123)
        p_raw = torch.tensor(y, dtype=torch.float32)
        r2_perm, nperm, pval_perm = gxe_core.adaptive_permutation_scan(
            Xt, Xinv, p_raw, Q_t, dof, r2_max, rng)
        results[name] = (r2_max, r2_perm, nperm, pval_perm)

    r2_max_int, r2_perm_int, nperm_int, p_int = results["int"]
    _, _, nperm_null, p_null = results["null"]
    # strong signal: no permutation exceeds -> floor p-value, full schedule
    assert p_int == pytest.approx(1.0 / (nperm_int + 1))
    assert nperm_int == 10000
    assert np.max(r2_perm_int) < r2_max_int
    # null: stops early (phat > 0.10 at the first block) or runs on with a
    # non-tiny empirical p; either way p_null >> p_int
    assert p_null > p_int
    if nperm_null == 100:
        assert p_null > 0.10


def test_beta_approx_wrapper():
    rng = np.random.default_rng(0)
    r2_perm = np.sort(rng.beta(0.5, 50, 1000))  # null-like max-r2 distribution
    out = gxe_core.beta_approx_pval(r2_perm, r2_nominal=0.9, dof=150)
    assert np.isfinite(out["pval_beta"])
    assert 0 <= out["pval_beta"] <= 1
    assert out["pval_beta"] < 0.01  # 0.9 is far in the tail
    # degenerate input (zero variance) -> NaN fallback, not an exception
    bad = gxe_core.beta_approx_pval(np.full(100, 0.01), 0.5, 150)
    assert not np.isfinite(bad["pval_beta"]) or 0 <= bad["pval_beta"] <= 1


def test_ivw_meta():
    # hand-computed: betas 0.5 (se 0.1), 0.7 (se 0.2)
    out = gxe_core.ivw_meta([0.5, 0.7], [0.1, 0.2])
    w1, w2 = 100.0, 25.0
    exp_beta = (w1 * 0.5 + w2 * 0.7) / 125.0
    assert out["beta"] == pytest.approx(exp_beta)
    assert out["se"] == pytest.approx(np.sqrt(1 / 125.0))
    assert out["k"] == 2
    assert 0 <= out["p_het"] <= 1
    # single study: no heterogeneity
    one = gxe_core.ivw_meta([0.5], [0.1])
    assert one["k"] == 1 and not np.isfinite(one["q_het"])
    # all-NaN input
    none = gxe_core.ivw_meta([np.nan, np.nan], [0.1, 0.2])
    assert none["k"] == 0 and not np.isfinite(none["beta"])


def test_tnt_table():
    # exhaustive mother x child genotype table
    gm = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2], dtype=float)
    gc = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2], dtype=float)
    T = gxe_core.transmitted_dosage(gm, gc)
    expected = [0, 0, np.nan,   # hom-ref mother: T=0; gc=2 is a Mendelian error
                0, 0.5, 1,      # het mother
                np.nan, 1, 1]   # hom-alt mother: T=1; gc=0 is an error
    assert np.allclose(T, expected, equal_nan=True)
    # missing propagates
    assert np.isnan(gxe_core.transmitted_dosage([np.nan], [1.0]))[0]


# ---------------------------------------------------------------------------
# End-to-end CLI test on tiny two-ancestry pgens
# ---------------------------------------------------------------------------

def _write_bed(pheno_df, path):
    """pheno_df: samples x phenotypes DataFrame -> bgzipped tabix BED."""
    bed = pd.DataFrame({
        "#chr": CHROM,
        "start": 100_000 + np.arange(pheno_df.shape[1]),
        "end": 100_000 + np.arange(pheno_df.shape[1]),
        "phenotype_id": pheno_df.columns,
    })
    bed = pd.concat([bed, pheno_df.T.reset_index(drop=True)], axis=1)
    bed.to_csv(path, sep="\t", index=False)
    subprocess.run(["bgzip", "-f", str(path)], check=True)
    subprocess.run(["tabix", "-p", "bed", "-f", str(path) + ".gz"], check=True)


@pytest.fixture(scope="module")
def cli_fixture(tmp_path_factory):
    """Two tiny ancestry pgens + pooled inputs; planted GxE at one variant."""
    tmp = tmp_path_factory.mktemp("gxe_cli")
    rng = np.random.default_rng(11)
    n_eas, n_eur = 60, 40
    G_eas = simulate_genotypes(n_eas, seed=21)
    G_eur = simulate_genotypes(n_eur, seed=22)
    eas_ids = [f"EAS{i}" for i in range(n_eas)]
    eur_ids = [f"EUR{i}" for i in range(n_eur)]
    for anc, G, ids in [("EAS", G_eas, eas_ids), ("EUR", G_eur, eur_ids)]:
        vcf = tmp / f"{anc}.vcf"
        write_vcf(G, VAR_IDS, ids, str(vcf))
        make_pgen(str(vcf), str(tmp / f"{anc}_qtl"))

    samples = eas_ids + eur_ids
    n = len(samples)
    e = np.concatenate([rng.normal(0, 1, n_eas), rng.normal(0.5, 1, n_eur)])
    causal = 60
    G_all = np.concatenate([G_eas, G_eur], 0)
    y1 = 0.4 * G_all[:, causal] + 1.0 * G_all[:, causal] * e + rng.normal(size=n)
    y2 = 0.6 * G_all[:, causal] + rng.normal(size=n)  # null interaction
    pheno = pd.DataFrame({"PHENO_INT": y1, "PHENO_NULL": y2},
                         index=samples)
    _write_bed(pheno, tmp / "pooled_test.bed")

    cov = pd.DataFrame(rng.normal(size=(n, 2)), index=samples,
                       columns=["PC1", "PC2"]).T
    cov.index.name = "covariate"
    cov.to_csv(tmp / "pooled_covariates_test.tsv", sep="\t")

    exp = pd.DataFrame([e], index=["GA"], columns=samples)
    exp.index.name = "exposure_id"
    exp.to_csv(tmp / "exposures.tsv", sep="\t")

    manifest = pd.DataFrame({
        "array_id": samples,
        "ancestry": ["EAS"] * n_eas + ["EUR"] * n_eur,
    })
    manifest.to_csv(tmp / "pooled_sample_manifest.tsv", sep="\t", index=False)
    return tmp, causal


def test_scanner_cli_end_to_end(cli_fixture):
    tmp, causal = cli_fixture
    out = tmp / "out.parquet"
    env = dict(os.environ)
    cmd = [
        sys.executable, str(HERE / "52_gxe_scan.py"),
        "--mode", "cis-perm",
        "--bed", str(tmp / "pooled_test.bed.gz"),
        "--pgens", f"EAS={tmp}/EAS_qtl", f"EUR={tmp}/EUR_qtl",
        "--manifest", str(tmp / "pooled_sample_manifest.tsv"),
        "--covariates", str(tmp / "pooled_covariates_test.tsv"),
        "--exposures", str(tmp / "exposures.tsv"),
        "--exposure", "GA",
        "--chrom", CHROM,
        "--maf-threshold", "0.01",
        "--maf-threshold-interaction", "0.05",
        "--perm-blocks", "50", "150",
        "--seed", "5",
        "--out", str(out),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert proc.returncode == 0, f"scanner failed:\n{proc.stderr}\n{proc.stdout}"
    df = pd.read_parquet(out)
    assert set(df.index) == {"PHENO_INT", "PHENO_NULL"}
    row = df.loc["PHENO_INT"]
    # planted interaction recovered at/near the causal variant
    assert abs(row["b_gi"] - 1.0) < 0.35
    assert row["pval_nominal"] < 1e-4
    assert row["pval_perm"] <= 0.05
    assert int(row["nperm"]) == 200  # 50 + 150 (both blocks: strong signal)
    # null phenotype: not significant
    assert df.loc["PHENO_NULL", "pval_nominal"] > 1e-3
    # top variant of the interaction phenotype is in tight LD with causal
    top_vid = row["variant_id"]
    top_pos = int(top_vid.split(":")[1])
    assert abs(top_pos - POSITIONS[causal]) <= 500
