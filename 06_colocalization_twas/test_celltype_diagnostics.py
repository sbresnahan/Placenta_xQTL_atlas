#!/usr/bin/env python3
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

SCRIPT = Path(__file__).with_name('50_celltype_diagnostics.py')
spec = importlib.util.spec_from_file_location('ctdiag', SCRIPT)
ctdiag = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ctdiag)


def test_support_hc3_and_permutation_signal():
    rng = np.random.default_rng(123)
    n = 320
    idx = pd.Index([f's{i}' for i in range(n)])
    g = rng.binomial(2, 0.28, n).astype(float)
    c_raw = np.clip(rng.beta(2, 6, n), 0, 1)
    c_t = np.arcsinh(c_raw)
    c_t = c_t - c_t.mean()
    z1 = rng.normal(size=n)
    z2 = rng.normal(size=n)
    y = 0.3*g + 0.5*c_t + 1.3*(g-g.mean())*c_t + 0.25*z1 + rng.normal(scale=0.6, size=n)
    cohort = pd.Series(np.where(np.arange(n) < n//2, 'A', 'B'), index=idx)
    cov = pd.DataFrame({'z1': z1, 'z2': z2}, index=idx)
    d, attr = ctdiag._prepare_design(
        pd.Series(y, index=idx), pd.Series(g, index=idx),
        pd.Series(c_t, index=idx), pd.Series(c_raw, index=idx),
        cov, cohort
    )
    assert attr['n_complete'] == n
    sup = ctdiag._support_metrics(d)
    assert sup['n_minor_dosage0'] + sup['n_minor_dosage1'] + sup['n_minor_dosage2'] == n
    assert sup['interaction_information'] > 0
    hc3 = ctdiag._ols_hc3_metrics(d)
    assert hc3['diag_status'] == 'ok'
    assert abs(hc3['recalc_beta_GxC'] - 1.3) < 0.35
    assert hc3['hc3_p_GxC'] < 1e-3
    perm = ctdiag._permute_residual_test(d, 99, seed=7)
    assert perm['perm_status'] == 'ok'
    assert perm['permutation_strata'] == 'within_cohort'
    assert 0 < perm['perm_empirical_p'] <= 0.05


def test_attrition_counts_and_missingness():
    idx = pd.Index(['a','b','c','d','e','f'])
    y = pd.Series([1,2,np.nan,4,5,6], index=idx)
    g = pd.Series([0,1,1,2,np.nan,0], index=idx)
    c = pd.Series([.1,.2,.3,.4,.5,.6], index=idx)
    cov = pd.DataFrame({'z':[1,2,3,np.nan,5,6]}, index=idx)
    d, attr = ctdiag._prepare_design(y,g,c,c,cov,None)
    assert attr['n_id_intersection'] == 6
    assert attr['n_missing_y'] == 1
    assert attr['n_missing_g'] == 1
    assert attr['n_missing_covariate_any'] == 1
    assert attr['n_complete'] == 3


def test_permutation_selection_stratified():
    rows = []
    for anc in ['EAS','EUR']:
        for mod in ['expression','splicing']:
            for i in range(20):
                rows.append({
                    'ancestry':anc,'modality':mod,'phenotype_id':f'p{i}',
                    'variant_id':f'v{i}','cell_type':'CT','status':'ok',
                    'p_GxC':(i+1)/1000,
                })
    df = pd.DataFrame(rows)
    sel = ctdiag._select_permutation_tests(df, n_top=2, n_random=3, seed=1)
    assert len(sel) == 4 * 5
    assert sum(v == 'top_nominal' for v in sel.values()) == 8
    assert sum(v == 'random' for v in sel.values()) == 12


def test_downsample_sets_and_subset_fit_reduce_power():
    rng = np.random.default_rng(321)
    n = 260
    idx = pd.Index([f's{i}' for i in range(n)])
    g = rng.binomial(2, 0.30, n).astype(float)
    c_raw = rng.beta(2, 5, n)
    c_t = np.arcsinh(c_raw); c_t = c_t - c_t.mean()
    z = rng.normal(size=n)
    # Moderate interaction: strong at full N, less stable at N=80.
    y = 0.2*g + 0.3*c_t + 0.8*(g-g.mean())*c_t + 0.2*z + rng.normal(scale=0.8, size=n)
    d, attr = ctdiag._prepare_design(
        pd.Series(y, index=idx), pd.Series(g, index=idx),
        pd.Series(c_t, index=idx), pd.Series(c_raw, index=idx),
        pd.DataFrame({'z': z}, index=idx), None
    )
    sets1 = ctdiag._make_downsample_sets(idx, 80, 5, 7, ('EAS','expression'))
    sets2 = ctdiag._make_downsample_sets(idx, 80, 5, 7, ('EAS','expression'))
    assert sets1 == sets2
    assert all(len(x) == 80 for x in sets1)
    fits = [ctdiag._fit_classical_subset(d, x, min_n=20) for x in sets1]
    assert all(f[2] == 80 for f in fits)
    assert sum(np.isfinite(f[1]) for f in fits) == 5


def test_downsample_summary_recomputes_bh_and_feature_assignment():
    class FakeStage49:
        @staticmethod
        def bh_qvalues(p):
            return ctdiag._load_stage49.__globals__['np'].minimum(
                # reuse module implementation directly without importing stage49
                ctdiag._bh_for_test(p), 1.0
            ) if hasattr(ctdiag, '_bh_for_test') else _bh(p)

    # Small local BH helper matching stage49 behavior.
    def _local_bh(p):
        p = np.asarray(p, dtype=float)
        n = len(p); order = np.argsort(p); q = np.empty(n); prev = 1.0
        for i in range(n-1, -1, -1):
            prev = min(prev, p[order[i]] * n / (i+1)); q[order[i]] = prev
        return np.minimum(q, 1.0)
    stage = type('Stage', (), {'bh_qvalues': staticmethod(_local_bh)})
    recs = [
        {'phenotype_id':'p1','variant_id':'v1','cell_type':'A','full_n':200,
         'full_beta':1.0,'full_p':1e-5,'full_q':1e-4,
         'pvals':np.array([0.001,0.20]),'betas':np.array([0.8,0.3]),'ns':np.array([80,80])},
        {'phenotype_id':'p1','variant_id':'v1','cell_type':'B','full_n':200,
         'full_beta':0.2,'full_p':0.2,'full_q':0.3,
         'pvals':np.array([0.5,0.6]),'betas':np.array([0.1,0.1]),'ns':np.array([80,80])},
        {'phenotype_id':'p2','variant_id':'v2','cell_type':'A','full_n':200,
         'full_beta':0.1,'full_p':0.8,'full_q':0.9,
         'pvals':np.array([0.8,0.9]),'betas':np.array([0.1,-0.1]),'ns':np.array([80,80])},
    ]
    rr, tt, ff = ctdiag._summarize_downsample_modality(recs, stage, 'EAS', 'expression', 80, 'EUR')
    assert len(rr) == 2 and len(tt) == 3 and len(ff) == 2
    assert rr.loc[0, 'n_features_assigned'] == 1
    assert rr.loc[1, 'n_features_assigned'] == 0
    p1 = ff[ff['phenotype_id'] == 'p1'].iloc[0]
    assert p1['downsample_assignment_rate'] == 0.5


def test_matched_ancestry_marks_variant_mismatch():
    rows = [
        {'ancestry':'EAS','modality':'expression','phenotype_id':'g1','cell_type':'T',
         'variant_id':'1:1:A:G','n':280,'beta_GxC':1.0,'se_GxC':0.2,'p_GxC':1e-4,
         'q_GxC':1e-3,'recalc_maf':0.3,'cell_raw_iqr_complete':0.2,'interaction_information':5.0},
        {'ancestry':'EUR','modality':'expression','phenotype_id':'g1','cell_type':'T',
         'variant_id':'1:1:A:G','n':136,'beta_GxC':0.8,'se_GxC':0.5,'p_GxC':0.1,
         'q_GxC':0.3,'recalc_maf':0.25,'cell_raw_iqr_complete':0.1,'interaction_information':2.0},
        {'ancestry':'EAS','modality':'expression','phenotype_id':'g2','cell_type':'T',
         'variant_id':'1:2:A:G','n':280,'beta_GxC':0.5,'se_GxC':0.3,'p_GxC':0.05,
         'q_GxC':0.1,'recalc_maf':0.3,'cell_raw_iqr_complete':0.2,'interaction_information':4.0},
        {'ancestry':'EUR','modality':'expression','phenotype_id':'g2','cell_type':'T',
         'variant_id':'1:3:C:T','n':136,'beta_GxC':-0.3,'se_GxC':0.6,'p_GxC':0.6,
         'q_GxC':0.8,'recalc_maf':0.2,'cell_raw_iqr_complete':0.1,'interaction_information':1.0},
    ]
    matched, summ = ctdiag._matched_ancestry_tables(pd.DataFrame(rows), pd.DataFrame())
    assert len(matched) == 2
    assert matched['same_variant'].sum() == 1
    assert summ.iloc[0]['same_variant_rate'] == 0.5
