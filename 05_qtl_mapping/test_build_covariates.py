#!/usr/bin/env python3
"""Tests for 25_build_covariates.py (prune-then-optimize schema).

Covers the round-4 behavior changes:
  - --exclude-covariates applied BEFORE correlation pruning (ct_Maternal
    cannot drag a correlated cell type out with it)
  - no covariate cap by default; explicit --max-covariates still caps
  - fail-fast on an unknown exclusion name
  - correlated HCPs (>0.9 vs a fixed covariate) are dropped, fixed kept
  - --hcp-k 0 still works with exclusion

Runs the script as a subprocess (main() is argparse-driven and sys.exits
on error) against synthetic per-ancestry fixtures.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      '25_build_covariates.py')

PASS, FAIL = 0, 0

def check(name, cond, detail=''):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  PASS: {name}")
    else:
        FAIL += 1
        print(f"  FAIL: {name} {detail}")


# ---------- fixtures ----------
N = 24
SAMPLES = [f'S{i:02d}' for i in range(N)]
rng = np.random.default_rng(42)

# Independent random values keep spurious |r| > 0.9 off the table
# (trig-harmonic fixtures alias at small N — do not use them)
pc_data = {f'PC{i+1}': rng.standard_normal(N) for i in range(5)}
hcp_data = {f'HCP_{i+1}': rng.standard_normal(N) for i in range(15)}

CELL_TYPES = ['Syncytiotrophoblast', 'Cytotrophoblast', 'EVT', 'Endothelial',
              'Fibroblast_Stromal', 'Hofbauer', 'Fetal_Immune', 'Maternal']


def write_fixtures(root, maternal_corr_with=None, hcp3_eq_pc1=False):
    """Write a minimal EAS fixture set into {root}/qtl and {root}/pcs.
    maternal_corr_with: if set, Maternal = 0.5 * <that column> + tiny noise.
    hcp3_eq_pc1: if True, HCP_3 = 2 * PC1 + tiny noise (|r| ~ 1).
    Returns (qtl_dir, pc_dir)."""
    qtl_dir = os.path.join(root, 'qtl')
    pc_dir = os.path.join(root, 'pcs')
    os.makedirs(qtl_dir)
    os.makedirs(pc_dir)

    # HCP factors: rows HCP_1..15, cols samples
    hcp = pd.DataFrame(hcp_data, index=SAMPLES).T
    if hcp3_eq_pc1:
        hcp.loc['HCP_3'] = 2 * pc_data['PC1'] + 1e-4 * rng.standard_normal(N)
    hcp.to_csv(os.path.join(qtl_dir, 'EAS_hcp_factors_harmonized.tsv'), sep='\t')

    # selected PCs + genotype PC table
    with open(os.path.join(qtl_dir, 'EAS_selected_pcs.txt'), 'w') as f:
        f.write('\n'.join(f'PC{i+1}' for i in range(5)) + '\n')
    pcs = pd.DataFrame({'sample_id': SAMPLES, **pc_data})
    pcs.to_csv(os.path.join(pc_dir, 'EAS_genotype_pcs.tsv'), sep='\t', index=False)

    # metadata
    meta = pd.DataFrame({
        'array_id': SAMPLES,
        'sex': ['M', 'F'] * (N // 2),
        'GA': 28 + 12 * rng.random(N),
    })
    meta.to_csv(os.path.join(qtl_dir, 'EAS_metadata.tsv'), sep='\t', index=False)

    # deconvolution: 8 collapsed types, Syncytiotrophoblast dominant
    props = {}
    for ct in CELL_TYPES:
        if ct == 'Syncytiotrophoblast':
            props[ct] = 0.55 + 0.05 * rng.random(N)
        else:
            props[ct] = 0.03 + 0.02 * rng.random(N)
    if maternal_corr_with:
        props['Maternal'] = 0.5 * props[maternal_corr_with] + 1e-4 * rng.standard_normal(N)
    deconv = pd.DataFrame({'sample_id': SAMPLES, 'cohort': 'cohort1', **props})
    deconv.to_csv(os.path.join(qtl_dir, 'EAS_deconvolution_harmonized.tsv'),
                  sep='\t', index=False)
    return qtl_dir, pc_dir


def run_25(qtl_dir, pc_dir, extra_args=()):
    cmd = [sys.executable, SCRIPT, '--qtl-dir', qtl_dir, '--pcair-dir', pc_dir,
           '--ancestries', 'EAS'] + list(extra_args)
    return subprocess.run(cmd, capture_output=True, text=True)


def read_covs(qtl_dir, suffix=''):
    path = os.path.join(qtl_dir, f'EAS_covariates{suffix}.tsv')
    return pd.read_csv(path, sep='\t', index_col=0)


def read_pruning(qtl_dir, suffix=''):
    path = os.path.join(qtl_dir, f'EAS_covariate_pruning{suffix}.tsv')
    return pd.read_csv(path, sep='\t')


# ---------- test 1: exclusion before correlation pruning ----------
tmp1 = tempfile.mkdtemp()
qtl1, pc1 = write_fixtures(tmp1, maternal_corr_with='Hofbauer')

r = run_25(qtl1, pc1, ['--exclude-covariates', 'ct_Maternal'])
check('exclusion run exits 0', r.returncode == 0, r.stderr[-500:])
covs = read_covs(qtl1)
pruning = read_pruning(qtl1)
check('ct_Maternal excluded from output', 'ct_Maternal' not in covs.index)
check('correlated ct_Hofbauer survives when ct_Maternal excluded first',
      'ct_Hofbauer' in covs.index)
excl = pruning[pruning['reason'] == 'excluded']
check('exclusion recorded in pruning TSV',
      len(excl) == 1 and excl.iloc[0]['covariate'] == 'ct_Maternal')

# control: without exclusion, the >0.9 pair triggers correlation pruning
r = run_25(qtl1, pc1, ['--out-suffix', '_noexcl'])
check('no-exclusion control exits 0', r.returncode == 0, r.stderr[-500:])
pruning_c = read_pruning(qtl1, '_noexcl')
check('control: Maternal/Hofbauer pair correlation-pruned',
      ((pruning_c['reason'] == 'correlated')
       & (pruning_c['covariate'].isin(['ct_Maternal', 'ct_Hofbauer']))).any())

# ---------- test 2: no cap by default ----------
tmp2 = tempfile.mkdtemp()
qtl2, pc2 = write_fixtures(tmp2)
r = run_25(qtl2, pc2, ['--exclude-covariates', 'ct_Maternal'])
check('no-cap run exits 0', r.returncode == 0, r.stderr[-500:])
covs2 = read_covs(qtl2)
# 15 HCP + 5 PC + sex + GA + 6 ct (8 minus Syncytiotrophoblast reference
# minus ct_Maternal) = 28; nothing should be cap-dropped
check('no cap: all 28 covariates retained', covs2.shape[0] == 28,
      f'got {covs2.shape[0]}')
check('no cap: HCP_15 retained', 'HCP_15' in covs2.index)

# ---------- test 3: explicit cap still works ----------
r = run_25(qtl2, pc2, ['--exclude-covariates', 'ct_Maternal',
                       '--max-covariates', '25', '--out-suffix', '_capped'])
check('capped run exits 0', r.returncode == 0, r.stderr[-500:])
covs3 = read_covs(qtl2, '_capped')
check('explicit cap: 25 covariates', covs3.shape[0] == 25,
      f'got {covs3.shape[0]}')
check('explicit cap: highest-index HCPs dropped first',
      'HCP_15' not in covs3.index and 'HCP_11' in covs3.index)

# ---------- test 4: fail-fast on unknown exclusion name ----------
r = run_25(qtl2, pc2, ['--exclude-covariates', 'ct_Nope'])
check('unknown exclusion name fails', r.returncode != 0)
check('failure message lists the bad name', 'ct_Nope' in r.stdout + r.stderr)

# ---------- test 5: correlated HCP dropped, fixed covariate kept ----------
tmp5 = tempfile.mkdtemp()
qtl5, pc5 = write_fixtures(tmp5, hcp3_eq_pc1=True)
r = run_25(qtl5, pc5, ['--exclude-covariates', 'ct_Maternal'])
check('correlated-HCP run exits 0', r.returncode == 0, r.stderr[-500:])
covs5 = read_covs(qtl5)
pruning5 = read_pruning(qtl5)
check('HCP_3 (|r|~1 with PC1) dropped', 'HCP_3' not in covs5.index)
check('PC1 (fixed) kept', 'PC1' in covs5.index)
h3 = pruning5[pruning5['covariate'] == 'HCP_3']
check('HCP_3 drop recorded as correlated with PC1',
      len(h3) == 1 and h3.iloc[0]['reason'] == 'correlated'
      and h3.iloc[0]['correlated_with'] == 'PC1')

# ---------- test 6: --hcp-k 0 with exclusion ----------
r = run_25(qtl2, pc2, ['--exclude-covariates', 'Maternal',  # bare name form
                       '--hcp-k', '0', '--out-suffix', '_k0'])
check('k=0 run exits 0', r.returncode == 0, r.stderr[-500:])
covs6 = read_covs(qtl2, '_k0')
check('k=0: no HCP covariates',
      not any(str(i).startswith('HCP_') for i in covs6.index))
check('k=0: bare name Maternal matched ct_Maternal',
      'ct_Maternal' not in covs6.index)
check('k=0: fixed covariates present',
      {'PC1', 'sex', 'GA', 'ct_Hofbauer'}.issubset(set(covs6.index)))

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
