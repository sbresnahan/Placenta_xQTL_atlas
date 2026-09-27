#!/usr/bin/env python3
"""Smoke tests for optimize_hcp_modalities.py (module 25b driver).

Tests the pure-python helpers with synthetic fixtures and runs optimize_one
end-to-end with a mocked subprocess runner (base.run); the full pipeline
requires tensorQTL + genotypes + Rhcpp and is exercised on seadragon.
"""
import argparse
import os
import re
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import optimize_hcp_modalities as m

PASS, FAIL = 0, 0

def check(name, cond, detail=''):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  PASS: {name}")
    else:
        FAIL += 1
        print(f"  FAIL: {name} {detail}")


def make_fixtures(tmp, anc='EAS', mod='splicing', n_chr1=400, n_other=100,
                  with_groups=True):
    """Synthetic qtl_dir: metadata (4 array_id samples), harmonized BED,
    optional groups file. Returns (qtl_dir, bed_path, meta_path)."""
    qtl_dir = os.path.join(tmp, f'qtl_{mod}_{n_chr1}')
    os.makedirs(qtl_dir, exist_ok=True)
    meta = pd.DataFrame({
        'rnaseq_id': ['RNA1', 'RNA2', 'RNA3', 'RNA4'],
        'array_id':  ['ARR1', 'ARR2', 'ARR3', 'ARR4'],
    })
    meta_path = os.path.join(qtl_dir, f'{anc}_metadata.tsv')
    meta.to_csv(meta_path, sep='\t', index=False)

    rng = np.random.default_rng(3)
    rows = []
    for i in range(n_chr1):
        rows.append(('chr1', 1000 + i * 10, 1001 + i * 10, f'chr1_g{i}'))
    for i in range(n_other):
        rows.append(('chr2', 1000 + i * 10, 1001 + i * 10, f'chr2_g{i}'))
    bed = pd.DataFrame(rows, columns=['#chr', 'start', 'end', 'phenotype_id'])
    for s in meta['array_id']:
        bed[s] = rng.normal(size=len(bed))
    bed_path = os.path.join(qtl_dir, f'{anc}_{mod}_harmonized.bed')
    bed.to_csv(bed_path, sep='\t', index=False)

    if with_groups:
        groups = pd.DataFrame({
            'phenotype_id': bed['phenotype_id'],
            'group_id': [p.rsplit('_', 1)[0] + '_gene' + p[-1] for p in bed['phenotype_id']],
        })
        groups.to_csv(os.path.join(qtl_dir, f'{anc}_{mod}.phenotype_groups.txt'),
                      sep='\t', index=False, header=False)
    return qtl_dir, bed_path, meta_path


def make_args(qtl_dir, tmp, **over):
    args = argparse.Namespace(
        qtl_dir=qtl_dir,
        qc_metrics=os.path.join(tmp, 'all_qc_metrics.tsv'),
        pcair_dir=tmp,
        scripts_dir=tmp,
        chr1_min_phenotypes=300,
        skip_existing=False,
        lambda1=0.5, lambda2=1.0, lambda3=1.0,
        qc_cor_threshold=0.9,
        max_hcp_phenotypes=40000,
        seed=1,
        exclude_covariates=['ct_Maternal'],
        cis_window=1000000,
        maf_threshold=0.01,
        fdr=0.05,
    )
    for k, v in over.items():
        setattr(args, k, v)
    # pooled QC metrics fixture (rnaseq_id rows)
    qc = pd.DataFrame(np.random.default_rng(4).normal(size=(4, 3)),
                      index=['RNA1', 'RNA2', 'RNA3', 'RNA4'],
                      columns=['GC_CONTENT', 'PCT_PF', 'MEDIAN_5PRIME_BIAS'])
    qc.index.name = 'sample_id'
    qc.to_csv(args.qc_metrics, sep='\t')
    return args


def mock_run_factory(egene_counts, calls):
    """Fake base.run: writes the per-k HCP file, covariates + pruning files,
    and a cisqtl parquet whose significant-gene count is egene_counts[k]."""
    def fake_run(cmd, desc, log_path=None):
        calls.append(desc)
        joined = ' '.join(cmd) if isinstance(cmd, list) else cmd
        if 'hcp_from_matrix.R' in joined:
            out = cmd[cmd.index('--output') + 1]
            k = int(cmd[cmd.index('--k') + 1])
            hdr = ['covariate', 'ARR1', 'ARR2', 'ARR3', 'ARR4']
            with open(out, 'w') as f:
                f.write('\t'.join(hdr) + '\n')
                for i in range(k):
                    f.write(f'HCP_{i+1}\t0.1\t0.2\t0.3\t0.4\n')
        elif '25_build_covariates.py' in joined:
            staging = cmd[cmd.index('--qtl-dir') + 1]
            anc = cmd[cmd.index('--ancestries') + 1]
            k = int(cmd[cmd.index('--hcp-k') + 1])
            with open(os.path.join(staging, f'{anc}_covariates_k{k}.tsv'), 'w') as f:
                f.write('covariate\tARR1\tARR2\tARR3\tARR4\n')
            pruned = pd.DataFrame({'covariate': ['ct_Maternal'] + (['HCP_9'] if k >= 10 else []),
                                   'reason': ['excluded'] + (['correlation'] if k >= 10 else [])})
            pruned.to_csv(os.path.join(staging, f'{anc}_covariate_pruning_k{k}.tsv'),
                          sep='\t', index=False)
        elif '27_run_tensorqtl.py' in joined:
            out_dir = cmd[cmd.index('--output-dir') + 1]
            anc = cmd[cmd.index('--ancestry') + 1]
            mod = cmd[cmd.index('--modality') + 1]
            cov = cmd[cmd.index('--covariates-file') + 1]
            k = int(re.search(r'_k(\d+)\.tsv$', cov).group(1))
            os.makedirs(out_dir, exist_ok=True)
            n_sig = egene_counts[k]
            df = pd.DataFrame({
                'phenotype_id': [f'g{i}' for i in range(n_sig)] + ['gN1', 'gN2'],
                'qval': [0.01] * n_sig + [0.5, 0.9],
            })
            df.to_parquet(os.path.join(out_dir, f'{anc}_{mod}_cisqtl.parquet'),
                          index=False)
        else:
            raise AssertionError(f'unexpected command: {joined}')
    return fake_run


tmp = tempfile.mkdtemp()

# ---------- 1. count_chr1_phenotypes ----------
print("\n[1] count_chr1_phenotypes")
qtl_dir, bed_path, _ = make_fixtures(tmp, n_chr1=400, n_other=100)
n_chr1, n_total = m.count_chr1_phenotypes(bed_path)
check("chr1 count", n_chr1 == 400, f"got {n_chr1}")
check("total count", n_total == 500, f"got {n_total}")

# ---------- 2. write_staging_bed scope ----------
print("\n[2] write_staging_bed")
out_bed = os.path.join(tmp, 'staging_chr1.bed')
ids = m.write_staging_bed(bed_path, out_bed, 'chr1')
check("chr1 scope keeps only chr1", len(ids) == 400 and all(i.startswith('chr1_') for i in ids))
out_bed2 = os.path.join(tmp, 'staging_gw.bed')
ids2 = m.write_staging_bed(bed_path, out_bed2, 'genome-wide')
check("genome-wide scope keeps all", len(ids2) == 500, f"got {len(ids2)}")

# ---------- 3. subset_groups_file ----------
print("\n[3] subset_groups_file")
groups_path = os.path.join(qtl_dir, 'EAS_splicing.phenotype_groups.txt')
sub_out = os.path.join(tmp, 'groups_sub.txt')
n_kept = m.subset_groups_file(groups_path, ids, sub_out)
sub = pd.read_csv(sub_out, sep='\t', header=None, names=['phenotype_id', 'group_id'])
check("groups subset to staging phenotypes", n_kept == 400 and sub.shape[0] == 400,
      f"got {n_kept}")
check("no chr2 phenotypes retained", not sub['phenotype_id'].str.startswith('chr2_').any())

# ---------- 4. optimize_one end-to-end (mocked run), chr1 scope ----------
print("\n[4] optimize_one: chr1 scope, k* tie -> smaller k, install")
work_root = os.path.join(tmp, 'work1')
os.makedirs(work_root, exist_ok=True)
args = make_args(qtl_dir, tmp)
egene_counts = {0: 10, 5: 42, 10: 42, 15: 30}
calls = []
orig_run = m.base.run
m.base.run = mock_run_factory(egene_counts, calls)
try:
    res = m.optimize_one('EAS', 'splicing', args, [0, 5, 10, 15], work_root, ['Rscript'])
finally:
    m.base.run = orig_run

k_star = int(res.loc[res['chosen'], 'k'].iloc[0])
check("k* = 5 (tie 5/10 -> smaller)", k_star == 5, f"got {k_star}")
check("mapping scope is chr1 (400 >= 300)",
      (res['mapping_scope'] == 'chr1').all())
check("k=0 wrote header-only HCP (no R call for k=0)",
      sum('HCP estimation (k=0)' in c for c in calls) == 0 and
      sum('HCP estimation' in c for c in calls) == 3,
      f"HCP calls: {[c for c in calls if 'HCP' in c]}")
k0_file = os.path.join(work_root, 'EAS', 'splicing', 'EAS_splicing_hcp_k0.tsv')
with open(k0_file) as f:
    check("k=0 HCP file is header-only", len(f.read().strip().split('\n')) == 1)
check("n_hcp_dropped counted at k>=10",
      int(res.loc[res['k'] == 10, 'n_hcp_dropped'].iloc[0]) == 1 and
      int(res.loc[res['k'] == 10, 'n_hcp_used'].iloc[0]) == 9)
installed = os.path.join(qtl_dir, 'EAS_splicing_hcp_factors_optimized.tsv')
check("optimized HCP installed", os.path.exists(installed))
k5_file = os.path.join(work_root, 'EAS', 'splicing', 'EAS_splicing_hcp_k5.tsv')
check("installed file matches k*=5 solution",
      open(installed).read() == open(k5_file).read())
check("optimal_hcp.tsv written",
      os.path.exists(os.path.join(work_root, 'EAS_splicing_optimal_hcp.tsv')))
png = os.path.join(work_root, 'EAS_splicing_optimal_hcp.png')
check("optimal_hcp.png written and non-trivial",
      os.path.exists(png) and os.path.getsize(png) > 5000)
# groups file subset to chr1 in staging
staged_groups = os.path.join(work_root, 'EAS', 'splicing', 'EAS_splicing.phenotype_groups.txt')
sg = pd.read_csv(staged_groups, sep='\t', header=None, names=['phenotype_id', 'group_id'])
check("staging groups file subset to chr1", sg.shape[0] == 400, f"got {sg.shape[0]}")

# ---------- 5. optimize_one genome-wide fallback ----------
print("\n[5] optimize_one: genome-wide fallback below chr1 minimum")
qtl_dir2, _, _ = make_fixtures(tmp, mod='RNA_editing', n_chr1=30, n_other=300)
work_root2 = os.path.join(tmp, 'work2')
os.makedirs(work_root2, exist_ok=True)
args2 = make_args(qtl_dir2, tmp)
calls2 = []
m.base.run = mock_run_factory({0: 3, 5: 8}, calls2)
try:
    res2 = m.optimize_one('EAS', 'RNA_editing', args2, [0, 5], work_root2, ['Rscript'])
finally:
    m.base.run = orig_run
check("mapping scope is genome-wide (30 < 300)",
      (res2['mapping_scope'] == 'genome-wide').all())
check("k* = 5", int(res2.loc[res2['chosen'], 'k'].iloc[0]) == 5)
check("optimized HCP installed (fallback arm)",
      os.path.exists(os.path.join(qtl_dir2, 'EAS_RNA_editing_hcp_factors_optimized.tsv')))

# ---------- 6. skip-existing resume ----------
print("\n[6] skip-existing resume")
args3 = make_args(qtl_dir2, tmp, skip_existing=True)
calls3 = []
m.base.run = mock_run_factory({0: 3, 5: 8}, calls3)
try:
    res3 = m.optimize_one('EAS', 'RNA_editing', args3, [0, 5], work_root2, ['Rscript'])
finally:
    m.base.run = orig_run
check("resume makes zero subprocess calls", len(calls3) == 0,
      f"got {len(calls3)}: {calls3}")
check("resume reproduces k*", int(res3.loc[res3['chosen'], 'k'].iloc[0]) == 5)

print(f"\n{'='*50}\nTOTAL: {PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
