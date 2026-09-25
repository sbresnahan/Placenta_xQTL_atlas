#!/usr/bin/env python3
"""Smoke tests for optimize_hcp_chr1.py helper functions (module 25a).

Tests the pure-python helpers with synthetic fixtures; the full end-to-end
run requires tensorQTL + genotypes and is exercised on seadragon.
"""
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import optimize_hcp_chr1 as m

PASS, FAIL = 0, 0

def check(name, cond, detail=''):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  PASS: {name}")
    else:
        FAIL += 1
        print(f"  FAIL: {name} {detail}")

tmp = tempfile.mkdtemp()

# ---------- fixtures ----------
# metadata: 4 post-outlier samples (of 5 sequenced; S5 dropped by outlier removal)
meta = pd.DataFrame({
    'rnaseq_id': ['RNA1', 'RNA2', 'RNA3', 'RNA4'],
    'array_id':  ['ARR1', 'ARR2', 'ARR3', 'ARR4'],
})
meta_path = os.path.join(tmp, 'EAS_metadata.tsv')
meta.to_csv(meta_path, sep='\t', index=False)

# HCP file in rnaseq_id space, including the dropped sample RNA5
hcp = pd.DataFrame(
    np.random.default_rng(1).normal(size=(3, 5)),
    index=['HCP1', 'HCP2', 'HCP3'],
    columns=['RNA1', 'RNA2', 'RNA3', 'RNA4', 'RNA5'])
hcp.index.name = 'covariate'
hcp_path = os.path.join(tmp, 'EAS_hcp_factors.tsv')
hcp.to_csv(hcp_path, sep='\t')

# ---------- 1. harmonize_hcp ----------
print("\n[1] harmonize_hcp")
out = os.path.join(tmp, 'harmonized.tsv')
h = m.harmonize_hcp(hcp_path, meta, out)
check("all 4 post-outlier samples retained", h.shape[1] == 4, f"got {h.shape[1]}")
check("dropped sample RNA5 excluded", 'ARR5' not in h.columns and 'RNA5' not in h.columns)
check("columns renamed to array_id", list(h.columns) == ['ARR1', 'ARR2', 'ARR3', 'ARR4'])
check("3 factors preserved", h.shape[0] == 3)
check("values follow the rename (RNA2->ARR2)",
      np.isclose(h.loc['HCP1', 'ARR2'], hcp.loc['HCP1', 'RNA2']))
# idempotent on already-array_id columns
h2 = m.harmonize_hcp(out, meta, os.path.join(tmp, 'harmonized2.tsv'))
check("idempotent on array_id columns", h2.shape == (3, 4))

# ---------- 2. write_empty_hcp ----------
print("\n[2] write_empty_hcp (k=0)")
empty_path = os.path.join(tmp, 'empty_hcp.tsv')
m.write_empty_hcp(meta, empty_path)
with open(empty_path) as f:
    lines = f.read().strip().split('\n')
check("header-only file (1 line)", len(lines) == 1, f"got {len(lines)} lines")
hdr = lines[0].split('\t')
check("header = covariate + sorted array_ids",
      hdr == ['covariate', 'ARR1', 'ARR2', 'ARR3', 'ARR4'], str(hdr))
# confirm 25_build_covariates.py can consume it (sample set from columns)
eh = pd.read_csv(empty_path, sep='\t', index_col=0)
check("readable with 0 rows, 4 sample cols", eh.shape == (0, 4), str(eh.shape))

# ---------- 3. count_egenes ----------
print("\n[3] count_egenes")
rng = np.random.default_rng(2)
n = 1000
pq = pd.DataFrame({
    'phenotype_id': np.repeat([f'g{i}' for i in range(200)], 5),  # 5 variants/gene
    'pval_nominal': rng.uniform(0, 1, n),
    'qval': np.concatenate([rng.uniform(0, 0.04, 37),      # 37 significant tests
                            rng.uniform(0.06, 1, n - 37)]),
})
# make two extra significant tests share one gene -> unique count < hit-row count
pq.loc[500, 'qval'] = 0.01
pq.loc[501, 'qval'] = 0.02
pq.loc[501, 'phenotype_id'] = pq.loc[500, 'phenotype_id']
pq_path = os.path.join(tmp, 'cisqtl.parquet')
pq.to_parquet(pq_path, index=False)
n_eg, n_tested = m.count_egenes(pq_path, 0.05)
n_sig_rows = int((pq['qval'] <= 0.05).sum())
n_sig_genes = int(pq.loc[pq['qval'] <= 0.05, 'phenotype_id'].nunique())
check("eGene count = unique significant genes", n_eg == n_sig_genes,
      f"{n_eg} vs {n_sig_genes}")
check("unique count deflated vs hit rows", n_eg < n_sig_rows,
      f"{n_eg} vs {n_sig_rows}")
# sharing the pair saves exactly 1 gene vs not sharing
pq2 = pq.copy()
pq2.loc[501, 'phenotype_id'] = 'g_unique_501'
pq2_path = os.path.join(tmp, 'cisqtl2.parquet')
pq2.to_parquet(pq2_path, index=False)
check("shared-gene pair saves exactly 1 eGene",
      m.count_egenes(pq2_path, 0.05)[0] == n_eg + 1)
check("n_tested = total rows", n_tested == n)
check("stricter FDR reduces count", m.count_egenes(pq_path, 0.01)[0] <= n_eg)

# ---------- 4. k* selection logic (argmax, ties -> smaller k) ----------
print("\n[4] k* selection (argmax, ties -> smaller k)")
def select_k(results):
    res_df = pd.DataFrame(results).sort_values('k').reset_index(drop=True)
    best = res_df['n_egenes'].max()
    return int(res_df.loc[res_df['n_egenes'] == best, 'k'].min())

check("clear maximum wins", select_k([
    {'k': 0, 'n_egenes': 10}, {'k': 5, 'n_egenes': 42},
    {'k': 10, 'n_egenes': 30}]) == 5)
check("tie -> smaller k", select_k([
    {'k': 0, 'n_egenes': 42}, {'k': 5, 'n_egenes': 42},
    {'k': 10, 'n_egenes': 42}]) == 0)
check("tie at nonzero -> smaller k", select_k([
    {'k': 0, 'n_egenes': 10}, {'k': 15, 'n_egenes': 50},
    {'k': 25, 'n_egenes': 50}]) == 15)
check("monotone increasing -> largest k", select_k([
    {'k': 0, 'n_egenes': 1}, {'k': 5, 'n_egenes': 2},
    {'k': 30, 'n_egenes': 3}]) == 30)

# ---------- 5. plot_optimization renders ----------
print("\n[5] plot_optimization")
res_df = pd.DataFrame({'k': [0, 5, 10, 15], 'n_egenes': [10, 42, 38, 35]})
png = os.path.join(tmp, 'opt.png')
m.plot_optimization(res_df, 'EAS', 5, 0.05, png)
check("png written and non-trivial", os.path.exists(png) and os.path.getsize(png) > 5000,
      f"size={os.path.getsize(png) if os.path.exists(png) else 0}")

# ---------- 6. bgzip_tabix round-trip ----------
print("\n[6] bgzip_tabix")
bed = pd.DataFrame({
    '#chr': ['chr1', 'chr1'], 'start': [100, 200], 'end': [200, 300],
    'phenotype_id': ['g1', 'g2'], 'ARR1': [0.1, -0.2], 'ARR2': [0.3, 0.05],
})
bed_path = os.path.join(tmp, 'expr.bed')
bed.to_csv(bed_path, sep='\t', index=False)
gz = os.path.join(tmp, 'expr.bed.gz')
m.bgzip_tabix(bed_path, gz)
check("bgzip + tbi written",
      os.path.exists(gz) and os.path.exists(gz + '.tbi'))
try:
    import pysam
    with pysam.TabixFile(gz) as tbx:
        rows = list(tbx.fetch('chr1', 150, 250))   # overlaps both 100-200 and 200-300
        rows1 = list(tbx.fetch('chr1', 210, 290))  # overlaps only 200-300
    check("tabix overlap query returns 2 rows", len(rows) == 2, f"got {len(rows)}")
    check("tabix narrow query returns 1 row", len(rows1) == 1, f"got {len(rows1)}")
except ImportError:
    check("tabix query (skipped: no pysam)", True)

# ---------- 6b. bgzip_tabix sorts unsorted input (tabix requires it) ----------
print("\n[6b] bgzip_tabix unsorted input")
unsorted_bed = pd.DataFrame({
    '#chr': ['chr1', 'chr2', 'chr1', 'chr1'],
    'start': [500, 100, 100, 300], 'end': [600, 200, 200, 400],
    'phenotype_id': ['g3', 'gx', 'g1', 'g2'],
    'ARR1': [0.1, 9.9, -0.2, 0.4], 'ARR2': [0.3, 9.9, 0.05, 0.6],
})
unsorted_path = os.path.join(tmp, 'expr_unsorted.bed')
unsorted_bed.to_csv(unsorted_path, sep='\t', index=False)
gz_u = os.path.join(tmp, 'expr_unsorted.bed.gz')
m.bgzip_tabix(unsorted_path, gz_u)
try:
    import pysam
    with pysam.TabixFile(gz_u) as tbx:
        starts = [int(r.split('\t')[1]) for r in tbx.fetch('chr1', 0, 10**9)]
    check("unsorted input indexed; chr1 starts come back sorted",
          starts == [100, 300, 500], f"got {starts}")
except ImportError:
    check("unsorted-input query (skipped: no pysam)", True)
check("no leftover sorted-temp file",
      not os.path.exists(gz_u.replace('.bed.gz', '') + '.sorted.tmp.bed'))

# ---------- 6c. bgzip_tabix htslib-binary fallback (pysam blocked) ----------
print("\n[6c] bgzip_tabix binary fallback")
import shutil
if shutil.which('bgzip') and shutil.which('tabix'):
    sys.modules['pysam'] = None  # force ImportError -> binary fallback
    try:
        gz_b = os.path.join(tmp, 'expr_bin.bed.gz')
        log_b = os.path.join(tmp, 'bgzip_tabix.log')
        m.bgzip_tabix(unsorted_path, gz_b, log_b)
        check("binary fallback wrote .bed.gz + .tbi",
              os.path.exists(gz_b) and os.path.exists(gz_b + '.tbi'))
        check("binary fallback accepted log_path", os.path.exists(log_b))
    finally:
        del sys.modules['pysam']
else:
    check("binary fallback (skipped: no bgzip/tabix on PATH)", True)

print(f"\n{'='*50}\nTOTAL: {PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
