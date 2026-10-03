#!/usr/bin/env python3
"""Tests for 23_prepare_intersection.py replicate handling.

Covers the keep-all back-mapping fix and the replicate-averaging guards:
  - backmap_rnaseq_ids keeps ALL runs whose array_id survives (the old
    dict inversion silently kept one hash-order-dependent survivor)
  - average_duplicate_columns averages technical-replicate sample columns
  - average_duplicate_rows averages replicate rows (numeric mean, other first)
  - rename_bed_columns averages duplicate columns after rnaseq_id->array_id
    rename (the guard the fix is designed to actually reach)

Runs standalone (python3 test_prepare_intersection.py) or under pytest.
"""
import importlib.util
import os
import sys
import tempfile

import numpy as np
import pandas as pd

MODULE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '23_prepare_intersection.py')
spec = importlib.util.spec_from_file_location('prepare_intersection', MODULE_PATH)
pi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pi)


def test_backmap_keeps_all_replicate_runs():
    rna_to_array = {'SRR1': 'ARRAY1', 'SRR2': 'ARRAY1', 'SRR3': 'ARRAY2'}
    out = pi.backmap_rnaseq_ids(rna_to_array, {'ARRAY1'})
    assert out == {'SRR1', 'SRR2'}, f"both runs must survive, got {out}"


def test_backmap_drops_non_intersection():
    rna_to_array = {'SRR1': 'ARRAY1', 'SRR2': 'ARRAY1', 'SRR3': 'ARRAY2'}
    out = pi.backmap_rnaseq_ids(rna_to_array, {'ARRAY2'})
    assert out == {'SRR3'}, f"only the intersecting run kept, got {out}"


def test_backmap_old_inversion_would_lose_data():
    # Document the regression: dict inversion keeps exactly one survivor.
    rna_to_array = {'SRR1': 'ARRAY1', 'SRR2': 'ARRAY1'}
    inverted = {v: k for k, v in rna_to_array.items()}
    assert len(inverted) == 1
    assert len(pi.backmap_rnaseq_ids(rna_to_array, {'ARRAY1'})) == 2


def test_average_duplicate_columns():
    df = pd.DataFrame({'A': [1.0, 3.0], 'B': [10.0, 20.0],
                       'A2': [5.0, 7.0]})
    df.columns = ['A', 'B', 'A']  # duplicate sample column
    out, dup = pi.average_duplicate_columns(df)
    assert dup == ['A']
    assert list(out.columns) == ['B', 'A'] or set(out.columns) == {'A', 'B'}
    assert np.isclose(out['A'].tolist(), [3.0, 5.0]).all()  # mean of 1,5 and 3,7
    assert np.isclose(out['B'].tolist(), [10.0, 20.0]).all()


def test_average_duplicate_columns_noop():
    df = pd.DataFrame({'A': [1.0], 'B': [2.0]})
    out, dup = pi.average_duplicate_columns(df)
    assert dup == [] and out.equals(df)


def test_average_duplicate_rows():
    df = pd.DataFrame({
        'sample_id': ['A', 'A', 'B'],
        'ct_Troph': [0.2, 0.4, 0.5],
        'cohort': ['NIGMS', 'NIGMS', 'GUSTO'],
    })
    out, n = pi.average_duplicate_rows(df, 'sample_id')
    assert n == 1
    assert list(out.columns) == ['sample_id', 'ct_Troph', 'cohort']
    row_a = out[out['sample_id'] == 'A'].iloc[0]
    assert np.isclose(row_a['ct_Troph'], 0.3)
    assert row_a['cohort'] == 'NIGMS'  # non-numeric takes first
    assert out['sample_id'].is_unique


def test_rename_bed_columns_averages_replicates():
    with tempfile.TemporaryDirectory() as td:
        bed = os.path.join(td, 'in.bed')
        out = os.path.join(td, 'out.bed')
        with open(bed, 'w') as f:
            f.write('#chr\tstart\tend\tphenotype_id\tSRR1\tSRR2\tSRR3\n')
            f.write('chr1\t100\t200\tgeneA\t1.0\t3.0\t10.0\n')
            f.write('chr1\t300\t400\tgeneB\t2.0\t4.0\t20.0\n')
        # SRR1 + SRR2 are technical replicates of ARRAY1; SRR3 -> ARRAY2
        renamed, dropped = pi.rename_bed_columns(
            bed, out, {'SRR1': 'ARRAY1', 'SRR2': 'ARRAY1', 'SRR3': 'ARRAY2'})
        assert renamed == 3 and dropped == 0
        with open(out) as f:
            header = f.readline().strip().split('\t')
        samples = [c for c in header
                   if c not in ('#chr', 'start', 'end', 'phenotype_id')]
        assert samples == ['ARRAY1', 'ARRAY2'], f"got {samples}"
        df = pd.read_csv(out, sep='\t')
        assert np.isclose(df['ARRAY1'].tolist(), [2.0, 3.0]).all()
        assert np.isclose(df['ARRAY2'].tolist(), [10.0, 20.0]).all()


def test_rename_bed_columns_drops_unmapped():
    with tempfile.TemporaryDirectory() as td:
        bed = os.path.join(td, 'in.bed')
        out = os.path.join(td, 'out.bed')
        with open(bed, 'w') as f:
            f.write('#chr\tstart\tend\tphenotype_id\tSRR1\tSRR9\n')
            f.write('chr1\t100\t200\tgeneA\t1.0\t9.9\n')
        renamed, dropped = pi.rename_bed_columns(bed, out, {'SRR1': 'ARRAY1'})
        assert renamed == 1 and dropped == 1
        df = pd.read_csv(out, sep='\t')
        assert list(df.columns) == ['#chr', 'start', 'end', 'phenotype_id',
                                    'ARRAY1']


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items())
           if k.startswith('test_') and callable(v)]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"  PASS {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL {fn.__name__}: {e}")
    print(f"\n{len(fns) - failed} passed, {failed} failed")
    sys.exit(1 if failed else 0)
