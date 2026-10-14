"""Nearest-TSS annotations preserve the original scalar selection rules."""

import numpy as np
import polars as pl

from epykit import annotate as annotation
from epykit.annotate import annotate_features


def scalar_lookup(chroms, mids, index):
    genes = np.full(len(mids), '', object)
    distances = np.full(len(mids), np.iinfo(np.int32).min, np.int64)
    for i, chrom in enumerate(chroms):
        if chrom not in index:
            continue
        pos, names, strands = index[chrom]
        center = int(mids[i])
        j = int(np.searchsorted(pos, center, side='left'))
        candidates = ([j-1] if j else []) + ([j] if j < len(pos) else [])
        if not candidates:
            continue
        best = min(candidates, key=lambda k: abs(int(pos[k])-center))
        genes[i] = str(names[best])
        distances[i] = (-1 if strands[best] == '-' else 1)*(center-int(pos[best]))
    return genes, distances


def test_vector_lookup_preserves_scalar_ties_strands_and_missing_chromosomes():
    index = {
        'chr1': (np.array([10, 10, 20]), np.array(['a', 'b', 'c']), np.array(['+', '-', '-'])),
        'chr2': (np.array([40]), np.array(['d']), np.array(['+'])),
        'empty': (np.array([], int), np.array([], str), np.array([], str)),
    }
    rng = np.random.default_rng(5511)
    mids = np.r_[10, 15, 10.5, -10, 2**33, rng.uniform(-100, 100, 1000)]
    chroms = np.r_[['chr1']*5, rng.choice(['chr1', 'chr2', 'empty', 'absent'], 1000)]
    expected = scalar_lookup(chroms, mids, index)
    got = annotation._nearest_tss_arrays(chroms, mids, index)
    for actual, want in zip(got, expected, strict=True):
        np.testing.assert_array_equal(actual, want)
    assert got[0][:3].tolist() == ['a', 'b', 'a']
    assert got[1][:3].tolist() == [0, -5, 0]
    empty = annotation._nearest_tss_arrays([], np.array([]), index)
    assert empty[0].shape == empty[1].shape == (0,)


def test_full_annotation_retains_site_order_dtype_and_sentinels(tmp_path, monkeypatch):
    gtf = tmp_path/'tiny.gtf'
    gtf.write_text('chr1\ttest\tgene\t101\t200\t.\t+\t.\tgene_id "g1"; gene_name "Gene1"; gene_type "protein_coding";\n'
                   'chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; gene_name "Gene1"; gene_type "protein_coding";\n')
    sites = pl.DataFrame({'chrom': ['chr1', 'missing', 'chr1'], 'pos': [110, 10, 90]})
    result = annotate_features(sites, str(gtf), multi_annotation=True)
    assert result['chrom'].to_list() == sites['chrom'].to_list()
    assert result['pos'].to_list() == sites['pos'].to_list()
    assert result['nearest_tss_gene'].to_list() == ['Gene1', '', 'Gene1']
    assert result['nearest_tss_distance'].dtype == pl.Int32
    assert result['nearest_tss_distance'].to_list() == [10, np.iinfo(np.int32).min, -10]

    def unexpected(*args):
        raise AssertionError('single annotation must not calculate nearest TSS')

    monkeypatch.setattr(annotation, '_nearest_tss_arrays', unexpected)
    single = annotate_features(sites, str(gtf), multi_annotation=False)
    assert 'nearest_tss_gene' not in single.columns
