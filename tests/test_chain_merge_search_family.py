"""Selection-aware correction for the default chain-merge search."""
import numpy as np
import polars as pl
from statsmodels.stats.multitest import multipletests

from epykit.dmr import apply_region_qfilter, call_dmr_chain_merge
from epykit._region_search import adjust_selected, count_intervals


def frame(chrom, positions, pvalue, effect):
    return pl.DataFrame({
        "chrom": [chrom] * len(positions), "pos": positions,
        "pvalue": [pvalue] * len(positions), "meth_diff": [effect] * len(positions),
    })


def test_background_without_selected_regions_counts_in_search_correction():
    signal = frame("chr1", np.arange(5) * 20 + 100, .02, .2)
    background = frame("chr2", np.arange(500) * 20 + 100, 1., 0.)
    kwargs = dict(min_cpgs=5, minlen_bp=50, dis_merge_bp=500)
    alone = call_dmr_chain_merge(signal, **kwargs)
    together = call_dmr_chain_merge(pl.concat([signal, background]), **kwargs)
    assert len(alone) == len(together) == 1
    assert alone["combined_pvalue"][0] == together["combined_pvalue"][0]
    assert len(apply_region_qfilter(alone, .05)) == 1
    assert len(apply_region_qfilter(together, .05)) == 0


def test_default_qvalues_match_by_with_unreported_intervals_set_to_one():
    positions = np.arange(9) * 20 + 100
    calls = call_dmr_chain_merge(frame("chr1", positions, .02, .2),
                               min_cpgs=3, minlen_bp=50, dis_merge_bp=500)
    # All contiguous intervals with >=3 CpGs and inclusive length >=50 bp.
    family_size = sum(j - i + 1 >= 3 and positions[j] - positions[i] + 1 >= 50
                      for i in range(9) for j in range(i, 9))
    padded = np.r_[calls["combined_pvalue"].to_numpy(), np.ones(family_size - len(calls))]
    expected = multipletests(padded, method="fdr_by")[1][:len(calls)]
    np.testing.assert_allclose(calls["combined_qvalue"].to_numpy(), expected)


def test_interval_count_matches_enumeration_for_gaps_and_length_boundaries():
    rng = np.random.default_rng(705321)
    for _ in range(100):
        positions = np.cumsum(rng.integers(1, 60, 25))
        for min_cpgs, minlen, max_gap in [(1, 1, 25), (3, 50, 35), (5, 100, 60)]:
            expected = sum(
                j - i + 1 >= min_cpgs
                and positions[j] - positions[i] + 1 >= minlen
                and np.all(np.diff(positions[i:j + 1]) <= max_gap)
                for i in range(len(positions)) for j in range(i, len(positions))
            )
            assert count_intervals(positions, min_cpgs, minlen, max_gap) == expected
    assert count_intervals([], 5, 50, 500) == 0


def test_sparse_adjustment_matches_full_by_with_ties_zeros_and_missing_values():
    pvalues = np.array([.003, np.nan, .0001, .003, 1., 0.])
    # Missing scores still occupy a hypothesis in the family (assigned p=1).
    padded = np.r_[np.nan_to_num(pvalues, nan=1.), np.ones(100 - len(pvalues))]
    expected = multipletests(padded, method="fdr_by")[1][:len(pvalues)]
    actual = adjust_selected(pvalues, 100)
    np.testing.assert_allclose(actual[np.isfinite(pvalues)], expected[np.isfinite(pvalues)])
    assert np.isnan(actual[1])


def test_cache_key_changes_for_corrected_search_family(tmp_path):
    from epykit import _cache
    from epykit._dmc_store import DMCStore
    from epykit.dmr import _dmr_chain_merge_cache_key
    store = DMCStore(tmp_path, "lr", {"input_sig": "same-input"})
    legacy = _cache.fingerprint([
        ("base", "same-input"), ("alpha", "0.05"), ("delta", "0.1"),
        ("dismerge", "500"), ("mincp", "5"), ("pctsig", "0.5"),
        ("minlen", "50"), ("pcol", "pvalue"),
    ])
    revised = _dmr_chain_merge_cache_key(store, .05, .1, 500, 5, .5, 50, "pvalue")
    assert revised != legacy
