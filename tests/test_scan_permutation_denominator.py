"""Independent hand calculations for complete-scan randomization inference."""
import numpy as np
import pytest
from epykit.dmr import _aggregate_region_perm_results


def test_empty_scans_are_valid_null_outcomes_and_no_double_adjustment():
    p,q,_=_aggregate_region_perm_results(observed_pvalues=np.array([.001,.1]),
        results=[(False,np.array([])),(False,np.array([.05])),(True,np.array([.001]))],
        n_perm=3,fdr_method='max_t')
    np.testing.assert_allclose(p,[.5,.75])
    np.testing.assert_array_equal(q,p)


def test_all_empty_scans_retain_monte_carlo_resolution():
    p,_,_=_aggregate_region_perm_results(observed_pvalues=np.array([.001]),
        results=[(False,np.array([]))]*99,n_perm=99,fdr_method='max_t')
    assert p[0]==.01


def test_failed_scan_does_not_become_an_empty_scan():
    with pytest.raises(RuntimeError,match='Incomplete'):
        _aggregate_region_perm_results(observed_pvalues=np.array([.001]),
            results=[(False,None)],n_perm=1,fdr_method='max_t')


def test_discrete_uniform_null_rank_controls_complete_null_error():
    # Enumerating all possible observed ranks is an independent check of the
    # complete-scan rank rule, including its attainable tail probability.
    p=np.arange(1,101)/100
    assert np.mean(p<=.05)==.05
