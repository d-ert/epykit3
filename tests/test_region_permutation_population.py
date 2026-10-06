"""Regression checks for the full label-assignment population, including failures."""
import numpy as np
import pytest
from epykit.dmr import _aggregate_region_perm_results


def aggregate(results, mode='region', observed=None):
    return _aggregate_region_perm_results(observed_pvalues=np.array([.001,.002]) if observed is None else observed,
        results=results,n_perm=len(results),fdr_method=mode)


def test_identity_and_reverse_labels_cannot_disappear_from_region_null():
    # Independent count: two observed regions recur in two of four draws,
    # so average null discoveries = 1, observed discoveries = 2, FDP = .5.
    scans=[(True,np.array([.001,.002])),(True,np.array([.001,.002])),
           (False,np.array([])),(False,np.array([]))]
    _,q,fdp=aggregate(scans)
    np.testing.assert_allclose(q,[.5,.5])
    assert fdp==.5


@pytest.mark.parametrize('mode',['region','max_t'])
@pytest.mark.parametrize('is_self',[False,True])
def test_any_failed_assignment_aborts_inference(mode,is_self):
    with pytest.raises(RuntimeError,match='Incomplete'):
        aggregate([(is_self,None),(False,np.array([]))],mode)


@pytest.mark.parametrize('mode',['region','max_t'])
@pytest.mark.parametrize('bad',[np.array([np.nan]),np.array([np.inf]),np.array([-1.]),np.array([1.01]),np.array([[.1]])])
def test_invalid_null_scores_cannot_produce_significance(mode,bad):
    with pytest.raises(ValueError,match='Invalid permutation'):
        aggregate([(True,bad)],mode)


@pytest.mark.parametrize('mode',['region','max_t'])
def test_incomplete_result_count_is_rejected(mode):
    with pytest.raises(RuntimeError,match='Incomplete'):
        _aggregate_region_perm_results(observed_pvalues=np.array([.01]),
            results=[(False,np.array([]))],n_perm=2,fdr_method=mode)


@pytest.mark.parametrize('mode',['region','max_t'])
def test_out_of_range_observed_score_rejected(mode):
    with pytest.raises(ValueError,match='Invalid observed'):
        aggregate([(False,np.array([]))],mode,np.array([-.1]))
