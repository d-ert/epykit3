"""False-discovery estimates use full scan exposure, not pooled region count."""

import numpy as np
import pytest

from epykit import _regional_gls


def test_empty_scans_change_expected_false_call_exposure():
    q = _regional_gls._permutation_fdp_qvalues([5.0], [[6.0]])
    more_exposure = _regional_gls._permutation_fdp_qvalues([5.0], [[6.0], [], [], []])
    assert q.tolist() == [1.0]
    assert more_exposure.tolist() == [.25]


def test_low_scoring_decoy_regions_do_not_change_high_threshold_fdp():
    a = _regional_gls._permutation_fdp_qvalues([5.0], [[6.0], [], []])
    b = _regional_gls._permutation_fdp_qvalues([5.0], [[6.0]+[.1]*1000, [], []])
    assert a == pytest.approx([1/3])
    assert b == pytest.approx(a)


def test_complete_scan_duplication_is_invariant_but_false_region_duplication_is_not():
    scans = [[6.0], []]
    a = _regional_gls._permutation_fdp_qvalues([5.0], scans)
    b = _regional_gls._permutation_fdp_qvalues([5.0], scans+scans)
    c = _regional_gls._permutation_fdp_qvalues([5.0], [[6.0, 6.0], []])
    assert a == pytest.approx([.5])
    assert b == pytest.approx(a)
    assert c == pytest.approx([1.0])


def test_score_ties_receive_same_monotone_fdp_estimate():
    q = _regional_gls._permutation_fdp_qvalues([5., 3., 5.], [[5.], [5.], []])
    assert q == pytest.approx([2/9, 2/9, 2/9])
    assert np.isfinite(q).all()


def test_zero_tail_is_an_estimate_and_failed_scores_cannot_be_dropped():
    assert _regional_gls._permutation_fdp_qvalues([5.], [[]]).tolist() == [0.]
    with pytest.raises(ValueError, match="Failed"):
        _regional_gls._permutation_fdp_qvalues([5.], [[float("nan")]])


def test_maximum_statistic_tail_includes_empty_scans_and_observed_mirror():
    p = _regional_gls._maximum_scan_pvalues([10., 5.], [[9.], [], [4.]], observed_copies=2)
    assert p == pytest.approx([.4, .6])


def test_three_versus_three_exact_maximum_tail_cannot_claim_five_percent():
    p = _regional_gls._maximum_scan_pvalues([100.], [[] for _ in range(18)], observed_copies=2)
    assert p == pytest.approx([.1])
    assert not np.any(p <= .05)
