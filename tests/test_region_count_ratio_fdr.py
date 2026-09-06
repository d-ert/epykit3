"""Unit tests for the shared count-ratio (target-decoy) region FDR helper.

Hand-computable fixtures. The helper is caller-agnostic: it takes observed
survivor p-values plus per-permutation null survivor pools and returns
``(empirical_pvalue, empirical_qvalue, fdr_set)``. ``_is_self_or_mirror_perm``
is the assignment filter the region mode applies before pooling.
"""

from __future__ import annotations

import numpy as np
import pytest

from epykit.dmr import _is_self_or_mirror_perm, _region_count_ratio_fdr


def test_self_perm_detected():
    """A shuffle reproducing the observed treatment set is the observed
    contrast itself, not a null draw."""
    assert _is_self_or_mirror_perm(
        perm_treatment=["b", "a", "c"],  # same set as observed treatment
        observed_treatment=["a", "b", "c"],
        observed_control=["d", "e", "f"],
    )


def test_mirror_perm_detected():
    """A shuffle assigning the observed control samples as treatment is the
    mirror swap: identical two-sided statistics."""
    assert _is_self_or_mirror_perm(
        perm_treatment=["e", "d", "f"],  # == observed control set
        observed_treatment=["a", "b", "c"],
        observed_control=["d", "e", "f"],
    )


def test_genuine_perm_not_excluded():
    assert not _is_self_or_mirror_perm(
        perm_treatment=["a", "d", "e"],  # mixed -> valid null draw
        observed_treatment=["a", "b", "c"],
        observed_control=["d", "e", "f"],
    )


def test_unequal_groups_have_no_mirror():
    """With 2 vs 3 samples the control set can never be a treatment draw, so
    only the exact self assignment is excluded."""
    assert not _is_self_or_mirror_perm(
        perm_treatment=["c", "d"],
        observed_treatment=["a", "b"],
        observed_control=["c", "d", "e"],
    )
    assert _is_self_or_mirror_perm(
        perm_treatment=["b", "a"],
        observed_treatment=["a", "b"],
        observed_control=["c", "d", "e"],
    )


def test_count_ratio_basic_hand_computed():
    """R=3 observed survivors, 2 permutations.

      observed p   = [0.001, 0.01, 0.02]           -> R = 3
      perm1 nulls  = [0.005, 0.5]   (2 survivors)
      perm2 nulls  = [0.9]          (1 survivor)
      pooled null  = [0.005, 0.5, 0.9]  -> N_null = 3

    set FDR = mean(2, 1) / 3 = 1.5/3 = 0.5

    empirical_pvalue = (#pooled <= p) / N_null
      0.001 -> 0/3 = 0
      0.01  -> 1/3
      0.02  -> 1/3

    count-ratio q: fdr(t) = (#pooled<=t / n_perm) / (#obs<=t), suffix-min
      t=0.001: V=0/2,  R=1 -> 0
      t=0.01 : V=1/2,  R=2 -> 0.25
      t=0.02 : V=1/2,  R=3 -> 0.16667
      suffix-min -> [0, 0.16667, 0.16667]
    """
    observed = np.array([0.001, 0.01, 0.02])
    null_pools = [np.array([0.005, 0.5]), np.array([0.9])]

    emp_p, emp_q, fdr_set = _region_count_ratio_fdr(
        observed_pvalues=observed,
        null_pools=null_pools,
        n_perm_used=2,
    )

    assert fdr_set == pytest.approx(0.5)
    assert emp_p == pytest.approx([0.0, 1 / 3, 1 / 3])
    assert emp_q == pytest.approx([0.0, 1 / 6, 1 / 6])


def test_count_ratio_preserves_input_order():
    """Output must align with the (unsorted) observed input order, not a
    sorted copy."""
    observed = np.array([0.02, 0.001, 0.01])  # deliberately unsorted
    null_pools = [np.array([0.005, 0.5]), np.array([0.9])]

    emp_p, emp_q, _ = _region_count_ratio_fdr(
        observed_pvalues=observed, null_pools=null_pools, n_perm_used=2
    )
    # same per-region values as the sorted case, reordered to match input
    assert emp_q == pytest.approx([1 / 6, 0.0, 1 / 6])
    assert emp_p == pytest.approx([1 / 3, 0.0, 1 / 3])


def test_count_ratio_q_capped_at_one_when_null_exceeds_observed():
    """When shuffles yield more survivors than observed, FDR saturates at 1."""
    observed = np.array([0.01, 0.02])  # R = 2
    null_pools = [np.array([0.001, 0.005, 0.008, 0.015])]  # 4 nulls, 1 perm
    _, emp_q, fdr_set = _region_count_ratio_fdr(
        observed_pvalues=observed, null_pools=null_pools, n_perm_used=1
    )
    assert fdr_set == pytest.approx(1.0)  # min(4/2, 1) = 1
    assert np.all(emp_q <= 1.0)
    assert emp_q[-1] == pytest.approx(1.0)  # loosest threshold = set FDR


def test_count_ratio_zero_survivor_perm_counts_as_zero():
    """A clean zero-survivor permutation is evidence of low noise: it stays in
    the divisor and contributes nothing to the pooled null.

      observed = [0.01, 0.02], pools = [[0.005], []] -> V(t) = 1/2 at both t
      fdr(0.01) = 0.5 / 1 = 0.5 ; fdr(0.02) = 0.5 / 2 = 0.25 ; suffix-min
      set FDR = mean(1, 0) / 2 = 0.25
    """
    observed = np.array([0.01, 0.02])
    emp_p, emp_q, fdr_set = _region_count_ratio_fdr(
        observed_pvalues=observed,
        null_pools=[np.array([0.005]), np.array([], dtype=np.float64)],
        n_perm_used=2,
    )
    assert fdr_set == pytest.approx(0.25)
    assert emp_q == pytest.approx([0.25, 0.25])
    assert emp_p == pytest.approx([1.0, 1.0])  # one pooled null <= both


def test_count_ratio_all_empty_pools_give_zero_fdr():
    """Only zero-survivor permutations: nothing in the pooled null, so every
    observed region has q = 0 and the set-level FDR is 0."""
    observed = np.array([0.3, 0.01])
    emp_p, emp_q, fdr_set = _region_count_ratio_fdr(
        observed_pvalues=observed,
        null_pools=[np.array([], dtype=np.float64)] * 3,
        n_perm_used=3,
    )
    assert fdr_set == 0.0
    assert emp_q == pytest.approx([0.0, 0.0])
    assert emp_p == pytest.approx([0.0, 0.0])


def test_count_ratio_empty_observed_returns_empty():
    emp_p, emp_q, fdr_set = _region_count_ratio_fdr(
        observed_pvalues=np.array([]), null_pools=[np.array([0.1])], n_perm_used=1
    )
    assert emp_p.size == 0 and emp_q.size == 0
    assert np.isnan(fdr_set)


def test_count_ratio_no_valid_perms_returns_nan_q():
    """Zero usable permutations -> q is NaN (cannot estimate), not 0."""
    observed = np.array([0.01, 0.02])
    emp_p, emp_q, fdr_set = _region_count_ratio_fdr(
        observed_pvalues=observed, null_pools=[], n_perm_used=0
    )
    assert np.all(np.isnan(emp_q))
    assert np.all(np.isnan(emp_p))
    assert np.isnan(fdr_set)


def test_count_ratio_q_is_monotone_nondecreasing_in_p():
    """A proper q-value never decreases as the p-value threshold loosens."""
    rng = np.random.default_rng(0)
    observed = np.sort(rng.uniform(0, 0.05, size=200))
    null_pools = [rng.uniform(0, 1, size=150) for _ in range(10)]
    _, emp_q, _ = _region_count_ratio_fdr(
        observed_pvalues=observed, null_pools=null_pools, n_perm_used=10
    )
    assert np.all(np.diff(emp_q) >= -1e-12)
