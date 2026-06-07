"""m-perm-1: failed permutations must NOT count in the empirical-p
denominator. Pre-fix: denominator was always n_perm+1, biasing empirical p
downward by the failure rate (which can be non-trivial under small-n DMR
calling where some permuted splits produce zero candidate regions)."""
import numpy as np
from epykit.dmr import _empirical_pvalues_from_null_pool


def test_failed_permutations_excluded_from_denominator():
    observed = np.array([1e-6, 1e-3, 0.05])
    null_pool = np.array([1e-4, 1e-2, 0.5])  # only 3 successful perms
    n_perm_successful = 3

    emp = _empirical_pvalues_from_null_pool(
        observed_pvalues=observed,
        null_pvalues_pool=null_pool,
        n_perm_used=n_perm_successful,
    )
    # For observed[0]=1e-6: 0 nulls <= it -> emp = (0+1)/(3+1) = 0.25.
    # The pre-fix code used (10+1) in the denominator -> 0.0909, biased low.
    assert abs(emp[0] - 0.25) < 1e-9, (
        f"Denominator must be n_perm_successful+1=4 -> emp=0.25; got {emp[0]}"
    )


def test_zero_successful_perms_returns_ones():
    """Defensive: if every permutation failed, empirical p must be 1.0
    (the most conservative possible value), not div-by-zero or undefined."""
    observed = np.array([0.01, 0.5])
    null_pool = np.array([], dtype=np.float64)
    emp = _empirical_pvalues_from_null_pool(
        observed_pvalues=observed,
        null_pvalues_pool=null_pool,
        n_perm_used=0,
    )
    np.testing.assert_array_equal(emp, np.array([1.0, 1.0]))


def test_searchsorted_counts_match_expected_p():
    """Sanity: every observed p gets the right (count<=obs + 1)/(n+1)."""
    observed = np.array([0.001, 0.05, 0.1, 0.5])
    null_pool = np.array([0.0005, 0.01, 0.02, 0.06, 0.07, 0.2, 0.3, 0.6])
    n_perm_successful = 4
    emp = _empirical_pvalues_from_null_pool(
        observed_pvalues=observed,
        null_pvalues_pool=null_pool,
        n_perm_used=n_perm_successful,
    )
    # counts: obs=0.001 -> 1 (only 0.0005); 0.05 -> 3; 0.1 -> 5; 0.5 -> 7.
    # denominator = 4+1 = 5
    expected = np.array([(1 + 1) / 5, (3 + 1) / 5, (5 + 1) / 5, (7 + 1) / 5])
    np.testing.assert_allclose(emp, expected, atol=1e-12)
