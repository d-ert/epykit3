"""Geometry-defined family counts and standard BY reference adjustment."""

import numpy as np
import pytest
from statsmodels.stats.multitest import multipletests

from epykit._region_search import adjust_selected, count_intervals


def test_interval_family_counts_all_endpoints_for_an_anchor():
    assert count_intervals(np.arange(10) * 20, 5, 50, 500) == 21


def test_no_admissible_interval_means_empty_family():
    assert count_intervals(np.arange(10) * 20, 5, 1000, 500) == 0


def test_interval_count_matches_brute_force_with_gaps():
    rng = np.random.default_rng(87)
    for _ in range(40):
        pos = np.cumsum(rng.integers(1, 300, 20))
        want = sum(b-a+1 >= 3 and pos[b]-pos[a]+1 >= 50
                   and np.all(np.diff(pos[a:b+1]) <= 200)
                   for a in range(len(pos)) for b in range(a, len(pos)))
        assert count_intervals(pos, 3, 50, 200) == want


def test_selected_adjustment_matches_complete_by_family():
    p = np.array([.005, .003, .2])
    full = np.r_[p, np.ones(18)]
    want = multipletests(full, method="fdr_by")[1][:3]
    assert adjust_selected(p, 21) == pytest.approx(want)
