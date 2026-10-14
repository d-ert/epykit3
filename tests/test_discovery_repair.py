"""Independent small-array oracles for regional discovery and fit scores."""

import numpy as np
import pytest

from epykit._regional_gls_native import (
    _trim_and_append_subregion,
    discover_prepared,
    fit_regions,
    prepare_discovery,
    profile_covariance,
)


def test_minimum_span_trimming_retains_best_admissible_prefix():
    pos = np.arange(7) * 15
    effect = np.array([.045, .015, .015, .015, .015, .045, .8])
    geometry = prepare_discovery(pos, 500, 3, 500)
    got = discover_prepared(pos, effect, *geometry, .1, 3, 500, 50)
    assert got.tolist() == [[2, 7]]


def test_constrained_trimming_matches_exhaustive_admissible_segments():
    rng = np.random.default_rng(895)
    for _ in range(100):
        pos = np.cumsum(rng.integers(1, 40, 12))
        effect = rng.uniform(-.2, .5, 12)
        min_count, min_span, tau = 3, 55, .08
        candidates = [(sum(effect[a:b] - tau), a, b)
                      for a in range(12) for b in range(a + min_count, 13)
                      if pos[b-1] - pos[a] + 1 >= min_span]
        out = []
        _trim_and_append_subregion.py_func(pos, effect, 0, 12, 1, tau, min_count, min_span, out)
        if candidates:
            best = max(candidates, key=lambda v: v[0])
            assert len(out) == 1
            a, b = out[0]
            assert sum(effect[a:b] - tau) == pytest.approx(best[0], abs=1e-12)
        else:
            assert out == []


def test_batch_returns_profiled_standard_error_without_length_boost():
    rng = np.random.default_rng(935)
    pos = np.cumsum(rng.integers(1, 80, 20)).astype(float)
    y = rng.normal(size=(6, 20))
    y -= y.mean(axis=0)
    h = np.ones(20)
    g = np.r_[np.zeros(3), np.ones(3)]
    row = fit_regions(pos, y, h, g, [[0, 20]], workers=1)[0]
    raw = profile_covariance(pos, y, h, g, row[0])
    assert row[3] == pytest.approx(raw[3], abs=1e-10)
