"""Independent small-array oracles for regional discovery and fit scores."""

import numpy as np
import pytest
from numba.typed import List

from epykit._regional_gls_native import (
    _refine_and_append_region,
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


def test_split_fragment_must_keep_its_own_smoothed_seed():
    from epykit._regional_gls_native import discover_prepared

    pos = np.arange(13) * 60
    raw = np.r_[np.full(5, .07), np.full(3, .01), np.full(5, .20)]
    smooth = np.r_[np.full(8, .06), np.full(5, .15)]
    geometry = prepare_discovery(pos, 500, 5, 500)

    positive = discover_prepared(
        pos, raw, *geometry, .1, 5, 500, 50, .05,
        smoothed=smooth, boundary_effects=raw,
    )
    negative = discover_prepared(
        pos, -raw, *geometry, .1, 5, 500, 50, .05,
        smoothed=-smooth, boundary_effects=-raw,
    )

    assert positive.tolist() == [[8, 13]]
    assert negative.tolist() == positive.tolist()


def test_discovery_accepts_readonly_smoothed_arrays_without_mutation():
    pos = np.arange(13) * 60
    raw = np.r_[np.full(5, .07), np.full(3, .01), np.full(5, .20)]
    smooth = np.r_[np.full(8, .06), np.full(5, .15)]
    smooth.flags.writeable = False
    geometry = prepare_discovery(pos, 500, 5, 500)
    got = discover_prepared(pos, raw, *geometry, .1, 5, 500, 50, .05,
                            smoothed=smooth, boundary_effects=raw)
    assert got.tolist() == [[8, 13]]
    np.testing.assert_array_equal(smooth, np.r_[np.full(8, .06), np.full(5, .15)])


def test_boundary_extension_does_not_lower_split_penalty():
    pos = np.arange(7) * 50
    effects = np.array([.2, .2, .052, .052, .052, .2, .2])
    for boundary in (0.0, .05, .08):
        typed_regions = List()
        typed_regions.append((np.int64(0), np.int64(0)))
        typed_regions.pop()
        _refine_and_append_region(
            pos, effects, 0, len(pos), 1, .1, 2, 50, typed_regions, boundary,
        )
        assert list(typed_regions) == [(0, 2), (5, 7)]


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


def test_seed_constrained_trimming_matches_exhaustive_intervals():
    rng = np.random.default_rng(7712)
    for direction in (-1, 1):
        for _ in range(100):
            pos = np.cumsum(rng.integers(1, 40, 12))
            effect = direction * rng.uniform(-.2, .5, 12)
            seeds = direction * rng.uniform(-.05, .2, 12)
            min_count, min_span, tau, min_seed = 3, 55, .08, .1
            admissible = [
                (sum(direction * effect[a:b] - tau), a, b)
                for a in range(12)
                for b in range(a + min_count, 13)
                if pos[b - 1] - pos[a] + 1 >= min_span
                and np.any(direction * seeds[a:b] >= min_seed)
            ]
            out = []
            _trim_and_append_subregion.py_func(
                pos, effect, 0, 12, direction, tau, min_count, min_span,
                out, seeds, min_seed,
            )
            if admissible:
                best_score = max(item[0] for item in admissible)
                got = out[0]
                got_score = sum(direction * effect[got[0]:got[1]] - tau)
                assert got_score == pytest.approx(best_score, abs=1e-12)
                assert np.any(direction * seeds[got[0]:got[1]] >= min_seed)
                assert pos[got[1] - 1] - pos[got[0]] + 1 >= min_span
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


def test_seed_threshold_equality_and_tied_optima_keep_first_interval():
    for direction in (-1, 1):
        out = []
        _trim_and_append_subregion.py_func(
            np.arange(5)*20, direction*np.full(5, .08), 0, 5,
            direction, .08, 3, 1, out,
            direction*np.array([0., 0., 0., .1, 0.]), .1,
        )
        assert out == [(0, 4)]
