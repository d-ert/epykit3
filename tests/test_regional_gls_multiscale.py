"""Precision-kernel discovery should resolve short and nearby opposing effects."""

import numpy as np


def test_multiscale_recovers_narrow_signal_diluted_by_five_hundred_bp_boxcar():
    from epykit._regional_gls_discovery import discover_multiscale
    from epykit._regional_gls_native import discover_prepared, prepare_discovery

    pos = np.arange(0, 1000, 20)
    effect = np.zeros(len(pos))
    effect[15:21] = .25
    geometry = prepare_discovery(pos, 500, 5, 500)
    old = discover_prepared(pos, effect, *geometry, .1, 5, 500, 50)
    assert old.tolist() == []
    new = discover_multiscale(pos, effect, np.ones(len(pos)), 500, .1, 5, 500, 50)
    assert new.tolist() == [[15, 21]]


def test_multiscale_retains_nearby_opposite_regions_without_scale_duplicates():
    from epykit._regional_gls_discovery import discover_multiscale

    pos = np.arange(0, 1200, 20)
    effect = np.zeros(len(pos))
    effect[15:21] = .25
    effect[25:31] = -.25
    bounds = discover_multiscale(pos, effect, np.ones(len(pos)), 500, .1, 5, 500, 50)
    assert bounds.tolist() == [[15, 21], [25, 31]]
    mirror = discover_multiscale(pos, -effect, np.ones(len(pos)), 500, .1, 5, 500, 50)
    assert mirror.tolist() == bounds.tolist()


def test_extreme_depth_cannot_turn_an_isolated_site_into_a_supported_region():
    from epykit._regional_gls_discovery import discover_multiscale

    pos = np.arange(0, 1000, 20)
    effect = np.zeros(len(pos))
    effect[20] = .2
    precision = np.ones(len(pos))
    precision[20] = 1e6
    assert discover_multiscale(pos, effect, precision, 500, .1, 5, 500, 50).tolist() == []
