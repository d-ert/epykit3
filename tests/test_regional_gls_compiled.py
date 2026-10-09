"""Independent numerical checks for native spatial REML acceleration."""

import numpy as np
import pytest
from scipy.optimize import minimize_scalar


def fixture(k=20, samples=6, seed=935):
    from epykit._regional_gls import transform_counts

    rng = np.random.default_rng(seed)
    p = np.cumsum(rng.integers(1, 80, k)).astype(float)
    n = rng.integers(5, 80, (k, samples))
    m = rng.binomial(n, rng.uniform(0.1, 0.9, (k, 1)))
    y = transform_counts(m, n).T
    h = np.sqrt(np.median(n, axis=1))
    return (
        p,
        (y - y.mean(axis=0)) * h,
        h,
        np.r_[np.zeros(samples // 2), np.ones(samples - samples // 2)],
    )


def dense_objective(p, yc, h, g, length, c):
    k = len(p)
    a = np.exp(-abs(p[:, None] - p[None, :]) / length) if length else np.eye(k)
    u = h / np.sqrt(np.median(h * h))
    r = (1 - c) * a + c * np.outer(u, u)
    q = np.linalg.inv(r)
    x = g - g.mean()
    gram = x @ x
    info = h @ q @ h
    effect = x @ yc @ q @ h / (gram * info)
    rss = np.sum(yc * (yc @ q)) - gram * info * effect**2
    df = (len(g) - 1) * k - 1
    value = df * np.log(rss / df) + (len(g) - 1) * np.linalg.slogdet(r)[1] + np.log(info)
    se = np.sqrt(rss / df / (gram * info))
    return value, effect, se


@pytest.mark.parametrize("k,samples", [(5, 6), (30, 10), (80, 6)])
def test_native_fixed_covariance_matches_dense_model(k, samples):
    from epykit._regional_gls_native import evaluate_covariance

    p, yc, h, g = fixture(k, samples)
    fast = evaluate_covariance(p, yc, h, g, 120.0, 0.3)
    dense = dense_objective(p, yc, h, g, 120.0, 0.3)
    assert fast == pytest.approx(dense, rel=1e-10, abs=1e-9)


@pytest.mark.parametrize("length", [0.0, 10.0, 250.0, 3000.0])
def test_closed_form_donor_variance_matches_scalar_numerical_reml(length):
    from epykit._regional_gls_native import profile_covariance

    p, yc, h, g = fixture()
    value, c, effect, se = profile_covariance(p, yc, h, g, length)
    numeric = minimize_scalar(
        lambda c: dense_objective(p, yc, h, g, length, c)[0],
        bounds=(0.0, 0.98),
        method="bounded",
        options={"xatol": 1e-12},
    )
    best = min(
        numeric.fun,
        dense_objective(p, yc, h, g, length, 0.0)[0],
        dense_objective(p, yc, h, g, length, 0.98)[0],
    )
    assert value == pytest.approx(best, abs=1e-7)
    assert (value, effect, se) == pytest.approx(dense_objective(p, yc, h, g, length, c), abs=1e-8)


def test_native_batch_parallel_matches_serial_and_mirror():
    from epykit._regional_gls_native import fit_regions

    p, yc, h, g = fixture(80)
    bounds = np.array([[0, 20], [20, 50], [50, 80]], np.int64)
    serial = fit_regions(p, yc, h, g, bounds, workers=1)
    parallel = fit_regions(p, yc, h, g, bounds, workers=2)
    mirror = fit_regions(p, yc, h, 1 - g, bounds, workers=1)
    assert serial == pytest.approx(parallel, abs=1e-10)
    assert serial[:, 2] == pytest.approx(-mirror[:, 2], abs=1e-10)
    assert serial[:, [0, 1, 3, 4]] == pytest.approx(mirror[:, [0, 1, 3, 4]], abs=1e-10)


@pytest.mark.parametrize("window", [100, 500, 1600])
def test_native_discovery_matches_segmented_numpy_geometry(window):
    from epykit._regional_gls import discover_regions
    from epykit._regional_gls_native import discover_prepared, prepare_discovery

    rng = np.random.default_rng(9021)
    pos = np.cumsum(rng.integers(1, 700, 10000))
    effects = rng.uniform(-0.4, 0.4, len(pos))
    geom = prepare_discovery(pos, window, 5, 500)
    got = discover_prepared(pos, effects, *geom, 0.1, 5, 500, 50)
    expected = discover_regions(pos, effects, window_bp=window)
    assert got.tolist() == [list(r) for r in expected]
