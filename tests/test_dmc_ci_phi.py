"""phi-aware Newcombe CI for the lr engine (M1).

The quasi-binomial ``lr`` p-value divides its statistic by a per-site
McCullagh-Nelder dispersion ``phi`` and references F(1, df) / chi^2(1). Before
M1 the ``meth_diff`` CI was a plain *binomial* Newcombe interval on the pooled
counts -- it ignored ``phi`` entirely and was anti-conservatively narrow next
to the dispersion-aware p-value (a CI could exclude 0 while q > 0.05, or the
reverse). These tests pin the phi-aware behaviour and the backward-compatible
``phi=None`` path (used by ``fisher``, which has no dispersion estimate).
"""
from __future__ import annotations

import numpy as np

from epykit._glm import newcombe_diff_ci


def _widths(lo, hi):
    return np.asarray(hi) - np.asarray(lo)


def test_phi_none_is_binomial_baseline():
    """phi=None reproduces the exact binomial Newcombe interval (fisher path)."""
    rng = np.random.default_rng(0)
    cov_a = rng.integers(20, 200, size=50).astype(float)
    cov_b = rng.integers(20, 200, size=50).astype(float)
    meth_a = np.round(cov_a * 0.6)
    meth_b = np.round(cov_b * 0.4)

    lo0, hi0 = newcombe_diff_ci(meth_a, cov_a, meth_b, cov_b)
    lo1, hi1 = newcombe_diff_ci(meth_a, cov_a, meth_b, cov_b, phi=None, df=None)
    assert np.allclose(lo0, lo1) and np.allclose(hi0, hi1)


def test_phi_clamped_to_one_matches_binomial():
    """Where phi == 1 (clamped, the quasi-binomial collapses to a binomial)
    the interval is identical to the phi-free one."""
    cov = np.full(20, 100.0)
    meth_a = np.full(20, 60.0)
    meth_b = np.full(20, 40.0)
    df = np.full(20, 1e6)  # large df -> t ~= z anyway

    lo_base, hi_base = newcombe_diff_ci(meth_a, cov, meth_b, cov)
    lo_phi, hi_phi = newcombe_diff_ci(
        meth_a, cov, meth_b, cov, phi=np.ones(20), df=df,
    )
    assert np.allclose(lo_base, lo_phi) and np.allclose(hi_base, hi_phi)


def test_overdispersion_widens_interval_toward_sqrt_phi():
    """At phi > 1 (with large df so t ~= z) each half-width grows toward
    sqrt(phi); the interval is strictly wider than the phi-free one."""
    n = 30
    cov = np.full(n, 100.0)
    meth_a = np.full(n, 60.0)
    meth_b = np.full(n, 40.0)
    df = np.full(n, 1e6)        # isolate the sqrt(phi) factor (t -> z)
    phi = np.full(n, 4.0)       # sqrt(phi) = 2

    lo0, hi0 = newcombe_diff_ci(meth_a, cov, meth_b, cov)             # phi-free
    lo1, hi1 = newcombe_diff_ci(meth_a, cov, meth_b, cov, phi=phi, df=df)

    w0, w1 = _widths(lo0, hi0), _widths(lo1, hi1)
    assert np.all(w1 > w0), "overdispersed CI must be wider than the binomial one"
    ratio = w1 / w0
    # Wilson half-widths are sub-linear in z (the 1 + z^2/n denominator
    # grows too), so the ratio lands below sqrt(phi)=2 but well above 1.
    assert np.all(ratio > 1.3) and np.all(ratio < 2.0)


def test_small_df_widens_further_than_normal():
    """A small df (e.g. n=3+3 -> df~=4, floored to 50 by the caller) uses a
    t tail, which is slightly wider than the normal z at the same phi."""
    n = 10
    cov = np.full(n, 80.0)
    meth_a = np.full(n, 50.0)
    meth_b = np.full(n, 30.0)
    phi = np.full(n, 3.0)

    _, hi_bigdf = newcombe_diff_ci(meth_a, cov, meth_b, cov, phi=phi, df=np.full(n, 1e6))
    _, hi_smalldf = newcombe_diff_ci(meth_a, cov, meth_b, cov, phi=phi, df=np.full(n, 50.0))
    # t(50) > z, so the floored-df interval is at least as wide.
    assert np.all(hi_smalldf >= hi_bigdf - 1e-9)


def test_ci_bracketed_and_bounded():
    """phi-inflated intervals stay bracketing and clamped to [-1, 1]."""
    rng = np.random.default_rng(7)
    cov_a = rng.integers(5, 300, size=200).astype(float)
    cov_b = rng.integers(5, 300, size=200).astype(float)
    meth_a = np.round(cov_a * rng.uniform(0.0, 1.0, size=200))
    meth_b = np.round(cov_b * rng.uniform(0.0, 1.0, size=200))
    phi = rng.uniform(0.5, 6.0, size=200)
    df = np.full(200, 50.0)

    lo, hi = newcombe_diff_ci(meth_a, cov_a, meth_b, cov_b, phi=phi, df=df)
    diff = meth_a / cov_a - meth_b / cov_b
    ok = np.isfinite(lo) & np.isfinite(hi)
    assert np.all(lo[ok] <= diff[ok] + 1e-9) and np.all(diff[ok] <= hi[ok] + 1e-9)
    assert np.all(lo[ok] >= -1.0 - 1e-9) and np.all(hi[ok] <= 1.0 + 1e-9)
