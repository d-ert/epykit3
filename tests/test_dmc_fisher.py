"""P1-1: vectorised Fisher two-sided p must match
scipy.stats.fisher_exact(alternative='two-sided') to machine precision."""
from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import fisher_exact

from epykit.dmc import fisher_exact_vectorized


def test_fisher_two_sided_matches_scipy_on_random_tables():
    rng = np.random.default_rng(0)
    n = 100
    a = rng.integers(0, 40, size=n).astype(np.int64)
    b = rng.integers(0, 40, size=n).astype(np.int64)
    c = rng.integers(0, 40, size=n).astype(np.int64)
    d = rng.integers(0, 40, size=n).astype(np.int64)
    # Avoid degenerate tables.
    keep = (a + b > 0) & (c + d > 0) & (a + c > 0) & (b + d > 0)
    a, b, c, d = a[keep], b[keep], c[keep], d[keep]

    epy_p, _epy_log2_or = fisher_exact_vectorized(a, b, c, d)
    ref_p = np.array([
        fisher_exact([[ai, bi], [ci, di]], alternative="two-sided")[1]
        for ai, bi, ci, di in zip(a, b, c, d, strict=True)
    ])
    np.testing.assert_allclose(
        epy_p, ref_p, atol=1e-12, rtol=1e-9,
        err_msg="vectorised Fisher two-sided p must match scipy reference",
    )
