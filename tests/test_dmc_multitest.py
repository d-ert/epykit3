"""P1-10: Storey pi0 clamped at 1/n."""
from __future__ import annotations

import numpy as np
import pytest

from epykit.dmc import _storey_pi0


def test_storey_pi0_clamped_at_one_over_n():
    """When all p-values are below lam=0.5, numerator = 0 -> pi0_hat = 0.
    Post-fix: clamp at 1/n (Storey's standard floor)."""
    rng = np.random.default_rng(0)
    n = 1000
    pvals = rng.uniform(0, 0.4, size=n)  # all p < 0.5 = lam
    pi0 = _storey_pi0(pvals)
    assert pi0 >= 1.0 / n - 1e-12, (
        f"_storey_pi0 returned {pi0}; expected >= 1/n = {1.0/n:.6f}"
    )
    assert pi0 <= 1.0, f"pi0={pi0} exceeds 1.0"


def test_storey_pi0_unclamped_when_above_floor():
    """On uniform p-values, pi0 should be close to 1 (not clamped)."""
    rng = np.random.default_rng(1)
    pvals = rng.uniform(0, 1, size=1000)
    pi0 = _storey_pi0(pvals)
    assert 0.9 <= pi0 <= 1.1, f"expected pi0 ~ 1 on uniform; got {pi0}"
