"""NaN (masked) p-values must be EXCLUDED from the FDR denominator, not
counted as p=1.0.

Regression test for M6. Sites masked to NaN (min-samples guard, degenerate
sites) are untested hypotheses, not p=1 observations. The old code filled
them with 1.0 and fed the full vector to ``multipletests`` / Storey, so they
inflated the BH denominator ``n`` and made every real q-value conservative
(a power loss). The in-code comment claimed they were "effectively excluded";
this pins that they actually are.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest
from statsmodels.stats.multitest import multipletests

from epykit.dmc import apply_multiple_testing_correction

REAL_P = np.array(
    [0.0001, 0.001, 0.008, 0.02, 0.03, 0.2, 0.4, 0.6, 0.8, 0.99],
    dtype=np.float64,
)


@pytest.mark.parametrize("method", ["fdr_bh", "fdr_by", "fdr_tsbh", "fdr_storey"])
def test_nan_pvalues_do_not_inflate_denominator(method):
    # Reference: the FDR procedure over ONLY the real (finite) p-values.
    if method == "fdr_storey":
        from epykit.dmc import _apply_storey_qvalues
        _, q_ref, _ = _apply_storey_qvalues(REAL_P)
    else:
        _, q_ref, _, _ = multipletests(REAL_P, method=method)

    # Same real p-values plus 990 masked (NaN) sites. The masked sites must
    # not change the q-values of the finite ones.
    p = np.concatenate([REAL_P, np.full(990, np.nan)])
    out = apply_multiple_testing_correction(
        pl.DataFrame({"pvalue": p}), method=method
    )
    q = out["qvalue"].to_numpy()

    # Finite q-values match the n=10 reference (denominator is 10, not 1000).
    assert np.allclose(q[: len(REAL_P)], q_ref, atol=1e-12, equal_nan=True)
    # Masked sites carry NaN q and are not rejected.
    assert np.all(np.isnan(q[len(REAL_P) :]))
    assert not out["reject"].to_numpy()[len(REAL_P) :].any()


def test_no_nan_is_byte_identical_to_plain_bh():
    # With zero masked sites (the intersect-mode benchmark case), the result
    # must be unchanged vs a direct statsmodels BH -- so paper numbers that
    # use intersect mode are unaffected by this fix.
    _, q_ref, _, _ = multipletests(REAL_P, method="fdr_bh")
    out = apply_multiple_testing_correction(
        pl.DataFrame({"pvalue": REAL_P}), method="fdr_bh"
    )
    assert np.allclose(out["qvalue"].to_numpy(), q_ref, atol=1e-12)
