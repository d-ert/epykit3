"""impute_knn_beta optional was_imputed mask (M-SEC3).

Imputed cells must be distinguishable from observed ones so downstream code can
exclude them from variance/dispersion analyses (where kNN imputation is
destructive).
"""

from __future__ import annotations

import numpy as np

from epykit.impute import impute_knn_beta


def test_impute_returns_mask_flagging_filled_cells():
    pos = np.array([100, 200, 300, 400, 500])
    beta = np.array([
        [0.1, np.nan, 0.3, np.nan, 0.5],
        [0.2, 0.25, np.nan, 0.45, 0.5],
    ])
    out, mask = impute_knn_beta(pos, beta, k=3, return_mask=True)

    assert mask.dtype == bool
    assert mask.shape == beta.shape
    filled = np.isnan(beta) & ~np.isnan(out)
    assert np.array_equal(mask, filled)
    # The three NaNs got filled and are flagged.
    assert mask[0, 1] and mask[0, 3] and mask[1, 2]
    # Observed cells are unchanged and not flagged.
    assert not mask[0, 0]
    assert out[0, 0] == 0.1


def test_impute_default_returns_array_only():
    pos = np.array([100, 200, 300])
    beta = np.array([[0.1, np.nan, 0.3]])
    out = impute_knn_beta(pos, beta)
    assert isinstance(out, np.ndarray)
