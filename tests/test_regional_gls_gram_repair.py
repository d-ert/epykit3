"""Matrix donor products agree with independent scalar quadratics."""

import numpy as np
import pytest
from test_regional_gls_masked import dense_masked_oracle

from epykit import _regional_gls_masked as masked
from epykit._regional_gls_native import evaluate_covariance


@pytest.mark.parametrize('samples', [6, 10])
@pytest.mark.parametrize('k', [5, 20, 100])
def test_donor_products_match_scalar_car_quadratics(samples, k):
    rng = np.random.default_rng(337+samples+k)
    diagonal = rng.uniform(1, 3, k)
    off = rng.uniform(-.3, 0, k-1)
    donor = rng.normal(size=(k, samples))
    qr = rng.normal(size=k)
    q = np.diag(diagonal) + np.diag(off, 1) + np.diag(off, -1)
    expected = np.array([[donor[:, a] @ q @ donor[:, b] for b in range(samples)] for a in range(samples)])
    gram, response = masked._donor_products(diagonal, off, donor, qr)
    np.testing.assert_allclose(gram, expected, rtol=2e-8, atol=2e-8)
    np.testing.assert_allclose(response, [np.dot(donor[:, a], qr) for a in range(samples)], rtol=2e-8, atol=2e-8)


@pytest.mark.parametrize('samples,k,dropout', [(6, 5, .01), (10, 20, .2), (6, 100, .2)])
@pytest.mark.parametrize('length,c', [(0., 0.), (10., .3), (120., .3), (10000., .98)])
def test_optimized_masked_model_retains_dense_observed_cell_oracle(samples, k, dropout, length, c):
    rng = np.random.default_rng(993+samples+k)
    pos = np.cumsum(rng.integers(5, 50, k)).astype(float)
    y = rng.normal(size=(samples, k))
    h = np.sqrt(rng.integers(5, 50, k)).astype(float)
    g = np.r_[np.zeros(samples//2), np.ones(samples-samples//2)]
    mask = rng.random((samples, k)) >= dropout
    mask[[0, samples//2]] = True
    mask[-1] = False  # An entirely absent donor contributes no cells.
    got = evaluate_covariance(pos, y, h, g, length, c, obs_mask=mask)
    expected = dense_masked_oracle(pos, y, h, g, mask, length, c)
    np.testing.assert_allclose(got, expected[:3], rtol=2e-8, atol=2e-8)
