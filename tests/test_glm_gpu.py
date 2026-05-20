"""Parity tests for the GPU IRLS backend.

The CPU implementation in :mod:`epykit._glm` is the source of truth; the
GPU implementation in :mod:`epykit._glm_gpu` must produce numerically
identical results (within float64 round-off) for the same input. These
tests run only when CuPy is importable -- they cleanly skip on CPU-only
runners so CI doesn't fail.

Tolerance rationale: both paths run the same algorithm with the same
linear algebra primitives at float64; differences come only from the
order of GPU floating-point reductions. 1e-6 is well above that noise
floor and well below scientific significance.
"""

from __future__ import annotations

import numpy as np
import pytest

# Skip the whole module when cupy is not installed. tests on a CPU-only
# CI machine will print a skip line and move on.
cp = pytest.importorskip("cupy", reason="CuPy required for GPU IRLS parity tests")

from epykit import _glm
from epykit._glm_gpu import irls_binomial_batch_gpu


def _make_glm_inputs(n_sites=512, n_samples=8, p=3, seed=42):
    """Build deterministic synthetic GLM inputs.

    A small (~few hundred sites, few samples) batch keeps the GPU
    parity test fast -- the goal is bit-equivalence, not throughput.
    """
    rng = np.random.default_rng(seed)
    # design: intercept + binary treatment + continuous covariate
    X = np.column_stack([
        np.ones(n_samples),
        rng.integers(0, 2, size=n_samples).astype(np.float64),
        rng.standard_normal(n_samples),
    ])
    beta_true = rng.standard_normal((n_sites, p)) * 0.5
    eta = beta_true @ X.T
    mu = 1.0 / (1.0 + np.exp(-eta))
    cov = rng.integers(8, 25, size=(n_sites, n_samples)).astype(np.int32)
    meth = rng.binomial(cov, mu).astype(np.int32)
    return meth, cov, X


def test_irls_gpu_matches_cpu_basic():
    """Coefficients, SE, deviance, Pearson agree to 1e-6 on a clean fixture."""
    meth, cov, X = _make_glm_inputs()
    beta_cpu, se_cpu, dev_cpu, pearson_cpu, n_eff_cpu = _glm.irls_binomial_batch(
        meth, cov, X,
    )
    beta_gpu, se_gpu, dev_gpu, pearson_gpu, n_eff_gpu = irls_binomial_batch_gpu(
        meth, cov, X,
    )

    np.testing.assert_array_equal(n_eff_cpu, n_eff_gpu)
    np.testing.assert_allclose(beta_cpu, beta_gpu, rtol=1e-6, atol=1e-6,
                               equal_nan=True, err_msg="beta diverges CPU vs GPU")
    np.testing.assert_allclose(se_cpu, se_gpu, rtol=1e-6, atol=1e-6,
                               equal_nan=True, err_msg="se_beta diverges")
    np.testing.assert_allclose(dev_cpu, dev_gpu, rtol=1e-6, atol=1e-6,
                               equal_nan=True, err_msg="deviance diverges")
    np.testing.assert_allclose(pearson_cpu, pearson_gpu, rtol=1e-6, atol=1e-6,
                               equal_nan=True, err_msg="pearson diverges")


def test_irls_gpu_return_cov_matches_cpu():
    """The optional cov_beta return matches CPU within tolerance."""
    meth, cov, X = _make_glm_inputs()
    out_cpu = _glm.irls_binomial_batch(meth, cov, X, return_cov=True)
    out_gpu = irls_binomial_batch_gpu(meth, cov, X, return_cov=True)
    assert len(out_cpu) == 6 and len(out_gpu) == 6
    np.testing.assert_allclose(
        out_cpu[5], out_gpu[5], rtol=1e-6, atol=1e-6, equal_nan=True,
        err_msg="cov_beta diverges CPU vs GPU",
    )


def test_irls_dispatch_routes_to_gpu():
    """irls_dispatch(backend='gpu') hits the GPU code path."""
    meth, cov, X = _make_glm_inputs(n_sites=64)
    out_gpu_direct = irls_binomial_batch_gpu(meth, cov, X)
    out_dispatch = _glm.irls_dispatch(meth, cov, X, backend="gpu")
    for a, b in zip(out_gpu_direct, out_dispatch):
        np.testing.assert_array_equal(a, b)


def test_irls_dispatch_unknown_backend_raises():
    meth, cov, X = _make_glm_inputs(n_sites=4)
    with pytest.raises(ValueError, match="Unknown glm_backend"):
        _glm.irls_dispatch(meth, cov, X, backend="quantum")
