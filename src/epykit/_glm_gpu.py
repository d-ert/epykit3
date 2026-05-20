"""GPU-accelerated batched binomial IRLS via CuPy.

This module mirrors :func:`epykit._glm.irls_binomial_batch` but runs the
matrix maths on a CUDA device through ``cupy``. The two implementations
must stay numerically identical -- :func:`tests.test_glm_gpu` enforces a
1e-6 tolerance on the synth fixture, and the same algorithmic structure
is preserved line-by-line below so reviewers can diff the two files.

Usage
-----
``cupy`` is optional. Install with ``pip install 'epykit[gpu]'``. The
module imports cupy lazily so a missing wheel surfaces only when the
GPU backend is explicitly requested.

The public entry point is :func:`irls_binomial_batch_gpu`. It accepts
the same arguments as the CPU version and returns numpy arrays (not
cupy arrays) so downstream code stays GPU-agnostic.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

# Constants kept in sync with _glm.py
_EPS = 1e-9
_PROP_CLIP = 1e-6


def _require_cupy():
    """Import cupy lazily with a clear install hint on ImportError."""
    try:
        import cupy as cp  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "cupy is required for backend='gpu'. "
            "Install with: pip install 'epykit[gpu]' "
            "(or 'epykit[gpu_jax]' for the JAX-based alternative). "
            "Heavy CUDA wheels; not included in 'epykit[all]'."
        ) from exc
    return cp


def _solve_weighted_lsq_gpu(cp, X_d, z_d, w_d):
    """GPU equivalent of :func:`epykit._glm._solve_weighted_lsq`.

    X_d : (n_samples, p)        cupy.float64
    z_d : (n_sites, n_samples)  cupy.float64
    w_d : (n_sites, n_samples)  cupy.float64
    Returns (n_sites, p) cupy.float64.
    """
    XtWX = cp.einsum("jp,ij,jq->ipq", X_d, w_d, X_d)
    XtWz = cp.einsum("jp,ij,ij->ip", X_d, w_d, z_d)
    p = X_d.shape[1]
    n_sites = z_d.shape[0]

    try:
        beta = cp.linalg.solve(XtWX, XtWz[..., None])[..., 0]
        return beta
    except cp.linalg.LinAlgError:
        pass

    # Per-site fallback with ridge. Loop on GPU is slow but only fires
    # when batched solve already failed (rare in practice).
    eye = cp.eye(p, dtype=cp.float64)
    beta = cp.full((n_sites, p), cp.nan, dtype=cp.float64)
    for i in range(n_sites):
        A = XtWX[i]
        try:
            beta[i] = cp.linalg.solve(A, XtWz[i])
        except cp.linalg.LinAlgError:
            try:
                beta[i] = cp.linalg.solve(A + 1e-8 * eye, XtWz[i])
            except cp.linalg.LinAlgError:
                pass
    return beta


def irls_binomial_batch_gpu(
    meth: np.ndarray,
    cov: np.ndarray,
    X: np.ndarray,
    max_iter: int = 25,
    tol: float = 1e-6,
    return_cov: bool = False,
):
    """GPU port of :func:`epykit._glm.irls_binomial_batch`.

    Same algorithm, same outputs (numpy arrays). All maths is done on
    device; numpy <-> cupy transfers happen at the boundary.

    Parameters
    ----------
    meth, cov : (n_sites, n_samples)
        Methylated / total read counts. Must be numpy arrays -- they get
        promoted to float64 and transferred to GPU.
    X : (n_samples, p)
        Shared design matrix.
    max_iter, tol, return_cov
        Same semantics as the CPU version.

    Returns
    -------
    Tuple of numpy arrays with the same shapes / dtypes as the CPU
    version. Always returns float64 host arrays for downstream code.
    """
    cp = _require_cupy()

    n_sites, n_samples = meth.shape
    p = X.shape[1]
    assert X.shape[0] == n_samples, "Design rows must match number of samples"

    # ---- Transfer inputs to device --------------------------------------
    cov_d = cp.asarray(cov, dtype=cp.float64)
    meth_d = cp.asarray(meth, dtype=cp.float64)
    X_d = cp.asarray(X, dtype=cp.float64)
    has_cov = cov_d > 0
    n_eff = has_cov.sum(axis=1).astype(cp.int32)

    # ---- Initialise beta from clipped logit proportions -----------------
    prop_raw = cp.where(has_cov, meth_d / cp.maximum(cov_d, 1.0), 0.5)
    prop_init = cp.clip(prop_raw, _PROP_CLIP, 1.0 - _PROP_CLIP)
    z_init = cp.log(prop_init / (1.0 - prop_init))
    w_init = cp.where(has_cov, cov_d, 0.0)

    beta = _solve_weighted_lsq_gpu(cp, X_d, z_init, w_init)
    converged = cp.zeros(n_sites, dtype=cp.bool_)

    # ---- IRLS iterations ------------------------------------------------
    for _ in range(max_iter):
        eta = beta @ X_d.T
        eta = cp.clip(eta, -30.0, 30.0)
        mu = 1.0 / (1.0 + cp.exp(-eta))
        mu = cp.clip(mu, _PROP_CLIP, 1.0 - _PROP_CLIP)
        var = mu * (1.0 - mu)

        w = cp.where(has_cov, cov_d * var, 0.0)
        resid_over_var = cp.where(
            has_cov, (meth_d / cp.maximum(cov_d, 1.0) - mu) / var, 0.0
        )
        z = eta + resid_over_var

        beta_new = _solve_weighted_lsq_gpu(cp, X_d, z, w)
        delta = cp.max(cp.abs(beta_new - beta), axis=1)
        beta = beta_new
        converged |= delta < tol
        if bool(converged.all().get()):
            break

    # ---- Final diagnostics ----------------------------------------------
    eta_unclipped = beta @ X_d.T
    SATURATION_THRESHOLD = 29.999
    separated_per_sample = (cp.abs(eta_unclipped) >= SATURATION_THRESHOLD) & has_cov
    site_separated = separated_per_sample.any(axis=1)

    eta = cp.clip(eta_unclipped, -30.0, 30.0)
    mu = 1.0 / (1.0 + cp.exp(-eta))
    mu = cp.clip(mu, _PROP_CLIP, 1.0 - _PROP_CLIP)
    var = mu * (1.0 - mu)

    y = meth_d
    n = cov_d
    n_mu = n * mu
    n_minus_y = cp.maximum(n - y, 0.0)
    n_one_minus_mu = n * (1.0 - mu)
    term_a = cp.where(
        (y > 0) & has_cov,
        y * cp.log(cp.maximum(y, _EPS) / cp.maximum(n_mu, _EPS)),
        0.0,
    )
    term_b = cp.where(
        (n_minus_y > 0) & has_cov,
        n_minus_y * cp.log(
            cp.maximum(n_minus_y, _EPS) / cp.maximum(n_one_minus_mu, _EPS)
        ),
        0.0,
    )
    deviance = 2.0 * (term_a + term_b).sum(axis=1)

    pearson_per = cp.where(
        has_cov & (var > 0),
        (y - n_mu) ** 2 / cp.maximum(n * var, _EPS),
        0.0,
    )
    pearson_chi2 = pearson_per.sum(axis=1)

    # Standard errors
    w_final = cp.where(has_cov, cov_d * var, 0.0)
    XtWX = cp.einsum("jp,ij,jq->ipq", X_d, w_final, X_d)
    se_beta = cp.full((n_sites, p), cp.nan, dtype=cp.float64)
    cov_beta = cp.full((n_sites, p, p), cp.nan, dtype=cp.float64)
    try:
        XtWX_inv = cp.linalg.inv(XtWX)
        cov_beta = XtWX_inv
        diag = cp.einsum("ipp->ip", XtWX_inv)
        se_beta = cp.sqrt(cp.where(diag > 0, diag, cp.nan))
    except cp.linalg.LinAlgError:
        n_fallback_failures = 0
        for i in range(n_sites):
            try:
                inv = cp.linalg.inv(XtWX[i])
                cov_beta[i] = inv
                d_i = cp.diag(inv)
                se_beta[i] = cp.sqrt(cp.where(d_i > 0, d_i, cp.nan))
            except cp.linalg.LinAlgError:
                n_fallback_failures += 1
        logger.warning(
            "GPU IRLS: GLM design matrix singular under batched inversion; "
            "fell back to per-site solve. %d / %d sites still failed and "
            "were NaN'd.", n_fallback_failures, n_sites,
        )

    degenerate = (n_eff < 2) | site_separated
    deviance = cp.where(degenerate, cp.nan, deviance)
    pearson_chi2 = cp.where(degenerate, cp.nan, pearson_chi2)
    beta = cp.where(degenerate[:, None], cp.nan, beta)
    se_beta = cp.where(degenerate[:, None], cp.nan, se_beta)
    cov_beta = cp.where(degenerate[:, None, None], cp.nan, cov_beta)

    n_separated = int(site_separated.sum().get())
    if n_separated > 0:
        sep_frac = n_separated / max(n_sites, 1)
        log_fn = logger.warning if sep_frac >= 0.05 else logger.info
        log_fn(
            "GPU IRLS: separation detected at %d / %d sites (%.1f%%); "
            "NaN'd for dispersion + p-value computation.",
            n_separated, n_sites, 100.0 * sep_frac,
        )

    # ---- Transfer back to host -----------------------------------------
    beta_h = cp.asnumpy(beta)
    se_beta_h = cp.asnumpy(se_beta)
    deviance_h = cp.asnumpy(deviance)
    pearson_chi2_h = cp.asnumpy(pearson_chi2)
    n_eff_h = cp.asnumpy(n_eff)

    if return_cov:
        cov_beta_h = cp.asnumpy(cov_beta)
        return beta_h, se_beta_h, deviance_h, pearson_chi2_h, n_eff_h, cov_beta_h
    return beta_h, se_beta_h, deviance_h, pearson_chi2_h, n_eff_h


__all__ = ["irls_binomial_batch_gpu"]
