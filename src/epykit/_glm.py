"""Vectorised binomial GLM for covariate-aware methylation testing.

This module provides the engine for ``test="glm"`` in the DMC/DMR pipeline.
It fits one binomial logistic GLM per site (or per tile) on a *shared*
design matrix derived from ``md.obs``, then exposes the per-site deviance
and Pearson chi-squared so the caller can run a deviance LR test of the
full design against a reduced design with the treatment column removed.

The implementation is a single-pass batched IRLS:

    eta_ij  = X_j * beta_i              # linear predictor at site i, sample j
    mu_ij   = sigmoid(eta_ij)
    w_ij    = n_ij * mu_ij * (1 - mu_ij) # binomial GLM weight (cov = trials)
    z_ij    = eta_ij + (y_ij/n_ij - mu_ij) / (mu_ij * (1 - mu_ij))
    beta_i  = (X' W_i X)^{-1} X' W_i z_i

At each iteration we build a (n_sites, p, p) batch of normal equations via
``np.einsum`` and solve them in a single batched ``np.linalg.solve`` call.
With p <= 6 and tile-level n_sites ~ 10^4 to 10^5, the whole IRLS converges
in well under a second.

The hot path deliberately avoids ``statsmodels`` so we don't pay per-site
Python overhead. ``patsy`` is used only once at design-matrix construction
time (a tiny up-front cost).
"""

from __future__ import annotations

import logging
from typing import Optional, Sequence

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)

_EPS = 1e-9
_PROP_CLIP = 1e-6


# Design matrix construction

def build_design(
    obs: pl.DataFrame,
    samples_ordered: Sequence[str],
    formula: Optional[str] = None,
    covariates: Optional[Sequence[str]] = None,
    treatment_col: str = "treatment",
    require_treatment_col: bool = True,
    return_design_info: bool = False,
) -> tuple[np.ndarray, np.ndarray, int, list[str], str]:
    """Build full + reduced model matrices from ``md.obs``.

    Parameters
    ----------
    obs
        ``md.obs`` polars DataFrame. Must contain ``sample_id`` and every
        column referenced by ``formula`` / ``covariates``.
    samples_ordered
        Sample ids in the exact order the DMC/DMR engine will load them
        (case ids first, then control ids). Rows of the returned design
        matrices follow this order.
    formula
        patsy formula string, e.g. ``"~ treatment + age + donor"``. The
        leading tilde is optional. If ``None`` it is synthesised from
        ``treatment_col`` and ``covariates``.
    covariates
        Convenience list of column names. Combined with ``formula`` by
        appending any name not already present.
    treatment_col
        Name of the binary 0/1 treatment column in ``obs``. The reduced
        design drops exactly this column from the full design; that single
        coefficient is the contrast tested by the deviance LR test.
    require_treatment_col
        When True (default, backwards-compatible) the function enforces
        that ``treatment_col`` is present and produces a reduced design
        with that column dropped. When False (used by the multi-group /
        contrast path in :func:`tl.dmc`), the reduced design is set to
        ``None`` and ``coef_idx`` becomes ``-1`` -- the caller is expected
        to use :func:`wald_test` against a contrast matrix instead of a
        full-vs-reduced deviance test.
    return_design_info
        When True, the function returns a 6-tuple with the patsy
        ``DesignInfo`` object appended. Allows callers to resolve
        contrast strings via :func:`resolve_contrast`.

    Returns
    -------
    X_full        (n_samples, p)        full design matrix
    X_reduced     (n_samples, p - 1)    X_full with the treatment column removed
    coef_idx      column index of treatment in X_full
    term_names    column names from patsy (length p)
    formula_used  the formula actually fitted, after merging covariates
    """
    if require_treatment_col and treatment_col not in obs.columns:
        raise ValueError(
            f"Treatment column '{treatment_col}' not found in md.obs. "
            f"Available: {obs.columns}"
        )

    # ---- Synthesise / normalise the formula --------------------------------
    terms: list[str] = []
    if formula is not None:
        rhs = formula.strip()
        if rhs.startswith("~"):
            rhs = rhs[1:].strip()
        # split on '+' but keep interactions ('a:b') and transforms intact
        for term in rhs.split("+"):
            t = term.strip()
            if t and t not in terms:
                terms.append(t)

    if covariates:
        for c in covariates:
            if c not in terms:
                terms.append(c)

    if require_treatment_col and treatment_col not in terms:
        terms.insert(0, treatment_col)

    if not terms:
        raise ValueError(
            "Empty formula. Pass a non-empty `formula=` or `covariates=` "
            "(or use the default binary path with a treatment column on md.obs)."
        )

    formula_used = "~ " + " + ".join(terms)

    # ---- Reorder obs rows to samples_ordered -------------------------------
    obs_pd = obs.to_pandas().set_index("sample_id")
    missing = [s for s in samples_ordered if s not in obs_pd.index]
    if missing:
        raise ValueError(
            f"Samples missing from md.obs: {missing}"
        )
    obs_pd = obs_pd.loc[list(samples_ordered)]

    # ---- Validate covariate columns are present and complete ----------------
    referenced_cols: set[str] = set()
    for term in terms:
        for tok in term.replace(":", " ").replace("*", " ").split():
            tok = tok.strip("()")
            if tok in obs_pd.columns:
                referenced_cols.add(tok)
    for col in referenced_cols:
        nulls = obs_pd[col].isna()
        if nulls.any():
            bad = obs_pd.index[nulls].tolist()
            raise ValueError(
                f"Column '{col}' has missing values for samples {bad}. "
                "Drop or impute these samples before calling the GLM."
            )

    # ---- Build design matrix via patsy --------------------------------------
    try:
        import patsy
    except ImportError as e:  # pragma: no cover - patsy ships with statsmodels
        raise ImportError(
            "patsy is required for covariate-aware DMR. It ships with "
            "statsmodels, which is already in epykit's dependencies."
        ) from e

    X_design = patsy.dmatrix(formula_used, data=obs_pd, return_type="matrix")
    X_full = np.asarray(X_design, dtype=np.float64)
    design_info = X_design.design_info
    term_names: list[str] = list(design_info.column_names)

    n_samples, p_full = X_full.shape
    if p_full >= n_samples:
        raise ValueError(
            f"Too many covariates: design has p={p_full} parameters but "
            f"only n_samples={n_samples}. Add more samples or drop "
            "covariates before fitting."
        )

    if require_treatment_col:
        if treatment_col not in term_names:
            raise ValueError(
                f"Treatment column '{treatment_col}' did not appear as a column "
                f"in the resulting design matrix (got {term_names}). It must be "
                "numeric (0/1) so patsy keeps its name verbatim."
            )
        coef_idx = term_names.index(treatment_col)
        if p_full < 2:
            raise ValueError(
                "Design must contain at least the intercept and the treatment "
                "column (p >= 2). Did you pass '~0 + ...'?"
            )
        # ---- Reduced design: drop the treatment column ----------------------
        X_reduced = np.delete(X_full, coef_idx, axis=1)
    else:
        coef_idx = -1
        X_reduced = None  # type: ignore[assignment]

    if return_design_info:
        return X_full, X_reduced, coef_idx, term_names, formula_used, design_info
    return X_full, X_reduced, coef_idx, term_names, formula_used


# Batched IRLS for the binomial GLM

def irls_dispatch(
    meth: np.ndarray,
    cov: np.ndarray,
    X: np.ndarray,
    *,
    backend: str = "cpu",
    **kwargs,
):
    """Dispatch IRLS to the CPU (default) or GPU backend.

    ``backend="cpu"`` calls :func:`irls_binomial_batch` directly -- the
    historical, default path. ``backend="gpu"`` imports
    :mod:`epykit._glm_gpu` lazily and routes through CuPy. Both return
    the same shapes and dtypes; downstream code is GPU-agnostic.

    Anything other than ``"cpu"`` / ``"gpu"`` raises ``ValueError``.
    """
    backend = (backend or "cpu").lower()
    if backend == "cpu":
        return irls_binomial_batch(meth, cov, X, **kwargs)
    if backend == "gpu":
        from ._glm_gpu import irls_binomial_batch_gpu
        return irls_binomial_batch_gpu(meth, cov, X, **kwargs)
    raise ValueError(
        f"Unknown glm_backend {backend!r}. Use 'cpu' (default) or 'gpu'."
    )


def irls_binomial_batch(
    meth: np.ndarray,
    cov: np.ndarray,
    X: np.ndarray,
    max_iter: int = 25,
    tol: float = 1e-6,
    return_cov: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray] | tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray
]:
    """Batched IRLS for binomial GLMs sharing a design matrix.

    Fits one independent GLM per row (= per site / tile) where the response
    is ``meth_ij`` successes out of ``cov_ij`` trials at sample j, modelled
    as ``logit(mu_ij) = X_j * beta_i``.

    Parameters
    ----------
    meth      (n_sites, n_samples) int32   methylated read counts
    cov       (n_sites, n_samples) int32   total read counts (= weights)
    X         (n_samples, p)       float   design matrix (shared across sites)
    max_iter, tol
        IRLS convergence controls. ``tol`` applies to the max-norm of
        ``beta - beta_prev``.
    return_cov
        When True, append the per-site covariance matrix
        ``cov_beta`` (shape ``(n_sites, p, p)``) to the return tuple as a
        sixth element. Used by :func:`wald_test` for joint / multi-row
        contrasts. Existing five-element callers stay unchanged.

    Returns
    -------
    beta            (n_sites, p)        fitted coefficients
    se_beta         (n_sites, p)        sqrt(diag((X' W X)^{-1})) at convergence
    deviance        (n_sites,)          -2 logL at the fitted mu (binomial)
    pearson_chi2    (n_sites,)          Sigma_j (y - n mu)^2 / (n mu (1-mu))
    n_eff           (n_sites,)          number of samples with cov > 0
    cov_beta        (n_sites, p, p)     OPTIONAL, only when return_cov=True

    Sites where the IRLS hit logistic-regression separation (any covered
    sample's eta reached the +/-30 clip bound) are NaN'd in ``deviance``,
    ``pearson_chi2``, ``beta`` and ``se_beta``. This is required for a
    sane Pearson dispersion: without it the Pearson denominator
    ``n * mu * (1 - mu)`` collapses to ``n * _PROP_CLIP`` at saturated
    samples and per-site chi-sq blows up by ~6 OOM.
    """
    n_sites, n_samples = meth.shape
    p = X.shape[1]
    assert X.shape[0] == n_samples, "Design rows must match number of samples"

    cov_f  = cov.astype(np.float64, copy=False)
    meth_f = meth.astype(np.float64, copy=False)
    has_cov = cov_f > 0
    n_eff = has_cov.sum(axis=1).astype(np.int32)

    # ---- Initialise beta from per-site clipped logit proportions -----------
    # Use a single round of weighted least squares on the logit of the
    # clipped sample proportions. Samples with zero coverage contribute zero
    # weight so they fall out of the init naturally.
    with np.errstate(invalid="ignore", divide="ignore"):
        prop_raw = np.where(has_cov, meth_f / np.maximum(cov_f, 1.0), 0.5)
    prop_init = np.clip(prop_raw, _PROP_CLIP, 1.0 - _PROP_CLIP)
    z_init = np.log(prop_init / (1.0 - prop_init))     # (n_sites, n_samples)
    # weights for init: coverage (mass) only where covered
    w_init = np.where(has_cov, cov_f, 0.0)

    beta = _solve_weighted_lsq(X, z_init, w_init)       # (n_sites, p)
    converged = np.zeros(n_sites, dtype=bool)

    # ---- IRLS iterations ---------------------------------------------------
    for it in range(max_iter):
        eta = beta @ X.T                                # (n_sites, n_samples)
        # clip eta to avoid sigmoid saturation breaking the weights
        eta = np.clip(eta, -30.0, 30.0)
        mu = 1.0 / (1.0 + np.exp(-eta))
        mu = np.clip(mu, _PROP_CLIP, 1.0 - _PROP_CLIP)
        var = mu * (1.0 - mu)

        # Working weights: w = n * mu * (1 - mu). Zero where no coverage.
        w = np.where(has_cov, cov_f * var, 0.0)
        # Working response: z = eta + (y/n - mu) / var.
        with np.errstate(invalid="ignore", divide="ignore"):
            resid_over_var = np.where(
                has_cov, (meth_f / np.maximum(cov_f, 1.0) - mu) / var, 0.0
            )
        z = eta + resid_over_var

        beta_new = _solve_weighted_lsq(X, z, w)

        delta = np.max(np.abs(beta_new - beta), axis=1)
        beta = beta_new
        newly_converged = delta < tol
        converged |= newly_converged
        if converged.all():
            break

    # ---- Final diagnostics (deviance, Pearson chi-sq, SE) ------------------
    eta_unclipped = beta @ X.T

    # Detect logistic-regression separation: a covered sample whose linear
    # predictor reached the eta clip bound means the IRLS could not fit a
    # finite mu there (one stratum is fully methylated / unmethylated). At
    # those samples the Pearson denominator collapses to n * _PROP_CLIP and
    # (y - n*mu)^2 / (n * _PROP_CLIP) blows up by ~6 orders of magnitude,
    # poisoning the chrom-pooled dispersion estimator. Mark separated sites
    # as degenerate so they (a) drop out of compute_dispersion_phi's
    # `usable` mask and (b) carry NaN p-values downstream instead of
    # spurious "significant" calls. We detect separation directly from
    # the saturated eta rather than relying on a slower convergence-
    # warning round-trip.
    SATURATION_THRESHOLD = 29.999  # tiny FP margin below the +/-30 clip
    separated_per_sample = (np.abs(eta_unclipped) >= SATURATION_THRESHOLD) & has_cov
    site_separated = separated_per_sample.any(axis=1)

    eta = np.clip(eta_unclipped, -30.0, 30.0)
    mu = 1.0 / (1.0 + np.exp(-eta))
    mu = np.clip(mu, _PROP_CLIP, 1.0 - _PROP_CLIP)
    var = mu * (1.0 - mu)

    # Binomial deviance: 2 * Sigma_j [y log(y/(n mu)) + (n-y) log((n-y)/(n(1-mu)))]
    # with the convention 0 log 0 = 0. Zero-coverage samples contribute 0.
    y = meth_f
    n = cov_f
    with np.errstate(invalid="ignore", divide="ignore"):
        n_mu = n * mu
        n_minus_y = np.maximum(n - y, 0.0)
        n_one_minus_mu = n * (1.0 - mu)
        term_a = np.where(
            (y > 0) & has_cov,
            y * np.log(np.maximum(y, _EPS) / np.maximum(n_mu, _EPS)),
            0.0,
        )
        term_b = np.where(
            (n_minus_y > 0) & has_cov,
            n_minus_y * np.log(
                np.maximum(n_minus_y, _EPS) / np.maximum(n_one_minus_mu, _EPS)
            ),
            0.0,
        )
    deviance = 2.0 * (term_a + term_b).sum(axis=1)

    with np.errstate(invalid="ignore", divide="ignore"):
        pearson_per = np.where(
            has_cov & (var > 0),
            (y - n_mu) ** 2 / np.maximum(n * var, _EPS),
            0.0,
        )
    pearson_chi2 = pearson_per.sum(axis=1)

    # Standard errors: sqrt(diag((X' W X)^{-1}))
    w_final = np.where(has_cov, cov_f * var, 0.0)
    XtWX = np.einsum("jp,ij,jq->ipq", X, w_final, X)
    se_beta = np.full((n_sites, p), np.nan, dtype=np.float64)
    cov_beta = np.full((n_sites, p, p), np.nan, dtype=np.float64)
    try:
        XtWX_inv = np.linalg.inv(XtWX)
        cov_beta = XtWX_inv
        diag = np.einsum("ipp->ip", XtWX_inv)
        se_beta = np.sqrt(np.where(diag > 0, diag, np.nan))
    except np.linalg.LinAlgError:
        # Batched inversion failed: at least one site has a singular
        # X'WX. Fall back per-site and surface how many sites needed the
        # rescue so users can spot a globally ill-conditioned design
        # (collinear covariates, perfect separation, etc).
        n_fallback_failures = 0
        for i in range(n_sites):
            try:
                inv = np.linalg.inv(XtWX[i])
                cov_beta[i] = inv
                se_beta[i] = np.sqrt(np.where(np.diag(inv) > 0, np.diag(inv), np.nan))
            except np.linalg.LinAlgError:
                n_fallback_failures += 1
                continue
        logger.warning(
            "GLM design matrix singular under batched inversion; fell back "
            "to per-site solve. %d / %d sites still failed and were NaN'd. "
            "Check for collinear covariates or perfect separation.",
            n_fallback_failures, n_sites,
        )

    # Sites with no usable data are degenerate. Sites where the GLM
    # separated (any covered sample's eta hit the clip bound) are also
    # marked degenerate: the fitted mu there is at the boundary, the
    # Pearson denominator collapses to n*_PROP_CLIP, and per-site chi-sq
    # blows up by ~6 OOM (driving chrom-pooled phi from O(1) to O(10^6)
    # and producing nonsensical p-values).
    degenerate = (n_eff < 2) | site_separated
    deviance = np.where(degenerate, np.nan, deviance)
    pearson_chi2 = np.where(degenerate, np.nan, pearson_chi2)
    beta = np.where(degenerate[:, None], np.nan, beta)
    se_beta = np.where(degenerate[:, None], np.nan, se_beta)
    cov_beta = np.where(degenerate[:, None, None], np.nan, cov_beta)

    n_separated = int(site_separated.sum())
    if n_separated > 0:
        # >5 % of sites separated is loud enough to warn the user about
        # a likely model-specification issue; below that, leave it as an
        # info-level breadcrumb that doesn't spam normal runs.
        sep_frac = n_separated / max(n_sites, 1)
        log_fn = logger.warning if sep_frac >= 0.05 else logger.info
        log_fn(
            "GLM separation detected at %d / %d sites (%.1f%%); "
            "NaN'd for dispersion + p-value computation.",
            n_separated, n_sites, 100.0 * sep_frac,
        )

    if return_cov:
        return beta, se_beta, deviance, pearson_chi2, n_eff, cov_beta
    return beta, se_beta, deviance, pearson_chi2, n_eff


def _solve_weighted_lsq(
    X: np.ndarray,
    z: np.ndarray,
    w: np.ndarray,
) -> np.ndarray:
    """Batched solve of (X' W_i X) beta_i = X' W_i z_i across sites i.

    X is (n_samples, p), z and w are (n_sites, n_samples).
    Returns (n_sites, p). Sites with a singular X' W X get NaN.
    """
    # X' W_i X : (n_sites, p, p)
    XtWX = np.einsum("jp,ij,jq->ipq", X, w, X)
    # X' W_i z_i : (n_sites, p)
    XtWz = np.einsum("jp,ij,ij->ip", X, w, z)

    # Ridge guard for sites that have only a couple of covered samples or
    # collinearity in their effective design. We add a tiny multiple of I
    # *only* where the matrix is otherwise singular, so that well-posed
    # sites stay numerically identical to a plain solve.
    p = X.shape[1]
    eye = np.eye(p, dtype=np.float64)

    # NumPy 2.x: np.linalg.solve no longer broadcasts a 2-D RHS against a
    # batched 3-D LHS, so we promote b to (n_sites, p, 1) and squeeze back.
    try:
        beta = np.linalg.solve(XtWX, XtWz[..., None])[..., 0]
        return beta
    except np.linalg.LinAlgError:
        pass

    # Per-site fallback with a ridge.
    n_sites = z.shape[0]
    beta = np.full((n_sites, p), np.nan, dtype=np.float64)
    for i in range(n_sites):
        A = XtWX[i]
        try:
            beta[i] = np.linalg.solve(A, XtWz[i])
        except np.linalg.LinAlgError:
            try:
                beta[i] = np.linalg.solve(A + 1e-8 * eye, XtWz[i])
            except np.linalg.LinAlgError:
                pass
    return beta


# Dispersion + reference-distribution helpers (shared with _score_finalize)

def compute_dispersion_phi(
    pearson_per_site: np.ndarray,
    df_per_site: np.ndarray,
    dispersion: str = "site",
    min_dispersion: float = 1.0,
    shrink_pseudo_df: float = 4.0,
    min_disp_sites: int = 100,
    chrom_name: str = "?",
) -> tuple[np.ndarray, float]:
    """McCullagh-Nelder dispersion in the three modes used elsewhere.

    Parameters
    ----------
    pearson_per_site
        Per-site Pearson chi-sq from the full-model fit (n_sites,).
    df_per_site
        Per-site residual degrees of freedom (n_obs_i - p_full), n_sites.
    dispersion
        ``"site"`` (default): per-site phi_i = pearson_i / df_i.
        ``"chrom"``: single chromosome-pooled phi.
        ``"shrink"``: James-Stein-style weighted average of per-site and
        chromosome estimates.
    min_dispersion
        Clamp on phi (default 1.0). Underdispersion (phi < 1) usually
        reflects model misspecification rather than truly less-than-
        binomial variability; clamping at 1 is the conservative choice.

    Returns
    -------
    phi_eff      (n_sites,) per-site dispersion used by the test
    phi_hat      scalar chromosome-pooled phi (logged / returned for audit)
    """
    if dispersion not in {"site", "chrom", "shrink", "eb"}:
        raise ValueError(
            f"dispersion must be 'site', 'chrom', 'shrink', or 'eb'; got {dispersion!r}"
        )

    usable = (df_per_site > 0) & np.isfinite(pearson_per_site)
    n_usable = int(usable.sum())

    if n_usable < min_disp_sites:
        phi_hat = float(min_dispersion)
        phi_raw = float(min_dispersion)
        if dispersion == "chrom":
            logger.warning(
                "%s: only %d sites usable for dispersion estimation; "
                "falling back to phi = %.2f.",
                chrom_name, n_usable, min_dispersion,
            )
    else:
        pearson_sum = float(pearson_per_site[usable].sum())
        df_sum = float(df_per_site[usable].sum())
        df_sum = max(df_sum, 1.0)
        phi_raw = pearson_sum / df_sum
        phi_hat = float(max(min_dispersion, phi_raw))
        logger.info(
            "%s: chrom-pooled phi = %.3f (raw %.3f, %d sites); dispersion='%s'",
            chrom_name, phi_hat, phi_raw, n_usable, dispersion,
        )

    if dispersion == "chrom":
        phi_eff = np.full(pearson_per_site.shape, phi_hat, dtype=np.float64)
        return phi_eff, phi_hat

    # site / shrink: per-site estimate
    df_safe = np.where(df_per_site > 0, df_per_site, 1.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        phi_site = pearson_per_site / df_safe
    phi_site = np.where(usable, phi_site, min_dispersion)
    phi_site = np.maximum(phi_site, min_dispersion)

    if dispersion == "site":
        return phi_site, phi_hat

    if dispersion == "eb":
        n_usable_eb = int(usable.sum())
        if n_usable_eb >= min_disp_sites:
            phi_obs = phi_site[usable]
            m = float(np.mean(phi_obs))
            v = float(np.var(phi_obs))
            if v > 1e-9 and m > 0:
                a_mom = m * m / v + 2.0
                w_eb = max(1.0, 2.0 * a_mom)
            else:
                w_eb = float(shrink_pseudo_df)
        else:
            w_eb = float(shrink_pseudo_df)
        num = df_safe * phi_site + w_eb * phi_hat
        den = df_safe + w_eb
        phi_eff = np.maximum(num / den, min_dispersion)
        phi_eff = np.where(usable, phi_eff, phi_hat)
        logger.info(
            "%s: bb_lr empirical-Bayes shrinkage w_eb=%.2f",
            chrom_name, w_eb,
        )
        return phi_eff, phi_hat

    # shrink
    w = float(shrink_pseudo_df)
    num = df_safe * phi_site + w * phi_hat
    den = df_safe + w
    phi_eff = np.maximum(num / den, min_dispersion)
    phi_eff = np.where(usable, phi_eff, phi_hat)
    return phi_eff, phi_hat


def reference_pvalues(
    stat: np.ndarray,
    phi_eff: np.ndarray,
    df_resid: np.ndarray,
    reference: str = "adaptive",
) -> np.ndarray:
    """Convert a (already dispersion-corrected) chi-sq statistic to p-values.

    ``reference="adaptive"`` (default) switches per-site between
    ``F(1, df_resid)`` where ``phi_eff > 1`` (real overdispersion signal)
    and ``chi2(1)`` where ``phi_eff`` was clamped to 1 -- the right
    behaviour for quasi-binomial GLMs whose dispersion estimate is noisy
    at small samples. ``"F"`` and ``"chi2"`` force a single reference
    distribution regardless of per-site dispersion.
    """
    if reference == "methylkit":
        raise ValueError(
            "reference='methylkit' was renamed to 'adaptive' in this release."
        )
    if reference not in {"adaptive", "F", "chi2"}:
        raise ValueError(
            f"reference must be 'adaptive', 'F', or 'chi2'; got {reference!r}"
        )
    from scipy import stats as sp_stats

    if reference == "adaptive":
        p_F = sp_stats.f.sf(stat, dfn=1, dfd=df_resid)
        p_chi2 = sp_stats.chi2.sf(stat, df=1)
        return np.where(phi_eff > 1.0, p_F, p_chi2)
    if reference == "F":
        return sp_stats.f.sf(stat, dfn=1, dfd=df_resid)
    return sp_stats.chi2.sf(stat, df=1)


# Contrast / Wald-test helpers -- used by tl.dmc(formula=..., contrast=...)

def resolve_contrast(
    contrast,
    term_names: Sequence[str],
    design_info=None,
) -> tuple[np.ndarray, str]:
    """Build a contrast matrix C of shape (k, p) from a flexible spec.

    Accepted forms:

    * ``contrast=str`` naming a single column in ``term_names`` -- produces a
      1xp contrast vector selecting that coefficient (e.g. ``"age"`` for a
      continuous covariate primary effect).
    * ``contrast=str`` naming a factor (no exact column match) -- every term
      whose name starts with ``"<factor>["`` is included (patsy treatment-
      coded dummies). Returns a kxp contrast for a joint F-test.
    * ``contrast=str`` containing ``"="`` or arithmetic operators (``+``,
      ``-``, ``*``) -- passed to patsy ``DesignInfo.linear_constraint`` for
      named linear contrasts like ``"group[T.KO] - group[T.WT]"``.
      Requires ``design_info`` to be supplied.
    * ``contrast=np.ndarray`` shape (k, p) -- used verbatim.

    Returns ``(C, label)`` where ``label`` describes the contrast for
    provenance.
    """
    p = len(term_names)
    if isinstance(contrast, np.ndarray):
        C = np.atleast_2d(contrast).astype(np.float64)
        if C.shape[1] != p:
            raise ValueError(
                f"Contrast matrix has {C.shape[1]} columns but design has p={p}"
            )
        return C, f"matrix[{C.shape[0]}x{C.shape[1]}]"

    if not isinstance(contrast, str):
        raise TypeError(
            f"contrast must be str or np.ndarray; got {type(contrast).__name__}"
        )

    # Exact column match -- single-row contrast selecting that coefficient
    if contrast in term_names:
        col = term_names.index(contrast)
        C = np.zeros((1, p), dtype=np.float64)
        C[0, col] = 1.0
        return C, contrast

    # Linear-constraint expression -- delegate to patsy
    expr_chars = set("=+-*")
    if any(c in contrast for c in expr_chars) and design_info is not None:
        try:
            constraint = design_info.linear_constraint(contrast)
        except Exception as exc:
            raise ValueError(
                f"Could not parse contrast {contrast!r} with patsy: {exc}"
            ) from exc
        C = np.asarray(constraint.coefs, dtype=np.float64)
        if C.ndim == 1:
            C = C[None, :]
        if C.shape[1] != p:
            raise ValueError(
                f"patsy returned contrast with {C.shape[1]} columns but design has p={p}"
            )
        return C, contrast

    # Factor name -- collect every term beginning with "<factor>["
    factor_terms = [
        (i, t) for i, t in enumerate(term_names)
        if t.startswith(f"{contrast}[") or t == f"C({contrast})"
        or t.startswith(f"C({contrast})[")
    ]
    if not factor_terms:
        raise ValueError(
            f"Could not resolve contrast {contrast!r} against design columns "
            f"{term_names}. Either use an exact column name, a factor name "
            "that appears as treatment-coded dummies, or a patsy "
            "linear-combination expression like 'group[T.A] - group[T.B]'."
        )
    C = np.zeros((len(factor_terms), p), dtype=np.float64)
    for row, (col, _) in enumerate(factor_terms):
        C[row, col] = 1.0
    return C, contrast


def wald_test(
    beta: np.ndarray,
    cov_beta: np.ndarray,
    C: np.ndarray,
    phi_eff: Optional[np.ndarray] = None,
    df_resid: Optional[np.ndarray] = None,
    reference: str = "adaptive",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-site Wald / joint-F test of H0: C*beta_i = 0.

    Parameters
    ----------
    beta            (n_sites, p)             fitted coefficients
    cov_beta        (n_sites, p, p)          parameter covariance matrix
    C               (k, p)                   contrast matrix
    phi_eff         (n_sites,) or None       per-site dispersion to scale
                                             cov_beta by (quasi-binomial).
                                             None -> no scaling (binomial).
    df_resid        (n_sites,) or None       per-site residual df for the F
                                             reference distribution. None ->
                                             chi^2 reference.
    reference       {"adaptive","F","chi2"}  How to convert stat -> p-value.
                                             ``"adaptive"`` switches per-site:
                                             F where phi_eff>1, chi^2 where
                                             phi_eff was clamped. ``"chi2"``
                                             forces chi^2. ``"F"`` forces F.

    Returns
    -------
    stat            (n_sites,)               F-statistic / k (or chi^2/k);
                                             reduces to Wald-z^2 at k=1.
    pvalue          (n_sites,)               two-sided p-value
    k               int                      contrast rank (df1)

    Implementation notes
    --------------------
    For a single-row contrast (k=1), Wald^2 = (Cbeta)^2/Var(Cbeta) is chi^2(1) under
    H0. For multi-row contrasts (k>1), the joint Wald statistic is
    (Cbeta)T[CSigmaCT]^-^1(Cbeta); divided by k it follows F(k, df_resid). When phi_eff
    is supplied, cov_beta is scaled by phi (quasi-binomial). When df_resid
    is None we use chi^2(k)/k as the reference.
    """
    from scipy import stats as sp_stats

    beta = np.asarray(beta, dtype=np.float64)
    cov_beta = np.asarray(cov_beta, dtype=np.float64)
    C = np.atleast_2d(np.asarray(C, dtype=np.float64))
    n_sites, p = beta.shape
    k = C.shape[0]
    if C.shape[1] != p:
        raise ValueError(
            f"Contrast has {C.shape[1]} columns but beta has p={p}"
        )
    if cov_beta.shape != (n_sites, p, p):
        raise ValueError(
            f"cov_beta shape {cov_beta.shape} != (n_sites={n_sites}, p={p}, p={p})"
        )

    # Cbeta -- shape (n_sites, k)
    Cb = beta @ C.T

    # C Sigma CT -- shape (n_sites, k, k). When phi_eff is supplied we scale here.
    if phi_eff is not None:
        scale = np.asarray(phi_eff, dtype=np.float64)
        cov_scaled = cov_beta * scale[:, None, None]
    else:
        cov_scaled = cov_beta
    CSCt = np.einsum("kp,ipq,lq->ikl", C, cov_scaled, C)

    stat = np.full(n_sites, np.nan, dtype=np.float64)
    finite = (
        np.isfinite(Cb).all(axis=1)
        & np.isfinite(CSCt).reshape(n_sites, -1).all(axis=1)
    )
    if k == 1:
        var = CSCt[:, 0, 0]
        good = finite & (var > 0)
        stat[good] = (Cb[good, 0] ** 2) / var[good]
    else:
        try:
            inv_block = np.linalg.inv(CSCt[finite])
            stat[finite] = np.einsum(
                "ij,ijk,ik->i", Cb[finite], inv_block, Cb[finite]
            )
        except np.linalg.LinAlgError:
            for i in np.where(finite)[0]:
                try:
                    inv_i = np.linalg.inv(CSCt[i])
                    stat[i] = Cb[i] @ inv_i @ Cb[i]
                except np.linalg.LinAlgError:
                    pass

    # Reference distribution -> p-value
    if df_resid is None or reference == "chi2":
        pvalue = sp_stats.chi2.sf(stat, df=k)
    elif reference == "F":
        f_stat = stat / k
        pvalue = sp_stats.f.sf(f_stat, dfn=k, dfd=df_resid)
    else:  # adaptive: per-site switch between F and chi2
        f_stat = stat / k
        p_F = sp_stats.f.sf(f_stat, dfn=k, dfd=df_resid)
        p_chi2 = sp_stats.chi2.sf(stat, df=k)
        if phi_eff is None:
            pvalue = p_chi2
        else:
            pvalue = np.where(np.asarray(phi_eff) > 1.0, p_F, p_chi2)

    return stat, pvalue, np.int32(k)


def delta_method_meth_diff_ci(
    coef: np.ndarray,
    coef_se: np.ndarray,
    ref_eta: Optional[np.ndarray] = None,
    alpha: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """Delta-method CI on the meth-scale difference Deltabeta from a logit-scale Wald CI.

    For a single binary-treatment coefficient, Deltaeta = coef and the meth-scale
    difference is Deltabeta ~= expit(eta_ref + coef) - expit(eta_ref). At the reference
    fitted mean eta_ref, the local Jacobian of the inverse link is
    ``dbeta/deta = expit(eta)*(1 - expit(eta))``. The Wald CI on coef is mapped to
    the beta scale by ``Deltabeta +/- z * |J(eta_ref)| * SE(coef)``, clamped to [-1, 1].

    Parameters
    ----------
    coef       (n_sites,) treatment coefficient on logit scale
    coef_se    (n_sites,) Wald SE of that coefficient
    ref_eta    (n_sites,) reference linear predictor (e.g. control-group
               fitted eta). When None we use ``eta_ref = 0``, the maximum-
               Jacobian point ``J(0) = 0.25``; this is the most conservative
               (widest) bound.
    alpha      Significance level (default 0.05 -> 95% CI).

    Returns
    -------
    (ci_lo, ci_hi) -- float64 arrays of shape (n_sites,), clamped to [-1, 1].
    """
    from scipy import stats as sp_stats
    z = float(sp_stats.norm.isf(alpha / 2.0))
    coef = np.asarray(coef, dtype=np.float64)
    coef_se = np.asarray(coef_se, dtype=np.float64)
    if ref_eta is None:
        # J(0) = 0.25
        jac = np.full_like(coef, 0.25)
    else:
        ref_eta = np.asarray(ref_eta, dtype=np.float64)
        # expit, bounded
        with np.errstate(over="ignore", under="ignore"):
            p_ref = 1.0 / (1.0 + np.exp(-np.clip(ref_eta, -30.0, 30.0)))
        jac = p_ref * (1.0 - p_ref)
    # meth-scale Deltabeta at the reference (single binary coefficient case): use
    # expit(eta+coef) - expit(eta) for the central estimate (more accurate than
    # the linearisation), but use jac * coef_se for the half-width.
    if ref_eta is None:
        eta_ref = np.zeros_like(coef)
    else:
        eta_ref = np.clip(np.asarray(ref_eta, dtype=np.float64), -30.0, 30.0)
    with np.errstate(over="ignore", under="ignore"):
        p_treat = 1.0 / (1.0 + np.exp(-np.clip(eta_ref + coef, -30.0, 30.0)))
        p_ref_arr = 1.0 / (1.0 + np.exp(-eta_ref))
    diff = p_treat - p_ref_arr
    half = z * jac * coef_se
    lo = np.clip(diff - half, -1.0, 1.0)
    hi = np.clip(diff + half, -1.0, 1.0)
    return lo, hi


def welch_meth_diff_ci(
    mean_case: np.ndarray,
    var_mean_case: np.ndarray,
    mean_ctrl: np.ndarray,
    var_mean_ctrl: np.ndarray,
    alpha: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """Normal CI on Deltabeta from per-group Welford accumulators.

    ``var_mean_*`` is the variance OF THE MEAN (= s^2 / n), so the SE of the
    difference is ``sqrt(var_mean_case + var_mean_ctrl)``. CI clamped to
    [-1, 1].
    """
    from scipy import stats as sp_stats
    z = float(sp_stats.norm.isf(alpha / 2.0))
    se = np.sqrt(np.maximum(var_mean_case + var_mean_ctrl, 0.0))
    diff = mean_case - mean_ctrl
    lo = np.clip(diff - z * se, -1.0, 1.0)
    hi = np.clip(diff + z * se, -1.0, 1.0)
    return lo, hi


def newcombe_diff_ci(
    meth_a: np.ndarray,
    cov_a: np.ndarray,
    meth_b: np.ndarray,
    cov_b: np.ndarray,
    alpha: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """Newcombe (1998) hybrid Wilson-score CI for pi_a - pi_b on POOLED counts.

    Used by the binomial-pool tests (lr, score, fisher, cmh) where no per-
    replicate variance is accumulated. Uses Wilson-score CIs on each
    pooled proportion, then combines them per Newcombe method 10.
    """
    from scipy import stats as sp_stats
    z = float(sp_stats.norm.isf(alpha / 2.0))

    def _wilson(m, n):
        m = m.astype(np.float64); n = n.astype(np.float64)
        with np.errstate(invalid="ignore", divide="ignore"):
            p_hat = np.where(n > 0, m / n, np.nan)
            denom = 1.0 + (z * z) / np.maximum(n, 1e-12)
            centre = (p_hat + (z * z) / (2.0 * np.maximum(n, 1e-12))) / denom
            half = (
                z * np.sqrt(
                    np.maximum(p_hat * (1.0 - p_hat) / np.maximum(n, 1e-12)
                               + (z * z) / (4.0 * np.maximum(n, 1e-12) ** 2),
                               0.0)
                )
                / denom
            )
            lo = np.clip(centre - half, 0.0, 1.0)
            hi = np.clip(centre + half, 0.0, 1.0)
        return p_hat, lo, hi

    p_a, l_a, u_a = _wilson(np.asarray(meth_a), np.asarray(cov_a))
    p_b, l_b, u_b = _wilson(np.asarray(meth_b), np.asarray(cov_b))
    diff = p_a - p_b
    lo = diff - np.sqrt(np.maximum((p_a - l_a) ** 2 + (u_b - p_b) ** 2, 0.0))
    hi = diff + np.sqrt(np.maximum((u_a - p_a) ** 2 + (p_b - l_b) ** 2, 0.0))
    lo = np.clip(lo, -1.0, 1.0)
    hi = np.clip(hi, -1.0, 1.0)
    return lo, hi
