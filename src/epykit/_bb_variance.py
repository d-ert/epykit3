"""Depth-aware variance for pooled beta-binomial group proportions.

The F reference is approximate; empirical calibration is required. This module
uses only observed counts, never simulation truth or external dispersion.
"""

import numpy as np


def moderate_depth_scale(
    sn_case,
    sn2_case,
    nv_case,
    chi_case,
    sn_ctrl,
    sn2_ctrl,
    nv_ctrl,
    chi_ctrl,
    min_dispersion=1.0,
    min_disp_sites=100,
):
    """Moderate depth-relative Pearson scale, then map rho to contrast variance.

    Constant biological rho does not imply constant quasi-binomial phi:
    phi grows with read depth. Learn the signed rho center before clipping,
    normalize residual sums by their depth-specific expectation, and fit the
    residual log-F prior. The resulting F degrees of freedom are approximate;
    count-model null calibration remains necessary at small sample sizes.
    """
    from .dmc import _smyth_eb_moderation

    a, b = np.maximum(sn_case, 1.0), np.maximum(sn_ctrl, 1.0)
    df = np.maximum(nv_case + nv_ctrl - 2.0, 1.0)
    exposure = (sn_case - nv_case - (sn2_case / a - 1.0)) + (
        sn_ctrl - nv_ctrl - (sn2_ctrl / b - 1.0)
    )
    chi = chi_case + chi_ctrl
    valid = (exposure > 0) & (nv_case >= 2) & (nv_ctrl >= 2) & np.isfinite(chi)
    if valid.sum() < min_disp_sites:
        rho = np.clip(np.divide(chi - df, exposure, out=np.zeros_like(df), where=valid), 0, 0.95)
        return (
            np.maximum(
                pooled_design_effect(sn_case, sn2_case, sn_ctrl, sn2_ctrl, rho), min_dispersion
            ),
            df,
            0.0,
        )
    rho_center = float(np.clip(np.sum(chi[valid] - df[valid]) / np.sum(exposure[valid]), 0, 0.95))
    expected_phi = 1.0 + rho_center * np.maximum(exposure, 0.0) / df
    relative_phi, df_mod, _, _ = _smyth_eb_moderation(
        chi / expected_phi,
        df,
        valid,
        min_dispersion=1e-6,
        min_disp_sites=min_disp_sites,
    )
    rho = np.clip(
        np.divide(
            (relative_phi * expected_phi - 1.0) * df,
            exposure,
            out=np.full_like(df, rho_center),
            where=valid,
        ),
        0,
        0.95,
    )
    phi = pooled_design_effect(sn_case, sn2_case, sn_ctrl, sn2_ctrl, rho)
    return np.maximum(phi, min_dispersion), df_mod, rho_center


def pooled_design_effect(sn_case, sn2_case, sn_ctrl, sn2_ctrl, rho):
    a = np.maximum(sn_case, 1.0)
    b = np.maximum(sn_ctrl, 1.0)
    leverage = ((sn2_case / a - 1) / a + (sn2_ctrl / b - 1) / b) / (1 / a + 1 / b)
    return 1 + rho * np.maximum(leverage, 0)


def estimate_scale(
    sn_case,
    sn2_case,
    nv_case,
    chi_case,
    sn_ctrl,
    sn2_ctrl,
    nv_ctrl,
    chi_ctrl,
    pooled_mu,
    moderate=False,
):
    """Estimate rho from the Pearson residual expectation and map to contrast variance.

    Empirical moderation uses the observed between-site variance after an
    approximate Pearson sampling-variance subtraction in pooled-mean strata.
    This is a moment shrinkage estimator, not an exact conjugate posterior.
    """
    a = np.maximum(sn_case, 1.0)
    b = np.maximum(sn_ctrl, 1.0)
    df = np.maximum(nv_case + nv_ctrl - 2, 1.0)
    exposure = (sn_case - nv_case - (sn2_case / a - 1)) + (sn_ctrl - nv_ctrl - (sn2_ctrl / b - 1))
    valid = (exposure > 0) & (nv_case >= 2) & (nv_ctrl >= 2) & (pooled_mu > 0) & (pooled_mu < 1)
    raw = np.divide(chi_case + chi_ctrl - df, exposure, out=np.zeros_like(a), where=valid)
    rho = np.clip(raw, 0, 0.95)
    effective_df = df.copy()
    if moderate and valid.sum() >= 100:
        # Fold methylation around 0.5: hyper/hypomethylation has symmetric
        # sampling behavior. Fit before clipping negative noise estimates.
        labels = np.digitize(np.minimum(pooled_mu, 1 - pooled_mu), [0.02, 0.05, 0.1, 0.2, 0.35])
        global_center = float(np.clip(np.mean(raw[valid]), 0, 0.95))
        for label in np.unique(labels):
            select = (labels == label) & valid
            fit = select if select.sum() >= 100 else valid
            center = float(np.clip(np.mean(raw[fit]), 0, 0.95)) if fit.sum() else global_center
            noise = 2 * (df[fit] + center * exposure[fit]) ** 2 / (df[fit] * exposure[fit] ** 2)
            between = max(float(np.var(raw[fit], ddof=1) - noise.mean()), 0.0)
            local_noise = (
                2 * (df + center * exposure) ** 2 / (df * np.maximum(exposure, 1e-12) ** 2)
            )
            weight = np.divide(
                between,
                between + local_noise,
                out=np.zeros_like(a),
                where=(between + local_noise) > 0,
            )
            rho[select] = np.clip(
                weight[select] * raw[select] + (1 - weight[select]) * center, 0, 0.95
            )
            # Keep actual residual df: the moment shrinkage approximation does
            # not justify adding arbitrary pseudo degrees of freedom.
    phi = pooled_design_effect(sn_case, sn2_case, sn_ctrl, sn2_ctrl, rho)
    return phi, effective_df, rho
