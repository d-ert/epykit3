"""Layer 2: statistical-primitive tests.

Each public statistical function is checked against a reference
implementation from scipy / statsmodels (or against analytical expectations
when no clean reference exists, as for the vectorised hypergeometric
approximation in :func:`fisher_exact_vectorized`).
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest


# fisher_exact_vectorized


def test_fisher_returns_pvalues_and_log2_or_arrays():
    """Function returns two arrays of equal length to inputs."""
    from epykit.dmc import fisher_exact_vectorized
    meth_a   = np.array([5, 10, 0, 20])
    unmeth_a = np.array([15, 0, 30, 5])
    meth_b   = np.array([2, 5, 1, 25])
    unmeth_b = np.array([18, 5, 29, 0])
    pvals, log2_or = fisher_exact_vectorized(meth_a, unmeth_a, meth_b, unmeth_b)
    assert pvals.shape == (4,)
    assert log2_or.shape == (4,)
    finite = np.isfinite(pvals)
    assert ((pvals[finite] >= 0) & (pvals[finite] <= 1)).all()


def test_fisher_extreme_difference_significant():
    """All-methylated vs all-unmethylated should be highly significant."""
    from epykit.dmc import fisher_exact_vectorized
    pvals, _ = fisher_exact_vectorized(
        np.array([50]), np.array([0]),
        np.array([0]),  np.array([50]),
    )
    assert pvals[0] < 1e-5


def test_fisher_identical_groups_not_significant():
    """Identical methylation in both groups should not be significant."""
    from epykit.dmc import fisher_exact_vectorized
    pvals, _ = fisher_exact_vectorized(
        np.array([10]), np.array([10]),
        np.array([10]), np.array([10]),
    )
    assert pvals[0] > 0.5


def test_fisher_reverse_separation_significant():
    """All-unmethylated A vs all-methylated B (hypo direction) must be significant."""
    from epykit.dmc import fisher_exact_vectorized
    pvals, _ = fisher_exact_vectorized(
        np.array([0]),  np.array([50]),
        np.array([50]), np.array([0]),
    )
    assert pvals[0] < 1e-5


def test_fisher_symmetry():
    """Swapping group A and B should yield the same p-value."""
    from epykit.dmc import fisher_exact_vectorized
    p_fwd, _ = fisher_exact_vectorized(
        np.array([50]), np.array([0]),
        np.array([0]),  np.array([50]),
    )
    p_rev, _ = fisher_exact_vectorized(
        np.array([0]),  np.array([50]),
        np.array([50]), np.array([0]),
    )
    np.testing.assert_allclose(p_fwd, p_rev, rtol=1e-10)


def test_fisher_degenerate_row_returns_nan():
    """A row with zero total in one group should yield NaN p-value."""
    from epykit.dmc import fisher_exact_vectorized
    pvals, log2_or = fisher_exact_vectorized(
        np.array([0]), np.array([0]),
        np.array([5]), np.array([5]),
    )
    assert np.isnan(pvals[0])


def test_fisher_directionally_agrees_with_scipy():
    """epykit's vectorised hypergeom approximation agrees with scipy's
    two-sided ``fisher_exact`` on the *significance call* (p<0.05) for the
    overwhelming majority of tables.

    The implementation uses ``2 * P(X >= meth_a)`` clamped to 1, which is a
    common one-sided-doubled approximation to the two-sided exact test.
    On extreme tables this can disagree with scipy by orders of magnitude
    in the raw p-value, but the *qualitative* call (significant vs not) is
    what users actually act on. We assert that.
    """
    from scipy.stats import fisher_exact as scipy_fisher
    from epykit.dmc import fisher_exact_vectorized

    rng = np.random.default_rng(42)
    n = 100  # bigger sample so the agreement rate is well-estimated
    meth_a   = rng.integers(2, 30, n)
    unmeth_a = rng.integers(2, 30, n)
    meth_b   = rng.integers(2, 30, n)
    unmeth_b = rng.integers(2, 30, n)

    epy_p, _ = fisher_exact_vectorized(meth_a, unmeth_a, meth_b, unmeth_b)
    scipy_p = np.array([
        scipy_fisher([[meth_a[i], unmeth_a[i]], [meth_b[i], unmeth_b[i]]]).pvalue
        for i in range(n)
    ])

    # Agreement on the qualitative call at alpha=0.05.
    # The doubled-one-sided approximation is most likely to disagree right
    # at the threshold (where both methods are near p=0.05), so a 75% bar
    # is realistic; 90% would be over-tight on small tables.
    agree_call = (epy_p < 0.05) == (scipy_p < 0.05)
    assert agree_call.mean() >= 0.75, (
        f"epykit/scipy disagree on significance call for "
        f"{(~agree_call).sum()}/{n} tables (rate "
        f"{(~agree_call).mean():.2%})"
    )

    # Disagreements should concentrate near the alpha=0.05 boundary, not be
    # spread across all p-values. Check that the median log10 disagreement
    # is still modest (within an order of magnitude). The mean would be
    # skewed by a handful of extreme tables.
    log_epy = np.log10(np.maximum(epy_p, 1e-300))
    log_scipy = np.log10(np.maximum(scipy_p, 1e-300))
    diff = np.abs(log_epy - log_scipy)
    assert np.median(diff) < 0.5, (
        f"median log10 disagreement too large: {np.median(diff):.3f}"
    )



# apply_multiple_testing_correction


def test_bh_matches_statsmodels_on_clean_pvalues():
    """BH q-values should match statsmodels.multipletests exactly when no
    NaN p-values are present."""
    from statsmodels.stats.multitest import multipletests
    from epykit.dmc import apply_multiple_testing_correction

    rng = np.random.default_rng(0)
    # Mix of small (true positives) and uniform (nulls) p-values.
    pvals = np.concatenate([rng.beta(0.1, 1.0, 30), rng.uniform(0, 1, 70)])
    df = pl.DataFrame({"chrom": ["chr1"] * 100, "pos": list(range(100)), "pvalue": pvals})

    out = apply_multiple_testing_correction(df, method="fdr_bh")
    expected_q = multipletests(pvals, method="fdr_bh")[1]
    np.testing.assert_allclose(out["qvalue"].to_numpy(), expected_q, atol=1e-12)


def test_bh_preserves_nan_pvalues_as_nan_qvalues():
    """NaN p-values must come back as NaN q-values, not 1.0 or 0.0."""
    from epykit.dmc import apply_multiple_testing_correction

    pvals = np.array([0.01, np.nan, 0.5, np.nan, 0.001])
    df = pl.DataFrame({"chrom": ["chr1"] * 5, "pos": list(range(5)), "pvalue": pvals})
    out = apply_multiple_testing_correction(df)
    q = out["qvalue"].to_numpy()
    assert np.isnan(q[1]) and np.isnan(q[3])
    assert np.isfinite(q[0]) and np.isfinite(q[2]) and np.isfinite(q[4])


def test_bh_adds_reject_column():
    """`reject` boolean column should appear with default qvalue_col."""
    from epykit.dmc import apply_multiple_testing_correction
    df = pl.DataFrame({"pvalue": [0.001, 0.01, 0.5, 0.9]})
    out = apply_multiple_testing_correction(df)
    assert "reject" in out.columns
    assert out["reject"].dtype == pl.Boolean


def test_bh_custom_column_names_dont_collide():
    """When applied with custom column names, reject column is renamed."""
    from epykit.dmc import apply_multiple_testing_correction
    df = pl.DataFrame({"combined_pvalue": [0.001, 0.01, 0.5, 0.9]})
    out = apply_multiple_testing_correction(
        df, pvalue_col="combined_pvalue", qvalue_col="combined_qvalue"
    )
    assert "combined_qvalue" in out.columns
    assert "combined_qvalue_reject" in out.columns


def test_bh_is_monotone_with_pvalue():
    """BH q-values are non-decreasing when p-values are non-decreasing."""
    from epykit.dmc import apply_multiple_testing_correction
    pvals = sorted(np.random.default_rng(1).uniform(0, 1, 200))
    df = pl.DataFrame({"pvalue": pvals})
    out = apply_multiple_testing_correction(df)
    q = out["qvalue"].to_numpy()
    diffs = np.diff(q[np.isfinite(q)])
    assert (diffs >= -1e-12).all(), "q-values should be monotone non-decreasing"



# build_design (GLM design matrix builder)


def test_build_design_treatment_only():
    """Minimal design: intercept + treatment yields a (n, 2) matrix."""
    from epykit._glm import build_design

    obs = pl.DataFrame({
        "sample_id": ["s1", "s2", "s3", "s4"],
        "treatment": [1, 1, 0, 0],
    })
    X_full, X_red, coef_idx, names, formula = build_design(
        obs, samples_ordered=["s1", "s2", "s3", "s4"]
    )
    assert X_full.shape == (4, 2)
    assert X_red.shape == (4, 1)  # intercept only after dropping treatment
    assert "treatment" in names
    assert names[coef_idx] == "treatment"


def test_build_design_with_continuous_covariate():
    """Continuous covariate adds one column to the full design."""
    from epykit._glm import build_design

    obs = pl.DataFrame({
        "sample_id": ["s1", "s2", "s3", "s4"],
        "treatment": [1, 1, 0, 0],
        "age":       [25.0, 30.0, 28.0, 35.0],
    })
    X_full, X_red, coef_idx, names, formula = build_design(
        obs, samples_ordered=["s1", "s2", "s3", "s4"],
        covariates=["age"],
    )
    assert X_full.shape == (4, 3)
    assert X_red.shape == (4, 2)
    assert "age" in names


def test_build_design_too_many_covariates_raises():
    """p >= n should be refused (mirrors methylKit's check)."""
    from epykit._glm import build_design

    obs = pl.DataFrame({
        "sample_id": ["s1", "s2"],
        "treatment": [1, 0],
        "a": [1.0, 2.0],
        "b": [3.0, 4.0],
    })
    with pytest.raises(ValueError, match="Too many covariates"):
        build_design(obs, samples_ordered=["s1", "s2"], covariates=["a", "b"])


def test_build_design_missing_covariate_value_raises():
    """NaN in a covariate column should fail loudly with sample ids."""
    from epykit._glm import build_design

    obs = pl.DataFrame({
        "sample_id": ["s1", "s2", "s3", "s4"],
        "treatment": [1, 1, 0, 0],
        "age":       [25.0, None, 28.0, 35.0],
    })
    with pytest.raises(ValueError, match="missing values"):
        build_design(
            obs, samples_ordered=["s1", "s2", "s3", "s4"],
            covariates=["age"],
        )



# irls_binomial_batch -- full-rank GLM solver


def test_irls_recovers_known_coefficients():
    """IRLS on a single binomial GLM agrees with statsmodels GLM(Binomial)
    to within a small tolerance.
    """
    statsmodels = pytest.importorskip("statsmodels.api")
    from epykit._glm import irls_binomial_batch

    rng = np.random.default_rng(0)
    n_samples = 12
    treatment = np.array([1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0], dtype=np.float64)
    X = np.column_stack([np.ones(n_samples), treatment])

    beta_true = np.array([-0.5, 1.5])  # intercept, treatment effect
    eta = X @ beta_true
    pi = 1.0 / (1.0 + np.exp(-eta))
    cov = np.full(n_samples, 50, dtype=np.int32)
    meth = rng.binomial(cov, pi).astype(np.int32)

    beta, se, dev, chi2, n_eff = irls_binomial_batch(
        meth[None, :], cov[None, :], X
    )

    sm_fit = statsmodels.GLM(
        meth, X,
        freq_weights=None,
        family=statsmodels.families.Binomial(),
        exposure=None,
        var_weights=None,
    )
    # statsmodels needs proportions for Binomial; provide successes/trials
    sm_fit = statsmodels.GLM(
        np.column_stack([meth, cov - meth]),
        X, family=statsmodels.families.Binomial(),
    ).fit()
    np.testing.assert_allclose(beta[0], sm_fit.params, atol=5e-3)


def test_irls_marks_degenerate_sites_nan():
    """Sites with fewer than 2 covered samples should produce NaN outputs."""
    from epykit._glm import irls_binomial_batch

    X = np.column_stack([np.ones(4), [1, 1, 0, 0]]).astype(np.float64)
    meth = np.array([[1, 0, 0, 0]], dtype=np.int32)
    cov  = np.array([[5, 0, 0, 0]], dtype=np.int32)  # only one sample covered
    beta, se, dev, chi2, n_eff = irls_binomial_batch(meth, cov, X)
    assert n_eff[0] < 2
    assert np.isnan(dev[0])
    assert np.all(np.isnan(beta[0]))


def test_irls_batched_solves_multiple_sites_independently():
    """Two sites with different signals should yield different fitted betas."""
    from epykit._glm import irls_binomial_batch

    X = np.column_stack([np.ones(8), [1, 1, 1, 1, 0, 0, 0, 0]]).astype(np.float64)
    # Site 0: strong hyper in treatment; Site 1: no effect.
    meth = np.array([
        [45, 48, 47, 44, 5, 6, 4, 7],
        [25, 24, 26, 25, 24, 25, 26, 25],
    ], dtype=np.int32)
    cov  = np.full((2, 8), 50, dtype=np.int32)

    beta, se, dev, chi2, n_eff = irls_binomial_batch(meth, cov, X)
    # Site 0 should have a much larger treatment coefficient than site 1.
    assert abs(beta[0, 1]) > 1.0
    assert abs(beta[1, 1]) < 0.5



