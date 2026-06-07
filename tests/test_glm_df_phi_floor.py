"""M4: the GLM p-value path must apply the same DF_PHI_FLOOR (=50) that the
`lr` path applies to the F-reference denominator df. Without the floor,
F(1, ~4) is ~250x more conservative than chi2(1) at typical statistics, so
`test="glm"` was systematically less powerful than `test="lr"` on identical
data despite both being quasi-binomial F/chi2 references."""
import numpy as np
from scipy import stats as sp_stats

from epykit import _glm
from epykit.dmc import DF_PHI_FLOOR


def test_reference_pvalues_df_floor_matches_floored_df():
    stat = np.array([3.84, 6.0, 10.0])
    phi_eff = np.array([2.0, 2.0, 2.0])  # > 1 -> F branch
    small_df = np.array([4.0, 4.0, 4.0])

    # Unfloored F(1, 4) is the over-conservative reference.
    p_unfloored = _glm.reference_pvalues(stat, phi_eff, small_df, reference="F")
    # With df_floor=50, the small df is lifted to 50.
    p_floored = _glm.reference_pvalues(
        stat, phi_eff, small_df, reference="F", df_floor=DF_PHI_FLOOR
    )
    p_at_50 = _glm.reference_pvalues(
        stat, phi_eff, np.full_like(small_df, 50.0), reference="F"
    )

    # The floored result equals evaluating directly at df=50 ...
    np.testing.assert_allclose(p_floored, p_at_50, rtol=1e-12)
    # ... and is strictly less conservative (smaller p) than F(1, 4).
    assert np.all(p_floored < p_unfloored)
    # At stat=3.84 the floored F(1,50) p is close to chi2(1)=0.05; F(1,4) is ~0.12.
    assert p_floored[0] < 0.07
    assert p_unfloored[0] > 0.10


def test_reference_pvalues_df_floor_no_op_when_df_already_large():
    stat = np.array([3.84])
    phi_eff = np.array([2.0])
    big_df = np.array([200.0])
    p_no_floor = _glm.reference_pvalues(stat, phi_eff, big_df, reference="F")
    p_floor = _glm.reference_pvalues(
        stat, phi_eff, big_df, reference="F", df_floor=DF_PHI_FLOOR
    )
    np.testing.assert_allclose(p_no_floor, p_floor, rtol=1e-12)


def test_reference_pvalues_default_df_floor_is_zero_backward_compatible():
    stat = np.array([3.84])
    phi_eff = np.array([2.0])
    small_df = np.array([4.0])
    # Default (no df_floor arg) must reproduce the raw-df behaviour.
    p_default = _glm.reference_pvalues(stat, phi_eff, small_df, reference="F")
    p_raw = sp_stats.f.sf(stat, dfn=1, dfd=small_df)
    np.testing.assert_allclose(p_default, p_raw, rtol=1e-12)


def test_wald_test_df_floor_lifts_small_df():
    # Single-coef contrast (k=1), one site, tiny residual df.
    beta = np.array([[0.0, 1.5]])          # (1 site, 2 coefs)
    cov_beta = np.array([[[1.0, 0.0],
                          [0.0, 0.25]]])    # SE of coef 1 = 0.5 -> Wald = 9.0
    C = np.array([[0.0, 1.0]])             # test coef index 1
    phi_eff = np.array([2.0])              # > 1 -> F branch
    df = np.array([4.0])

    _, p_unfloored, _ = _glm.wald_test(
        beta, cov_beta, C, phi_eff=phi_eff, df_resid=df, reference="F"
    )
    _, p_floored, _ = _glm.wald_test(
        beta, cov_beta, C, phi_eff=phi_eff, df_resid=df, reference="F",
        df_floor=DF_PHI_FLOOR,
    )
    assert p_floored[0] < p_unfloored[0]
