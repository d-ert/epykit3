"""Pin the empirical justification of ``DF_PHI_FLOOR`` in dmc.py.

The floor is used in ``_score_finalize`` to keep the F(1, df) tail
close to chi^2(1) when ``reference="adaptive"`` or ``"F"``. The value
50 was introduced in epykit 0.7.2 to fix EB-mode FPR inflation; at the
small df_phi (~ 4) that EB shrinkage can produce, F(1, 4) is ~5x
inflated vs chi^2(1) at the p = 0.05 critical region, causing the
adaptive branch to over-report extreme p-values. Lifting the floor to
50 reduces the F-vs-chi^2 disagreement at p = 0.05 to ~11% relative,
which is below the per-CpG calibration noise floor observed on real
data.

This test makes that justification traceable and breaks any future
change to ``DF_PHI_FLOOR`` that would re-introduce the
pre-0.7.2 calibration problem. The choice is ultimately validated
empirically by the null-calibration figure produced by the Linux
re-run; this test ensures the analytical story matches what the
code does.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from epykit.dmc import (
    DF_PHI_FLOOR,
    _DF_PHI_FLOOR_F_VS_CHI2_TOL_AT_P05,
)


# Critical statistic where chi^2(1) returns p = 0.05.
_CHI2_CRIT = stats.chi2.ppf(0.95, df=1)


def test_floor_value_unchanged():
    """The 50 was hand-calibrated in 0.7.2. Any change requires a
    re-run of the null calibration suite and a methods appendix
    update; assert here so the change cannot land silently.
    """
    assert DF_PHI_FLOOR == 50.0, (
        f"DF_PHI_FLOOR has been changed from its calibrated default of "
        f"50.0 to {DF_PHI_FLOOR}. Run "
        f"benchmark/scripts/run_null_calibration.py and update the "
        f"methods appendix with the new calibration figure before "
        f"committing this change."
    )


def test_floor_drives_f_within_documented_agreement():
    """At ``DF_PHI_FLOOR``, F(1, df) and chi^2(1) agree at p = 0.05
    to within the documented tolerance.
    """
    p_chi2 = stats.chi2.sf(_CHI2_CRIT, df=1)
    p_F = stats.f.sf(_CHI2_CRIT, dfn=1, dfd=DF_PHI_FLOOR)
    relative_excess = (p_F - p_chi2) / p_chi2
    assert 0.0 <= relative_excess <= _DF_PHI_FLOOR_F_VS_CHI2_TOL_AT_P05, (
        f"F(1, {DF_PHI_FLOOR}).sf({_CHI2_CRIT:.3f}) vs chi^2(1).sf(...) "
        f"relative excess = {relative_excess:.4f}; documented bound is "
        f"{_DF_PHI_FLOOR_F_VS_CHI2_TOL_AT_P05}. Either tighten the floor "
        f"or update the bound and the comment in dmc.py."
    )


def test_pathological_low_df_would_be_unacceptable():
    """At df_phi = 4 (the EB-shrunk pathological case) the F-vs-chi^2
    disagreement is multiples of the floor-corrected value.

    This is the failure mode the floor was introduced to prevent.
    """
    p_chi2 = stats.chi2.sf(_CHI2_CRIT, df=1)
    p_F_at_4 = stats.f.sf(_CHI2_CRIT, dfn=1, dfd=4)
    relative_excess = (p_F_at_4 - p_chi2) / p_chi2
    # Pathological excess is order-of-magnitude worse than the floored bound.
    assert relative_excess > 1.0, (
        f"F(1, 4) was expected to be at least 100% inflated vs chi^2(1) "
        f"at p = 0.05 (the bug DF_PHI_FLOOR fixed). Got "
        f"{relative_excess:.4f}. If scipy's F/chi^2 semantics changed, "
        f"the EB-mode story needs to be re-derived."
    )


def test_chi2_limit_is_exact():
    """At infinity, F(1, inf) == chi^2(1). Sanity check that scipy
    behaves as the derivation assumes."""
    p_chi2 = stats.chi2.sf(_CHI2_CRIT, df=1)
    p_F_inf = stats.f.sf(_CHI2_CRIT, dfn=1, dfd=1e8)
    assert p_F_inf == pytest.approx(p_chi2, rel=1e-3)
