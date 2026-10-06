"""Estimated dispersion degrees of freedom must retain small-sample uncertainty."""
import numpy as np
from scipy import stats
from epykit.dmc import DF_PHI_FLOOR


def test_numerical_guard_preserves_valid_small_df():
    df=np.array([2.,4.,8.,100.])
    np.testing.assert_array_equal(np.maximum(df,DF_PHI_FLOOR),df)


def test_small_sample_f_tail_is_not_replaced_by_asymptotic_tail():
    # A heavier finite-df tail is uncertainty, not a calibration defect.
    assert stats.f.sf(25,1,8)>100*stats.f.sf(25,1,50)
