"""Numerical dispersion guards must not invent residual information."""
import numpy as np
from scipy import stats
from epykit import _glm
from epykit.dmc import DF_PHI_FLOOR


def test_reference_preserves_estimated_degrees_of_freedom():
    statistic=np.array([3.84,6.,10.]);df=np.array([4.,8.,20.]);phi=np.full(3,2.)
    actual=_glm.reference_pvalues(statistic,phi,df,reference='F',df_floor=DF_PHI_FLOOR)
    np.testing.assert_allclose(actual,stats.f.sf(statistic,1,df),rtol=1e-12)


def test_wald_test_small_sample_guard_is_noop():
    beta=np.array([[0.,1.5]]);cov=np.array([[[1.,0.],[0.,.25]]]);contrast=np.array([[0.,1.]])
    args=dict(phi_eff=np.array([2.]),df_resid=np.array([4.]),reference='F')
    _,first,_=_glm.wald_test(beta,cov,contrast,**args)
    _,second,_=_glm.wald_test(beta,cov,contrast,df_floor=DF_PHI_FLOOR,**args)
    np.testing.assert_allclose(first,second)
