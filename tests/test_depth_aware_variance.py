import numpy as np
from scipy import stats
from epykit._bb_variance import pooled_design_effect
from epykit.dmc import _score_finalize


def test_design_effect_equals_direct_variance_of_weighted_counts():
    a=np.array([5.,10.,15.,30.,40.]);b=np.array([20.,20.,20.,20.,20.]);rho=.2
    scale=pooled_design_effect(np.array([a.sum()]),np.array([(a*a).sum()]),
                              np.array([b.sum()]),np.array([(b*b).sum()]),rho)[0]
    direct=(np.sum(a*(1+(a-1)*rho))/a.sum()**2+np.sum(b*(1+(b-1)*rho))/b.sum()**2)/(1/a.sum()+1/b.sum())
    assert np.isclose(scale,direct)
    equal=pooled_design_effect(np.array([100.]),np.array([2000.]),np.array([100.]),np.array([2000.]),rho)
    assert np.isclose(equal[0],1+19*rho)


def test_adaptive_keeps_f_reference_when_estimated_scale_is_clamped():
    n=np.array([[20]*5]*3,dtype=float)
    m=np.array([[1]*5,[10]*5,[15]*5],dtype=float)
    def accum(m): return n.sum(1),m.sum(1),(m*m/n).sum(1),np.full(3,5)
    inp=(*accum(m),*accum(m+1))
    adaptive=_score_finalize(*inp,dispersion='site',reference='adaptive')[0]
    explicit=_score_finalize(*inp,dispersion='site',reference='F')[0]
    chi=_score_finalize(*inp,dispersion='site',reference='chi2')[0]
    np.testing.assert_array_equal(adaptive,explicit)
    assert np.all(adaptive>chi)


def test_depth_specific_confidence_intervals_are_group_symmetric():
    from epykit._glm import newcombe_diff_ci
    m1,n1,m2,n2=map(np.array,([40.],[100.],[10.],[20.]))
    first=newcombe_diff_ci(m1,n1,m2,n2,phi_a=np.array([8.]),phi_b=np.array([2.]),df=np.array([8.]))
    reversed_ci=newcombe_diff_ci(m2,n2,m1,n1,phi_a=np.array([2.]),phi_b=np.array([8.]),df=np.array([8.]))
    np.testing.assert_allclose(first[0],-reversed_ci[1])
    np.testing.assert_allclose(first[1],-reversed_ci[0])
