"""Independent count-likelihood oracles for the new engine.

These catch replacing BB with quasi-binomial deviance, pooling replicates,
using the pooled fraction as the BB mean MLE, or mishandling boundaries.
"""
import importlib.util

import numpy as np
import pytest
from scipy.optimize import minimize_scalar
from scipy.stats import betabinom, binom


def _core():
    spec = importlib.util.find_spec('epykit._beta_binomial')
    assert spec is not None, 'The genuine beta-binomial likelihood core is missing'
    from epykit import _beta_binomial
    return _beta_binomial


@pytest.mark.parametrize('mu,rho', [(.01,.001),(.5,.15),(.95,.35),(.3,1e-7)])
def test_logpmf_matches_scipy(mu, rho):
    bb = _core()
    n = np.array([1,4,20,80,1000])
    m = np.array([0,1,8,79,900])
    k = (1-rho)/rho
    expected = betabinom.logpmf(m,n,mu*k,(1-mu)*k)
    np.testing.assert_allclose(bb.logpmf(m,n,mu,rho), expected, atol=8e-8, rtol=1e-8)


def test_binomial_limit_and_point_mass_means():
    bb = _core()
    n = np.array([4,20,40])
    m = np.array([0,9,40])
    np.testing.assert_allclose(bb.logpmf(m,n,.3,0), binom.logpmf(m,n,.3), atol=1e-12)
    assert bb.logpmf(0,20,0,.15) == 0
    assert bb.logpmf(20,20,1,.15) == 0
    assert np.isneginf(bb.logpmf(1,20,0,.15))
    assert np.isneginf(bb.logpmf(19,20,1,.15))
    assert bb.logpmf(0,0,.3,.15) == 0


def test_uneven_coverage_mean_is_beta_binomial_mle():
    bb = _core()
    m = np.array([0,0,0,20,40])
    n = np.array([4,7,13,23,40])
    rho = .15
    k = (1-rho)/rho
    oracle = minimize_scalar(lambda mu: -betabinom.logpmf(m,n,mu*k,(1-mu)*k).sum(),
                             bounds=(1e-10,1-1e-10),method='bounded',options={'xatol':1e-12})
    mu,ll,information = bb.fit_mean(m,n,rho)
    assert abs(mu-oracle.x)<2e-7
    assert abs(ll+oracle.fun)<1e-9
    assert information>0
    assert abs(mu-m.sum()/n.sum())>.1


def test_true_two_group_lr_matches_independent_fits():
    bb = _core()
    m = np.array([0,1,5,17,39,1,1,1,3,5])
    n = np.array([4,7,13,23,40,4,7,13,23,40])
    rho = .15
    k = (1-rho)/rho
    def fit(x,y):
        return minimize_scalar(lambda mu: -betabinom.logpmf(x,y,mu*k,(1-mu)*k).sum(),
            bounds=(1e-10,1-1e-10),method='bounded',options={'xatol':1e-12})
    a,b,h0 = fit(m[:5],n[:5]),fit(m[5:],n[5:]),fit(m,n)
    result = bb.fit_fixed_dispersion(m[None,:],n[None,:],5,np.array([rho]))
    assert abs(result['lr'][0]-2*(h0.fun-a.fun-b.fun))<1e-8
    np.testing.assert_allclose([result['mu_case'][0],result['mu_control'][0],result['mu_null'][0]],
        [a.x,b.x,h0.x],atol=2e-7)
    swap = bb.fit_fixed_dispersion(m[None,::-1],n[None,::-1],5,np.array([rho]))
    np.testing.assert_allclose(swap['lr'],result['lr'],atol=1e-10)
    np.testing.assert_allclose(swap['mu_case'],result['mu_control'],atol=1e-10)


def test_zero_coverage_does_not_become_zero_methylation():
    bb = _core()
    a = bb.fit_mean(np.array([0,1,3]),np.array([0,4,10]),.1)
    b = bb.fit_mean(np.array([1,3]),np.array([4,10]),.1)
    np.testing.assert_allclose(a,b,atol=1e-12)
    assert bb.fit_mean(np.array([0,0]),np.array([4,20]),.1)[0] == 0
    assert bb.fit_mean(np.array([4,20]),np.array([4,20]),.1)[0] == 1
    assert np.isnan(bb.fit_mean(np.array([0]),np.array([0]),.1)[0])


@pytest.mark.parametrize('rho',[.0003,.01,.35])
def test_high_coverage_mean_and_information_match_likelihood_oracle(rho):
    bb = _core()
    m = np.array([1,40,100,4000])
    n = np.array([32,300,5000,10000])
    k = (1-rho)/rho
    def ll(mu):
        return betabinom.logpmf(m,n,mu*k,(1-mu)*k).sum()
    oracle = minimize_scalar(lambda p:-ll(p),bounds=(1e-9,1-1e-9),method='bounded',
                             options={'xatol':1e-12})
    mu,actual,information = bb.fit_mean(m,n,rho)
    assert abs(mu-oracle.x)<3e-7
    assert abs(actual+oracle.fun)<1e-7
    h = 1e-5
    independent_information = -(ll(mu+h)-2*ll(mu)+ll(mu-h))/(h*h)
    assert abs(information-independent_information)/information<.002


@pytest.mark.parametrize('m,mu', [([0,0],0.),([300,300],1.)])
def test_high_coverage_endpoint_fits_have_finite_information(m,mu):
    bb = _core()
    fitted,ll,information = bb.fit_mean(np.array(m),np.array([300,300]),.1)
    assert fitted == mu
    assert ll == 0
    assert np.isfinite(information) and information>0


@pytest.mark.parametrize('m,n', [([2],[1]),([-1],[10]),([1.5],[10]),([1],[10.5])])
def test_invalid_counts_are_refused(m,n):
    bb = _core()
    with pytest.raises(ValueError):
        bb.fit_mean(np.array(m),np.array(n),.1)
