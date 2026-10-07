"""Count-based prior and conditional-LR oracles, independent of benchmark truth."""
import numpy as np
from scipy.special import logsumexp
from scipy.stats import betabinom, f
from scipy.stats import norm
from scipy.optimize import brentq, minimize_scalar
from epykit import _beta_binomial as bb


def test_subsampled_prior_and_pvalues_are_invariant_to_group_relabeling():
    rng = np.random.default_rng(4710921)
    sites,nrep = 6000,3
    mean = rng.uniform(.01,.99,sites)
    rho = rng.choice([.002,.04,.25],sites,p=[.4,.4,.2])
    n = np.full((sites,2*nrep),20,dtype=np.int64)
    k = (1-rho)/rho
    latent = rng.beta((mean*k)[:,None],((1-mean)*k)[:,None],size=n.shape)
    m = rng.binomial(n,latent)
    order = np.r_[np.arange(nrep,2*nrep),np.arange(nrep)]
    prior = bb.estimate_prior(m,n,nrep)
    swapped = bb.estimate_prior(m[:,order],n[:,order],nrep)
    np.testing.assert_allclose([prior['log_rho_mean'],prior['log_rho_sd']],
                              [swapped['log_rho_mean'],swapped['log_rho_sd']],atol=1e-9,rtol=0)
    a = bb.test_lognormal(m,n,nrep,prior=prior)
    b = bb.test_lognormal(m[:,order],n[:,order],nrep,prior=swapped)
    # Near LR=0, the F tail's square-root slope magnifies likelihood roundoff.
    np.testing.assert_allclose(a['pvalue'],b['pvalue'],atol=1e-6,rtol=0)
    np.testing.assert_equal(a['pvalue']<.05,b['pvalue']<.05)
    np.testing.assert_allclose(a['mu_case']-a['mu_control'],
                              -(b['mu_case']-b['mu_control']),atol=1e-7,rtol=0)


def test_factorized_count_prior_em_matches_direct_likelihood():
    assert hasattr(bb,'fit_prior'), 'Count-likelihood empirical prior is missing'
    rng = np.random.default_rng(7710821)
    n = np.full((2000,5),40)
    mu = rng.choice([.15,.5,.85],len(n),p=[.25,.5,.25])
    rho = rng.choice([.005,.2],len(n),p=[.7,.3])
    k = (1-rho)/rho
    m = rng.binomial(n,rng.beta((mu*k)[:,None],((1-mu)*k)[:,None],size=n.shape))
    prior = bb.fit_prior(m,n,mean_grid=np.array([.15,.5,.85]),rho_grid=np.array([.005,.2]))
    np.testing.assert_allclose(prior['rho_weights'],[.7,.3],atol=.06)
    np.testing.assert_allclose(prior['mean_weights'],[.25,.5,.25],atol=.055)
    logL = np.empty((len(m),3,2))
    for a,p in enumerate(prior['mean_grid']):
        for b,r in enumerate(prior['rho_grid']):
            c = (1-r)/r
            logL[:,a,b] = betabinom.logpmf(m,n,p*c,(1-p)*c).sum(1)
    oracle = logsumexp(logL+np.log(prior['mean_weights'])[None,:,None]
                       +np.log(prior['rho_weights'])[None,None,:],axis=(1,2)).sum()
    assert abs(prior['loglik']-oracle)<1e-7
    assert np.min(np.diff(prior['loglik_trace']))>=-1e-7
    perm = rng.permutation(len(m))
    again = bb.fit_prior(m[perm],n[perm],mean_grid=prior['mean_grid'],rho_grid=prior['rho_grid'])
    np.testing.assert_allclose(again['rho_weights'],prior['rho_weights'],atol=1e-9)


def test_profile_dispersion_map_and_f_reference_match_direct_objective():
    assert hasattr(bb,'test_counts'), 'Moderated count-likelihood testing is missing'
    m = np.array([[0,1,5,17,39,1,1,1,3,5], [1,2,4,5,9,0,1,3,8,19]])
    n = np.tile([4,7,13,23,40,4,7,13,23,40],(2,1))
    prior = dict(rho_grid=np.array([.001,.03,.15,.4]),rho_weights=np.array([.1,.3,.4,.2]))
    objectives = []
    fits = []
    for rho,weight in zip(prior['rho_grid'],prior['rho_weights']):
        fit = bb.fit_fixed_dispersion(m,n,5,rho)
        fits.append(fit)
        objectives.append(fit['loglik_full']+np.log(weight))
    best = np.argmax(objectives,axis=0)
    result = bb.test_counts(m,n,5,prior=prior,reference='F',intervals=False)
    np.testing.assert_equal(result['rho'],prior['rho_grid'][best])
    lr = np.array([fits[best[i]]['lr'][i] for i in range(2)])
    np.testing.assert_allclose(result['lr'],lr,atol=1e-10)
    np.testing.assert_allclose(result['pvalue'],f.sf(lr,1,8),atol=1e-12)
    np.testing.assert_equal(result['df'],[8,8])


def test_missing_replicates_mask_inference_and_set_local_reference_df():
    assert hasattr(bb,'test_counts'), 'Count-likelihood testing is missing'
    m = np.array([[0,1,3,3,0,0,1,2], [0,0,0,1,0,0,1,2]])
    n = np.array([[0,4,8,10,0,4,8,10], [0,0,0,10,0,4,8,10]])
    prior = dict(rho_grid=np.array([.01,.1]),rho_weights=np.array([.5,.5]))
    result = bb.test_counts(m,n,4,prior=prior,intervals=False)
    assert result['df'][0] == 4
    assert np.isfinite(result['pvalue'][0])
    assert np.isnan(result['pvalue'][1])


def test_difference_profile_interval_inverts_the_same_conditional_lr():
    m = np.array([[0,1,5,17,39,1,1,1,3,5], [0,0,0,0,0,0,0,0,0,0]])
    n = np.tile([4,7,13,23,40,4,7,13,23,40],(2,1))
    prior = dict(rho_grid=np.array([.15]),rho_weights=np.array([1.]))
    result = bb.test_counts(m,n,5,prior=prior,reference='F',intervals=True)
    lo,hi = result['ci_lo'],result['ci_hi']
    assert np.all(lo <= result['mu_case']-result['mu_control'])
    assert np.all(hi >= result['mu_case']-result['mu_control'])
    assert lo[1]<0<hi[1]  # identical zero-count groups still have uncertainty
    critical = f.ppf(.95,1,8)
    k = (1-.15)/.15
    def conditional_loglik(delta):
        a,b = max(0.,-delta),min(1.,1-delta)
        def ll(mu):
            return (betabinom.logpmf(m[0,:5],n[0,:5],(mu+delta)*k,(1-mu-delta)*k).sum()
                    +betabinom.logpmf(m[0,5:],n[0,5:],mu*k,(1-mu)*k).sum())
        opt = minimize_scalar(lambda p:-ll(p),bounds=(a+1e-10,b-1e-10),method='bounded',
                              options={'xatol':1e-12})
        return -opt.fun
    def root(delta):
        return 2*(result['loglik_full'][0]-conditional_loglik(delta))-critical
    point = result['mu_case'][0]-result['mu_control'][0]
    oracle_lo = brentq(root,-.99,point)
    oracle_hi = brentq(root,point,.99)
    np.testing.assert_allclose([lo[0],hi[0]],[oracle_lo,oracle_hi],atol=2e-7)
    assert ((lo[0]>0 or hi[0]<0) == (result['pvalue'][0]<.05))


def test_continuous_log_dispersion_shrinkage_matches_numerical_profile_oracle():
    assert hasattr(bb,'test_lognormal'), 'Continuous dispersion shrinkage is missing'
    m = np.array([[0,1,5,17,39,1,1,1,3,5]])
    n = np.array([[4,7,13,23,40,4,7,13,23,40]])
    prior = dict(log_rho_mean=np.log(.04),log_rho_sd=1.2)
    def objective(z):
        r = np.exp(z)
        fit = bb.fit_fixed_dispersion(m,n,5,r)
        return -(fit['loglik_full'][0]+norm.logpdf(z,prior['log_rho_mean'],prior['log_rho_sd']))
    # Dense bracketing avoids assuming that the profile objective is unimodal.
    grid = np.linspace(np.log(1e-6),np.log(.95),101)
    best = np.argmin([objective(z) for z in grid])
    lo,hi = grid[max(0,best-1)],grid[min(len(grid)-1,best+1)]
    oracle = minimize_scalar(objective,bounds=(lo,hi),method='bounded',options={'xatol':1e-10})
    result = bb.test_lognormal(m,n,5,prior=prior,intervals=False)
    assert abs(np.log(result['rho'][0])-oracle.x)<2e-4
    fixed = bb.fit_fixed_dispersion(m,n,5,result['rho'])
    np.testing.assert_allclose(result['lr'],fixed['lr'],atol=1e-9)
    np.testing.assert_allclose(result['pvalue'],f.sf(fixed['lr'],1,8),atol=1e-12)


def test_log_prior_projection_uses_latent_dispersion_mixture_moments():
    assert hasattr(bb,'project_log_prior'), 'Continuous prior projection is missing'
    prior = dict(rho_grid=np.array([0.,.01,.2]),rho_weights=np.array([.1,.6,.3]))
    result = bb.project_log_prior(prior)
    log_rho = np.log([.005,.01,.2])
    center = .1*log_rho[0]+.6*log_rho[1]+.3*log_rho[2]
    variance = np.dot([.1,.6,.3],(log_rho-center)**2)
    assert abs(result['log_rho_mean']-center)<1e-12
    assert abs(result['log_rho_sd']-np.sqrt(variance))<1e-12
    assert result['binomial_representative'] == .005


def test_high_depth_zero_group_runs_through_shrinkage_and_intervals():
    result = bb.test_lognormal(np.array([[0,0,30,30]]),np.full((1,4),300),2,
        prior={'log_rho_mean':np.log(.05),'log_rho_sd':1.},intervals=True)
    assert np.isfinite(result['pvalue']).all()
    assert result['mu_case'][0] == 0
    assert result['ci_lo'][0]<result['ci_hi'][0]


def test_sparse_union_training_samples_eligible_groups_before_subsampling():
    m = np.zeros((100000,4),dtype=int)
    n = np.zeros_like(m)
    n[:,0] = 20
    n[-100:] = 20
    m[-100:] = np.tile([2,8,3,7],(100,1))
    prior = bb.estimate_prior(m,n,2)
    assert prior['n_groups'] == 200
    assert np.isfinite(prior['log_rho_mean'])
