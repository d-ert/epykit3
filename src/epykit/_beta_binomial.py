"""Replicate-count beta-binomial likelihood and numerical fitting.

rho is biological overdispersion: Var(M)=N*mu*(1-mu)*(1+(N-1)*rho).
Coverage is conditioned on; missing coverage is not an observed zero.
This module owns numerical inference, not store I/O or region selection.
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit
from scipy import special, stats

ALGORITHM_REVISION = 'count-likelihood-eb-bb-lognormal-F-v3-experimental'


@njit(cache=True)
def _rising_correction(z: float, count: int) -> float:
    """log[(z)_count / z**count], without large-concentration cancellation."""
    if count <= 1:
        return 0.0
    if count <= 256:
        result = 0.0
        for j in range(1, count):
            result += math.log1p(j / z)
        return result
    if z > 1e6 and count / z < 1e-4:
        c = float(count)
        s1 = c * (c-1) / 2
        s2 = c * (c-1) * (2*c-1) / 6
        s3 = s1*s1
        return s1/z - s2/(2*z*z) + s3/(3*z*z*z)
    return math.lgamma(z+count) - math.lgamma(z) - count*math.log(z)


@njit(cache=True)
def _logpmf(m: int, n: int, mu: float, rho: float) -> float:
    if n == 0:
        return 0.0
    if mu <= 0.0:
        return 0.0 if m == 0 else -math.inf
    if mu >= 1.0:
        return 0.0 if m == n else -math.inf
    result = (math.lgamma(n+1) - math.lgamma(m+1) - math.lgamma(n-m+1)
              + m*math.log(mu) + (n-m)*math.log1p(-mu))
    if rho == 0.0:
        return result
    k = (1-rho)/rho
    return (result + _rising_correction(mu*k,m)
            + _rising_correction((1-mu)*k,n-m) - _rising_correction(k,n))


@njit(cache=True)
def _digamma_trigamma(x: float) -> tuple[float,float]:
    d = 0.0
    t = 0.0
    while x < 10:
        d -= 1/x
        t += 1/(x*x)
        x += 1
    r = 1/x
    r2 = r*r
    d += math.log(x)-r/2-r2*(1/12-r2*(1/120-r2*(1/252-r2*(1/240-r2/132))))
    t += r+r2/2+r*r2*(1/6-r2*(1/30-r2*(1/42-r2*(1/30-r2*5/66))))
    return d,t


@njit(cache=True)
def _score_information(m, n, mu: float, rho: float, start: int, stop: int):
    score = 0.0
    information = 0.0
    if rho == 0:
        for j in range(start,stop):
            if n[j] == 0:
                continue
            if m[j] > 0:
                score += m[j]/mu
                information += m[j]/(mu*mu)
            if n[j] > m[j]:
                score -= (n[j]-m[j])/(1-mu)
                information += (n[j]-m[j])/((1-mu)*(1-mu))
        return score,information
    k = (1-rho)/rho
    for j in range(start,stop):
        if n[j] <= 256:
            # k/(mu*k+h)=1/(mu+h/k): stable even for rho near zero.
            for h in range(m[j]):
                v = 1/(mu+h/k)
                score += v
                information += v*v
            for h in range(n[j]-m[j]):
                v = 1/(1-mu+h/k)
                score -= v
                information += v*v
        else:
            a,b = mu*k,(1-mu)*k
            if m[j]>0:
                da,ta = _digamma_trigamma(a)
                dam,tam = _digamma_trigamma(a+m[j])
                score += k*(dam-da)
                information += k*k*(ta-tam)
            if n[j]>m[j]:
                db,tb = _digamma_trigamma(b)
                dbu,tbu = _digamma_trigamma(b+n[j]-m[j])
                score -= k*(dbu-db)
                information += k*k*(tb-tbu)
    return score,information


@njit(cache=True)
def _group_loglik(m,n,mu,rho,start,stop):
    total = 0.0
    if 0<mu<1 and rho>=1e-5:
        k = (1-rho)/rho
        a,b = mu*k,(1-mu)*k
        base = math.lgamma(k)-math.lgamma(a)-math.lgamma(b)
        for j in range(start,stop):
            if n[j]>0:
                total += (math.lgamma(n[j]+1)-math.lgamma(m[j]+1)-math.lgamma(n[j]-m[j]+1)
                          +math.lgamma(m[j]+a)+math.lgamma(n[j]-m[j]+b)
                          -math.lgamma(n[j]+k)+base)
        return total
    for j in range(start,stop):
        total += _logpmf(m[j],n[j],mu,rho)
    return total


@njit(cache=True)
def _fit_mean(m,n,rho,start,stop):
    total_m = 0
    total_n = 0
    wm = 0.0
    wn = 0.0
    for j in range(start,stop):
        total_m += m[j]
        total_n += n[j]
        if n[j] > 0:
            weight = 1/(1+(n[j]-1)*rho)
            wm += weight*m[j]
            wn += weight*n[j]
    if total_n == 0:
        return math.nan,0.0,0.0,True
    if total_m == 0:
        _,info = _score_information(m,n,0.0,rho,start,stop)
        return 0.0,0.0,info,True
    if total_m == total_n:
        _,info = _score_information(m,n,1.0,rho,start,stop)
        return 1.0,0.0,info,True
    mu = wm/wn
    lo,hi = 0.0,1.0
    converged = False
    for _ in range(60):
        score,information = _score_information(m,n,mu,rho,start,stop)
        if abs(score/information) < 5e-12:
            converged = True
            break
        if score > 0:
            lo = mu
        else:
            hi = mu
        proposal = mu+score/information
        if proposal <= lo or proposal >= hi:
            proposal = (lo+hi)/2
        if abs(proposal-mu) < 5e-12:
            mu = proposal
            converged = True
            break
        mu = proposal
    _,information = _score_information(m,n,mu,rho,start,stop)
    return mu,_group_loglik(m,n,mu,rho,start,stop),information,converged


@njit(cache=True)
def _fixed_batch(m,n,n_case,rho):
    output = np.empty((len(m),7),dtype=np.float64)
    convergence = np.ones(len(m),dtype=np.bool_)
    for i in range(len(m)):
        a,la,ia,ca = _fit_mean(m[i],n[i],rho[i],0,n_case)
        b,lb,ib,cb = _fit_mean(m[i],n[i],rho[i],n_case,m.shape[1])
        c,lc,_,cc = _fit_mean(m[i],n[i],rho[i],0,m.shape[1])
        d = max(0.0,2*(la+lb-lc)) if math.isfinite(a) and math.isfinite(b) else math.nan
        output[i] = a,b,c,d,ia,ib,la+lb
        convergence[i] = ca and cb and cc
    return output,convergence


def _counts(m,n,ndim):
    m,n = np.asarray(m),np.asarray(n)
    if m.shape != n.shape or m.ndim != ndim:
        raise ValueError(f'Count arrays must have the same {ndim}-dimensional shape')
    if (not np.all(np.isfinite(m)) or not np.all(np.isfinite(n))
        or np.any(m<0) or np.any(n<0) or np.any(m>n)
        or np.any(m!=np.floor(m)) or np.any(n!=np.floor(n))
        or np.any(n>np.iinfo(np.int64).max)):
        raise ValueError('Require finite integer counts with 0 <= methylated <= coverage')
    return np.ascontiguousarray(m,dtype=np.int64),np.ascontiguousarray(n,dtype=np.int64)


def _check_rho(rho):
    value = np.asarray(rho,dtype=float)
    if np.any(~np.isfinite(value)) or np.any(value<0) or np.any(value>=1):
        raise ValueError('rho must be finite in [0,1)')
    return value


def logpmf(m,n,mu,rho):
    m,n,mu,rho = np.broadcast_arrays(m,n,mu,rho)
    shape = m.shape
    m,n = _counts(m.reshape(-1),n.reshape(-1),1)
    rho = _check_rho(rho).reshape(-1)
    mu = np.asarray(mu,dtype=float).reshape(-1)
    if np.any(~np.isfinite(mu)) or np.any(mu<0) or np.any(mu>1):
        raise ValueError('mu must be finite in [0,1]')
    result = np.array([_logpmf(m[i],n[i],mu[i],rho[i]) for i in range(len(m))]).reshape(shape)
    return result.item() if shape == () else result


def fit_mean(m,n,rho):
    m,n = _counts(m,n,1)
    rho = float(_check_rho(rho))
    mu,ll,information,converged = _fit_mean(m,n,rho,0,len(m))
    if not converged:
        raise RuntimeError('Beta-binomial mean optimization did not converge')
    return mu,ll,information


def fit_fixed_dispersion(m,n,n_case,rho):
    m,n = _counts(m,n,2)
    if not isinstance(n_case,(int,np.integer)) or not 0<n_case<m.shape[1]:
        raise ValueError('n_case must split the columns into two nonempty groups')
    rho = np.broadcast_to(_check_rho(rho),(len(m),)).copy()
    output,convergence = _fixed_batch(m,n,n_case,rho)
    if not np.all(convergence):
        raise RuntimeError(f'Beta-binomial mean optimization failed at {(~convergence).sum()} sites')
    names = ('mu_case','mu_control','mu_null','lr','information_case','information_control','loglik_full')
    return {name:output[:,j] for j,name in enumerate(names)}


def _mixture_loglik(m,n,mean_grid,rho_grid):
    """Training likelihood, bounded in group count; no fitted site variances."""
    choose = special.gammaln(n+1)-special.gammaln(m+1)-special.gammaln(n-m+1)
    result = np.empty((len(m),len(mean_grid),len(rho_grid)))
    for a,mu in enumerate(mean_grid):
        if mu == 0 or mu == 1:
            compatible = (m==0).all(1) if mu == 0 else (m==n).all(1)
            result[:,a,:] = np.where(compatible,0.,-np.inf)[:,None]
            continue
        for b,rho in enumerate(rho_grid):
            if rho == 0:
                logp = choose+m*np.log(mu)+(n-m)*np.log1p(-mu)
            else:
                k = (1-rho)/rho
                logp = (choose+special.betaln(m+mu*k,n-m+(1-mu)*k)
                        -special.betaln(mu*k,(1-mu)*k))
            result[:,a,b] = np.where(n>0,logp,0.).sum(1)
    return result


def _fit_mixture(log_likelihood,mean_initial,rho_initial,tol,max_iter):
    shift = np.max(log_likelihood,axis=(1,2))
    if not np.all(np.isfinite(shift)):
        raise ValueError('Prior grid gives zero probability to an observed group')
    likelihood = np.exp(log_likelihood-shift[:,None,None])
    mean_weight = np.asarray(mean_initial,dtype=float).copy()
    rho_weight = np.asarray(rho_initial,dtype=float).copy()
    trace = []
    converged = False
    for iteration in range(max_iter+1):
        posterior = likelihood*mean_weight[None,:,None]*rho_weight[None,None,:]
        marginal = posterior.sum(axis=(1,2))
        loglik = float((np.log(marginal)+shift).sum())
        trace.append(loglik)
        if len(trace)>1 and abs(trace[-1]-trace[-2]) <= tol*len(likelihood):
            converged = True
            break
        if iteration == max_iter:
            break
        posterior /= marginal[:,None,None]
        mean_weight = np.maximum(posterior.sum(axis=(0,2))/len(likelihood),1e-300)
        rho_weight = np.maximum(posterior.sum(axis=(0,1))/len(likelihood),1e-300)
        mean_weight /= mean_weight.sum()
        rho_weight /= rho_weight.sum()
    return dict(mean_weights=mean_weight,rho_weights=rho_weight,
                loglik=loglik,loglik_trace=np.asarray(trace),converged=converged)


def fit_prior(m,n,*,mean_grid=None,rho_grid=None,max_groups=4096,tol=1e-5,max_iter=1000):
    """Likelihood EB on a factorized mean × dispersion mixing distribution.

    Every row is one replicate group, not one pooled count. Observed count
    likelihood accounts for sampling noise. Mean–rho independence and the
    discretization are working assumptions, recorded in the result.
    """
    m,n = _counts(m,n,2)
    valid = (n>0).sum(1)>=2
    m,n = m[valid],n[valid]
    if len(m)==0:
        raise ValueError('Dispersion prior requires groups with at least two covered replicates')
    if max_groups<1:
        raise ValueError('max_groups must be positive')
    if len(m)>max_groups:
        indices = np.random.default_rng(20261007).choice(len(m),max_groups,replace=False)
        m,n = m[indices],n[indices]
    if mean_grid is None:
        mean_grid = (1-np.cos(np.linspace(0,np.pi,41)))/2
    if rho_grid is None:
        rho_grid = np.r_[0.,np.geomspace(1e-4,.8,25)]
    mean_grid = np.asarray(mean_grid,dtype=float)
    rho_grid = _check_rho(rho_grid)
    if (mean_grid.ndim!=1 or rho_grid.ndim!=1 or len(mean_grid)<1 or len(rho_grid)<1
        or np.any(~np.isfinite(mean_grid)) or np.any(mean_grid<0) or np.any(mean_grid>1)
        or np.any(np.diff(mean_grid)<=0) or np.any(np.diff(rho_grid)<=0)):
        raise ValueError('Prior grids must be finite, increasing, and in the model parameter domain')
    likelihood = _mixture_loglik(m,n,mean_grid,rho_grid)
    uniform_mean = np.full(len(mean_grid),1/len(mean_grid))
    uniform_rho = np.full(len(rho_grid),1/len(rho_grid))
    # Two deterministic starts expose local optima of the structured mixture.
    second_rho = np.exp(-.5*((np.log(np.maximum(rho_grid,1e-5))-np.log(.05))/2)**2)
    second_rho /= second_rho.sum()
    starts = [_fit_mixture(likelihood,uniform_mean,initial,tol,max_iter)
              for initial in (uniform_rho,second_rho)]
    best = max(starts,key=lambda result: result['loglik'])
    result = dict(best)
    result.update(mean_grid=mean_grid,rho_grid=rho_grid,n_groups=len(m),
                  algorithm_revision=ALGORITHM_REVISION,
                  prior_assumption='factorized mean and rho count mixture',
                  start_logliks=[s['loglik'] for s in starts],
                  start_rho_weights=[s['rho_weights'] for s in starts],
                  informative_groups=int(((m.sum(1)>0)&(m.sum(1)<n.sum(1))).sum()))
    return result


@njit(cache=True)
def _map_batch(m,n,n_case,rho_grid,log_weights):
    output = np.empty((len(m),8),dtype=np.float64)
    converged = np.ones(len(m),dtype=np.bool_)
    for i in range(len(m)):
        best = -math.inf
        best_a,best_b,best_ia,best_ib,best_ll,best_rho = 0.,0.,0.,0.,0.,0.
        best_converged = True
        for j in range(len(rho_grid)):
            rho = rho_grid[j]
            a,la,ia,ca = _fit_mean(m[i],n[i],rho,0,n_case)
            b,lb,ib,cb = _fit_mean(m[i],n[i],rho,n_case,m.shape[1])
            objective = la+lb+log_weights[j]
            if objective>best:
                best = objective
                best_a,best_b,best_ia,best_ib = a,b,ia,ib
                best_ll,best_rho = la+lb,rho
                best_converged = ca and cb
        c,lc,_,cc = _fit_mean(m[i],n[i],best_rho,0,m.shape[1])
        lr = max(0.,2*(best_ll-lc)) if math.isfinite(best_a) and math.isfinite(best_b) else math.nan
        output[i] = best_a,best_b,c,lr,best_ia,best_ib,best_ll,best_rho
        converged[i] = best_converged and cc
    return output,converged


@njit(cache=True)
def _difference_loglik(m,n,n_case,rho,delta,initial):
    """Profile the common nuisance mean with mu_case-mu_control=delta."""
    lo,hi = max(0.,-delta),min(1.,1-delta)
    if hi<=lo:
        return (_group_loglik(m,n,lo+delta,rho,0,n_case)
                +_group_loglik(m,n,lo,rho,n_case,len(m)))
    left,right = lo,hi
    epsilon = min(1e-10,(hi-lo)*1e-5)
    mu = max(lo+epsilon,min(hi-epsilon,initial))
    for _ in range(60):
        sa,ia = _score_information(m,n,mu+delta,rho,0,n_case)
        sb,ib = _score_information(m,n,mu,rho,n_case,len(m))
        step = (sa+sb)/(ia+ib)
        if abs(step)<5e-12:
            break
        if sa+sb>0:
            left = mu
        else:
            right = mu
        candidate = mu+step
        if candidate<=left or candidate>=right:
            candidate = (left+right)/2
        if candidate==mu:
            break
        mu = candidate
    best = (_group_loglik(m,n,mu+delta,rho,0,n_case)
            +_group_loglik(m,n,mu,rho,n_case,len(m)))
    # Finite endpoint likelihoods matter when a group has only zeros or ones.
    for endpoint in (lo,hi):
        ll = (_group_loglik(m,n,endpoint+delta,rho,0,n_case)
              +_group_loglik(m,n,endpoint,rho,n_case,len(m)))
        best = max(best,ll)
    return best


@njit(cache=True)
def _difference_intervals(m,n,n_case,fits,critical):
    result = np.full((len(m),2),np.nan)
    for i in range(len(m)):
        a,b,rho,ll = fits[i,0],fits[i,1],fits[i,7],fits[i,6]
        if not math.isfinite(a) or not math.isfinite(b):
            continue
        point = a-b
        for side in range(2):
            outside = -1. if side==0 else 1.
            inside = point
            endpoint_ll = _difference_loglik(m[i],n[i],n_case,rho,outside,b)
            if 2*(ll-endpoint_ll)<=critical[i]:
                result[i,side] = outside
                continue
            for _ in range(28):
                delta = (outside+inside)/2
                constrained = _difference_loglik(m[i],n[i],n_case,rho,delta,b)
                if 2*(ll-constrained)>critical[i]:
                    outside = delta
                else:
                    inside = delta
            result[i,side] = (outside+inside)/2
    return result


def test_counts(m,n,n_case,*,prior=None,reference='adaptive',intervals=False):
    """Conditional BB LR with count-likelihood EB profile-MAP dispersion.

    Adaptive uses the methylSig-inspired F(1,n_valid-2) approximation.
    This reference is not an exact null law or posterior-uncertainty integral.
    Explicit chi2 is available for calibration diagnostics.
    """
    m,n = _counts(m,n,2)
    if not isinstance(n_case,(int,np.integer)) or not 0<n_case<m.shape[1]:
        raise ValueError('n_case must split the columns into two nonempty groups')
    if reference not in {'adaptive','F','chi2'}:
        raise ValueError('reference must be adaptive, F or chi2')
    if prior is None:
        count = min(len(m),2048)
        indices = np.random.default_rng(20261007).choice(len(m),count,replace=False)
        cols = max(n_case,m.shape[1]-n_case)
        training_m = np.zeros((2*count,cols),dtype=np.int64)
        training_n = np.zeros_like(training_m)
        training_m[:count,:n_case] = m[indices,:n_case]
        training_n[:count,:n_case] = n[indices,:n_case]
        training_m[count:,:m.shape[1]-n_case] = m[indices,n_case:]
        training_n[count:,:m.shape[1]-n_case] = n[indices,n_case:]
        prior = fit_prior(training_m,training_n)
    rho_grid = _check_rho(prior['rho_grid'])
    weights = np.asarray(prior['rho_weights'],dtype=float)
    if (rho_grid.ndim!=1 or weights.shape!=rho_grid.shape or len(weights)==0
        or np.any(~np.isfinite(weights)) or np.any(weights<0) or weights.sum()<=0):
        raise ValueError('Dispersion prior must have nonnegative finite weights on its rho grid')
    log_weights = np.log(np.maximum(weights/weights.sum(),1e-300))
    output,converged = _map_batch(m,n,n_case,rho_grid,log_weights)
    if not np.all(converged):
        raise RuntimeError(f'Beta-binomial likelihood optimization failed at {(~converged).sum()} sites')
    names = ('mu_case','mu_control','mu_null','lr','information_case','information_control','loglik_full','rho')
    result = {name:output[:,j] for j,name in enumerate(names)}
    n_case_valid = (n[:,:n_case]>0).sum(1)
    n_ctrl_valid = (n[:,n_case:]>0).sum(1)
    df = np.maximum(n_case_valid+n_ctrl_valid-2,1)
    p = stats.chi2.sf(result['lr'],1) if reference=='chi2' else stats.f.sf(result['lr'],1,df)
    p[(n_case_valid<2)|(n_ctrl_valid<2)] = np.nan
    result.update(pvalue=p,df=df,prior=prior,reference=reference)
    if intervals:
        critical = np.full(len(m),stats.chi2.ppf(.95,1)) if reference=='chi2' else stats.f.ppf(.95,1,df)
        ci = _difference_intervals(m,n,n_case,output,critical)
        ci[(n_case_valid<2)|(n_ctrl_valid<2)] = np.nan
        result.update(ci_lo=ci[:,0],ci_hi=ci[:,1])
    return result


@njit(cache=True)
def _log_dispersion_objective(m,n,n_case,z,center,sd):
    rho = math.exp(z)
    a,la,ia,ca = _fit_mean(m,n,rho,0,n_case)
    b,lb,ib,cb = _fit_mean(m,n,rho,n_case,len(m))
    value = la+lb-.5*((z-center)/sd)**2
    return value,a,b,ia,ib,la+lb,ca and cb


@njit(cache=True)
def _lognormal_batch(m,n,n_case,center,sd):
    output = np.empty((len(m),8))
    convergence = np.ones(len(m),dtype=np.bool_)
    lower,upper = math.log(1e-6),math.log(.95)
    grid = np.linspace(lower,upper,17)
    ratio = (math.sqrt(5)-1)/2
    for i in range(len(m)):
        ma,na = np.sum(m[i,:n_case]),np.sum(n[i,:n_case])
        mb,nb = np.sum(m[i,n_case:]),np.sum(n[i,n_case:])
        if (ma==0 or ma==na) and (mb==0 or mb==nb):
            # Both mean fits are exact endpoints. Their likelihood is one
            # at every rho, so log-rho's mode is exactly the prior center.
            z = min(upper,max(lower,center))
            _,a,b,ia,ib,ll,ca = _log_dispersion_objective(m[i],n[i],n_case,z,center,sd)
            rho = math.exp(z)
            c,lc,_,cc = _fit_mean(m[i],n[i],rho,0,m.shape[1])
            lr = max(0.,2*(ll-lc)) if math.isfinite(a) and math.isfinite(b) else math.nan
            output[i] = a,b,c,lr,ia,ib,ll,rho
            convergence[i] = ca and cc
            continue
        values = np.empty(len(grid))
        best_index = 0
        for j in range(len(grid)):
            values[j] = _log_dispersion_objective(m[i],n[i],n_case,grid[j],center,sd)[0]
            if values[j]>values[best_index]:
                best_index = j
        lo = grid[max(0,best_index-1)]
        hi = grid[min(len(grid)-1,best_index+1)]
        x1,x2 = hi-ratio*(hi-lo),lo+ratio*(hi-lo)
        f1 = _log_dispersion_objective(m[i],n[i],n_case,x1,center,sd)[0]
        f2 = _log_dispersion_objective(m[i],n[i],n_case,x2,center,sd)[0]
        for _ in range(28):
            if hi-lo<1e-5:
                break
            if f1>f2:
                hi,x2,f2 = x2,x1,f1
                x1 = hi-ratio*(hi-lo)
                f1 = _log_dispersion_objective(m[i],n[i],n_case,x1,center,sd)[0]
            else:
                lo,x1,f1 = x1,x2,f2
                x2 = lo+ratio*(hi-lo)
                f2 = _log_dispersion_objective(m[i],n[i],n_case,x2,center,sd)[0]
        z = (lo+hi)/2
        if values[best_index]>_log_dispersion_objective(m[i],n[i],n_case,z,center,sd)[0]:
            z = grid[best_index]
        _,a,b,ia,ib,ll,ca = _log_dispersion_objective(m[i],n[i],n_case,z,center,sd)
        rho = math.exp(z)
        c,lc,_,cc = _fit_mean(m[i],n[i],rho,0,m.shape[1])
        lr = max(0.,2*(ll-lc)) if math.isfinite(a) and math.isfinite(b) else math.nan
        output[i] = a,b,c,lr,ia,ib,ll,rho
        convergence[i] = ca and cc
    return output,convergence


def test_lognormal(m,n,n_case,*,prior,reference='adaptive',intervals=False):
    """True BB profile likelihood with a continuous normal prior on log rho.

    The mode is in log-rho coordinates, as in log-dispersion shrinkage.
    F remains a small-sample approximation; intervals condition on fitted rho.
    """
    m,n = _counts(m,n,2)
    if not isinstance(n_case,(int,np.integer)) or not 0<n_case<m.shape[1]:
        raise ValueError('n_case must split columns into two nonempty groups')
    center,sd = float(prior['log_rho_mean']),float(prior['log_rho_sd'])
    if not np.isfinite(center) or not np.isfinite(sd) or sd<=0:
        raise ValueError('Log-dispersion prior requires a finite center and positive sd')
    if reference not in {'adaptive','F','chi2'}:
        raise ValueError('reference must be adaptive, F or chi2')
    output,converged = _lognormal_batch(m,n,n_case,center,sd)
    if not np.all(converged):
        raise RuntimeError(f'Beta-binomial optimization failed at {(~converged).sum()} sites')
    names = ('mu_case','mu_control','mu_null','lr','information_case','information_control','loglik_full','rho')
    result = {name:output[:,j] for j,name in enumerate(names)}
    na,nb = (n[:,:n_case]>0).sum(1),(n[:,n_case:]>0).sum(1)
    df = np.maximum(na+nb-2,1)
    p = stats.chi2.sf(result['lr'],1) if reference=='chi2' else stats.f.sf(result['lr'],1,df)
    p[(na<2)|(nb<2)] = np.nan
    result.update(pvalue=p,df=df,prior=prior,reference=reference)
    if intervals:
        critical = np.full(len(m),stats.chi2.ppf(.95,1)) if reference=='chi2' else stats.f.ppf(.95,1,df)
        ci = _difference_intervals(m,n,n_case,output,critical)
        ci[(na<2)|(nb<2)] = np.nan
        result.update(ci_lo=ci[:,0],ci_hi=ci[:,1])
    return result


def project_log_prior(prior):
    """Smooth the count-learned latent rho distribution into a lognormal prior.

    Grid weights can be weakly identified individually. Log-scale moments
    give a continuous shrinkage density rather than selecting a grid spike.
    The binomial grid state represents the lowest unresolved rho interval;
    map it to half the smallest positive node and record that approximation.
    """
    grid = _check_rho(prior['rho_grid'])
    weights = np.asarray(prior['rho_weights'],dtype=float)
    if (grid.ndim!=1 or weights.shape!=grid.shape or len(grid)==0
        or np.any(~np.isfinite(weights)) or np.any(weights<0) or weights.sum()<=0):
        raise ValueError('Invalid count-mixture dispersion prior')
    weights = weights/weights.sum()
    positive = grid[grid>0]
    representative = float(positive.min()/2) if len(positive) else 1e-6
    log_rho = np.log(np.where(grid>0,grid,representative))
    center = float(weights@log_rho)
    sd = float(max(np.sqrt(weights@((log_rho-center)**2)),.05))
    return dict(prior,log_rho_mean=center,log_rho_sd=sd,
                binomial_representative=representative,
                scoring_prior='normal distribution on log rho; moment projection of latent count mixture')


def estimate_prior(m,n,n_case,max_groups=4096):
    """Learn one continuous shrinkage prior from a bounded replicate-group sample."""
    m,n = _counts(m,n,2)
    if not 0<n_case<m.shape[1] or len(m)==0:
        raise ValueError('Prior estimation requires sites and two nonempty replicate groups')
    if max_groups<2:
        raise ValueError('max_groups must allow both replicate groups')
    case_indices = np.flatnonzero((n[:,:n_case]>0).sum(1)>=2)
    ctrl_indices = np.flatnonzero((n[:,n_case:]>0).sum(1)>=2)
    per_group = max_groups//2
    if len(case_indices)>per_group:
        case_indices = np.random.default_rng(20261007).choice(case_indices,per_group,replace=False)
    if len(ctrl_indices)>per_group:
        # Each biological group uses the same role-independent selection.
        # Consuming one RNG case-first would change priors on label reversal.
        ctrl_indices = np.random.default_rng(20261007).choice(ctrl_indices,per_group,replace=False)
    ncase,nctrl = len(case_indices),len(ctrl_indices)
    if ncase+nctrl==0:
        raise ValueError('No groups have two covered biological replicates for dispersion estimation')
    columns = max(n_case,m.shape[1]-n_case)
    training_m = np.zeros((ncase+nctrl,columns),dtype=np.int64)
    training_n = np.zeros_like(training_m)
    training_m[:ncase,:n_case] = m[case_indices,:n_case]
    training_n[:ncase,:n_case] = n[case_indices,:n_case]
    training_m[ncase:,:m.shape[1]-n_case] = m[ctrl_indices,n_case:]
    training_n[ncase:,:m.shape[1]-n_case] = n[ctrl_indices,n_case:]
    return project_log_prior(fit_prior(training_m,training_n,max_groups=max_groups))
