"""Depth-aware variance for pooled beta-binomial group proportions.

The F reference is approximate; empirical calibration is required. This module
uses only observed counts, never simulation truth or external dispersion.
"""
import numpy as np


def pooled_design_effect(sn_case, sn2_case, sn_ctrl, sn2_ctrl, rho):
    a=np.maximum(sn_case,1.);b=np.maximum(sn_ctrl,1.)
    leverage=((sn2_case/a-1)/a+(sn2_ctrl/b-1)/b)/(1/a+1/b)
    return 1+rho*np.maximum(leverage,0)


def estimate_scale(sn_case,sn2_case,nv_case,chi_case,
                   sn_ctrl,sn2_ctrl,nv_ctrl,chi_ctrl,pooled_mu,moderate=False):
    """Estimate rho from the Pearson residual expectation and map to contrast variance.

    Empirical moderation uses the observed between-site variance after an
    approximate Pearson sampling-variance subtraction in pooled-mean strata.
    This is a moment shrinkage estimator, not an exact conjugate posterior.
    """
    a=np.maximum(sn_case,1.);b=np.maximum(sn_ctrl,1.)
    df=np.maximum(nv_case+nv_ctrl-2,1.)
    exposure=(sn_case-nv_case-(sn2_case/a-1))+(sn_ctrl-nv_ctrl-(sn2_ctrl/b-1))
    valid=(exposure>0)&(nv_case>=2)&(nv_ctrl>=2)&(pooled_mu>0)&(pooled_mu<1)
    raw=np.divide(chi_case+chi_ctrl-df,exposure,out=np.zeros_like(a),where=valid)
    rho=np.clip(raw,0,.95)
    effective_df=df.copy()
    if moderate and valid.sum()>=100:
        # Fold methylation around 0.5: hyper/hypomethylation has symmetric
        # sampling behavior. Fit before clipping negative noise estimates.
        labels=np.digitize(np.minimum(pooled_mu,1-pooled_mu),[.02,.05,.1,.2,.35])
        global_center=float(np.clip(np.mean(raw[valid]),0,.95))
        for label in np.unique(labels):
            select=(labels==label)&valid
            fit=select if select.sum()>=100 else valid
            center=float(np.clip(np.mean(raw[fit]),0,.95)) if fit.sum() else global_center
            noise=2*(df[fit]+center*exposure[fit])**2/(df[fit]*exposure[fit]**2)
            between=max(float(np.var(raw[fit],ddof=1)-noise.mean()),0.)
            local_noise=2*(df+center*exposure)**2/(df*np.maximum(exposure,1e-12)**2)
            weight=np.divide(between,between+local_noise,out=np.zeros_like(a),where=(between+local_noise)>0)
            rho[select]=np.clip(weight[select]*raw[select]+(1-weight[select])*center,0,.95)
            # Keep actual residual df: the moment shrinkage approximation does
            # not justify adding arbitrary pseudo degrees of freedom.
    phi=pooled_design_effect(sn_case,sn2_case,sn_ctrl,sn2_ctrl,rho)
    return phi,effective_df,rho
