"""Hand-rolled HMM segmentation for methylation signals.

Powers three downstream callers in 0.6:
  - :mod:`epykit.pmd`  (partially methylated domains, 2-state long-range)
  - :mod:`epykit.hmr`  (hypo / low-methylated regions, 2-state short-range)
  - :mod:`epykit.dmr_hmm` (group-contrast DMR, 3-state on meth_diff)

Why a hand roll instead of pomegranate / hmmlearn?
  - pomegranate is ~50 MB and has frequently changed its API between
    minor versions; pinning it would create a heavy + brittle dep.
  - hmmlearn is well-maintained but its discrete-output HMM only
    supports multinomial emissions; we want Bernoulli per CpG with
    an explicit beta prior, which is cleaner to express directly.
  - The math here is ~150 LoC of numpy and runs in milliseconds on
    typical chrom-sized inputs.

The exposed API is a single :func:`segment` function. It runs forward-
backward (log-space) to compute posterior state probabilities, then
Viterbi to produce a hard MAP state assignment. Inputs are float
observations (beta) and integer state counts; outputs are int state
labels (and optionally the posterior matrix).
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


_LOG_EPS = -1e300  # log(0) sentinel


def _logsumexp(a: np.ndarray, axis: int = -1) -> np.ndarray:
    """Numerically stable log-sum-exp."""
    a = np.asarray(a)
    a_max = np.max(a, axis=axis, keepdims=True)
    # Avoid -inf - -inf in the subtraction
    a_max = np.where(np.isfinite(a_max), a_max, 0.0)
    out = np.log(np.sum(np.exp(a - a_max), axis=axis))
    return np.squeeze(a_max, axis=axis) + out


def _gaussian_emission_logprob(
    observations: np.ndarray,
    state_means: np.ndarray,
    state_sd: float,
) -> np.ndarray:
    """Gaussian log-likelihood with a shared per-state standard deviation.

    Use this for continuous-valued emissions like ``meth_diff`` where
    the underlying distribution isn't a Bernoulli. NaN observations
    fall back to uniform log-prob across states (no bias).
    """
    obs = np.asarray(observations, dtype=np.float64)
    sm = np.asarray(state_means, dtype=np.float64)
    sd = max(float(state_sd), 1e-6)
    # -0.5 * ((y - mu)/sd)^2  - 0.5*log(2pi) - log(sd)
    diff = (obs[:, None] - sm[None, :]) / sd
    logp = -0.5 * diff * diff - 0.5 * np.log(2.0 * np.pi) - np.log(sd)
    nan_mask = ~np.isfinite(obs)
    if nan_mask.any():
        logp[nan_mask, :] = -np.log(len(sm))
    return logp


def _bernoulli_emission_logprob(
    observations: np.ndarray,
    state_means: np.ndarray,
    state_clip: float = 1e-6,
) -> np.ndarray:
    """Bernoulli log-likelihood matrix.

    observations : (n_sites,) float in [0, 1]   (beta, possibly NaN)
    state_means  : (n_states,) float in (0, 1)  (one beta per state)
    Returns (n_sites, n_states) log-likelihoods. Sites with NaN
    observations get equal log-prob across states (= log(1/n_states))
    so they don't bias the path.
    """
    obs = np.clip(observations, 0.0, 1.0)
    sm = np.clip(state_means, state_clip, 1.0 - state_clip)
    # log(p^y * (1-p)^(1-y))
    logp = (
        obs[:, None] * np.log(sm[None, :])
        + (1.0 - obs)[:, None] * np.log(1.0 - sm[None, :])
    )
    # NaN observations -> uniform.
    nan_mask = ~np.isfinite(observations)
    if nan_mask.any():
        logp[nan_mask, :] = -np.log(len(state_means))
    return logp


def _build_transition(
    n_states: int,
    transition_priors: Optional[np.ndarray] = None,
    self_loop: float = 0.99,
) -> np.ndarray:
    """Return an (n_states, n_states) row-stochastic transition matrix.

    Defaults to a sticky chain: each state has ``self_loop`` self-prob
    and (1 - self_loop) / (n_states - 1) on every other state. Override
    by passing a fully-specified ``transition_priors`` matrix.
    """
    if transition_priors is not None:
        priors = np.asarray(transition_priors, dtype=np.float64)
        if priors.shape != (n_states, n_states):
            raise ValueError(
                f"transition_priors must be shape ({n_states}, {n_states}); "
                f"got {priors.shape}"
            )
        # Row-normalise defensively.
        return priors / priors.sum(axis=1, keepdims=True)

    A = np.full((n_states, n_states), (1.0 - self_loop) / max(n_states - 1, 1))
    np.fill_diagonal(A, self_loop)
    return A


def segment(
    observations: np.ndarray,
    n_states: int = 2,
    *,
    state_means: Optional[np.ndarray] = None,
    transition_priors: Optional[np.ndarray] = None,
    self_loop: float = 0.99,
    return_posteriors: bool = False,
    emission: str = "bernoulli",
    emission_sd: float = 0.1,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Forward-backward + Viterbi segmentation of a per-site beta signal.

    Parameters
    ----------
    observations
        Per-site beta in ``[0, 1]``. NaNs are allowed (treated as uniform
        emission).
    n_states
        Number of HMM states. ``state_means`` (if given) must have the
        same length.
    state_means
        Per-state Bernoulli means. Default: evenly spaced from 0.05 to
        0.95.
    transition_priors
        Full (n_states, n_states) row-stochastic matrix. Default: a
        sticky chain controlled by ``self_loop``.
    self_loop
        Diagonal of the default sticky-chain transition matrix.
    return_posteriors
        When True, return ``(viterbi, posteriors)`` where
        ``posteriors`` is the (n_sites, n_states) posterior probability
        matrix from forward-backward.

    Returns
    -------
    np.ndarray
        Integer state labels (length n_sites) -- the Viterbi MAP path.
    """
    obs = np.asarray(observations, dtype=np.float64)
    n_sites = len(obs)
    if n_sites == 0:
        empty = np.zeros(0, dtype=np.int32)
        if return_posteriors:
            return empty, np.zeros((0, n_states))
        return empty

    if state_means is None:
        # Evenly spaced means, avoiding 0 and 1.
        state_means = np.linspace(0.05, 0.95, n_states)
    state_means = np.asarray(state_means, dtype=np.float64)
    if len(state_means) != n_states:
        raise ValueError(
            f"state_means must have length n_states={n_states}; "
            f"got {len(state_means)}"
        )

    A = _build_transition(n_states, transition_priors, self_loop)
    logA = np.log(np.maximum(A, 1e-300))
    log_pi = np.full(n_states, -np.log(n_states))  # uniform start
    if emission == "bernoulli":
        logB = _bernoulli_emission_logprob(obs, state_means)
    elif emission == "gaussian":
        logB = _gaussian_emission_logprob(obs, state_means, state_sd=emission_sd)
    else:
        raise ValueError(
            f"emission must be 'bernoulli' (default) or 'gaussian'; got {emission!r}"
        )

    # ---- Forward (log space) ----
    log_alpha = np.full((n_sites, n_states), _LOG_EPS, dtype=np.float64)
    log_alpha[0] = log_pi + logB[0]
    for t in range(1, n_sites):
        # log_alpha[t, j] = logB[t,j] + logsumexp_i(log_alpha[t-1,i] + logA[i,j])
        log_alpha[t] = logB[t] + _logsumexp(
            log_alpha[t - 1, :, None] + logA, axis=0,
        )

    # ---- Backward (log space) ----
    log_beta = np.full((n_sites, n_states), _LOG_EPS, dtype=np.float64)
    log_beta[-1] = 0.0
    for t in range(n_sites - 2, -1, -1):
        log_beta[t] = _logsumexp(
            logA + logB[t + 1, None, :] + log_beta[t + 1, None, :], axis=1,
        )

    log_post = log_alpha + log_beta
    log_post -= _logsumexp(log_post, axis=1, keepdims=True) if False else \
                _logsumexp(log_post, axis=1)[:, None]
    posteriors = np.exp(log_post)

    # ---- Viterbi ----
    log_delta = np.full((n_sites, n_states), _LOG_EPS, dtype=np.float64)
    psi = np.zeros((n_sites, n_states), dtype=np.int32)
    log_delta[0] = log_pi + logB[0]
    for t in range(1, n_sites):
        scores = log_delta[t - 1, :, None] + logA  # (i, j)
        psi[t] = np.argmax(scores, axis=0)
        log_delta[t] = logB[t] + scores[psi[t], np.arange(n_states)]

    viterbi = np.empty(n_sites, dtype=np.int32)
    viterbi[-1] = int(np.argmax(log_delta[-1]))
    for t in range(n_sites - 2, -1, -1):
        viterbi[t] = psi[t + 1, viterbi[t + 1]]

    if return_posteriors:
        return viterbi, posteriors
    return viterbi


def runs_of_state(
    viterbi: np.ndarray,
    target_state: int,
    positions: Optional[np.ndarray] = None,
) -> list[tuple[int, int, int]]:
    """Extract contiguous runs of ``target_state`` from a Viterbi path.

    Returns a list of ``(start_idx, end_idx_exclusive, length)`` tuples.
    When ``positions`` is supplied, the indices are translated to
    positions and the returned tuples are ``(start_pos, end_pos, length_in_sites)``
    instead.
    """
    if len(viterbi) == 0:
        return []
    runs: list[tuple[int, int, int]] = []
    in_run = False
    start = 0
    for i, s in enumerate(viterbi):
        if s == target_state and not in_run:
            in_run = True
            start = i
        elif s != target_state and in_run:
            runs.append((start, i, i - start))
            in_run = False
    if in_run:
        runs.append((start, len(viterbi), len(viterbi) - start))

    if positions is not None:
        translated: list[tuple[int, int, int]] = []
        for s, e, l in runs:
            # +1 on end so range is inclusive of the last site
            translated.append((int(positions[s]), int(positions[e - 1]) + 1, l))
        return translated
    return runs


__all__ = ["segment", "runs_of_state"]
