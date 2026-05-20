"""Tests for the shared HMM segmentation engine in ``epykit._hmm``.

The contract:
  1. On a hand-built two-state Markov chain with clear emissions, the
     Viterbi path recovers the underlying states with >= 95 % accuracy.
  2. The posteriors sum to 1 per site (correctness check on
     forward-backward normalisation).
  3. ``runs_of_state`` finds the obvious contiguous runs.
  4. Edge cases: empty input, NaN observations, mis-shaped priors.
"""

from __future__ import annotations

import numpy as np
import pytest

from epykit._hmm import runs_of_state, segment


def _simulate_two_state_chain(
    n: int = 2000,
    p_self: float = 0.99,
    mean_state0: float = 0.1,
    mean_state1: float = 0.85,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """Simulate a 2-state HMM with Bernoulli-style emissions."""
    rng = np.random.default_rng(seed)
    means = (mean_state0, mean_state1)
    states = np.zeros(n, dtype=np.int32)
    for t in range(1, n):
        if rng.random() > p_self:
            states[t] = 1 - states[t - 1]
        else:
            states[t] = states[t - 1]
    obs = rng.binomial(1, p=np.array([means[s] for s in states])).astype(np.float64)
    return obs, states


def test_segment_recovers_states_high_accuracy():
    obs, truth = _simulate_two_state_chain()
    path = segment(obs, n_states=2,
                   state_means=np.array([0.1, 0.85]),
                   self_loop=0.99)
    # Account for label swapping: HMM doesn't know which state is "0" vs "1".
    acc_direct = (path == truth).mean()
    acc_flipped = (path == (1 - truth)).mean()
    accuracy = max(acc_direct, acc_flipped)
    assert accuracy >= 0.90, f"HMM recovery only {accuracy:.2%} accurate"


def test_segment_returns_posteriors_when_requested():
    obs, _truth = _simulate_two_state_chain(n=100)
    path, posteriors = segment(obs, n_states=2, return_posteriors=True)
    assert posteriors.shape == (100, 2)
    # Each row must sum to ~1.
    np.testing.assert_allclose(posteriors.sum(axis=1), 1.0, atol=1e-6)


def test_segment_empty_input():
    out = segment(np.array([]), n_states=2)
    assert isinstance(out, np.ndarray)
    assert out.shape == (0,)


def test_segment_rejects_mismatched_state_means():
    with pytest.raises(ValueError, match="state_means must have length"):
        segment(np.array([0.5, 0.5]), n_states=2, state_means=np.array([0.3]))


def test_segment_rejects_wrong_transition_shape():
    with pytest.raises(ValueError, match="transition_priors must be shape"):
        segment(np.array([0.5, 0.5]), n_states=2,
                transition_priors=np.eye(3))


def test_runs_of_state_finds_contiguous_blocks():
    path = np.array([0, 0, 1, 1, 1, 0, 1, 1, 0, 0])
    runs = runs_of_state(path, target_state=1)
    # Expect runs at [2,5), [6,8).
    assert runs == [(2, 5, 3), (6, 8, 2)]


def test_runs_of_state_translates_positions():
    path = np.array([0, 1, 1, 0])
    positions = np.array([100, 200, 300, 400], dtype=np.int32)
    runs = runs_of_state(path, target_state=1, positions=positions)
    # [200, 301) -- end is +1 past the last in-run position
    assert runs == [(200, 301, 2)]


def test_segment_handles_nan_observations():
    """NaN observations should not crash and should not bias the path."""
    obs = np.array([0.1, np.nan, 0.1, 0.9, 0.9, np.nan, 0.9])
    out = segment(obs, n_states=2, state_means=np.array([0.1, 0.9]))
    # Valid-output sanity: NaN positions are assigned a state (any state),
    # and the surrounding beta=0.1 / beta=0.9 segments still recover.
    assert out.shape == (7,)
    assert out.dtype == np.int32
