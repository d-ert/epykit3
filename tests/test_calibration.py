"""End-to-end calibration scaffolding test (CI-gated).

This test verifies that the null-calibration pipeline correctly
classifies a calibrated test (uniform p-values) and a non-calibrated
test (Beta(2, 1) p-values, stochastically larger than uniform). It
does NOT run the full ``ep.tl.dmc`` engine on synthetic data -- that
is the job of the slow benchmark tests in
``benchmark/scripts/tests/test_null_engines.py``. The point here is
to catch regressions in the calibration *scaffolding*: the
``(pvalues, qvalues)`` contract from ``_null_engines.py``, the K-S
test wiring, and the median + IQR aggregator.

The real Linux re-run (Track 1 Layer B) drives ``ep.tl.dmc`` through
the same scaffolding on a 10k+ CpG slice of GSE263850 and writes the
calibration figures the methods appendix cites. Both layers share
this scaffolding -- a bug here breaks both.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest

# The calibration scripts live under benchmark/scripts. Add to path
# explicitly so this test is invokable from the project root.
_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "benchmark" / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))


def _uniform_engine(samples_treatment, samples_control, seed=0, n_sites=20_000):
    """Calibrated test stand-in: p-values truly uniform on [0, 1]."""
    rng = np.random.default_rng(seed)
    p = rng.uniform(0, 1, size=n_sites)
    # q-values from a simple BH on uniform p's are also approximately
    # uniform; for the calibration scaffolding test we just pass the
    # p-values through as q-values (no real engine here).
    return p, p


def _conservative_engine(samples_treatment, samples_control, seed=0, n_sites=20_000):
    """Conservative test stand-in: p-values from Beta(2, 1) -- CDF
    F(x) = x^2 -- stochastically larger than uniform."""
    rng = np.random.default_rng(seed)
    p = rng.beta(2.0, 1.0, size=n_sites)
    return p, p


@pytest.mark.slow
def test_uniform_engine_passes_ks_calibration():
    """A truly uniform-p engine: K-S vs Uniform(0, 1) should not reject
    even at moderate alpha. This pins the "calibrated" path of the
    calibration scaffolding -- if a refactor breaks the KS wiring or
    the p-value sampling, this test fails."""
    from run_null_calibration import (
        ks_against_uniform,
        run_null_calibration,
    )

    samples = [f"s{i}" for i in range(1, 7)]
    per_shuffle, pvalue_samples = run_null_calibration(
        engine_fn=_uniform_engine,
        engine_name="uniform_mock",
        scenario_name="cov10_3v3",
        samples=samples,
        n_per_group=3,
        k_shuffles=50,
        q_thresh=0.05,
        seed=42,
        qq_shuffles=5,
        qq_samples_per_shuffle=10_000,
    )
    assert per_shuffle.height == 50
    assert not pvalue_samples.is_empty()

    ks = ks_against_uniform(pvalue_samples["pvalue"].to_numpy())
    # n >> 1e3, calibrated test: K-S should NOT reject at alpha = 0.01.
    assert ks["pvalue"] > 0.01, (
        f"calibrated mock engine rejected by K-S (p={ks['pvalue']:.4g}); "
        f"calibration scaffolding may be miswired"
    )


@pytest.mark.slow
def test_conservative_engine_fails_ks_calibration():
    """A conservative test (Beta(2, 1) p-values) must be detected as
    non-calibrated by the same scaffolding. If THIS test passes, the
    scaffolding has lost its discriminative power and the
    ``test_uniform_engine_passes_ks_calibration`` positive result above
    is meaningless."""
    from run_null_calibration import (
        ks_against_uniform,
        run_null_calibration,
    )

    samples = [f"s{i}" for i in range(1, 7)]
    _, pvalue_samples = run_null_calibration(
        engine_fn=_conservative_engine,
        engine_name="conservative_mock",
        scenario_name="cov10_3v3",
        samples=samples,
        n_per_group=3,
        k_shuffles=50,
        q_thresh=0.05,
        seed=42,
        qq_shuffles=5,
        qq_samples_per_shuffle=10_000,
    )
    ks = ks_against_uniform(pvalue_samples["pvalue"].to_numpy())
    # Beta(2, 1) at n=50_000 is detected with overwhelming power.
    assert ks["pvalue"] < 1e-6, (
        f"conservative mock engine not detected by K-S "
        f"(p={ks['pvalue']:.4g}); discriminative power lost"
    )


@pytest.mark.slow
def test_summarize_observed_fdr_within_expected_ci():
    """For a calibrated test at q < 0.05, median observed FDR over
    enough shuffles should be close to 0.05. This is the headline
    calibration claim the methods appendix makes -- pin it here."""
    from run_null_calibration import (
        run_null_calibration,
        summarize_observed_fdr,
    )

    samples = [f"s{i}" for i in range(1, 7)]
    per_shuffle, _ = run_null_calibration(
        engine_fn=_uniform_engine,
        engine_name="uniform_mock",
        scenario_name="cov10_3v3",
        samples=samples,
        n_per_group=3,
        k_shuffles=200,
        q_thresh=0.05,
        seed=0,
    )
    summary = summarize_observed_fdr(per_shuffle)
    median = summary["median_observed_fdr"][0]
    # For truly uniform q's the empirical fraction below 0.05 is ~ 0.05
    # by definition; allow generous slack so this test does not flake.
    assert 0.04 < median < 0.06, (
        f"calibrated mock median observed FDR = {median:.4f}; expected "
        f"close to 0.05 (slack 0.04-0.06)"
    )
