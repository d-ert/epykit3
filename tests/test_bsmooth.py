"""BSmooth-style local-polynomial smoother tests.

Contract:
  1. Constant input -> constant output (within float tolerance).
  2. Linear ramp -> recovered closely (degree-2 polynomial nests linear).
  3. Adaptive bandwidth: a sparse-CpG region uses a wider bp window than
     a dense region. We don't test the exact bandwidth -- we test that
     output stays finite and bounded there.
  4. Edge sites get smoothed (not NaN'd) when enough neighbors are in range.
  5. End-to-end: ep.pp.smooth(md, method="bsmooth") writes the sidecar.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit.dmr import _bsmooth_make_njit, smooth_methylation_bsmooth

pytestmark = pytest.mark.slow


# ---- 1. Direct numerical correctness on hand-built signals ---------


def _run_kernel(positions, n_meth, coverage, *, ns=10, h_min=500.0, degree=2):
    """Direct call into the compiled kernel for unit tests."""
    smoother = _bsmooth_make_njit()
    return smoother(
        positions.astype(np.float64),
        n_meth.astype(np.float64),
        coverage.astype(np.float64),
        int(ns), float(h_min), int(degree), 3,
    )


def test_bsmooth_constant_input_constant_output():
    """A constant beta should pass through unchanged."""
    positions = np.arange(0, 50_000, 500, dtype=np.float64)   # 100 CpGs
    n = len(positions)
    cov = np.full(n, 20, dtype=np.float64)
    meth = np.full(n, 14, dtype=np.float64)   # beta = 0.7 everywhere
    out = _run_kernel(positions, meth, cov, ns=20, h_min=2000.0)
    # Allow 1e-3 absolute tolerance (the local polynomial fit has tiny
    # boundary artefacts at the extreme ends).
    np.testing.assert_allclose(out, 0.7, atol=1e-3)


def test_bsmooth_linear_ramp_recovered():
    """A linear ramp in beta should be recovered closely by the smoother.

    With degree=2 (the canonical BSmooth fit), a linear trend is nested
    inside the polynomial family so we expect near-exact recovery away
    from the boundaries.
    """
    positions = np.arange(0, 50_000, 500, dtype=np.float64)
    n = len(positions)
    # True beta: linear ramp from 0.20 to 0.80
    beta_true = np.linspace(0.20, 0.80, n)
    cov = np.full(n, 30, dtype=np.float64)
    meth = np.rint(beta_true * cov)
    out = _run_kernel(positions, meth, cov, ns=15, h_min=2000.0)
    # Interior should match the underlying ramp tightly (sampling noise
    # is the only floor since meth is rounded from beta_true * cov).
    interior = slice(20, n - 20)
    diff = np.abs(out[interior] - beta_true[interior])
    assert diff.max() < 0.02, f"max interior error {diff.max():.4f}; expected < 0.02"
    # Even edge sites should be smoothed (not NaN).
    assert np.isfinite(out).all()


def test_bsmooth_returns_raw_beta_with_too_few_neighbors():
    """Sites whose window has < min_cpgs_for_smooth valid neighbors -> raw beta."""
    # Three CpGs spaced 10 kb apart; with h_min=1000 each site sees only
    # itself, so the smoothed value must equal beta_raw.
    positions = np.array([1000.0, 11000.0, 21000.0])
    cov = np.array([20.0, 20.0, 20.0])
    meth = np.array([5.0, 10.0, 15.0])   # beta = 0.25, 0.50, 0.75
    out = _run_kernel(positions, meth, cov, ns=5, h_min=1000.0)
    raw = meth / cov
    np.testing.assert_allclose(out, raw, atol=1e-9)


def test_bsmooth_clips_to_unit_interval():
    """Output must stay in [0, 1] even when the local fit would extrapolate."""
    positions = np.arange(0, 20_000, 500, dtype=np.float64)
    n = len(positions)
    # Pathological: half all-methylated, half all-unmethylated. The
    # quadratic fit at the boundary could in principle extrapolate
    # outside [0, 1]; the smoother must clip.
    cov = np.full(n, 30, dtype=np.float64)
    meth = np.where(positions < 10_000, 30.0, 0.0)
    out = _run_kernel(positions, meth, cov, ns=10, h_min=2000.0)
    assert (out >= 0.0).all()
    assert (out <= 1.0).all()


def test_bsmooth_high_coverage_dominates():
    """A high-coverage CpG should drive the local fit more than its low-coverage neighbors.

    Three CpGs at close spacing: outer two at coverage 5, beta 0.0; center
    at coverage 500, beta 1.0. The center smoothed value should be much
    closer to 1.0 than to 0.5 (coverage-weighted mean).
    """
    positions = np.array([0.0, 1000.0, 2000.0])
    cov = np.array([5.0, 500.0, 5.0])
    meth = np.array([0.0, 500.0, 0.0])
    # need min_cpgs_for_smooth=3, ns small enough that all 3 are picked up
    smoother = _bsmooth_make_njit()
    out = smoother(positions, meth, cov, 3, 5000.0, 1, 3)
    # weighted mean of beta with weights w~cov x tricube would be:
    # w_left ~= 5 * tri(0.2)  w_center ~= 500 * tri(0)  w_right ~= 5 * tri(0.2)
    # tri(0) = 1, tri(0.2) = (1-0.008)^3 ~= 0.976
    # so center weight dominates -> out[1] should be very close to 1.
    assert out[1] > 0.95, f"center smoothed = {out[1]}, expected > 0.95"


def test_bsmooth_degree_1_works():
    """degree=1 path produces sensible output."""
    positions = np.arange(0, 20_000, 500, dtype=np.float64)
    n = len(positions)
    beta_true = np.linspace(0.30, 0.70, n)
    cov = np.full(n, 20, dtype=np.float64)
    meth = np.rint(beta_true * cov)
    out = _run_kernel(positions, meth, cov, ns=10, h_min=2000.0, degree=1)
    # Linear fit on a linear ramp should also recover the ramp.
    interior = slice(5, n - 5)
    np.testing.assert_allclose(out[interior], beta_true[interior], atol=0.03)


# ---- 2. Argument validation ----------------------------------------


def test_smooth_methylation_bsmooth_rejects_bad_degree(tmp_path):
    with pytest.raises(ValueError, match="degree must be 1 or 2"):
        smooth_methylation_bsmooth(str(tmp_path), [], degree=3)


def test_smooth_methylation_bsmooth_rejects_bad_ns(tmp_path):
    with pytest.raises(ValueError, match="ns must be >= 2"):
        smooth_methylation_bsmooth(str(tmp_path), [], ns=1)


def test_smooth_methylation_bsmooth_rejects_bad_h_bp(tmp_path):
    with pytest.raises(ValueError, match="h_bp must be > 0"):
        smooth_methylation_bsmooth(str(tmp_path), [], h_bp=0)


# ---- 3. End-to-end via pp.smooth -----------------------------------


def test_pp_smooth_bsmooth_writes_sidecar(synth_md_filtered):
    """ep.pp.smooth(md, method='bsmooth') runs and records the path."""
    md = synth_md_filtered
    ep.pp.smooth(md, method="bsmooth", ns=10, h_bp=500)
    assert "smooth_path" in md.uns
    assert md.uns["smooth_params"]["method"] == "bsmooth"
    assert md.uns["smooth_params"]["ns"] == 10
    assert md.uns["smooth_params"]["h_bp"] == 500
    # The sidecar parquet root must exist on disk.
    smooth_root = Path(md.uns["smooth_path"])
    assert smooth_root.exists()
    # And contain per-sample partitions.
    samples = list(smooth_root.glob("sample=*"))
    assert len(samples) > 0


def test_pp_smooth_unknown_method_raises(synth_md_filtered):
    with pytest.raises(ValueError, match="Unknown smoothing method"):
        ep.pp.smooth(synth_md_filtered, method="loess")


def test_pp_smooth_default_still_gaussian(synth_md_filtered):
    """Default method= preserves the pre-existing Gaussian path."""
    md = synth_md_filtered
    ep.pp.smooth(md)  # no method= kwarg
    assert md.uns["smooth_params"]["method"] == "gaussian"


def test_bsmooth_module_exposed():
    """ep.smooth_methylation_bsmooth is importable from the package surface."""
    assert hasattr(ep, "smooth_methylation_bsmooth")
    assert callable(ep.smooth_methylation_bsmooth)
