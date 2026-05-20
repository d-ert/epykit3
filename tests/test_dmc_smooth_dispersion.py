"""Tests for DSS-style per-sample count smoothing in the quasi-binomial LR test.

Replicates DSS::DMLfit.multiFactor(smoothing=TRUE): each sample's raw
counts (meth, cov) are replaced by a uniform-box moving average over a
+/-smoothing_span_bp/2 neighborhood before they hit the score
accumulators. Dispersion is per-CpG (unchanged), matching DSS exactly.

The previous attempt smoothed the *dispersion* across positions, which
DSS does not do and which empirically killed power (89 sig CpGs vs
~22 K). This file covers the corrected primitive.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit.dmc import (
    _dmc_input_signature,
    _resolve_dmc_store_dir,
    _smooth_sample_counts_box,
)

pytestmark = pytest.mark.slow


# 1. Unit tests for `_smooth_sample_counts_box`


def test_box_smooth_constant_input_returns_constant():
    """Constant counts -> smoothed counts equal that constant."""
    n = 50
    positions = np.arange(0, n * 10, 10, dtype=np.int64)
    meth = np.full(n, 5, dtype=np.int32)
    cov  = np.full(n, 20, dtype=np.int32)
    meth_sm, cov_sm = _smooth_sample_counts_box(meth, cov, positions, window_bp=200)
    np.testing.assert_allclose(meth_sm, 5.0, atol=1e-9)
    np.testing.assert_allclose(cov_sm,  20.0, atol=1e-9)


def test_box_smooth_single_cpg():
    """A single CpG smoothed against itself returns its own count."""
    positions = np.array([1000], dtype=np.int64)
    meth = np.array([3], dtype=np.int32)
    cov  = np.array([10], dtype=np.int32)
    meth_sm, cov_sm = _smooth_sample_counts_box(meth, cov, positions, window_bp=500)
    assert meth_sm[0] == pytest.approx(3.0)
    assert cov_sm[0]  == pytest.approx(10.0)


def test_box_smooth_isolated_cpg_unchanged():
    """A CpG with no neighbors in window keeps its own count."""
    positions = np.array([100, 10_000_000], dtype=np.int64)  # 10 Mb apart
    meth = np.array([3, 7], dtype=np.int32)
    cov  = np.array([10, 20], dtype=np.int32)
    meth_sm, cov_sm = _smooth_sample_counts_box(meth, cov, positions, window_bp=500)
    np.testing.assert_array_equal(meth_sm, np.array([3.0, 7.0]))
    np.testing.assert_array_equal(cov_sm,  np.array([10.0, 20.0]))


def test_box_smooth_window_zero_is_passthrough():
    """window_bp=0 -> smoothing is a no-op (return inputs cast to float)."""
    positions = np.array([0, 100, 200], dtype=np.int64)
    meth = np.array([1, 2, 3], dtype=np.int32)
    cov  = np.array([10, 20, 30], dtype=np.int32)
    meth_sm, cov_sm = _smooth_sample_counts_box(meth, cov, positions, window_bp=0)
    np.testing.assert_array_equal(meth_sm, np.array([1.0, 2.0, 3.0]))
    np.testing.assert_array_equal(cov_sm,  np.array([10.0, 20.0, 30.0]))


def test_box_smooth_matches_brute_force_on_random_fixture():
    """Smoothed value equals (sum of x[j] in window) / window_size -- brute force.

    Mirrors DSS's smooth.chr(..., method="avg") + nitem_bin geometry:
    for each CpG i, window is the set of j with |pos_j - pos_i| <= half.
    """
    rng = np.random.default_rng(7)
    n = 80
    # Random-ish spaced positions on chr1.
    spacings = rng.integers(20, 300, size=n).astype(np.int64)
    positions = np.cumsum(spacings)
    meth = rng.integers(0, 30, size=n).astype(np.int32)
    cov  = rng.integers(5, 60, size=n).astype(np.int32)
    window_bp = 500
    half = window_bp // 2

    meth_sm, cov_sm = _smooth_sample_counts_box(meth, cov, positions, window_bp)

    # Brute force: O(n^2) reference implementation.
    for i in range(n):
        mask = np.abs(positions - positions[i]) <= half
        expected_meth = meth[mask].mean()
        expected_cov  = cov[mask].mean()
        assert meth_sm[i] == pytest.approx(expected_meth, rel=1e-9, abs=1e-9), (
            f"site {i} meth: got {meth_sm[i]}, expected {expected_meth}"
        )
        assert cov_sm[i] == pytest.approx(expected_cov, rel=1e-9, abs=1e-9), (
            f"site {i} cov: got {cov_sm[i]}, expected {expected_cov}"
        )


def test_box_smooth_first_and_last_use_truncated_window():
    """At chromosome ends, the window is truncated to in-range neighbors."""
    positions = np.array([0, 100, 200, 300, 400], dtype=np.int64)
    meth = np.array([1, 1, 1, 1, 1], dtype=np.int32)
    cov  = np.array([10, 20, 30, 40, 50], dtype=np.int32)
    # window=200 -> half=100; site 0 sees self + site 100 (gap 100 <= 100).
    meth_sm, cov_sm = _smooth_sample_counts_box(meth, cov, positions, window_bp=200)
    # Site 0 window: {0, 100} -> cov avg = 15. Site 4 window: {300, 400} -> cov avg = 45.
    assert cov_sm[0] == pytest.approx(15.0)
    assert cov_sm[4] == pytest.approx(45.0)
    # Middle site 2 (pos=200): window = {100, 200, 300} -> cov avg = 30.
    assert cov_sm[2] == pytest.approx(30.0)


def test_box_smooth_empty_input():
    """Length-0 inputs return empty arrays without crashing."""
    meth_sm, cov_sm = _smooth_sample_counts_box(
        np.array([], dtype=np.int32),
        np.array([], dtype=np.int32),
        np.array([], dtype=np.int64),
        window_bp=500,
    )
    assert meth_sm.shape == (0,)
    assert cov_sm.shape  == (0,)


# 2. Cache routing


def test_smoothing_signature_changes_with_span(tmp_path):
    """Different smoothing_span_bp -> different cache signatures."""
    common = dict(
        methylstore_path=tmp_path,
        samples_case=["A", "B"],
        samples_control=["C", "D"],
        test="lr",
        chromosomes=["chr1"],
        unite=True,
        min_samples_case=0,
        min_samples_control=0,
        dispersion="site",
        reference="adaptive",
        samples_all_ordered=None,
        group_labels_per_sample=None,
        contrast_label=None,
    )
    sig_500  = _dmc_input_signature(**common, smoothing=True, smoothing_span_bp=500)
    sig_1000 = _dmc_input_signature(**common, smoothing=True, smoothing_span_bp=1000)
    sig_off  = _dmc_input_signature(**common, smoothing=False, smoothing_span_bp=500)
    assert sig_500 != sig_1000
    # Flipping the smoothing bit invalidates the cache even at the same span.
    assert sig_500 != sig_off


def test_smoothing_signature_ignores_span_when_smoothing_off(tmp_path):
    """When smoothing=False, smoothing_span_bp doesn't enter the hash.

    Otherwise touching the default would invalidate every existing
    un-smoothed cache.
    """
    common = dict(
        methylstore_path=tmp_path,
        samples_case=["A", "B"],
        samples_control=["C", "D"],
        test="lr",
        chromosomes=["chr1"],
        unite=True,
        min_samples_case=0,
        min_samples_control=0,
        dispersion="site",
        reference="adaptive",
        samples_all_ordered=None,
        group_labels_per_sample=None,
        contrast_label=None,
        smoothing=False,
    )
    sig_a = _dmc_input_signature(**common, smoothing_span_bp=500)
    sig_b = _dmc_input_signature(**common, smoothing_span_bp=999)
    assert sig_a == sig_b


def test_smoothing_store_dir_is_suffixed(tmp_path):
    """smoothing=True routes to a separate '<test>_smooth' bucket."""
    store = tmp_path / "store"
    store.mkdir()
    plain  = _resolve_dmc_store_dir(store, test="lr", out_dir=None, smoothing=False)
    smooth = _resolve_dmc_store_dir(store, test="lr", out_dir=None, smoothing=True)
    assert plain != smooth
    assert plain.name  == "lr"
    assert smooth.name == "lr_smooth"


# 3. End-to-end through `tl.dmc`


def test_dmc_smoothing_runs_and_records_metadata(synth_md_filtered):
    """tl.dmc(..., smoothing=True) runs end-to-end and records params."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr", smoothing=True, smoothing_span_bp=500)
    assert "dmc_lr" in md.varm
    meta = md.uns["dmc"]
    assert meta["smoothing"] is True
    assert meta["smoothing_span_bp"] == 500
    # Store path should live under a '_smooth' bucket.
    assert "lr_smooth" in str(meta["store_path"])


def test_dmc_smoothing_off_records_metadata(synth_md_filtered):
    """When smoothing=False, span isn't surfaced (set to None)."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr", smoothing=False)
    meta = md.uns["dmc"]
    assert meta["smoothing"] is False
    assert meta["smoothing_span_bp"] is None


def test_dmc_smoothing_does_not_collide_with_unsmoothed_cache(synth_md_filtered):
    """Running un-smoothed then smoothed produces two distinct caches; neither poisons the other."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr", smoothing=False)
    plain_store_path = md.uns["dmc"]["store_path"]
    plain_pvals = md.varm["dmc_lr"]["pvalue"].to_numpy()

    ep.tl.dmc(md, test="lr", smoothing=True, smoothing_span_bp=500)
    smooth_store_path = md.uns["dmc"]["store_path"]
    smooth_pvals = md.varm["dmc_lr"]["pvalue"].to_numpy()

    assert plain_store_path != smooth_store_path
    assert Path(plain_store_path).exists()
    assert Path(smooth_store_path).exists()
    # And the smoothing actually changed p-values vs the un-smoothed mode.
    diff = np.abs(plain_pvals - smooth_pvals)
    diff = diff[np.isfinite(diff)]
    assert (diff > 0).any()


def test_dmc_smoothing_output_schema_matches_unsmoothed(synth_md_filtered):
    """Smoothed DMC table has the same columns + dtypes as the un-smoothed run."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr", smoothing=False)
    plain_schema = md.varm["dmc_lr"].schema

    ep.tl.dmc(md, test="lr", smoothing=True)
    smooth_schema = md.varm["dmc_lr"].schema

    assert dict(plain_schema) == dict(smooth_schema)


def test_use_smoothed_emits_deprecation_warning(synth_md_filtered):
    """The legacy pseudo-count path now emits a DeprecationWarning pointing at smoothing=True."""
    md = synth_md_filtered
    ep.pp.smooth(md, method="bsmooth", ns=20, h_bp=2000)
    with pytest.warns(DeprecationWarning, match="smoothing=True"):
        ep.tl.dmc(md, test="lr", use_smoothed=True)
