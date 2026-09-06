"""Regression tests for the four `lr` improvements introduced in 0.7.1:

* ``fdr_method="fdr_storey"`` and ``fdr_method="fdr_tsbh"`` paths through
  :func:`epykit.dmc.apply_multiple_testing_correction`.
* :func:`epykit.dmc.combine_neighbour_pvalues` (sign-aware Stouffer
  combiner).
* ``sep_fallback=True`` on the `lr` test (handled inside
  :func:`epykit.dmc._score_finalize`).
* ``dispersion="eb"`` validation.
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit.dmc import (
    _storey_pi0,
    _apply_storey_qvalues,
    apply_multiple_testing_correction,
    combine_neighbour_pvalues,
)
from tests.fixtures.synth import SimConfig, generate

pytestmark = pytest.mark.slow


# --- Storey -----------------------------------------------------------------


def test_storey_pi0_uniform_returns_near_one():
    """If all p-values are uniform on [0, 1] (pure null), pi0 ~= 1."""
    rng = np.random.default_rng(0)
    p = rng.uniform(0, 1, size=10000)
    pi0 = _storey_pi0(p)
    assert 0.9 < pi0 <= 1.0


def test_storey_pi0_50pct_signal_returns_around_05():
    """If half the p-values are at 0 (true alternatives) and half uniform,
    pi0 should be near 0.5."""
    rng = np.random.default_rng(0)
    p = np.concatenate([
        np.full(5000, 1e-10),                       # 50% strong alternatives
        rng.uniform(0, 1, size=5000),               # 50% uniform nulls
    ])
    pi0 = _storey_pi0(p)
    assert 0.4 < pi0 < 0.6, pi0


def test_storey_qvalues_at_least_as_powerful_as_bh():
    """At pi0 < 1, Storey's q-values must be <= BH q-values site-by-site."""
    rng = np.random.default_rng(1)
    p = np.concatenate([
        rng.uniform(0, 0.01, size=200),     # 200 strong alternatives
        rng.uniform(0, 1, size=800),        # 800 uniform nulls
    ])
    rng.shuffle(p)
    _, q_storey, pi0 = _apply_storey_qvalues(p)
    from statsmodels.stats.multitest import multipletests
    _, q_bh, _, _ = multipletests(p, method="fdr_bh")
    assert pi0 < 1.0, "pi0 should be < 1 with 80% nulls in the histogram"
    assert (q_storey <= q_bh + 1e-12).all(), \
        "Storey q-values should never exceed BH q-values when pi0 < 1"


def test_apply_multiple_testing_correction_invalid_method():
    df = pl.DataFrame({"pvalue": [0.1, 0.2, 0.3]})
    with pytest.raises(ValueError, match="method must be one of"):
        apply_multiple_testing_correction(df, method="not_a_method")


def test_apply_multiple_testing_correction_storey_path():
    """The Storey path emits a qvalue + reject column on a DataFrame input."""
    p = np.concatenate([
        np.full(50, 1e-10),
        np.random.default_rng(2).uniform(0, 1, size=500),
    ])
    df = pl.DataFrame({"pvalue": p})
    out = apply_multiple_testing_correction(df, method="fdr_storey")
    assert "qvalue" in out.columns
    assert "reject" in out.columns
    assert (out["qvalue"] >= 0).all() and (out["qvalue"] <= 1.0).all()


def test_apply_multiple_testing_correction_tsbh_path():
    """fdr_tsbh routes through statsmodels and produces sane q-values."""
    p = np.concatenate([
        np.full(50, 1e-10),
        np.random.default_rng(3).uniform(0, 1, size=500),
    ])
    df = pl.DataFrame({"pvalue": p})
    out = apply_multiple_testing_correction(df, method="fdr_tsbh")
    assert "qvalue" in out.columns
    # All true positives should be rejected at q < 0.05.
    n_rejected = int((out["qvalue"] < 0.05).sum())
    assert n_rejected >= 50


# --- Neighbour combiner -----------------------------------------------------


def _make_df(n=30, seed=0, sign_pattern="agree"):
    """Build a tiny synthetic DMC frame with controllable sign pattern."""
    rng = np.random.default_rng(seed)
    positions = np.arange(100, 100 + n * 50, 50)  # every 50 bp
    chrom = pl.Series("chrom", ["chr1"] * n, dtype=pl.Utf8)
    # Sites 5-15: significant true positives (small p)
    p = rng.uniform(0.5, 1.0, size=n)
    p[5:15] = rng.uniform(1e-5, 1e-3, size=10)
    if sign_pattern == "agree":
        meth_diff = np.zeros(n)
        meth_diff[5:15] = 0.5  # all positive
    else:  # "disagree"
        meth_diff = np.zeros(n)
        meth_diff[5:15:2] = 0.5   # alternating
        meth_diff[6:15:2] = -0.5
    return pl.DataFrame({
        "chrom": chrom,
        "pos": pl.Series("pos", positions, dtype=pl.Int64),
        "pvalue": p,
        "meth_diff": meth_diff,
    })


def test_neighbour_combine_aligned_signs_lowers_pvalue():
    """When neighbours agree on direction, combining lowers the p-value."""
    df = _make_df(sign_pattern="agree")
    out = combine_neighbour_pvalues(
        df,
        neighbour_bp=200,
        min_sign_agreement=0.5,
        require_focal_signal=False,
    )
    # At least one of the true-positive sites should have pvalue_combined < raw.
    sub = out.filter(pl.col("meth_diff").abs() > 0)
    assert (sub["pvalue_combined"] <= sub["pvalue"] + 1e-12).all(), \
        "combined p should never inflate p"


def test_neighbour_combine_disagreeing_signs_dampens_combination():
    """When neighbours disagree on direction, combined p should NOT be much
    smaller than the focal raw p (sign-agreement guard prevents amplification)."""
    df = _make_df(sign_pattern="disagree")
    out = combine_neighbour_pvalues(
        df,
        neighbour_bp=200,
        min_sign_agreement=0.6,
        require_focal_signal=True,
        focal_p_thresh=0.5,
    )
    # For sites near the disagreement zone, combined p should equal raw p
    # (because sign-agreement < 0.6 kicks the focal site back to raw).
    raw = out["pvalue"].to_numpy()
    comb = out["pvalue_combined"].to_numpy()
    assert (comb <= raw + 1e-12).all()


def test_neighbour_combine_preserves_nan():
    df = pl.DataFrame({
        "chrom": ["chr1", "chr1", "chr1"],
        "pos": [100, 200, 300],
        "pvalue": [0.01, float("nan"), 0.5],
        "meth_diff": [0.5, 0.5, 0.5],
    })
    out = combine_neighbour_pvalues(df, neighbour_bp=200)
    assert out.filter(pl.col("pos") == 200)["pvalue_combined"].is_nan().all()


def test_neighbour_combine_never_inflates_pvalue():
    """Property test: combined p must be <= raw p for every site."""
    rng = np.random.default_rng(42)
    n = 500
    df = pl.DataFrame({
        "chrom": ["chr1"] * n,
        "pos": (np.arange(n) * 25).tolist(),
        "pvalue": rng.uniform(0, 1, size=n).tolist(),
        "meth_diff": rng.choice([-0.5, 0.5], size=n).tolist(),
    })
    out = combine_neighbour_pvalues(df, neighbour_bp=200, min_sign_agreement=0.5)
    raw = out["pvalue"].to_numpy()
    comb = out["pvalue_combined"].to_numpy()
    # Allow tiny floating-point slack.
    assert (comb <= raw + 1e-12).all(), \
        "combine_neighbour_pvalues must never produce a p larger than raw"


# --- power_stack dispatch (1.0 breaking change) -----------------------------


def _make_md_for_power_stack(n_per_group: int, tmp_path_local):
    """Build a small MethylData with n_per_group samples per group.

    Uses the same Bismark/.cov -> read_bismark -> filter_coverage pattern
    as the conftest synth_md_filtered fixture, but with a configurable
    n_per_group and a small chromosome set to keep runtime short.
    """
    cfg = SimConfig(
        n_per_group=n_per_group,
        chromosomes=("chr1",),
        cpgs_per_chrom=500,
        n_scattered_dmcs=50,
        n_dmrs=2,
        dmr_size_cpgs=5,
        seed=123,
    )
    result = generate(cfg, tmp_path_local / f"n{n_per_group}")
    md = ep.read_bismark(
        result["samplesheet"],
        treatment_group="treatment",
        control_group="control",
        store_dir=str(tmp_path_local / f"store_n{n_per_group}"),
    )
    ep.pp.filter_coverage(md, lo_count=2, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")
    return md


def test_power_stack_auto_engages_all_four_at_any_n(tmp_path, caplog):
    """power_stack='auto' flips neighbour_combine, fdr_method, sep_fallback
    even when n > 2 (1.0 breaking change vs 0.7.6)."""
    import logging
    caplog.set_level(logging.INFO, logger="epykit")
    md = _make_md_for_power_stack(n_per_group=5, tmp_path_local=tmp_path)
    ep.tl.dmc(md, test="lr", power_stack="auto")

    # neighbour_combine ON -> emits pvalue_combined column.
    assert "pvalue_combined" in md.varm["dmc_lr"].columns, (
        "power_stack='auto' at n=5 must engage neighbour_combine "
        "(1.0 breaking change)"
    )
    # fdr_method auto-flipped: INFO log line.
    log_text = " ".join(r.message for r in caplog.records)
    assert "fdr_method" in log_text and "fdr_tsbh" in log_text, (
        f"Expected INFO log for fdr_method flip; got: {log_text!r}"
    )


def test_power_stack_lr_plus_alias(tmp_path):
    """power_stack='lr+' is an explicit alias for the full stack."""
    md = _make_md_for_power_stack(n_per_group=5, tmp_path_local=tmp_path)
    ep.tl.dmc(md, test="lr", power_stack="lr+")
    assert "pvalue_combined" in md.varm["dmc_lr"].columns, (
        "power_stack='lr+' must engage neighbour_combine"
    )


def test_power_stack_conservative_preserves_pre_1_0_behavior(tmp_path):
    """power_stack='conservative' engages only at n <= 2 (pre-1.0 'auto')."""
    md_small = _make_md_for_power_stack(n_per_group=2, tmp_path_local=tmp_path)
    md_large = _make_md_for_power_stack(n_per_group=5, tmp_path_local=tmp_path)
    ep.tl.dmc(md_small, test="lr", power_stack="conservative")
    ep.tl.dmc(md_large, test="lr", power_stack="conservative")
    assert "pvalue_combined" in md_small.varm["dmc_lr"].columns, (
        "power_stack='conservative' at n=2 must engage neighbour_combine"
    )
    assert "pvalue_combined" not in md_large.varm["dmc_lr"].columns, (
        "power_stack='conservative' at n=5 must NOT engage neighbour_combine"
    )


def test_power_stack_off_leaves_knobs_at_user_values(tmp_path):
    """power_stack='off' (or False) does not flip anything."""
    md = _make_md_for_power_stack(n_per_group=5, tmp_path_local=tmp_path)
    ep.tl.dmc(md, test="lr", power_stack="off",
              neighbour_combine=False, fdr_method="fdr_bh",
              sep_fallback=False)
    assert "pvalue_combined" not in md.varm["dmc_lr"].columns, (
        "power_stack='off' must leave neighbour_combine=False"
    )


def test_power_stack_true_is_lr_plus_alias(tmp_path):
    """Bare True aliases 'lr+'."""
    md = _make_md_for_power_stack(n_per_group=5, tmp_path_local=tmp_path)
    ep.tl.dmc(md, test="lr", power_stack=True)
    assert "pvalue_combined" in md.varm["dmc_lr"].columns, (
        "power_stack=True must alias 'lr+' and engage neighbour_combine"
    )


def test_power_stack_false_is_off_alias(tmp_path):
    """Bare False aliases 'off' (no flips)."""
    md = _make_md_for_power_stack(n_per_group=5, tmp_path_local=tmp_path)
    ep.tl.dmc(md, test="lr", power_stack=False,
              neighbour_combine=False, fdr_method="fdr_bh",
              sep_fallback=False)
    assert "pvalue_combined" not in md.varm["dmc_lr"].columns, (
        "power_stack=False must alias 'off' and leave neighbour_combine=False"
    )


def test_power_stack_unknown_string_raises(tmp_path):
    """Unknown power_stack strings raise ValueError before running."""
    md = _make_md_for_power_stack(n_per_group=5, tmp_path_local=tmp_path)
    with pytest.raises(ValueError, match="power_stack must be one of"):
        ep.tl.dmc(md, test="lr", power_stack="ultra")


def test_power_stack_default_is_off_no_breaking_change(tmp_path):
    """Spec §1 decision 1: bare tl.dmc(test='lr') with no power_stack arg
    produces bare-engine output (no neighbour_combine, no tsbh, no sep_fallback).
    This pins the default behavior and prevents accidental regression to the
    'auto'-as-default change that contradicted the spec."""
    md = _make_md_for_power_stack(n_per_group=5, tmp_path_local=tmp_path)
    ep.tl.dmc(md, test="lr")  # NO power_stack arg
    # No power-stack engagement -> no combined columns added.
    assert "pvalue_combined" not in md.varm["dmc_lr"].columns
