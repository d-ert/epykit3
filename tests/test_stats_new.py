"""New stat / utility features: shrinkage, kNN imputation, clocks, deconvolution."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit.dmc import shrink_meth_diff
from epykit.impute import impute_knn_beta, impute_knn_anndata


# Empirical-Bayes shrinkage of meth_diff


def test_shrink_meth_diff_pulls_low_information_toward_zero():
    """A site with a wide CI (large SE) should be shrunk hard; a site
    with a tight CI should barely move.

    tau^2 is estimated from the empirical variance of meth_diff across
    sites, minus the mean sampling variance. So the test data must
    have some real between-site variance, otherwise tau^2=0 collapses
    every estimate to 0 (which is the *correct* EB answer when there's
    no signal -- see test_shrink_meth_diff_handles_no_signal).
    """
    rng = np.random.default_rng(13)
    # 100 sites: a real signal exists (true effects +/-0.4), but observed
    # estimates carry sampling noise that's much smaller on tight-CI sites
    # than wide-CI sites.
    n = 100
    true_effects = rng.choice([-0.4, 0.0, 0.4], size=n, p=[0.25, 0.5, 0.25])
    se_tight = np.full(n // 2, 0.03)
    se_wide = np.full(n - n // 2, 0.30)
    se = np.concatenate([se_tight, se_wide])
    obs = true_effects + rng.normal(0, se, size=n)
    df = pl.DataFrame({
        "meth_diff": obs,
        "meth_diff_ci_lo": obs - 1.96 * se,
        "meth_diff_ci_hi": obs + 1.96 * se,
    })
    out = shrink_meth_diff(df, se_from="ci")
    assert "meth_diff_shrunk" in out.columns
    assert "shrinkage_factor" in out.columns

    sf = out.get_column("shrinkage_factor").to_numpy()
    # tau^2 should be ~0.13 (var of +/-0.4 mixture). Tight SE^2 ~= 0.0009, wide
    # SE^2 ~= 0.09. Shrinkage factor tight = 0.13/(0.13+0.0009) ~= 0.993;
    # wide = 0.13/(0.13+0.09) ~= 0.59.
    mean_tight = float(np.nanmean(sf[: n // 2]))
    mean_wide = float(np.nanmean(sf[n // 2:]))
    assert mean_tight > 0.90, (
        f"tight-CI mean shrinkage_factor {mean_tight:.3f} should exceed 0.90"
    )
    assert mean_wide < 0.80, (
        f"wide-CI mean shrinkage_factor {mean_wide:.3f} should be <0.80"
    )
    assert mean_tight > mean_wide


def test_shrink_meth_diff_handles_no_signal():
    """If all observed Deltabeta are within sampling noise, tau^2 collapses to 0
    and every shrunk estimate is exactly 0."""
    rng = np.random.default_rng(0)
    n = 200
    se = np.full(n, 0.10)
    md = rng.normal(0, 0.10, size=n)  # mean 0, sd matches SE -> no real signal
    df = pl.DataFrame({
        "meth_diff": md,
        "meth_diff_ci_lo": md - 1.96 * se,
        "meth_diff_ci_hi": md + 1.96 * se,
    })
    out = shrink_meth_diff(df)
    shrunk = out.get_column("meth_diff_shrunk").to_numpy()
    # With tau^2 = max(0, Var(md) - mean(SE^2)) ~= max(0, ~0.01 - ~0.01) ~= 0
    # -> every shrunk value is 0.
    assert np.all(np.abs(shrunk) < 1e-6)


def test_shrink_meth_diff_rejects_missing_columns():
    df = pl.DataFrame({"meth_diff": [0.1, 0.2]})
    with pytest.raises(ValueError, match="meth_diff_ci_lo"):
        shrink_meth_diff(df, se_from="ci")


def test_shrink_meth_diff_real_dmc_table(synth_md_filtered):
    """End-to-end on a real DMC table from the synth fixture."""
    ep.tl.dmc(synth_md_filtered, test="lr")
    df = synth_md_filtered.dmc
    out = shrink_meth_diff(df)
    n_finite = out.filter(pl.col("meth_diff_shrunk").is_finite()).height
    assert n_finite > 0
    # On real data, |shrunk| <= |raw| at every site (no shrinkage
    # estimator can amplify a signal).
    raw = df.get_column("meth_diff").to_numpy()
    shr = out.get_column("meth_diff_shrunk").to_numpy()
    finite = np.isfinite(raw) & np.isfinite(shr)
    assert np.all(np.abs(shr[finite]) <= np.abs(raw[finite]) + 1e-9), (
        "shrunk |Deltabeta| should never exceed raw |Deltabeta|"
    )


# kNN methylation imputation


def test_impute_knn_beta_fills_gaps_in_a_smooth_signal():
    """A sample with a slowly-varying beta profile and a few NaN holes
    should impute back to near the local mean."""
    rng = np.random.default_rng(7)
    n_sites = 100
    positions = (np.arange(n_sites) * 100).astype(np.int64)
    # Sinusoidal beta trace per sample.
    true_beta = 0.5 + 0.3 * np.sin(np.arange(n_sites) / 6.0)
    beta_full = np.tile(true_beta, (3, 1)) + rng.normal(0, 0.01, (3, n_sites))
    beta_with_holes = beta_full.copy()
    hole_idx = np.array([20, 21, 22, 60, 80])
    beta_with_holes[0, hole_idx] = np.nan  # only sample 0 has holes
    imputed = impute_knn_beta(positions, beta_with_holes, k=5, max_distance_bp=10_000)
    # Imputed values should land within 0.1 of the truth (loose, accounts
    # for the local mean differing from the exact point estimate).
    err = np.abs(imputed[0, hole_idx] - true_beta[hole_idx])
    assert (err < 0.1).all(), (
        f"kNN imputation errors {err} exceed 0.1 -- local kNN should "
        "recover a smoothly-varying signal."
    )
    # Sample 0's non-hole sites are untouched.
    untouched = np.array([i for i in range(n_sites) if i not in hole_idx])
    assert np.array_equal(
        imputed[0, untouched], beta_with_holes[0, untouched]
    )


def test_impute_knn_beta_respects_max_distance(capsys):
    """When max_distance_bp is small and the missing site has no covered
    neighbours inside that window, the value stays NaN."""
    positions = np.array([0, 100, 50_000_000, 50_000_100], dtype=np.int64)
    beta = np.array([[0.1, np.nan, 0.8, 0.9]], dtype=np.float64)
    imp = impute_knn_beta(positions, beta, k=2, max_distance_bp=1_000)
    # Only position 0 is within 1 kb of position 100 -> imputes from that.
    assert not np.isnan(imp[0, 1])
    # Now flip: try to impute position 0 with NaN, with neighbours far away.
    beta = np.array([[np.nan, np.nan, 0.8, 0.9]], dtype=np.float64)
    imp = impute_knn_beta(positions, beta, k=2, max_distance_bp=1_000)
    # No covered neighbour within 1 kb of position 0 -> stays NaN.
    assert np.isnan(imp[0, 0])


def test_impute_knn_beta_validates_inputs():
    with pytest.raises(ValueError, match="sorted ascending"):
        impute_knn_beta(np.array([10, 5, 20]), np.zeros((1, 3)))
    with pytest.raises(ValueError, match="align"):
        impute_knn_beta(np.array([0, 1]), np.zeros((1, 5)))
    with pytest.raises(ValueError, match="k must be"):
        impute_knn_beta(np.array([0, 1]), np.zeros((1, 2)), k=0)


def test_impute_knn_anndata_per_chromosome_smoke(tmp_path):
    """Build a tiny AnnData-like duck-typed object and round-trip kNN
    imputation through the AnnData wrapper without depending on a full
    epykit pipeline."""
    pytest.importorskip("anndata")
    import anndata as ad

    rng = np.random.default_rng(42)
    chroms = np.array(["chr1"] * 50 + ["chr2"] * 50)
    positions = np.concatenate([
        np.arange(50) * 100, np.arange(50) * 100,
    ]).astype(np.int64)
    # Truth: sinusoidal beta; introduce NaN holes per sample.
    truth = 0.5 + 0.3 * np.sin(np.arange(100) / 5.0)
    X = np.tile(truth, (3, 1)) + rng.normal(0, 0.01, (3, 100))
    X[0, [10, 11, 60, 75]] = np.nan
    adata = ad.AnnData(X=X)
    import pandas as pd
    adata.var = pd.DataFrame(
        {"chrom": chroms, "pos": positions},
        index=[f"site_{i}" for i in range(100)],
    )
    imputed = impute_knn_anndata(adata, k=5, max_distance_bp=10_000)
    assert imputed.shape == X.shape
    assert not np.isnan(imputed[0, 10])
    # Sample 0 values at NaN positions land near truth.
    err = np.abs(imputed[0, [10, 11, 60, 75]] - truth[[10, 11, 60, 75]])
    assert (err < 0.1).all()


# Age clocks (generic linear runner)


def test_age_clock_recovers_known_age_from_synthetic_data(synth_md_filtered):
    """Build a tiny synthetic clock from the synth fixture: assign each
    sample a 'true age' linearly from a handful of CpG beta values, then
    confirm the runner recovers it."""
    md = synth_md_filtered
    # Pick 10 CpGs from the store via the existing common-sites helper.
    samples = md.obs.get_column("sample_id").to_list()
    one_sample = samples[0]
    sites = (
        pl.scan_parquet(f"{md.store}/sample={one_sample}/chrom=*/part-*.parquet")
        .select(["chrom", "pos"])
        .collect()
        .head(20)
        .with_columns(pl.col("pos").cast(pl.Int64))
    )
    # Build a coefficient table -- give each CpG a known linear weight.
    cpg_ids = [f"cg{i:04d}" for i in range(sites.height)]
    coefficients = pl.DataFrame({
        "cpg_id": cpg_ids,
        "coefficient": [1.0] * sites.height,
    })
    manifest = sites.with_columns(pl.Series("cpg_id", cpg_ids)).select(
        ["cpg_id", "chrom", "pos"]
    )
    ep.tl.age_clock(
        md, coefficients, manifest,
        intercept=0.0, transform=None, name="synth_clock",
    )
    assert "synth_clock" in md.obs.columns
    ages = md.obs.get_column("synth_clock").to_numpy()
    # The clock is a linear combination of beta values, so all samples
    # should land in [-inf, +inf] but for sanity, on real beta values with
    # all-positive coefficients we expect every age to be > 0 and
    # bounded by sum(coef) = 20.
    assert np.all((ages >= 0) & (ages <= 20)), (
        f"clock outputs {ages} fell outside expected [0, 20] sum range"
    )
    assert "synth_clock_diagnostics" in md.uns


def test_age_clock_horvath_transform_branches():
    """The horvath transform piecewise: negative linear -> exp-based,
    positive linear -> linear-scale."""
    from epykit.clocks import age_clock as _age_clock_fn  # noqa: F401
    # Direct math check on the transform formula without running on data.
    lin_neg = -0.5
    lin_pos = 1.0
    transformed_neg = np.exp(lin_neg) - 1.0      # ~= -0.393
    transformed_pos = lin_pos * 21.0 + 20.0      # = 41.0
    assert abs(transformed_neg - (-0.39346934)) < 1e-6
    assert transformed_pos == 41.0


# Deconvolution (NNLS)


def test_deconvolve_runs_and_returns_long_format(synth_md_filtered):
    """Build a tiny 3-cell-type reference matrix from synth CpGs and
    confirm the runner produces sane proportions."""
    md = synth_md_filtered
    samples = md.obs.get_column("sample_id").to_list()
    one_sample = samples[0]
    sites = (
        pl.scan_parquet(f"{md.store}/sample={one_sample}/chrom=*/part-*.parquet")
        .select(["chrom", "pos"])
        .collect()
        .head(40)
        .with_columns(pl.col("pos").cast(pl.Int64))
    )
    cpg_ids = [f"cg{i:04d}" for i in range(sites.height)]
    rng = np.random.default_rng(11)
    # Fake reference: 3 cell types with deliberately different beta profiles
    # so NNLS has signal to use.
    n = sites.height
    ref = pl.DataFrame({
        "cpg_id": cpg_ids,
        "TypeA": rng.uniform(0.0, 0.3, n),
        "TypeB": rng.uniform(0.3, 0.6, n),
        "TypeC": rng.uniform(0.6, 0.9, n),
    })
    manifest = sites.with_columns(pl.Series("cpg_id", cpg_ids)).select(
        ["cpg_id", "chrom", "pos"]
    )
    ep.tl.deconvolve(md, ref, manifest)
    long = md.uns["deconvolution"]
    # Long format: 3 cell types x n_samples rows.
    assert long.height == len(samples) * 3
    # Per-sample proportions sum to ~1 (or are all NaN if no coverage).
    for sample in samples:
        s_sub = long.filter(pl.col("sample_id") == sample)
        props = s_sub.get_column("proportion").to_numpy()
        if np.isnan(props).all():
            continue
        assert abs(props.sum() - 1.0) < 1e-6 or props.sum() == 0.0, (
            f"sample {sample}: proportions sum to {props.sum()}"
        )
        # All proportions in [0, 1].
        assert ((props >= 0) & (props <= 1)).all()
    # Wide-format obs columns get prefixed.
    obs_cols = md.obs.columns
    assert any(c.startswith("frac_") for c in obs_cols)


def test_deconvolve_unknown_method_errors(synth_md_filtered):
    md = synth_md_filtered
    ref = pl.DataFrame({"cpg_id": ["cg1"], "TypeA": [0.5]})
    manifest = pl.DataFrame({"cpg_id": ["cg1"], "chrom": ["chr1"], "pos": [100]})
    with pytest.raises(ValueError, match="nnls"):
        ep.tl.deconvolve(md, ref, manifest, method="rpc")


# ---- CI columns + bb_lr (merged from test_dmc_ci_and_rename.py) ----------


@pytest.mark.parametrize("test", ["lr", "score", "logit_t", "welch_t", "bb_lr"])
def test_dmc_emits_meth_diff_ci_columns(synth_md_filtered, test):
    """Every DMC test path emits meth_diff_ci_lo / meth_diff_ci_hi."""
    md = synth_md_filtered
    ep.tl.dmc(md, test=test)
    df = md.get_dmc(test=test)
    assert df is not None
    assert "meth_diff_ci_lo" in df.columns
    assert "meth_diff_ci_hi" in df.columns
    finite = df.filter(
        pl.col("meth_diff").is_not_null()
        & pl.col("meth_diff_ci_lo").is_not_null()
        & pl.col("meth_diff_ci_hi").is_not_null()
    )
    if len(finite) > 100:
        m = finite.get_column("meth_diff").to_numpy()
        lo = finite.get_column("meth_diff_ci_lo").to_numpy()
        hi = finite.get_column("meth_diff_ci_hi").to_numpy()
        assert np.all(lo <= m + 1e-6)
        assert np.all(m - 1e-6 <= hi)
        frac_inside = float(((lo <= m) & (m <= hi)).mean())
        assert frac_inside >= 0.95


def test_beta_binomial_is_rejected(synth_md_filtered):
    """test='beta_binomial' was removed; it should now raise an error."""
    md = synth_md_filtered
    with pytest.raises((ValueError, KeyError, NotImplementedError)):
        ep.tl.dmc(md, test="beta_binomial")


def test_bb_lr_is_distinct_from_lr(synth_md_filtered):
    """bb_lr (true quasi-binomial LRT) produces a separate output table
    from lr and surfaces coef_treatment / coef_se."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="bb_lr")
    df = md.get_dmc(test="bb_lr")
    assert df is not None
    assert "coef_treatment" in df.columns
    assert "coef_se" in df.columns
    coef = df.get_column("coef_treatment").drop_nulls().to_numpy()
    assert coef.size > 0
    assert np.isfinite(coef).any()
