"""Numerical oracles for the experimental CpG-level regional backend."""

import numpy as np
import pytest
from statsmodels.regression.linear_model import GLS


def test_fixed_covariance_matches_dense_gls_with_locus_intercepts():
    from epykit._regional_gls import fit_region, transform_counts

    rng = np.random.default_rng(981)
    pos = np.array([0, 12, 35, 90, 150])
    n = rng.integers(10, 60, (5, 6))
    m = rng.binomial(n, np.array([0.2, 0.4, 0.6, 0.8, 0.3])[:, None])
    group = np.array([0, 0, 0, 1, 1, 1])
    fit = fit_region(pos, m, n, group, covariance=(120.0, 0.35))
    y = transform_counts(m, n).T.ravel()
    depth = np.median(n, axis=1)
    h = np.sqrt(depth)
    u = h / np.sqrt(np.median(depth))
    r = 0.65 * np.exp(-np.abs(pos[:, None] - pos[None, :]) / 120.0) + 0.35 * np.outer(u, u)
    sigma = np.kron(np.eye(6), r / h[:, None] / h[None, :])
    design = np.column_stack([np.tile(np.eye(5), (6, 1)), np.repeat(group, 5)])
    oracle = GLS(y, design, sigma=sigma).fit()
    assert fit["effect_transformed"] == pytest.approx(oracle.params[-1], abs=1e-10)
    assert fit["statistic"] == pytest.approx(oracle.tvalues[-1], abs=1e-9)
    assert fit["residual_df"] == oracle.df_resid


def test_donor_offset_information_saturates_with_more_cpgs():
    from epykit._regional_gls import covariance_precision

    information = []
    for k in (20, 200):
        depth = np.linspace(10, 80, k)
        h = np.sqrt(depth)
        apply, _ = covariance_precision(np.arange(k) * 20, depth, 30.0, 0.5)
        information.append(h @ apply(h))
    # Finite-K information approaches median(depth)/c = 90 rather than
    # growing tenfold with tenfold CpGs.
    assert information[1] / information[0] < 1.2
    assert max(information) < 90.0


def test_public_backend_runs_without_dmc_and_records_full_discovery(tmp_path):
    import polars as pl

    import epykit as ep

    rng = np.random.default_rng(781)
    rows = []
    for j in range(6):
        sample = f"s{j}"
        rows.append(dict(sample_id=sample, treatment=int(j >= 3)))
        part = tmp_path / "counts" / f"sample={sample}" / "chrom=chr1"
        part.mkdir(parents=True)
        m = rng.binomial(40, 0.65 if j >= 3 else 0.25, 30)
        pl.DataFrame(dict(pos=np.arange(30) * 20, N_meth=m, coverage=[40] * 30)).write_parquet(
            part / "part-0.parquet"
        )
    md = ep.MethylData(obs=pl.DataFrame(rows), store=str(tmp_path / "counts"))
    ep.tl.dmr(md, method="regional_gls", n_perm=20, min_cpgs=5, tsv=False)
    assert md.uns["dmr_params"]["null_scans"] == 18
    assert md.uns["dmr_params"]["excluded_self_or_mirror"] == 2
    assert md.uns["dmr_params"]["experimental"] is True
    assert len(md.uns["dmr_candidates"]) > 0
    assert md.uns["dmr_candidates"]["start"].min() == 0
    assert md.uns["dmr_candidates"]["end"].max() == 581
    with pytest.raises(NotImplementedError, match="exchangeable"):
        ep.tl.dmr(md, method="regional_gls", covariates=["batch"], tsv=False)
    md.store = str(tmp_path / "missing")
    with pytest.raises(ValueError, match="count store"):
        ep.tl.dmr(md, method="regional_gls", tsv=False)


def test_discovery_recovers_moderate_sites_and_never_crosses_large_gaps():
    from epykit._regional_gls import discover_regions

    pos = np.r_[np.arange(0, 120, 20), np.arange(900, 1020, 20)]
    effect = np.full(len(pos), 0.15)
    regions = discover_regions(pos, effect, window_bp=100, min_cpgs=5, max_gap=500, min_length=50)
    assert regions == [(0, 6), (6, 12)]


def test_balanced_permutations_include_observed_and_empty_scans_in_denominator():
    from epykit._regional_gls import calibrate_regions, label_assignments

    assignments = label_assignments(np.array([0, 0, 0, 1, 1, 1]), 100, 42)
    assert len(assignments) == 20
    assert len({tuple(a) for a in assignments}) == 20
    result, meta = calibrate_regions([dict(statistic=5.0)], [[5.0], [-5.0], []])
    assert result[0]["combined_pvalue"] == 1.0
    assert meta["null_scans"] == 3
    assert meta["empty_null_scans"] == 1


def test_self_mirror_pool_removal_matches_direct_conditional_calibration():
    from epykit._regional_gls import calibrate_regions, remove_self_mirror_scores

    observed = [dict(statistic=9.0), dict(statistic=-5.0), dict(statistic=2.0)]
    conditional = [[4.0, -3.0], [], [7.0, 2.0, -1.0]]
    all_scans = [[9.0, -5.0, 2.0], [-9.0, 5.0, -2.0], *conditional]
    old, meta = calibrate_regions([r.copy() for r in observed], all_scans)
    updated = remove_self_mirror_scores(old, meta["null_regions"], duplicates=2)
    direct, _ = calibrate_regions([r.copy() for r in observed], conditional)
    assert [r["combined_pvalue"] for r in updated] == pytest.approx(
        [r["combined_pvalue"] for r in direct]
    )
    assert [r["combined_qvalue"] for r in updated] == pytest.approx(
        [r["combined_qvalue"] for r in direct]
    )


def test_fast_regional_gls_execution_and_calibration(tmp_path):
    import polars as pl

    import epykit as ep

    rng = np.random.default_rng(42)
    rows = []
    for j in range(6):
        sample = f"s{j}"
        rows.append(dict(sample_id=sample, treatment=int(j >= 3)))
        part = tmp_path / "counts" / f"sample={sample}" / "chrom=chr1"
        part.mkdir(parents=True)
        # 100 CpGs, treatment difference in middle
        m = rng.binomial(40, 0.70 if (j >= 3) else 0.25, 100)
        pl.DataFrame(dict(pos=np.arange(100) * 20, N_meth=m, coverage=[40] * 100)).write_parquet(
            part / "part-0.parquet"
        )
    md = ep.MethylData(obs=pl.DataFrame(rows), store=str(tmp_path / "counts"))
    ep.tl.dmr(md, method="regional_gls", n_perm=20, min_cpgs=5, tsv=False)
    assert "dmr" in md.uns
    assert "dmr_candidates" in md.uns
    assert len(md.uns["dmr_candidates"]) > 0
    assert md.uns["dmr_params"]["null_scans"] == 18
    assert md.uns["dmr_params"]["decoy_covariance"] == "per-region-reml"
    assert md.uns["dmr_params"]["accelerator"] == "numba-native-analytical-donor-reml"
    assert md.uns["dmr_params"]["unique_null_fits"] == 9


def test_regional_cache_reuses_valid_fit_and_invalidates_changed_counts(tmp_path, monkeypatch):
    import polars as pl

    import epykit as ep
    import epykit._regional_gls_native as native

    rng = np.random.default_rng(775)
    rows = []
    for j in range(6):
        sample = f"s{j}"
        rows.append(dict(sample_id=sample, treatment=int(j >= 3)))
        part = tmp_path / "counts" / f"sample={sample}" / "chrom=chr1"
        part.mkdir(parents=True)
        m = rng.binomial(40, 0.7 if j >= 3 else 0.2, 40)
        pl.DataFrame(dict(pos=np.arange(40) * 20, N_meth=m, coverage=[40] * 40)).write_parquet(
            part / "part-0.parquet"
        )
    md = ep.MethylData(
        obs=pl.DataFrame(rows),
        store=str(tmp_path / "counts"),
        analysis_root=str(tmp_path / "analysis"),
    )
    ep.tl.dmr(md, method="regional_gls", n_perm=20, tsv=False)
    before = md.uns["dmr_candidates"].clone()

    def no_refit(*args, **kwargs):
        raise AssertionError("Unexpected regional refit")

    monkeypatch.setattr(native, "fit_regions", no_refit)
    ep.tl.dmr(md, method="regional_gls", n_perm=20, tsv=False)
    assert md.uns["dmr_params"]["cache_hit"] is True
    assert before.equals(md.uns["dmr_candidates"])
    assert not (tmp_path / "counts" / ".cache").exists()
    target = tmp_path / "counts" / "sample=s0" / "chrom=chr1" / "part-0.parquet"
    pl.read_parquet(target).with_columns((pl.col("N_meth") + 1).alias("N_meth")).write_parquet(
        target
    )
    with pytest.raises(AssertionError, match="Unexpected regional refit"):
        ep.tl.dmr(md, method="regional_gls", n_perm=20, tsv=False)
    with pytest.raises(ValueError, match="geometry"):
        ep.tl.dmr(md, method="regional_gls", window_bp=0, tsv=False)
