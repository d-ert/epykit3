"""Permutation empirical FDR for DMC."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

import epykit as ep
from tests.fixtures.synth import SimConfig, generate

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def small_md(tmp_path_factory):
    """A small fixture for the permutation test -- n_perm is expensive, so
    keep n_sites modest. 2 chromosomes x 400 CpGs = 800 candidate sites is
    enough to exercise the empirical FDR pipeline without blowing the test
    budget."""
    cfg = SimConfig(
        n_per_group=3,
        cpgs_per_chrom=400,
        chromosomes=("chr1", "chr2"),
        n_dmrs=2,
        dmr_size_cpgs=5,
        n_scattered_dmcs=80,
        seed=7,
    )
    out_dir = tmp_path_factory.mktemp("empfdr")
    result = generate(cfg, out_dir)
    md = ep.read_bismark(
        result["samplesheet"],
        treatment_group="treatment",
        control_group="control",
        store_dir=str(out_dir / "store"),
    )
    ep.pp.filter_coverage(md, lo_count=3, hi_perc=99.9)
    ep.pp.unite(md, type="intersect")
    truth = pl.read_parquet(result["truth"])
    return md, truth


def test_empirical_fdr_columns_present_and_in_range(small_md):
    """With empirical_fdr=True, tl.dmc must add finite empirical_pvalue
    and empirical_qvalue columns and write the params to md.uns."""
    md, _truth = small_md
    ep.tl.dmc(md, test="lr", empirical_fdr=True, n_perm=10, perm_seed=11)
    df = md.dmc

    assert "empirical_pvalue" in df.columns, "empirical_pvalue not added"
    assert "empirical_qvalue" in df.columns, "empirical_qvalue not added"

    emp_p = df.get_column("empirical_pvalue").to_numpy()
    emp_q = df.get_column("empirical_qvalue").to_numpy()

    # +1/+1 correction means the lowest possible empirical p at n_perm=10
    # is 1/(N_null + 1) > 0, and the upper bound is 1.0.
    emp_p_finite = emp_p[np.isfinite(emp_p)]
    assert emp_p_finite.size > 0, "all empirical p-values are NaN"
    assert (emp_p_finite > 0).all(), "empirical p-values collapsed to 0"
    assert (emp_p_finite <= 1.0).all(), "empirical p-values exceed 1.0"

    emp_q_finite = emp_q[np.isfinite(emp_q)]
    assert emp_q_finite.size > 0, "all empirical q-values are NaN"
    assert (emp_q_finite >= 0).all() and (emp_q_finite <= 1.0).all()

    # uns should record the permutation settings so the analysis is
    # reproducible.
    uns = md.uns.get("dmc", {})
    assert uns.get("empirical_fdr") is True
    assert uns.get("n_perm") == 10
    assert uns.get("perm_seed") == 11


def test_empirical_fdr_refuses_contrast_path(small_md):
    """The formula= / contrast= GLM path should refuse empirical_fdr with
    a ValueError pointing at the right design issue."""
    md, _truth = small_md
    with pytest.raises(ValueError, match="empirical_fdr"):
        ep.tl.dmc(
            md, formula="~ group", contrast="group",
            empirical_fdr=True, n_perm=5,
        )


def test_empirical_fdr_tracks_asymptotic_on_null_sites(small_md):
    """Sanity check: on null (non-DMC) sites, the empirical and asymptotic
    p-value distributions should agree roughly in the bulk. We check the
    median of each is non-trivially > 0.1 (null sites should be spread
    over the full unit interval)."""
    md, truth = small_md
    ep.tl.dmc(md, test="lr", empirical_fdr=True, n_perm=10, perm_seed=21)
    df = md.dmc
    joined = (
        truth.with_columns(pl.col("pos").cast(pl.Int64))
        .join(
            df.with_columns(pl.col("pos").cast(pl.Int64)),
            on=["chrom", "pos"], how="left",
        )
    )
    null_sites = joined.filter(
        ~pl.col("is_dmc") & ~pl.col("in_dmr")
    )
    asym = null_sites.get_column("pvalue").drop_nulls().to_numpy()
    emp = null_sites.get_column("empirical_pvalue").drop_nulls().to_numpy()
    assert asym.size > 50 and emp.size > 50, (
        f"too few null sites to evaluate: asym={asym.size}, emp={emp.size}"
    )
    # Null p-value medians should both clear 0.1; tight equality would
    # require many more perms than we run here.
    assert np.median(asym) > 0.1, (
        f"asymptotic null p-value median collapsed: {np.median(asym):.3f}"
    )
    assert np.median(emp) > 0.1, (
        f"empirical null p-value median collapsed: {np.median(emp):.3f}"
    )


# ---- DMR permutation (merged from test_dmr_permutation.py) ---------------


def test_dmr_empirical_fdr_columns(synth_md_filtered):
    """tl.dmr(empirical_fdr=True) appends empirical_pvalue / qvalue columns."""
    md = synth_md_filtered
    ep.tl.dmr(
        md,
        method="tile",
        empirical_fdr=True,
        n_perm=5,
        perm_seed=42,
        chromosomes=["chr1"],
    )
    dmr = md.uns["dmr"]
    assert isinstance(dmr, pl.DataFrame)
    if len(dmr) > 0:
        assert "empirical_pvalue" in dmr.columns
        assert "empirical_qvalue" in dmr.columns
        emp = dmr.get_column("empirical_pvalue").drop_nulls().to_numpy()
        assert (emp >= 0).all() and (emp <= 1).all()
    params = md.uns["dmr_params"]
    assert params.get("empirical_fdr") is True
    assert params.get("n_perm") == 5
