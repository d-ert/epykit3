"""Multi-group / continuous-covariate contrasts via tl.dmc.

These tests pin down statistical recovery on the multi-group and continuous-
covariate fixture extensions -- not just column presence -- so a regression in
the joint F-test or in the Wald-on-continuous-coef path fails loudly.
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

import epykit as ep
from tests.fixtures.synth import SimConfig, generate

pytestmark = pytest.mark.slow


# Calibrated against the fixture. The multi-group joint F-test (3 groups,
# 4 reps each, effect step 0.20 => Deltabeta_max~=0.40 across levels) is well-
# powered; the continuous-covariate fixture is weaker (4 vs 4 samples, age
# slope 0.5pp Deltabeta/year, ages drawn U(20,80), so the per-sample effective
# Deltabeta varies with the realised age sample -- 8 samples leaves a noisy slope
# estimate). The continuous test therefore checks structural correctness
# (engine runs, finite p-values, FDR not catastrophic) rather than power.
MULTIGROUP_POWER_MIN = 0.30
MULTIGROUP_FDR_MAX   = 0.15


@pytest.fixture(scope="module")
def multigroup_md(tmp_path_factory):
    cfg = SimConfig(
        n_groups=3,
        n_per_group=4,
        cpgs_per_chrom=600,
        chromosomes=("chr1", "chr2"),
        seed=4242,
    )
    out_dir = tmp_path_factory.mktemp("multigroup")
    result = generate(cfg, out_dir)
    md = ep.read_bismark(
        result["samplesheet"],
        groups=list(result["group_ids"].keys()),
        store_dir=str(out_dir / "store"),
    )
    ep.pp.filter_coverage(md, lo_count=3, hi_perc=99.9)
    ep.pp.unite(md, type="intersect")
    return md, pl.read_parquet(result["truth"]), cfg


@pytest.fixture(scope="module")
def continuous_md(tmp_path_factory):
    cfg = SimConfig(
        n_per_group=4,
        cpgs_per_chrom=600,
        chromosomes=("chr1", "chr2"),
        continuous_covariate=True,
        n_groups=2,
        seed=4243,
    )
    out_dir = tmp_path_factory.mktemp("continuous")
    result = generate(cfg, out_dir)
    md = ep.read_bismark(
        result["samplesheet"],
        treatment_group="treatment",
        control_group="control",
        store_dir=str(out_dir / "store"),
    )
    md.obs = md.obs.with_columns(
        pl.col("age").cast(pl.Float64, strict=False)
    )
    ep.pp.filter_coverage(md, lo_count=3, hi_perc=99.9)
    ep.pp.unite(md, type="intersect")
    return md, pl.read_parquet(result["truth"]), cfg


def _join_truth(df: pl.DataFrame, truth: pl.DataFrame) -> pl.DataFrame:
    return (
        truth.with_columns(pl.col("pos").cast(pl.Int64))
        .join(
            df.with_columns(pl.col("pos").cast(pl.Int64)),
            on=["chrom", "pos"],
            how="left",
        )
    )


def test_multigroup_factor_joint_test(multigroup_md):
    """A 3-group joint F-test recovers seeded multi-group DMCs and keeps
    empirical FDR near the nominal alpha.
    """
    md, truth, _cfg = multigroup_md
    ep.tl.dmc(md, formula="~ group", contrast="group")
    df = md.varm["dmc_glm_contrast"]
    assert df is not None
    assert "f_stat" in df.columns
    assert "df1" in df.columns
    assert df.get_column("df1")[0] >= 2  # 3 levels -> df1 == 2
    level_cols = [c for c in df.columns if c.startswith("mean_beta_")]
    assert len(level_cols) >= 3
    assert "qvalue" in df.columns

    joined = _join_truth(df, truth)
    mg_truth = joined.filter(pl.col("is_multigroup_dmc"))
    assert mg_truth.height > 0, (
        "Fixture should produce >0 multi-group seed sites; check "
        "SimConfig(n_groups=3, ...) wiring."
    )

    # Power: fraction of seeded multi-group DMCs called at q<0.05.
    sig_mg = mg_truth.filter(pl.col("qvalue") < 0.05).height
    power = sig_mg / mg_truth.height
    assert power >= MULTIGROUP_POWER_MIN, (
        f"Multi-group power too low: {sig_mg}/{mg_truth.height} = "
        f"{power:.3f} < {MULTIGROUP_POWER_MIN}"
    )

    # FDR: among all q<0.05 calls, how many are not real signal? "Real" =
    # any of the three seeded effect kinds (binary DMC, scattered DMC, or
    # multi-group DMC). Anything else flagged is a false positive.
    called = joined.filter(pl.col("qvalue") < 0.05)
    if called.height > 0:
        is_truly_diff = (
            pl.col("is_dmc")
            | pl.col("is_multigroup_dmc")
        )
        fp = called.filter(~is_truly_diff).height
        fdr = fp / called.height
        assert fdr <= MULTIGROUP_FDR_MAX, (
            f"Multi-group FDR too high: {fp}/{called.height} = "
            f"{fdr:.3f} > {MULTIGROUP_FDR_MAX}"
        )


def test_continuous_covariate_primary(continuous_md):
    """Continuous covariate as primary effect -- structural correctness check
    on the Wald-on-age path.

    The fixture (4 vs 4 samples, 200 age-DMCs at 0.5pp Deltabeta/year, ages drawn
    U(20,80)) is under-powered for an absolute power floor: a small but
    real slope at n=8 leaves wide SEs, and ~30 % of sites separate on the
    binomial GLM because low-coverage sites hit beta=0/1 boundaries. So we
    test what *should* hold structurally -- engine ran, finite p-values,
    FDR not catastrophic on the q<0.05 calls -- rather than recovery rate.
    Strengthening the fixture (or testing power at a tighter alpha) is on
    the 0.3 roadmap.
    """
    md, truth, _cfg = continuous_md
    ep.tl.dmc(md, formula="~ age", contrast="age")
    df = md.varm["dmc_glm_contrast"]
    assert df is not None
    assert "coef_treatment" in df.columns or "f_stat" in df.columns
    assert "qvalue" in df.columns

    # Engine produced finite p-values (didn't NaN out wholesale).
    finite_p = df.filter(pl.col("pvalue").is_finite()).height
    assert finite_p > 0.3 * df.height, (
        f"Continuous Wald produced too few finite p-values: "
        f"{finite_p}/{df.height}"
    )

    joined = _join_truth(df, truth)
    age_truth = joined.filter(pl.col("is_age_dmc"))
    assert age_truth.height > 0, (
        "Fixture should produce >0 age-DMC seed sites; check "
        "SimConfig(continuous_covariate=True, ...) wiring."
    )

    # FDR: among called sites, any that's not a seeded effect of any kind
    # is a false positive. With ~5% nominal alpha and BH correction the
    # observed FDR should stay below ~0.5 even with weak signal -- anything
    # at or above 0.8 would indicate a calibration bug.
    called = joined.filter(pl.col("qvalue") < 0.05)
    if called.height > 0:
        is_truly_diff = (
            pl.col("is_age_dmc") | pl.col("is_dmc")
            | pl.col("is_multigroup_dmc")
        )
        fp = called.filter(~is_truly_diff).height
        fdr = fp / called.height
        assert fdr < 0.80, (
            f"Continuous-covariate FDR catastrophically high: "
            f"{fp}/{called.height} = {fdr:.3f} (>=0.80 suggests broken "
            "q-value calibration, not just low power)"
        )


def test_named_contrast_single_row(multigroup_md):
    """A patsy-style linear contrast either runs successfully (produces a
    p-value column) or fails with a clean ValueError. Anything else (crash,
    silent NaN-only output) is a regression.

    Note: this is the "tough resolution path" -- we don't pin the result,
    just guarantee no unexpected exception type slips through.
    """
    md, _truth, _cfg = multigroup_md
    levels = sorted(md.obs.get_column("group").unique().to_list())
    a, b = levels[0], levels[1]
    contrast = f"group[T.{a}] - group[T.{b}]"
    try:
        ep.tl.dmc(md, formula="~ group", contrast=contrast)
    except ValueError as exc:
        # Acceptable: the contrast string failed to resolve. The message
        # should mention the contrast or the formula so users can debug.
        msg = str(exc).lower()
        assert any(token in msg for token in ("contrast", "formula", "group")), (
            f"ValueError raised but message is unhelpful: {exc!r}"
        )
        return
    df = md.varm["dmc_glm_contrast"]
    assert df is not None, "contrast resolved but no varm entry written"
    assert "pvalue" in df.columns, "contrast resolved but pvalue column missing"
    # Should have produced *some* finite p-values on a fixture this size.
    finite_p = df.filter(pl.col("pvalue").is_finite()).height
    assert finite_p > 0, "all p-values NaN -- contrast silently failed"
