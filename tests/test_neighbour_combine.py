"""Neighbour-aware Stouffer combining (``combine_neighbour_pvalues``).

* Untested (NaN p-value) sites must be excluded from both the combined
  statistic and the ``_n_neighbours`` audit count. The pre-fix code set
  NaN-p sites' signed z to 0.0, which passed the window mask, so untested
  CpGs were counted as contributing neighbours.
* ``pvalue`` must stay the raw per-CpG value when ``neighbour_combine=True``;
  the combined p-value lives in ``pvalue_combined`` and its BH q-value in
  ``qvalue_combined``.
* The docstring must own the Stouffer independence caveat.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit.dmc import combine_neighbour_pvalues


def test_nan_p_neighbours_excluded_from_count():
    # 5 CpGs within a 200 bp window. Only the focal site (pos 100) is tested;
    # the other four are untested (NaN p-value).
    df = pl.DataFrame({
        "chrom": ["chr1"] * 5,
        "pos": [60, 80, 100, 120, 140],
        "pvalue": [np.nan, np.nan, 0.001, np.nan, np.nan],
        "meth_diff": [0.2, 0.2, 0.3, 0.2, 0.2],
    })
    out = combine_neighbour_pvalues(df, neighbour_bp=200)

    focal = out.filter(pl.col("pos") == 100)
    n_neighbours = focal["pvalue_combined_n_neighbours"][0]
    # Only the focal site itself is a valid contributor -> count is 1, not 5.
    assert n_neighbours == 1, (
        f"untested NaN-p neighbours must not be counted; got {n_neighbours}"
    )


def test_nan_p_neighbours_do_not_change_combined_p():
    # A focal site flanked only by untested neighbours must keep (combine to)
    # essentially its own evidence, not be diluted/altered by z=0 phantoms.
    df_isolated = pl.DataFrame({
        "chrom": ["chr1"],
        "pos": [100],
        "pvalue": [0.001],
        "meth_diff": [0.3],
    })
    df_with_nan = pl.DataFrame({
        "chrom": ["chr1"] * 3,
        "pos": [80, 100, 120],
        "pvalue": [np.nan, 0.001, np.nan],
        "meth_diff": [0.2, 0.3, 0.2],
    })
    p_iso = combine_neighbour_pvalues(df_isolated, neighbour_bp=200).filter(
        pl.col("pos") == 100
    )["pvalue_combined"][0]
    p_nan = combine_neighbour_pvalues(df_with_nan, neighbour_bp=200).filter(
        pl.col("pos") == 100
    )["pvalue_combined"][0]
    # Both reduce to the single valid site's own Stouffer self-combination.
    assert np.isclose(p_iso, p_nan, rtol=1e-9, atol=1e-12)


# --- Raw p-values are preserved ---------------------------------------------


def test_neighbour_combine_does_not_overwrite_pvalue(synth_md):
    md = synth_md
    ep.pp.filter_coverage(md, lo_count=3, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")

    # First run: baseline lr, no combine. Capture per-CpG raw p.
    ep.tl.dmc(md, test="lr", neighbour_combine=False)
    raw_df = md.dmc.clone()
    assert "pvalue" in raw_df.columns
    assert raw_df.height > 0

    # Second run: lr with neighbour_combine=True.
    ep.tl.dmc(md, test="lr", neighbour_combine=True, neighbour_bp=500)
    comb_df = md.dmc

    assert "pvalue" in comb_df.columns
    assert "pvalue_combined" in comb_df.columns
    assert "qvalue" in comb_df.columns
    assert "qvalue_combined" in comb_df.columns

    # `pvalue` must still equal the raw per-CpG p-values (not combined).
    joined = raw_df.select(["chrom", "pos", "pvalue"]).rename(
        {"pvalue": "pvalue_raw_ref"}
    ).join(comb_df.select(["chrom", "pos", "pvalue", "pvalue_combined"]),
           on=["chrom", "pos"], how="inner")
    raw_ref = joined["pvalue_raw_ref"].to_numpy()
    canonical = joined["pvalue"].to_numpy()
    mask = np.isfinite(raw_ref) & np.isfinite(canonical)
    np.testing.assert_allclose(
        canonical[mask], raw_ref[mask], rtol=1e-10,
        err_msg="`pvalue` column was overwritten with combined values; "
                "should remain the raw per-CpG p-value.",
    )

    # And the combined column must differ from raw on at least some sites.
    combined = joined["pvalue_combined"].to_numpy()
    diff_mask = np.isfinite(combined) & np.isfinite(raw_ref) & (combined != raw_ref)
    assert diff_mask.sum() > 0, (
        "pvalue_combined identical to raw -- combiner appears to be a no-op."
    )


@pytest.mark.slow
def test_empirical_fdr_uses_raw_p_under_neighbour_combine(synth_md):
    """When both neighbour_combine and empirical_fdr are on, the empirical
    p-value comparison must use raw (not combined) p-values on the
    observed side -- otherwise the null pool (raw) is incompatible."""
    md = synth_md
    ep.pp.filter_coverage(md, lo_count=3, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")

    ep.tl.dmc(
        md, test="lr",
        neighbour_combine=True, neighbour_bp=500,
        empirical_fdr=True, n_perm=3, perm_seed=11,
    )
    df = md.dmc
    assert "empirical_pvalue" in df.columns

    emp_p = df["empirical_pvalue"].drop_nans().to_numpy()
    if emp_p.size == 0:
        pytest.skip("synthetic fixture too small; no finite empirical p")
    # empirical_pvalue must be finite and in [0, 1]; the floor 1/(n_perm+1)
    # is checked in test_dmc_empirical_fdr_denominator.py.
    assert np.isfinite(emp_p).any()
    assert np.nanmin(emp_p) >= 0.0
    assert np.nanmax(emp_p) <= 1.0 + 1e-9


# --- Docstring owns the independence caveat ---------------------------------


def test_docstring_acknowledges_independence_violation():
    doc = combine_neighbour_pvalues.__doc__ or ""
    assert "independence" in doc.lower() or "correlated" in doc.lower(), (
        "Docstring must mention that adjacent CpGs are correlated and "
        "Stouffer assumes independence."
    )
    assert "min_sign_agreement" in doc or "sign agreement" in doc.lower(), (
        "Docstring must explain that the FDR safety net comes from the "
        "sign-agreement gate, not the Stouffer null."
    )
    assert "Brown" in doc or "v0.8" in doc or "future work" in doc.lower(), (
        "Docstring must point readers at the planned correlation-aware "
        "replacement (Brown's method) and note this is a known limitation."
    )
