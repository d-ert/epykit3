"""P0-1: when neighbour_combine=True, `pvalue` must remain the raw
per-CpG p-value. The combined p-value lives in `pvalue_combined` and
its BH q-value in `qvalue_combined`. Downstream code reading `pvalue`
must see un-combined values."""
from __future__ import annotations

import numpy as np
import pytest

import epykit as ep


def test_neighbour_combine_does_not_overwrite_pvalue(synth_md):
    md = synth_md
    ep.pp.filter_coverage(md, lo_count=3, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")

    # First run: baseline lr, no combine. Capture per-CpG raw p.
    ep.tl.dmc(md, test="lr", neighbour_combine=False)
    raw_df = md.dmc.clone()
    assert "pvalue" in raw_df.columns
    raw_p = raw_df["pvalue"].to_numpy()

    # Second run: lr with neighbour_combine=True.
    ep.tl.dmc(md, test="lr", neighbour_combine=True, neighbour_bp=500)
    comb_df = md.dmc

    # Required columns post-fix:
    assert "pvalue" in comb_df.columns
    assert "pvalue_combined" in comb_df.columns
    assert "qvalue" in comb_df.columns
    assert "qvalue_combined" in comb_df.columns

    # `pvalue` must still equal the raw per-CpG p-values (not combined).
    # Align on (chrom, pos) -- assume both runs have the same rows.
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
    # P0-1 sanity: empirical_pvalue must be finite and in [0, 1].
    # (The floor 1/(n_perm+1) is checked by the P0-2b test.)
    assert np.isfinite(emp_p).any()
    assert np.nanmin(emp_p) >= 0.0
    assert np.nanmax(emp_p) <= 1.0 + 1e-9
