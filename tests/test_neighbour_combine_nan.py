"""D9: the neighbour Stouffer combine must exclude untested (NaN p-value)
sites from both the combined statistic and the `_n_neighbours` audit count.
The pre-fix code set NaN-p sites' signed z to 0.0, which passed the
window mask, so untested CpGs were counted as contributing neighbours."""
import numpy as np
import polars as pl

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
