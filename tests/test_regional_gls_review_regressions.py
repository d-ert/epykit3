"""Independent final-review failures: exact fits and incompatible null statistics."""

import numpy as np
import pytest

from epykit._regional_gls import fit_region, run_regional_gls
from epykit._regional_gls_native import evaluate_covariance, fit_regions


@pytest.mark.parametrize("k", [3, 5, 8, 20])
@pytest.mark.parametrize("missing", [False, True])
@pytest.mark.parametrize("complement", [False, True])
def test_exact_count_fit_is_rejected_symmetrically(k, missing, complement):
    p = np.arange(k)*20.
    g = np.r_[np.zeros(3), np.ones(3)]
    n = np.full((k, 6), 30.)
    m = np.tile([10., 10., 10., 20., 20., 20.], (k, 1))
    if missing:
        m[0, 0] = n[0, 0] = 0.
    if complement:
        g = 1-g
    with pytest.raises(ValueError, match="Degenerate"):
        fit_region(p, m, n, g)
    with pytest.raises(ValueError, match="Degenerate"):
        fit_region(p, m, n, g, covariance=(100., .3))


def test_batch_cannot_silently_keep_a_numerically_exact_fit():
    p = np.arange(5)*20.
    g = np.r_[np.zeros(3), np.ones(3)]
    y = np.tile(g-.5, (5, 1)).T
    with pytest.raises(ValueError, match="Degenerate"):
        fit_regions(p, y, np.ones(5), g, [[0, 5]], workers=1)


def test_relative_degeneracy_check_preserves_small_nonzero_variance():
    rng = np.random.default_rng(918)
    p = np.arange(8)*20.
    g = np.r_[np.zeros(3), np.ones(3)]
    mask = np.ones((6, 8), bool)
    mask[0, 0] = False
    y = np.tile(g-.5, (8, 1)).T + 1e-8*rng.normal(size=(6, 8))
    y[~mask] = 0
    for scale in [1e-8, 1., 1e8]:
        a = evaluate_covariance(p, y*scale, np.ones(8), g, 100., .3, obs_mask=mask)
        b = evaluate_covariance(p, y*scale, np.ones(8), 1-g, 100., .3, obs_mask=mask)
        assert np.isfinite(a).all()
        assert a[2] > 0
        assert a[1]/scale == pytest.approx(.999999997, abs=1e-7)
        assert a[1:] == pytest.approx((-b[1], b[2]), rel=1e-6, abs=1e-20*scale)


def test_maximum_selector_rejects_observed_pooled_decoy_covariance():
    with pytest.raises(ValueError, match="covariance replay"):
        run_regional_gls(None, selection_method="max_t", pooled_decoy_covariance=True)


def test_pooled_covariance_invalidates_maximum_diagnostic(tmp_path):
    import polars as pl

    import epykit as ep

    rng = np.random.default_rng(8978)
    rows = []
    for s in range(6):
        sample = f"s{s}"
        rows.append(dict(sample_id=sample, treatment=int(s >= 3)))
        directory = tmp_path / "counts" / f"sample={sample}" / "chrom=chr1"
        directory.mkdir(parents=True)
        pl.DataFrame(dict(pos=np.arange(15)*20, N_meth=rng.binomial(40, .75 if s >= 3 else .25, 15),
                          coverage=[40]*15)).write_parquet(directory / "part-0.parquet")
    md = ep.MethylData(obs=pl.DataFrame(rows), store=str(tmp_path / "counts"))
    _, candidates, meta = run_regional_gls(md, n_perm=20, pooled_decoy_covariance=True,
                                         q_cutoff=None, cache=False, workers=1)
    assert len(candidates) > 0
    assert not meta["maxT_available"]
    assert np.isnan(candidates["maxT_pvalue"].to_numpy()).all()
    assert "unavailable" in meta["maxT_scope"]
