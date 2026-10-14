"""Whole-scan coverage eligibility, union partitions and reported effects."""

import numpy as np
import polars as pl

import epykit as ep
from epykit import _regional_gls


def test_group_eligibility_is_recomputed_for_each_assignment():
    mask = np.array([[1, 1, 0, 1, 1, 0], [1, 1, 1, 1, 0, 0]], dtype=float)
    beta = mask * np.array([.2, .2, .2, .6, .6, .6])
    original = np.array([0, 0, 0, 1, 1, 1])
    permuted = np.array([0, 1, 1, 0, 0, 1])
    raw, shrunk, eligible = _regional_gls._scan_group_effects(beta, mask, original, 2)
    assert eligible.tolist() == [True, False]
    assert raw[0] == .4 or np.isclose(raw[0], .4)
    assert np.isclose(shrunk[0], .4*np.sqrt(2/3))
    _, _, eligible_permuted = _regional_gls._scan_group_effects(beta, mask, permuted, 2)
    assert eligible_permuted.tolist() == [False, True]


def test_union_scan_preserves_chromosome_absent_from_first_sample(tmp_path):
    rng = np.random.default_rng(978)
    rows, observed_m = [], []
    for s in range(6):
        sample = f"s{s}"
        rows.append(dict(sample_id=sample, treatment=int(s >= 3)))
        directory = tmp_path / "counts" / f"sample={sample}"
        directory.mkdir(parents=True)
        if s == 0:
            continue
        chrom = directory / "chrom=chr1"
        chrom.mkdir()
        m = rng.binomial(40, .75 if s >= 3 else .25, 15)
        observed_m.append(m)
        pl.DataFrame(dict(pos=np.arange(15)*20, N_meth=m, coverage=[40]*15)).write_parquet(
            chrom / "part-0.parquet"
        )
    md = ep.MethylData(obs=pl.DataFrame(rows), store=str(tmp_path / "counts"))
    _, candidates, meta = _regional_gls.run_regional_gls(
        md, n_perm=20, min_cpgs=5, min_effect=.1, q_cutoff=None, cache=False, workers=1,
        record_scan_scores=True,
    )
    assert meta["common_cpgs"]["chr1"] == 15
    assert meta["union_cpgs"]["chr1"] == 15
    assert len(candidates) > 0
    assert len(meta["scan_scores"]["observed"]) == len(candidates)
    assert len(meta["scan_scores"]["unique_null"]) == meta["unique_null_fits"]
    counts = np.array(observed_m)
    raw_effect = counts[2:].mean(axis=0)/40-counts[:2].mean(axis=0)/40
    for row in candidates.to_dicts():
        indices = (np.arange(15)*20 >= row["start"]) & (np.arange(15)*20 < row["end"])
        assert np.isclose(row["mean_meth_diff"], raw_effect[indices].mean())
        assert row["residual_df"] == 4*row["n_cpgs"]-1
