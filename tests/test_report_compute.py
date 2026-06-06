"""Tests for the new report compute helpers and Plotly twins."""

from __future__ import annotations

import numpy as np
import pytest


pytest.importorskip("plotly")


def test_pvalue_histogram(synth_md_filtered):
    import epykit as ep
    from epykit.pl._compute import compute_pvalue_histogram

    ep.tl.dmc(synth_md_filtered, test="lr")
    counts, edges = compute_pvalue_histogram(synth_md_filtered, bins=20)
    assert counts.sum() > 0
    assert len(edges) == len(counts) + 1
    assert edges[0] == 0.0 and abs(edges[-1] - 1.0) < 1e-9


def test_dmr_size_distribution(synth_md_filtered):
    import epykit as ep
    from epykit.pl._compute import compute_dmr_size_distribution

    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    ep.tl.dmr(md, method="tile", tile_size_bp=500, min_cpgs_per_tile=3,
              min_mean_qvalue=1.0)
    sizes = compute_dmr_size_distribution(md)
    assert sizes.ndim == 1


def test_global_methylation(synth_md_filtered):
    import epykit as ep
    from epykit.pl._compute import compute_global_methylation

    ep.tl.qc(synth_md_filtered)
    samples, values, groups = compute_global_methylation(synth_md_filtered)
    assert len(samples) == len(values) == synth_md_filtered.n_samples
    assert all(0.0 <= v <= 1.0 for v in values if v == v)


def test_correlation_matrix(synth_md_filtered):
    import epykit as ep
    from epykit.pl._compute import compute_sample_correlation_matrix

    ep.tl.qc(synth_md_filtered, run_sample_correlation=True)
    mat, labels = compute_sample_correlation_matrix(synth_md_filtered)
    assert mat.shape == (len(labels), len(labels))
    assert np.allclose(np.diag(mat), 1.0, atol=0.05)


def test_scree(synth_md_filtered):
    from epykit.pl._compute import compute_scree

    ev = compute_scree(synth_md_filtered, n_sites=2000, max_components=4)
    assert ev.ndim == 1 and len(ev) >= 1
    assert (ev >= 0).all()


def test_scatter_subsample_keeps_all_significant(synth_md_filtered):
    """max_points caps total points but never drops a significant CpG, so
    hyper/hypo counts and Manhattan peaks stay exact on huge tables."""
    import epykit as ep
    from epykit.pl._compute import (
        compute_volcano_data, compute_ma_data, compute_manhattan_data,
    )
    ep.tl.dmc(synth_md_filtered, test="lr")

    full = compute_volcano_data(synth_md_filtered, alpha=0.05, min_abs_diff=0.1)
    n_total = full.meth_diff.size
    n_sig = int(full.sig.sum())
    assert n_sig > 0 and n_total > n_sig  # fixture has both sig and ns

    cap = n_sig + 50
    v = compute_volcano_data(synth_md_filtered, alpha=0.05, min_abs_diff=0.1, max_points=cap)
    assert v.meth_diff.size <= cap
    assert int(v.sig.sum()) == n_sig  # every significant point retained
    assert int(v.hyper.sum()) + int(v.hypo.sum()) == n_sig

    m = compute_ma_data(synth_md_filtered, alpha=0.05, min_abs_diff=0.1, max_points=cap)
    assert m.mean_beta.size <= cap and int(m.sig.sum()) == n_sig

    # Manhattan: total rendered points (summed over chrom blocks) is capped.
    man = compute_manhattan_data(synth_md_filtered, alpha=0.05, max_points=cap)
    rendered = sum(b["n"] for b in man.chrom_blocks)
    assert rendered <= cap


def test_plotly_twins_smoke(synth_md_filtered):
    import epykit as ep
    from epykit.pl import _plotly as P

    md = synth_md_filtered
    ep.tl.qc(md, run_sample_correlation=True)
    ep.tl.dmc(md, test="lr")
    ep.tl.dmr(md, method="tile", tile_size_bp=500, min_cpgs_per_tile=3,
              min_mean_qvalue=1.0)
    for fn in (P.pvalue_histogram_plotly, P.global_methylation_bar_plotly,
               P.sample_correlation_plotly, P.scree_plotly,
               P.dmr_size_hist_plotly):
        fig = fn(md)
        assert fig is not None and len(fig.data) >= 1
