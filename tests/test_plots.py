"""Layer 6b: matplotlib plot smoke tests.

Renders every public ``ep.pl.*`` function once against the fixture using
the Agg backend (no display required). We don't assert on pixel content --
just that the function returns a Figure/Axes without raising.

A failure here typically means an API drift (e.g. a column rename in the
DMC output that the plot consumes by name).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import polars as pl
import pytest


# Belt-and-suspenders: even though conftest already forces Agg, ensure it
# here for tests that may be run in isolation.
matplotlib.use("Agg", force=True)


@pytest.fixture
def synth_md_with_dmc(synth_md_filtered):
    """A filtered MethylData with a DMC and DMR table populated."""
    import epykit as ep
    ep.tl.dmc(synth_md_filtered, test="lr")
    ep.tl.dmr(
        synth_md_filtered,
        method="tile",
        tile_size_bp=500,
        min_cpgs_per_tile=3,
    )
    return synth_md_filtered


def _close_all():
    """Close every open figure between tests so memory doesn't snowball."""
    plt.close("all")



# Differential analysis plots


def test_pl_volcano(synth_md_with_dmc, tmp_path):
    import epykit as ep
    ep.pl.volcano(synth_md_with_dmc, save=str(tmp_path / "volcano"))
    _close_all()


def test_pl_ma_plot(synth_md_with_dmc, tmp_path):
    import epykit as ep
    ep.pl.ma_plot(synth_md_with_dmc, save=str(tmp_path / "ma"))
    _close_all()


def test_pl_manhattan(synth_md_with_dmc, tmp_path):
    import epykit as ep
    ep.pl.manhattan(synth_md_with_dmc, save=str(tmp_path / "manhattan"))
    _close_all()



# QC plots


def test_pl_coverage_histogram(synth_md_with_dmc, tmp_path):
    import epykit as ep
    ep.pl.coverage_histogram(synth_md_with_dmc, save=str(tmp_path / "cov_hist"))
    _close_all()


def test_pl_methylation_heatmap(synth_md_with_dmc, tmp_path):
    import epykit as ep
    # n_top kept small to stay fast.
    ep.pl.methylation_heatmap(
        synth_md_with_dmc, n_top=50,
        save=str(tmp_path / "heatmap"),
    )
    _close_all()



# Clustering


def test_pl_pca(synth_md_with_dmc, tmp_path):
    import epykit as ep
    ep.pl.pca(synth_md_with_dmc, save=str(tmp_path / "pca"))
    _close_all()



# Genomic context plots (need annotated DMCs)


@pytest.fixture
def synth_md_with_annotation(synth_md_with_dmc):
    """Inject a synthetic ``feature_type`` / ``cpg_context`` column into the
    DMC table so the genomic-context plots have data to consume. This skips
    the full ``annotate_features`` path (which needs a real GTF)."""
    import numpy as np

    md = synth_md_with_dmc
    df = md.varm["dmc_lr"]
    rng = np.random.default_rng(0)
    features  = rng.choice(["promoter", "exon", "intron", "intergenic"], len(df))
    contexts  = rng.choice(["island", "shore", "shelf", "open_sea"], len(df))
    annotated = df.with_columns([
        pl.Series("feature_type", features.tolist()),
        pl.Series("cpg_context",  contexts.tolist()),
    ])
    md.varm["dmc_lr_annotated"] = annotated
    return md


def test_pl_genomic_context_bar(synth_md_with_annotation, tmp_path):
    import epykit as ep
    ep.pl.genomic_context_bar(
        synth_md_with_annotation, save=str(tmp_path / "ctx_bar")
    )
    _close_all()


def test_pl_cpg_island_pie(synth_md_with_annotation, tmp_path):
    import epykit as ep
    ep.pl.cpg_island_pie(
        synth_md_with_annotation, save=str(tmp_path / "cpg_pie")
    )
    _close_all()



# Lazy-load contract: pl namespace still works as an attribute, not a
# fragile lazy proxy.


def test_pl_namespace_is_real_module():
    """``epykit.pl`` should be the actual submodule, not a custom proxy
    (we removed the lazy-load shim in this session)."""
    import epykit
    import types
    assert isinstance(epykit.pl, types.ModuleType)
    for fn in ("volcano", "manhattan", "pca", "coverage_histogram"):
        assert hasattr(epykit.pl, fn), f"epykit.pl missing {fn}"
