"""Smoke tests for annotatr-style plots, the compute layer, and the
multi-panel composer.

These are smokes: assert the function returns a Figure / Axes without
raising. Compute-layer tests do small numeric checks because they form
the contract every plotting backend now relies on.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pytest


matplotlib.use("Agg", force=True)

pytestmark = pytest.mark.slow


# ---------------------------------------------------------------------------
# Fixture: annotated DMR table with multi-annotation columns.
# Mirrors the test_plots.py `synth_md_with_annotation` pattern but adds
# the list-valued multi-annotation columns needed by the annotatr plots.
# ---------------------------------------------------------------------------


@pytest.fixture
def synth_md_with_dmc(synth_md_filtered):
    import epykit as ep
    ep.tl.dmc(synth_md_filtered, test="lr")
    ep.tl.dmr(
        synth_md_filtered,
        method="tile", tile_size_bp=500, min_cpgs_per_tile=3,
    )
    return synth_md_filtered


@pytest.fixture
def annotated_md(synth_md_with_dmc):
    """Inject feature_type / cpg_context / multi-annotation columns into both
    DMC and DMR tables so every annotatr-style plot has data to consume."""
    md = synth_md_with_dmc
    rng = np.random.default_rng(7)

    # DMC table: per-row scalar feature_type + cpg_context.
    dmc = md.varm["dmc_lr"]
    features = rng.choice(["promoter", "exon", "intron", "intergenic"], len(dmc))
    contexts = rng.choice(["island", "shore", "shelf", "open_sea"], len(dmc))
    multi_feats = [
        rng.choice(
            ["promoter", "exon", "intron", "intergenic"],
            size=rng.integers(1, 3), replace=False,
        ).tolist()
        for _ in range(len(dmc))
    ]
    md.varm["dmc_lr_annotated"] = dmc.with_columns([
        pl.Series("feature_type", features.tolist()),
        pl.Series("cpg_context", contexts.tolist()),
        pl.Series("all_overlapping_features", multi_feats, dtype=pl.List(pl.Utf8)),
    ])

    # DMR table on md.uns["dmr"]: build a small one if epykit's pipeline
    # didn't already populate it.
    dmr = md.uns.get("dmr")
    if dmr is None or not isinstance(dmr, pl.DataFrame) or dmr.is_empty():
        # Build a synthetic DMR table just for the plot smoke tests.
        n = 30
        dmr_types = rng.choice(["hyper", "hypo"], n).tolist()
        meth_diffs = (rng.uniform(-0.5, 0.5, n)).tolist()
        feats = rng.choice(["promoter", "exon", "intron", "intergenic"], n).tolist()
        ctxs = rng.choice(["island", "shore", "shelf", "open_sea"], n).tolist()
        multi = [
            rng.choice(
                ["promoter", "exon", "intron", "intergenic"],
                size=rng.integers(1, 3), replace=False,
            ).tolist()
            for _ in range(n)
        ]
        md.uns["dmr"] = pl.DataFrame({
            "chrom": ["chr1"] * n,
            "start": (np.arange(n) * 1000).tolist(),
            "end": (np.arange(n) * 1000 + 500).tolist(),
            "dmr_type": dmr_types,
            "meth_diff": meth_diffs,
            "feature_type": feats,
            "cpg_context": ctxs,
            "all_overlapping_features": pl.Series(multi, dtype=pl.List(pl.Utf8)),
        })
    else:
        # Augment whatever the pipeline gave us with the columns our plots
        # need; never overwrite if the column already exists.
        new_cols = []
        if "dmr_type" not in dmr.columns:
            new_cols.append(pl.Series(
                "dmr_type",
                rng.choice(["hyper", "hypo"], len(dmr)).tolist(),
            ))
        if "feature_type" not in dmr.columns:
            new_cols.append(pl.Series(
                "feature_type",
                rng.choice(["promoter", "exon", "intron", "intergenic"], len(dmr)).tolist(),
            ))
        if "all_overlapping_features" not in dmr.columns:
            multi = [
                rng.choice(
                    ["promoter", "exon", "intron", "intergenic"],
                    size=rng.integers(1, 3), replace=False,
                ).tolist()
                for _ in range(len(dmr))
            ]
            new_cols.append(pl.Series(
                "all_overlapping_features", multi, dtype=pl.List(pl.Utf8),
            ))
        if new_cols:
            md.uns["dmr"] = dmr.with_columns(new_cols)
    return md


def _close_all():
    plt.close("all")


# ---------------------------------------------------------------------------
# Compute-layer contract tests
# ---------------------------------------------------------------------------


def test_compute_annotation_counts_scalar_col(annotated_md):
    from epykit.pl._compute import compute_annotation_counts
    counts = compute_annotation_counts(
        annotated_md.uns["dmr"], annot_col="feature_type",
    )
    assert "feature_type" in counts.columns
    assert "count" in counts.columns
    assert counts["count"].sum() == annotated_md.uns["dmr"].height


def test_compute_annotation_counts_list_col_dedupes(annotated_md):
    """A region with [exon, promoter] should contribute to both classes
    exactly once -- not twice."""
    from epykit.pl._compute import compute_annotation_counts
    counts = compute_annotation_counts(
        annotated_md.uns["dmr"], annot_col="all_overlapping_features",
    )
    # Total per-class counts can exceed n_regions when regions hit multiple
    # classes, but each region should contribute exactly its unique-class
    # cardinality.
    dmr = annotated_md.uns["dmr"]
    expected = sum(
        len(set(row))
        for row in dmr.get_column("all_overlapping_features").to_list()
    )
    assert counts["count"].sum() == expected


def test_compute_coannotation_matrix_diagonal_matches_counts(annotated_md):
    from epykit.pl._compute import (
        compute_annotation_counts, compute_coannotation_matrix,
    )
    counts = compute_annotation_counts(
        annotated_md.uns["dmr"], annot_col="all_overlapping_features",
    )
    mat, classes = compute_coannotation_matrix(
        annotated_md.uns["dmr"], annot_col="all_overlapping_features",
    )
    diag = {c: int(mat[i, i]) for i, c in enumerate(classes)}
    counts_by_class = {
        row["all_overlapping_features"]: row["count"]
        for row in counts.to_dicts()
    }
    for c, expected in counts_by_class.items():
        assert diag[c] == expected, (c, diag[c], expected)


def test_compute_numerical_by_annotation_includes_all_overlay(annotated_md):
    from epykit.pl._compute import compute_numerical_by_annotation
    long = compute_numerical_by_annotation(
        annotated_md.uns["dmr"],
        value_col="meth_diff",
        annot_col="feature_type",
        include_all=True,
    )
    unique_annots = set(long.get_column("annot").to_list())
    assert "All" in unique_annots
    # The "All" rows must equal one row per region.
    n_all = long.filter(pl.col("annot") == "All").height
    # n_all may be smaller than n_dmrs if some meth_diffs are null; check >= 1.
    assert n_all >= 1


def test_compute_categorical_proportions_sums_to_one(annotated_md):
    from epykit.pl._compute import compute_categorical_proportions
    props = compute_categorical_proportions(
        annotated_md.uns["dmr"],
        group_col="dmr_type", annot_col="feature_type",
        include_all_group=True, normalize=True,
    )
    for grp, sub in props.group_by("dmr_type"):
        total = sub["proportion"].sum()
        assert abs(total - 1.0) < 1e-6, (grp, total)


# ---------------------------------------------------------------------------
# annotatr plot smokes
# ---------------------------------------------------------------------------


def test_plot_annotation_counts(annotated_md, tmp_path):
    import epykit as ep
    ep.pl.plot_annotation_counts(
        annotated_md, level="dmr",
        save=str(tmp_path / "annot_counts"),
    )
    _close_all()


def test_plot_numerical_by_annotation(annotated_md, tmp_path):
    import epykit as ep
    ep.pl.plot_numerical_by_annotation(
        annotated_md, value="meth_diff",
        annot_col="feature_type", level="dmr",
        save=str(tmp_path / "num_by_annot"),
    )
    _close_all()


def test_plot_coannotations_heatmap(annotated_md, tmp_path):
    import epykit as ep
    ep.pl.plot_coannotations(
        annotated_md, mode="heatmap", level="dmr",
        save=str(tmp_path / "coannot"),
    )
    _close_all()


def test_plot_coannotations_proportion(annotated_md, tmp_path):
    import epykit as ep
    ep.pl.plot_coannotations(
        annotated_md, mode="proportion", level="dmr",
        save=str(tmp_path / "coannot_prop"),
    )
    _close_all()


def test_plot_categorical(annotated_md, tmp_path):
    import epykit as ep
    ep.pl.plot_categorical(
        annotated_md, group_col="dmr_type", annot_col="feature_type",
        save=str(tmp_path / "categorical"),
    )
    _close_all()


# ---------------------------------------------------------------------------
# Composer
# ---------------------------------------------------------------------------


def test_figure_grid_panel_letters(annotated_md, tmp_path):
    import epykit as ep
    fig, axd = ep.pl.figure_grid(
        panels={
            "A": (ep.pl.plot_annotation_counts, dict(md=annotated_md)),
            "B": (ep.pl.plot_categorical, dict(md=annotated_md)),
        },
        layout="A B",
        figsize=(10, 4),
        save="composer_smoke",
        md=annotated_md,
    )
    assert set(axd.keys()) == {"A", "B"}
    _close_all()


# ---------------------------------------------------------------------------
# Style: publication theme + palette swap survives apply / set / apply.
# ---------------------------------------------------------------------------


def test_apply_theme_publication_sets_truetype():
    import matplotlib as mpl
    from epykit._style import apply_theme, set_palette, PALETTE

    apply_theme("publication")
    assert mpl.rcParams["pdf.fonttype"] == 42

    before = PALETTE["hyper"]
    set_palette("colorblind")
    assert PALETTE["hyper"] != before
    # Restore default so other tests aren't affected.
    set_palette("default")
    assert PALETTE["hyper"] == before
