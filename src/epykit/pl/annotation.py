"""Annotatr-style annotation plots.

Five plot types mirroring the R ``annotatr`` package's
``plot_annotation`` / ``plot_numerical`` / ``plot_coannotations`` /
``plot_categorical`` family. All consume already-annotated DMC or DMR
tables on ``md.dmc`` / ``md.uns['dmr']`` (run :func:`epykit.tl.annotate`
first; pass ``multi_annotation=True`` if you want the co-annotation /
proportion plots).

Each function follows the existing epykit signature pattern:
``(md, *, ..., ax=None, figsize=..., save=None) -> (fig, axes)``.
"""

from __future__ import annotations

from math import ceil
from typing import Optional, Sequence

import numpy as np
import polars as pl

from .._style import PALETTE
from ._compute import (
    _resolve_annotated_table,
    compute_annotation_counts,
    compute_categorical_proportions,
    compute_coannotation_matrix,
    compute_numerical_by_annotation,
)
from ._utils import _get_ax, _save_fig
from ..methyldata import MethylData


# A modest colour-blind safe cycle for feature / context classes. Falls
# back through PALETTE for known names (promoter / island / etc.) and
# uses this cycle for anything else.
_ANNOT_CYCLE = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#F0E442", "#56B4E9", "#E69F00", "#999999",
]


def _color_for(label: str, fallback_idx: int) -> str:
    if label in PALETTE:
        return PALETTE[label]
    return _ANNOT_CYCLE[fallback_idx % len(_ANNOT_CYCLE)]


def plot_annotation_counts(
    md: MethylData,
    *,
    level: str = "dmr",
    annot_col: str = "feature_type",
    kind: str = "bar",
    autopct: str | None = "%1.1f%%",
    title: str | None = None,
    ax=None,
    figsize: tuple | None = None,
    save: str | None = None,
):
    """Bar or pie chart of region counts per annotation class.

    Equivalent to annotatr's ``plot_annotation``. Pass
    ``annot_col="all_overlapping_features"`` to count the multi-annotation
    explode (each region contributes once per overlapping class).

    Parameters
    ----------
    kind : {"bar", "pie"}
        Render as a vertical bar chart (default) or a pie chart. Pies are
        the right pick when proportions matter and there are <=8 classes;
        bars when absolute counts are the message.
    autopct : str or None
        Pie-only: format string for slice percentages. Pass ``None`` to
        hide the labels.
    """
    if kind not in ("bar", "pie"):
        raise ValueError(f"kind must be 'bar' or 'pie'; got {kind!r}")

    df = _resolve_annotated_table(md, level)
    counts = compute_annotation_counts(df, annot_col=annot_col)

    labels = counts[annot_col].to_list()
    values = counts["count"].to_list()
    colors = [_color_for(lbl, i) for i, lbl in enumerate(labels)]

    if figsize is None:
        figsize = (6, 5) if kind == "pie" else (7, 4)
    fig, ax = _get_ax(ax, figsize)

    if kind == "bar":
        ax.bar(labels, values, color=colors)
        ax.set_xlabel(annot_col)
        ax.set_ylabel(f"# {level.upper()}s")
        ax.tick_params(axis="x", rotation=30)
    else:
        # Don't crowd small slices with text -- only label slices >=3%.
        total = sum(values) or 1
        slice_labels = [
            lbl if (v / total) >= 0.03 else ""
            for lbl, v in zip(labels, values)
        ]
        wedges, texts, autotexts = ax.pie(
            values, labels=slice_labels, colors=colors,
            autopct=autopct, startangle=90,
            wedgeprops=dict(edgecolor="white", linewidth=1.0),
            textprops=dict(fontsize=9),
        )
        # Place a full legend off to the right so even small slices are
        # readable in print.
        ax.legend(
            wedges,
            [f"{lbl} ({v:,})" for lbl, v in zip(labels, values)],
            title=annot_col, loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False, fontsize=8,
        )

    ax.set_title(title or f"Annotation counts ({level})")
    if kind == "pie":
        fig.tight_layout()

    if save:
        _save_fig(md, fig, save)
    return fig, ax


def plot_numerical_by_annotation(
    md: MethylData,
    *,
    value: str = "meth_diff",
    annot_col: str = "feature_type",
    level: str = "dmr",
    bins: int = 30,
    background: bool = True,
    density: bool = True,
    ncols: int = 3,
    figsize: Optional[tuple] = None,
    sharex: bool = True,
    save: str | None = None,
):
    """Faceted histogram of a numeric column by annotation class.

    Reproduces annotatr's ``plot_numerical``: one panel per annotation
    class, drawn as grey filled bars, with a red-outlined "All"
    background overlay (when ``background=True``) so the per-class shape
    can be compared to the genome-wide distribution.

    The reference figure 1A in the design doc is this plot with
    ``value="meth_diff"`` and ``annot_col`` running across CpG-island
    contexts.
    """
    df = _resolve_annotated_table(md, level)
    long = compute_numerical_by_annotation(
        df, value_col=value, annot_col=annot_col, include_all=background,
    )
    classes = [c for c in long.get_column("annot").unique().to_list() if c != "All"]
    classes = sorted(classes)

    if background:
        bg_vals = (
            long.filter(pl.col("annot") == "All").get_column(value).to_numpy()
        )
    else:
        bg_vals = None

    # Bin edges shared across panels so the "All" overlay aligns.
    flat = long.filter(pl.col("annot") != "All").get_column(value).to_numpy()
    if flat.size == 0:
        raise ValueError(
            f"No values for {value!r} across any class; check the annotated table."
        )
    edges = np.linspace(float(np.nanmin(flat)), float(np.nanmax(flat)), bins + 1)

    nrows = ceil(len(classes) / ncols)
    if figsize is None:
        figsize = (3.4 * ncols, 2.4 * nrows)

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(
        nrows, ncols, figsize=figsize, sharex=sharex,
        squeeze=False,
    )
    axes_flat = axes.flat

    for i, cls in enumerate(classes):
        ax = axes_flat[i]
        vals = (
            long.filter(pl.col("annot") == cls).get_column(value).to_numpy()
        )
        ax.hist(
            vals, bins=edges, density=density,
            color=PALETTE.get("neutral", "#888888"),
            alpha=0.7, label="Data",
        )
        if background and bg_vals is not None and bg_vals.size:
            ax.hist(
                bg_vals, bins=edges, density=density,
                histtype="step", color=PALETTE.get("hyper", "#e05263"),
                linewidth=1.4, label="Background",
            )
        ax.set_title(str(cls), fontsize=10)
        ax.set_xlabel(value)
        if density:
            ax.set_ylabel("density")
        else:
            ax.set_ylabel("count")

    # Hide empty cells
    for j in range(len(classes), nrows * ncols):
        axes_flat[j].axis("off")

    # Single shared legend on the first axes.
    if background:
        axes_flat[0].legend(loc="best", fontsize=8, frameon=False)
    fig.tight_layout()

    if save:
        _save_fig(md, fig, save)
    return fig, axes


def plot_coannotations(
    md: MethylData,
    *,
    annot_col: str = "all_overlapping_features",
    level: str = "dmr",
    mode: str = "heatmap",
    annot_log: bool = False,
    cmap: str = "viridis",
    ax=None,
    figsize=(6, 5),
    save: str | None = None,
):
    """Pairwise co-occurrence of annotation classes.

    Two modes:

    * ``mode="heatmap"`` (default) -- count matrix of regions that hit
      every pair of annotations. Diagonal counts regions that hit that
      single class. Cells are labelled with the count.
    * ``mode="proportion"`` -- divide each cell by the row-class total so
      cell (i, j) reads "fraction of class-i regions that *also* hit j".

    Set ``annot_log=True`` to colour-map ``log10(count + 1)`` -- helpful
    when one cell dwarfs the others. Requires ``multi_annotation=True``
    on :func:`epykit.tl.annotate`.
    """
    df = _resolve_annotated_table(md, level)
    matrix, classes = compute_coannotation_matrix(df, annot_col=annot_col)

    display = matrix.astype(np.float64)
    if mode == "proportion":
        diag = np.diag(matrix).astype(np.float64)
        diag[diag == 0] = np.nan
        display = display / diag[:, None]
        cbar_label = "proportion"
        fmt = "{:.2f}"
    elif mode == "heatmap":
        cbar_label = "count"
        fmt = "{:,.0f}"
    else:
        raise ValueError("mode must be 'heatmap' or 'proportion'")

    if annot_log and mode == "heatmap":
        display = np.log10(display + 1.0)
        cbar_label = "log10(count + 1)"

    fig, ax = _get_ax(ax, figsize)
    im = ax.imshow(display, cmap=cmap, aspect="auto")
    ax.set_xticks(np.arange(len(classes)))
    ax.set_xticklabels(classes, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(classes)))
    ax.set_yticklabels(classes)
    # Cell labels: print raw counts even when log-colour is applied.
    for i in range(len(classes)):
        for j in range(len(classes)):
            if mode == "proportion":
                val = display[i, j]
                text = "--" if np.isnan(val) else fmt.format(val)
            else:
                text = fmt.format(matrix[i, j])
            ax.text(
                j, i, text, ha="center", va="center",
                color="white" if (display[i, j] if not np.isnan(display[i, j]) else 0)
                                > np.nanmean(display) else "black",
                fontsize=8,
            )
    fig.colorbar(im, ax=ax, label=cbar_label, fraction=0.04, pad=0.02)
    ax.set_title(f"Co-annotation ({level}, mode={mode})")

    if save:
        _save_fig(md, fig, save)
    return fig, ax


def plot_categorical(
    md: MethylData,
    *,
    group_col: str = "dmr_type",
    annot_col: str = "feature_type",
    level: str = "dmr",
    position: str = "fill",
    include_all_group: bool = True,
    annotation_order: Optional[Sequence[str]] = None,
    ax=None,
    figsize=(6, 5),
    save: str | None = None,
):
    """Stacked / dodged bar chart of annotation classes within categories.

    Equivalent to annotatr's ``plot_categorical``. Defaults match the
    reference figure panel C: x-axis is ``all/hyper/hypo``, fill is
    feature class, stacked as proportions of each bar.

    Parameters
    ----------
    group_col : str
        Column on the annotated table that drives the x-axis. Defaults
        to epykit's ``dmr_type`` (hyper/hypo); auto-falls-back to
        ``DM_status`` if that exists instead.
    annot_col : str
        Column that drives the fill. Accepts list-valued multi-annotation
        columns -- a region with two annotations contributes to both
        slices.
    position : {"fill", "stack"}
        ``"fill"`` (default) normalises each bar to 100 percent. ``"stack"``
        leaves raw counts.
    include_all_group : bool
        Prepend an ``"all"`` bar that aggregates every region regardless
        of ``group_col``. Matches the reference figure 1C layout.
    annotation_order : sequence of str, optional
        Order the stacks from bottom to top. Defaults to alphabetical;
        pass an explicit order for journal-style consistency across
        figures.
    """
    df = _resolve_annotated_table(md, level)
    props = compute_categorical_proportions(
        df, group_col=group_col, annot_col=annot_col,
        include_all_group=include_all_group,
        normalize=(position == "fill"),
    )
    # If autodetection in the compute layer swapped to DM_status / dmr_type,
    # honour it on the resulting columns.
    resolved_group = group_col if group_col in props.columns else (
        "DM_status" if "DM_status" in props.columns else "dmr_type"
    )
    value_col = "proportion" if position == "fill" else "count"
    groups_present = props.get_column(resolved_group).unique().to_list()
    # Put "all" first when included, then sorted remainder.
    if "all" in groups_present:
        ordered_groups = ["all"] + sorted(g for g in groups_present if g != "all")
    else:
        ordered_groups = sorted(groups_present, key=str)

    classes_present = props.get_column(annot_col).unique().to_list()
    if annotation_order is not None:
        ordered_classes = [c for c in annotation_order if c in classes_present]
        ordered_classes += [c for c in classes_present if c not in ordered_classes]
    else:
        ordered_classes = sorted(classes_present, key=str)

    # Build (len(groups) x len(classes)) matrix.
    mat = np.zeros((len(ordered_groups), len(ordered_classes)), dtype=np.float64)
    for row in props.to_dicts():
        gi = ordered_groups.index(row[resolved_group])
        ci = ordered_classes.index(row[annot_col])
        mat[gi, ci] = float(row[value_col] or 0.0)

    fig, ax = _get_ax(ax, figsize)
    bottom = np.zeros(len(ordered_groups))
    for ci, cls in enumerate(ordered_classes):
        ax.bar(
            ordered_groups, mat[:, ci],
            bottom=bottom, label=str(cls),
            color=_color_for(cls, ci),
            edgecolor="white", linewidth=0.4,
        )
        bottom += mat[:, ci]

    ax.set_xlabel(resolved_group)
    ax.set_ylabel("Proportion" if position == "fill" else "Count")
    ax.set_title(f"{annot_col} by {resolved_group}")
    if position == "fill":
        ax.set_ylim(0, 1.0)
    ax.legend(
        title=annot_col, bbox_to_anchor=(1.02, 1.0),
        loc="upper left", frameon=False, fontsize=9,
    )
    fig.tight_layout()

    if save:
        _save_fig(md, fig, save)
    return fig, ax


__all__ = [
    "plot_annotation_counts",
    "plot_numerical_by_annotation",
    "plot_coannotations",
    "plot_categorical",
]
