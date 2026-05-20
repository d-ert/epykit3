"""DMR / DMC set-overlap visualisations.

UpSet-style plot for comparing 2-6 DMR (or DMC) sets across contrasts or
methods. The implementation is matplotlib-only -- no upsetplot dependency --
because the typical methylation use case is 2-5 sets where a hand-rolled
UpSet matrix is small and avoids dragging in a transitive pandas-version
constraint from upsetplot's lower bounds.

For 2 sets the same function falls back to a clean 2-circle Venn diagram
(no matplotlib_venn dependency; we draw the circles directly).
"""

from __future__ import annotations

from itertools import chain, combinations
from typing import Iterable, Mapping

import numpy as np
import polars as pl

from ._utils import _get_ax, _save_fig


def _dmr_to_key_set(df: pl.DataFrame, key_cols: tuple[str, ...]) -> set[tuple]:
    """Materialise a DMR table as a set of identifier tuples."""
    if df is None or df.is_empty():
        return set()
    missing = [c for c in key_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"DMR table missing key columns {missing}; "
            f"available: {sorted(df.columns)[:15]}..."
        )
    return set(
        tuple(row[c] for c in key_cols)
        for row in df.select(list(key_cols)).iter_rows(named=True)
    )


def _powerset_nonempty(items: Iterable[str]) -> list[tuple[str, ...]]:
    items = list(items)
    return [
        combo for r in range(1, len(items) + 1)
        for combo in combinations(items, r)
    ]


def dmr_overlap(
    sets: Mapping[str, pl.DataFrame],
    *,
    key_cols: tuple[str, ...] = ("chrom", "start", "end"),
    min_size: int = 1,
    sort_by: str = "size",
    figsize=(8, 5),
    ax=None,
    save: str | None = None,
    md=None,
):
    """Plot the intersection structure of 2-6 DMR sets.

    For 2 sets, draws a 2-circle Venn (overlap proportional to actual
    intersection size). For 3+ sets, draws an UpSet plot:

    * Top bar chart: size of each non-empty intersection, descending.
    * Bottom dot matrix: which sets participate in that intersection
      (filled dots connected by a line).
    * Right-side bar: total size of each input set.

    ``key_cols`` controls how rows are compared. For DMR tables the
    default ``(chrom, start, end)`` matches identical tile coordinates;
    pass ``("chrom", "pos")`` to compare DMC tables instead. For "any
    DMR that overlaps another by >=1 bp" semantics, pre-merge your tables
    on intervals (e.g. via bioframe) and pass a unique merged-region ID.

    Parameters
    ----------
    sets : dict[str, pl.DataFrame]
        Label -> DMR / DMC table. 2 <= len(sets) <= 6.
    key_cols : tuple[str, ...]
        Columns that together identify "the same DMR" across tables.
    min_size : int
        Drop intersections with fewer than this many rows (after binning
        into combos). Default 1.
    sort_by : {"size", "degree"}
        UpSet only: order combos by intersection size (default) or by
        the number of sets in the combo.
    md : MethylData, optional
        Forwarded to ``_save_fig`` for path resolution when ``save`` is
        set; never read otherwise.

    Returns
    -------
    (Figure, Axes | tuple[Axes, ...])
    """
    if not isinstance(sets, Mapping) or len(sets) < 2:
        raise ValueError("Provide at least 2 named DMR sets as a dict.")
    if len(sets) > 6:
        raise ValueError(
            f"dmr_overlap supports up to 6 sets; got {len(sets)}. "
            "For larger comparisons aggregate sets first or use upsetplot."
        )

    labels = list(sets.keys())
    set_keys = {lbl: _dmr_to_key_set(sets[lbl], key_cols) for lbl in labels}

    if len(labels) == 2:
        return _venn2(
            labels, set_keys, ax=ax, figsize=figsize, save=save, md=md,
        )
    return _upset(
        labels, set_keys, min_size=min_size, sort_by=sort_by,
        ax=ax, figsize=figsize, save=save, md=md,
    )


# 2-set Venn (no external dep)

def _venn2(labels, set_keys, *, ax, figsize, save, md):
    a, b = labels
    A, B = set_keys[a], set_keys[b]
    only_a = len(A - B)
    only_b = len(B - A)
    both = len(A & B)

    fig, ax = _get_ax(ax, figsize)
    # Two circles, fixed radii, fixed centres.
    from matplotlib.patches import Circle
    r = 1.0
    cx_a, cx_b = -0.7, 0.7
    ax.add_patch(Circle((cx_a, 0), r, alpha=0.4, color="tab:blue", label=a))
    ax.add_patch(Circle((cx_b, 0), r, alpha=0.4, color="tab:orange", label=b))
    ax.text(cx_a - 0.55, 0, f"{only_a:,}", ha="center", va="center", fontsize=12)
    ax.text(cx_b + 0.55, 0, f"{only_b:,}", ha="center", va="center", fontsize=12)
    ax.text(0, 0, f"{both:,}", ha="center", va="center", fontsize=12, weight="bold")
    ax.text(cx_a, 1.1, a, ha="center", va="bottom", fontsize=10, weight="bold")
    ax.text(cx_b, 1.1, b, ha="center", va="bottom", fontsize=10, weight="bold")
    ax.set_xlim(-2, 2)
    ax.set_ylim(-1.4, 1.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"Overlap: {a} vs {b}  (|A & B|={both:,})")

    if save and md is not None:
        _save_fig(md, fig, save)
    return fig, ax


# UpSet for 3+ sets

def _upset(labels, set_keys, *, min_size, sort_by, ax, figsize, save, md):
    if ax is not None:
        # An UpSet needs three coupled axes; if the caller passed `ax`
        # we honour it for the bar chart only and skip the matrix /
        # totals to stay composable.
        return _upset_single_ax(
            labels, set_keys, min_size=min_size, sort_by=sort_by,
            ax=ax, save=save, md=md,
        )

    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    # Build the combo -> size table.
    combos = _powerset_nonempty(labels)
    combo_sizes: dict[tuple[str, ...], int] = {}
    for combo in combos:
        inside = set.intersection(*(set_keys[lbl] for lbl in combo))
        outside = set.union(*(set_keys[lbl] for lbl in labels if lbl not in combo)) \
            if len(combo) < len(labels) else set()
        size = len(inside - outside)
        if size >= min_size:
            combo_sizes[combo] = size

    if not combo_sizes:
        raise ValueError(
            f"No intersection has >={min_size} elements; nothing to plot."
        )

    if sort_by == "size":
        ordered = sorted(combo_sizes.items(), key=lambda kv: -kv[1])
    elif sort_by == "degree":
        ordered = sorted(combo_sizes.items(), key=lambda kv: (len(kv[0]), -kv[1]))
    else:
        raise ValueError(f"sort_by must be 'size' or 'degree'; got {sort_by!r}")

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(
        2, 2, figure=fig,
        width_ratios=[len(ordered), max(1.5, len(labels) * 0.6)],
        height_ratios=[3, max(1.0, len(labels) * 0.4)],
        hspace=0.05, wspace=0.05,
    )
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_mat = fig.add_subplot(gs[1, 0], sharex=ax_bar)
    ax_tot = fig.add_subplot(gs[1, 1], sharey=ax_mat)
    ax_corner = fig.add_subplot(gs[0, 1])
    ax_corner.axis("off")

    sizes = [s for _, s in ordered]
    xs = np.arange(len(ordered))
    ax_bar.bar(xs, sizes, color="tab:blue", edgecolor="black", linewidth=0.5)
    ax_bar.set_ylabel("Intersection size")
    ax_bar.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    for x, s in zip(xs, sizes):
        ax_bar.text(x, s, f"{s:,}", ha="center", va="bottom", fontsize=8)

    # Matrix: y is set, x is combo. Filled dot = set is in combo.
    label_to_y = {lbl: i for i, lbl in enumerate(labels)}
    ax_mat.set_xlim(-0.5, len(ordered) - 0.5)
    ax_mat.set_ylim(-0.5, len(labels) - 0.5)
    ax_mat.invert_yaxis()
    ax_mat.set_yticks(np.arange(len(labels)))
    ax_mat.set_yticklabels(labels)
    ax_mat.set_xticks([])
    for spine in ("top", "right", "bottom"):
        ax_mat.spines[spine].set_visible(False)
    for x, (combo, _) in enumerate(ordered):
        for y, lbl in enumerate(labels):
            colour = "black" if lbl in combo else "lightgrey"
            ax_mat.scatter([x], [y], s=80, color=colour, zorder=2)
        ys = sorted(label_to_y[lbl] for lbl in combo)
        if len(ys) >= 2:
            ax_mat.plot([x, x], [ys[0], ys[-1]], color="black", lw=1.5, zorder=1)
    for y in range(len(labels)):
        ax_mat.axhline(y, color="lightgrey", lw=0.5, alpha=0.5, zorder=0)

    # Totals on the right.
    totals = [len(set_keys[lbl]) for lbl in labels]
    ax_tot.barh(
        np.arange(len(labels)), totals,
        color="tab:gray", edgecolor="black", linewidth=0.5,
    )
    ax_tot.invert_xaxis()
    ax_tot.tick_params(axis="y", which="both", left=False, labelleft=False)
    ax_tot.set_xlabel("Set size")
    for spine in ("top", "right"):
        ax_tot.spines[spine].set_visible(False)
    for y, t in enumerate(totals):
        ax_tot.text(t, y, f"{t:,}", ha="right", va="center", fontsize=8)

    fig.suptitle(f"DMR overlap UpSet ({len(labels)} sets)", y=0.98)

    if save and md is not None:
        _save_fig(md, fig, save)
    return fig, (ax_bar, ax_mat, ax_tot)


def _upset_single_ax(labels, set_keys, *, min_size, sort_by, ax, save, md):
    """Compact fallback when the caller passes a pre-made `ax`."""
    combos = _powerset_nonempty(labels)
    combo_sizes = {}
    for combo in combos:
        inside = set.intersection(*(set_keys[lbl] for lbl in combo))
        outside = set.union(*(set_keys[lbl] for lbl in labels if lbl not in combo)) \
            if len(combo) < len(labels) else set()
        size = len(inside - outside)
        if size >= min_size:
            combo_sizes[combo] = size
    if not combo_sizes:
        raise ValueError(f"No intersection has >={min_size} elements.")
    if sort_by == "size":
        ordered = sorted(combo_sizes.items(), key=lambda kv: -kv[1])
    else:
        ordered = sorted(combo_sizes.items(), key=lambda kv: (len(kv[0]), -kv[1]))
    xs = np.arange(len(ordered))
    sizes = [s for _, s in ordered]
    ax.bar(xs, sizes, color="tab:blue", edgecolor="black", linewidth=0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels(
        [" & ".join(c) for c, _ in ordered], rotation=45, ha="right", fontsize=8,
    )
    ax.set_ylabel("Intersection size")
    ax.set_title("DMR overlap")
    if save and md is not None:
        _save_fig(md, ax.figure, save)
    return ax.figure, ax


__all__ = ["dmr_overlap"]
