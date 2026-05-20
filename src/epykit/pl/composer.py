"""Multi-panel figure composer for publication-ready layouts.

Wraps matplotlib's ``subplot_mosaic`` so callers can describe an A/B/C/D
layout as a mosaic string and drop pre-built plots into each slot.
Letters of the layout become uppercase panel labels in the corner of
each axes, like the figure reference in the design doc.

Example
-------
.. code-block:: python

    from epykit.pl.composer import figure_grid
    import epykit as ep

    figure_grid(
        panels={
            "A": (ep.pl.volcano, dict(md=md)),
            "B": (ep.pl.tss_metaplot, dict(md=md, gtf_path="genes.gtf.gz")),
            "C": (ep.pl.plot_categorical, dict(md=md)),
            "D": (ep.pl.plot_coannotations, dict(md=md)),
        },
        layout="A B / C D",
        figsize=(11, 8),
        save="figure1",
    )
"""

from __future__ import annotations

from typing import Callable, Optional, Sequence

import numpy as np

from .._style import apply_theme


def _render_table_axes(ax, df, *, max_rows: int = 20, max_cols: int = 8) -> None:
    """Render a polars / pandas / list-of-dicts frame as a matplotlib table."""
    try:
        import polars as pl
        if isinstance(df, pl.DataFrame):
            df = df.head(max_rows).to_pandas()
    except Exception:
        pass
    try:
        import pandas as pd
        if isinstance(df, pd.DataFrame):
            df = df.head(max_rows).iloc[:, :max_cols]
            cell_text = df.astype(str).values.tolist()
            col_labels = list(df.columns)
        else:
            raise TypeError
    except Exception:
        # Best-effort: assume list of dicts
        rows = list(df)[:max_rows]
        col_labels = list(rows[0].keys()) if rows else []
        cell_text = [[str(r.get(c, "")) for c in col_labels] for r in rows]

    ax.axis("off")
    tbl = ax.table(
        cellText=cell_text, colLabels=col_labels,
        loc="center", cellLoc="left",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1.0, 1.2)


def figure_grid(
    panels: dict,
    *,
    layout: str | Sequence[Sequence[str]] = "A B / C D",
    figsize: tuple = (11, 7),
    label_fontsize: float = 14,
    label_weight: str = "bold",
    label_offset: tuple = (-0.08, 1.02),
    theme: Optional[str] = "publication",
    width_ratios: Optional[Sequence[float]] = None,
    height_ratios: Optional[Sequence[float]] = None,
    save: str | None = None,
    md=None,
):
    """Assemble a multi-panel matplotlib figure.

    Parameters
    ----------
    panels : dict
        Maps layout letters (``"A"`` / ``"B"`` / ...) to one of:

        * ``(plot_fn, kwargs)`` -- ``plot_fn`` is called with the axes
          injected as ``ax=<axes>`` plus the keyword args. The function
          should *not* call ``plt.show`` / ``plt.savefig`` itself.
        * ``("table", df)`` -- render ``df`` as a static table in that
          slot.
        * a ``Callable`` -- treated as ``(callable, {})``.
        * pre-built ``(fig, ax)`` -- only the axes content is reused;
          the donor figure is closed automatically.
    layout : str | sequence[sequence[str]]
        Mosaic spec for :meth:`matplotlib.figure.Figure.subplot_mosaic`.
        String form uses ``"A B / C D"`` (rows separated by ``/``, cells
        by spaces; repeat a letter to merge cells). Sequence form is
        ``[["A", "B"], ["C", "C"]]``.
    label_fontsize, label_weight, label_offset
        Panel-label styling. Defaults match Nature-family conventions.
    theme : str or None
        Theme to apply before drawing. ``None`` skips re-theming. Use
        ``"publication"`` for the 300dpi / TrueType-embedded preset.
    width_ratios, height_ratios
        Forwarded to mosaic's ``gridspec_kw``.
    save : str, optional
        Saves the composite via :func:`epykit.pl._utils._save_fig` when
        ``md`` is provided (falls back to ``fig.savefig`` otherwise).
    md : MethylData, optional
        Required when ``save`` is set and you want the figure path tracked
        on ``md.uns['figures']``.

    Returns
    -------
    (fig, axes_dict)
        ``axes_dict`` maps each letter to its axes object for further
        tweaking before saving.
    """
    if theme:
        apply_theme(theme)

    import matplotlib.pyplot as plt

    if isinstance(layout, str):
        mosaic = [row.split() for row in layout.split("/") if row.strip()]
    else:
        mosaic = [list(row) for row in layout]

    gridspec_kw = {}
    if width_ratios is not None:
        gridspec_kw["width_ratios"] = list(width_ratios)
    if height_ratios is not None:
        gridspec_kw["height_ratios"] = list(height_ratios)

    fig, axd = plt.subplot_mosaic(
        mosaic, figsize=figsize, gridspec_kw=gridspec_kw,
    )

    for letter, ax in axd.items():
        entry = panels.get(letter)
        if entry is None:
            ax.axis("off")
            continue

        # 1. Pre-built (fig, ax) handoff: copy its content over.
        if (
            isinstance(entry, tuple) and len(entry) == 2
            and hasattr(entry[0], "axes") and hasattr(entry[1], "figure")
        ):
            src_fig, src_ax = entry
            _adopt_axes_content(src_ax, ax)
            plt.close(src_fig)

        # 2. Table.
        elif isinstance(entry, tuple) and len(entry) == 2 and entry[0] == "table":
            _render_table_axes(ax, entry[1])

        # 3. (callable, kwargs) -- the most common path.
        elif isinstance(entry, tuple) and len(entry) == 2 and callable(entry[0]):
            fn, kwargs = entry
            kwargs = dict(kwargs)
            kwargs.setdefault("ax", ax)
            fn(**kwargs)

        # 4. Bare callable.
        elif callable(entry):
            entry(ax=ax)

        else:
            ax.text(
                0.5, 0.5, f"unrecognised panel for {letter!r}",
                ha="center", va="center", transform=ax.transAxes,
            )
            ax.axis("off")

        # Panel label (A / B / C / ...) in the top-left corner.
        ax.text(
            label_offset[0], label_offset[1], letter,
            transform=ax.transAxes,
            fontsize=label_fontsize, fontweight=label_weight,
            va="top", ha="right",
        )

    fig.tight_layout()

    if save:
        from ._utils import _save_fig
        if md is not None:
            _save_fig(md, fig, save)
        else:
            from pathlib import Path
            path = Path(f"figures/{save}.png")
            path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(path, dpi=300, bbox_inches="tight")
            plt.close(fig)

    return fig, axd


def _adopt_axes_content(src_ax, dst_ax) -> None:
    """Replay common artists from ``src_ax`` onto ``dst_ax``.

    This is a best-effort port for compositing pre-built figures. It
    covers lines, scatter (PathCollection), bars (Rectangles), images
    and text -- which spans the artist set the rest of epykit's plots
    actually use. Anything fancier should be rebuilt by passing
    ``(fn, kwargs)`` instead so the renderer draws directly into the
    target axes.
    """
    from matplotlib.collections import PathCollection
    from matplotlib.patches import Rectangle
    from matplotlib.image import AxesImage

    for ln in src_ax.get_lines():
        dst_ax.plot(
            ln.get_xdata(), ln.get_ydata(),
            color=ln.get_color(), linestyle=ln.get_linestyle(),
            linewidth=ln.get_linewidth(), label=ln.get_label(),
        )
    for col in src_ax.collections:
        if isinstance(col, PathCollection):
            offs = col.get_offsets()
            if len(offs):
                dst_ax.scatter(
                    offs[:, 0], offs[:, 1],
                    s=col.get_sizes() if col.get_sizes().size else 12,
                    c=col.get_facecolors() if col.get_facecolors().size else None,
                    alpha=col.get_alpha(),
                )
    for patch in src_ax.patches:
        if isinstance(patch, Rectangle):
            dst_ax.add_patch(Rectangle(
                patch.get_xy(), patch.get_width(), patch.get_height(),
                facecolor=patch.get_facecolor(), edgecolor=patch.get_edgecolor(),
                linewidth=patch.get_linewidth(), alpha=patch.get_alpha(),
            ))
    for img in src_ax.get_images():
        if isinstance(img, AxesImage):
            arr = np.asarray(img.get_array())
            dst_ax.imshow(arr, cmap=img.get_cmap(), aspect="auto",
                          vmin=img.norm.vmin, vmax=img.norm.vmax)
    for txt in src_ax.texts:
        dst_ax.text(*txt.get_position(), txt.get_text(),
                    fontsize=txt.get_fontsize(), color=txt.get_color())

    dst_ax.set_xlabel(src_ax.get_xlabel())
    dst_ax.set_ylabel(src_ax.get_ylabel())
    dst_ax.set_title(src_ax.get_title())
    dst_ax.set_xlim(src_ax.get_xlim())
    dst_ax.set_ylim(src_ax.get_ylim())
    handles, labels = src_ax.get_legend_handles_labels()
    if handles:
        dst_ax.legend(handles, labels, frameon=False, fontsize=8)


__all__ = ["figure_grid"]
