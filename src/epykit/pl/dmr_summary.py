"""DMR-level summary plots: aggregate violin and clustered heatmap.

These complement the per-DMR strip plots in :mod:`epykit.pl.dmr_boxplot`
by summarising methylation level across entire DMR sets in a single
panel each. Designed to match the publication-grade Panel A
(hyper/hypo violins) and Panel B (clustered DMR heatmap with gene
labels) of multi-panel methylation figures.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import polars as pl

from .._style import PALETTE
from ._utils import _get_ax, _save_fig
from ..methyldata import MethylData


def _resolve_dmr(md: MethylData) -> pl.DataFrame:
    dmr = md.uns.get("dmr")
    if not isinstance(dmr, pl.DataFrame) or dmr.is_empty():
        raise ValueError(
            "md.uns['dmr'] is empty. Run ep.tl.dmr(md) first."
        )
    return dmr


def _group_label(md: MethylData, *, case_label: str | None,
                 control_label: str | None) -> tuple[str, str]:
    """Pick human-readable labels for the case/control violins/columns."""
    if control_label and case_label:
        return control_label, case_label
    if "group" in md.obs.columns:
        groups = list(dict.fromkeys(md.obs.get_column("group").to_list()))
        if len(groups) == 2:
            # md.control_ids / treatment_ids come from samplesheet groups;
            # fall back to alphabetical otherwise.
            ctrl = md.obs.filter(
                pl.col("sample_id").is_in(md.control_ids)
            ).get_column("group").unique().to_list()
            case = md.obs.filter(
                pl.col("sample_id").is_in(md.treatment_ids)
            ).get_column("group").unique().to_list()
            return (
                control_label or (ctrl[0] if ctrl else groups[0]),
                case_label or (case[0] if case else groups[1]),
            )
    return control_label or "control", case_label or "treatment"


def dmr_violin(
    md: MethylData,
    *,
    split: str = "dmr_type",
    value_case: str = "mean_beta_case",
    value_control: str = "mean_beta_control",
    case_label: str | None = None,
    control_label: str | None = None,
    show_box: bool = True,
    figsize: tuple | None = None,
    save: str | None = None,
):
    """Aggregate methylation-level violins per DMR direction.

    Faceted violin layout: one panel per value of ``split`` (default
    ``dmr_type`` -> ``hyper`` / ``hypo``), two violins per panel
    (control beta distribution vs case beta distribution across the
    DMRs in that direction).

    This is the "Panel A" of a typical multi-panel methylation figure
    (e.g. WT vs Het-KO across hyper-DMRs and hypo-DMRs).

    Parameters
    ----------
    split : str
        DMR-table column whose unique values produce the facets. Most
        common: ``"dmr_type"``. Pass any other categorical column
        (``"feature_type"``, ``"cpg_context"``) to slice by that
        instead -- one violin pair per category.
    value_case, value_control : str
        Columns on the DMR table holding the case-group and
        control-group per-DMR methylation level. Defaults match epykit's
        DMR output.
    case_label, control_label : str, optional
        Display names for the two violins. Default infers from
        ``md.obs["group"]`` mapped through ``md.control_ids`` /
        ``md.treatment_ids``.
    show_box : bool
        Overlay matplotlib's box-and-whisker on each violin (median +
        IQR + whiskers). Common in publication figures.
    """
    import matplotlib.pyplot as plt

    dmr = _resolve_dmr(md)
    if split not in dmr.columns:
        raise ValueError(
            f"split column {split!r} not on DMR table. "
            f"Try one of: {[c for c in dmr.columns if c in ('dmr_type', 'feature_type', 'cpg_context')]}"
        )
    if value_case not in dmr.columns or value_control not in dmr.columns:
        raise ValueError(
            f"need both {value_case!r} and {value_control!r} on DMR table; "
            "epykit's tl.dmr produces them by default."
        )
    ctrl_lbl, case_lbl = _group_label(
        md, case_label=case_label, control_label=control_label,
    )

    categories = (
        dmr.get_column(split).unique().to_list()
    )
    # Hyper first, then hypo, then anything else alphabetical -- matches the
    # typical figure ordering.
    rank = {"hyper": 0, "hypo": 1}
    categories = sorted(categories, key=lambda c: (rank.get(str(c), 99), str(c)))

    n = len(categories)
    if figsize is None:
        figsize = (3.4 * n, 4.2)
    fig, axes = plt.subplots(1, n, figsize=figsize, sharey=True, squeeze=False)
    axes = axes.flat

    color_ctrl = PALETTE.get("control", "#4a90d9")
    color_case = PALETTE.get("treatment", "#e05263")

    for i, cat in enumerate(categories):
        ax = axes[i]
        sub = dmr.filter(pl.col(split) == cat)
        ctrl = sub.get_column(value_control).drop_nulls().drop_nans().to_numpy()
        case = sub.get_column(value_case).drop_nulls().drop_nans().to_numpy()

        # Drop NaN that snuck through (drop_nans may not catch float64 nans
        # depending on polars version).
        ctrl = ctrl[np.isfinite(ctrl)]
        case = case[np.isfinite(case)]
        if ctrl.size == 0 and case.size == 0:
            ax.text(0.5, 0.5, "(no data)", ha="center", va="center",
                    transform=ax.transAxes)
            ax.axis("off")
            continue

        positions = [1, 2]
        data = [ctrl, case]
        parts = ax.violinplot(
            data, positions=positions, showmeans=False,
            showmedians=False, showextrema=False, widths=0.8,
        )
        for body, fill in zip(parts["bodies"], [color_ctrl, color_case]):
            body.set_facecolor(fill)
            body.set_edgecolor("black")
            body.set_alpha(0.85)
            body.set_linewidth(0.8)

        if show_box:
            ax.boxplot(
                data, positions=positions, widths=0.18, showfliers=False,
                patch_artist=True,
                medianprops=dict(color="black", linewidth=1.4),
                boxprops=dict(facecolor="white", edgecolor="black", linewidth=0.8),
                whiskerprops=dict(color="black", linewidth=0.8),
                capprops=dict(color="black", linewidth=0.8),
            )

        ax.set_xticks(positions)
        ax.set_xticklabels([ctrl_lbl, case_lbl])
        ax.set_ylim(0.0, 1.0)
        title_n = sub.height
        pretty_cat = str(cat).title() if str(cat) in ("hyper", "hypo") else str(cat)
        ax.set_title(f"{pretty_cat} DMRs (n={title_n})")
        if i == 0:
            ax.set_ylabel("Methylation level")

    fig.tight_layout()
    if save:
        _save_fig(md, fig, save)
    return fig, list(axes)


def _build_dmr_matrix(
    md: MethylData,
    dmr: pl.DataFrame,
    *,
    samples: list[str],
) -> tuple[np.ndarray, list[dict]]:
    """Return an (n_dmrs, n_samples) per-DMR mean beta matrix plus
    parallel metadata rows."""
    rows_meta: list[dict] = []
    rows_vals: list[np.ndarray] = []
    sample_to_col = {s: i for i, s in enumerate(samples)}
    for r in dmr.iter_rows(named=True):
        chrom = r["chrom"]
        start = int(r.get("start", r.get("pos", 0)))
        end = int(r.get("end", start + 1))
        beta_df = md.region_beta(chrom, start, end)
        vec = np.full(len(samples), np.nan, dtype=np.float64)
        for br in beta_df.iter_rows(named=True):
            ci = sample_to_col.get(br["sample"])
            if ci is not None:
                v = br.get("mean_beta")
                vec[ci] = float(v) if v is not None else np.nan
        rows_vals.append(vec)
        rows_meta.append(r)
    matrix = np.vstack(rows_vals) if rows_vals else np.zeros((0, len(samples)))
    return matrix, rows_meta


def dmr_heatmap(
    md: MethylData,
    *,
    n_top: int = 200,
    by: str = "qvalue",
    cluster_rows: bool = True,
    cluster_method: str = "average",
    cluster_metric: str = "euclidean",
    group_by: Optional[str] = "group",
    label_genes: int | bool = 20,
    cmap: str = "RdYlBu_r",
    vmin: float = 0.0,
    vmax: float = 1.0,
    figsize: tuple | None = None,
    save: str | None = None,
):
    """Per-DMR methylation heatmap with hierarchical row clustering.

    Rows are DMRs, columns are samples. Columns are reordered so samples
    sharing ``group_by`` end up adjacent (and a labelled bracket runs
    above each group). Rows are clustered via
    :func:`scipy.cluster.hierarchy.linkage` and a dendrogram is drawn on
    the left. The right margin gets a sparse set of gene-name labels so
    a viewer can read off the most-extreme DMRs without crowding.

    Parameters
    ----------
    n_top : int
        Number of DMRs to draw. Picked by ``by`` ascending (lowest qvalue
        first). Cap at a sane figure size with figsize=.
    by : str
        DMR-table column to rank by. Falls back gracefully when absent.
    cluster_rows : bool
        Reorder rows by hierarchical clustering. Off -> rows stay in
        ``by`` order.
    cluster_method, cluster_metric : str
        Forwarded to :func:`scipy.cluster.hierarchy.linkage`.
    group_by : str or None
        ``md.obs`` column used to (a) sort columns so a group's samples
        sit adjacent and (b) draw a labelled bracket / line at the top.
        ``None`` skips both.
    label_genes : int | bool
        Number of gene names to print along the right margin. Picks the
        DMRs with the largest |meth_diff|. ``True`` -> 20; ``False`` /
        ``0`` skips the column.
    cmap, vmin, vmax : str / float
        Heatmap colormap and clipping range. Defaults match the figure
        (RdYlBu_r, 0..1).
    """
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    dmr = _resolve_dmr(md)
    rank = by if by in dmr.columns else (
        "pvalue" if "pvalue" in dmr.columns else dmr.columns[0]
    )
    n = min(int(n_top), len(dmr))
    top = dmr.sort(rank).head(n)

    samples = md.obs.get_column("sample_id").to_list()

    # Sort columns by group so samples in the same condition are adjacent.
    if group_by and group_by in md.obs.columns:
        obs_sorted = md.obs.sort(group_by)
        samples_ordered = obs_sorted.get_column("sample_id").to_list()
        groups_ordered = obs_sorted.get_column(group_by).to_list()
    else:
        samples_ordered = samples
        groups_ordered = ["all"] * len(samples)

    matrix, rows_meta = _build_dmr_matrix(md, top, samples=samples_ordered)
    if matrix.size == 0:
        raise ValueError("DMR matrix is empty -- no rows to plot")

    # Row clustering on rows that have no NaN; rows with NaN are appended.
    valid_row_mask = ~np.isnan(matrix).any(axis=1)
    row_order = np.arange(matrix.shape[0])
    if cluster_rows and valid_row_mask.sum() > 2:
        from scipy.cluster.hierarchy import linkage, leaves_list
        # Cluster valid rows; keep NaN rows at the bottom in their existing
        # order so the dendrogram doesn't choke.
        valid_rows = np.where(valid_row_mask)[0]
        invalid_rows = np.where(~valid_row_mask)[0]
        Z = linkage(matrix[valid_rows], method=cluster_method, metric=cluster_metric)
        leaves = leaves_list(Z)
        row_order = np.concatenate([valid_rows[leaves], invalid_rows])
    matrix_ord = matrix[row_order]
    rows_meta_ord = [rows_meta[i] for i in row_order]

    # Figure layout: dendrogram + heatmap + (optional) gene labels + cbar.
    # Four columns -- the colorbar gets its own axis on the far right so it
    # never collides with the gene-label column.
    if figsize is None:
        figsize = (1.6 + 0.45 * len(samples_ordered) + (1.8 if label_genes else 0.4),
                   1.8 + 0.025 * matrix_ord.shape[0])
        figsize = (min(figsize[0], 17), min(max(figsize[1], 4.5), 14))
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(
        2, 4, figure=fig,
        width_ratios=[0.7, 4.0, (1.4 if label_genes else 0.05), 0.15],
        height_ratios=[0.35, 8.0],
        hspace=0.06, wspace=0.05,
    )
    ax_group = fig.add_subplot(gs[0, 1])
    ax_dend = fig.add_subplot(gs[1, 0])
    ax_heat = fig.add_subplot(gs[1, 1])
    ax_lbl = fig.add_subplot(gs[1, 2])
    ax_cbar = fig.add_subplot(gs[1, 3])

    # Dendrogram
    if cluster_rows and valid_row_mask.sum() > 2:
        from scipy.cluster.hierarchy import dendrogram
        dendrogram(
            Z, orientation="left", ax=ax_dend, no_labels=True,
            color_threshold=0, above_threshold_color="black",
        )
        ax_dend.invert_yaxis()  # match imshow row order
    ax_dend.set_xticks([])
    ax_dend.set_yticks([])
    for spine in ax_dend.spines.values():
        spine.set_visible(False)

    # Heatmap
    im = ax_heat.imshow(
        matrix_ord, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
        interpolation="nearest",
    )
    ax_heat.set_xticks(np.arange(len(samples_ordered)))
    ax_heat.set_xticklabels(samples_ordered, rotation=45, ha="right", fontsize=7)
    ax_heat.set_yticks([])
    ax_heat.set_xlabel("")
    ax_heat.set_ylabel("")
    # Row count goes in the title rather than the y-axis label, where it
    # would otherwise collide with the dendrogram column.
    ax_heat.set_title(f"{matrix_ord.shape[0]:,} DMRs (top by {rank})", fontsize=10)

    # Group bracket at the top
    if group_by and group_by in md.obs.columns:
        ax_group.set_xlim(-0.5, len(samples_ordered) - 0.5)
        ax_group.set_ylim(0, 1)
        ax_group.axis("off")
        # Find runs of identical groups
        i0 = 0
        for j in range(1, len(groups_ordered) + 1):
            if j == len(groups_ordered) or groups_ordered[j] != groups_ordered[i0]:
                mid = (i0 + j - 1) / 2.0
                ax_group.plot(
                    [i0 - 0.4, j - 1 + 0.4], [0.5, 0.5],
                    color="black", lw=1.4,
                )
                ax_group.text(
                    mid, 0.7, str(groups_ordered[i0]),
                    ha="center", va="bottom", fontsize=10, fontweight="bold",
                )
                i0 = j
    else:
        ax_group.axis("off")

    # Gene labels in the right margin
    if label_genes:
        n_labels = 20 if label_genes is True else int(label_genes)
        gene_col = (
            "gene_name" if "gene_name" in top.columns
            else "gene_id" if "gene_id" in top.columns
            else None
        )
        if gene_col is None or n_labels <= 0:
            ax_lbl.axis("off")
        else:
            # Pick by largest |meth_diff| (or any signed effect col).
            score_col = (
                "meth_diff" if "meth_diff" in top.columns
                else "mean_meth_diff" if "mean_meth_diff" in top.columns
                else None
            )
            if score_col:
                abs_diff = np.abs(np.fromiter(
                    (r.get(score_col) or 0.0 for r in rows_meta_ord),
                    count=len(rows_meta_ord), dtype=np.float64,
                ))
            else:
                abs_diff = np.arange(len(rows_meta_ord))[::-1]
            keep_idx = np.argsort(-abs_diff)[:n_labels]
            keep_idx_sorted = sorted(keep_idx.tolist())

            ax_lbl.set_xlim(0, 1)
            ax_lbl.set_ylim(matrix_ord.shape[0] - 0.5, -0.5)
            ax_lbl.axis("off")

            # Spread labels vertically so they don't overlap. Each label
            # occupies roughly `min_spacing` rows of vertical space; an
            # iterative two-pass sweep pushes the closer labels apart and
            # then re-centres so the displayed positions don't leak off
            # either edge of the axes.
            n_rows = matrix_ord.shape[0]
            # Reserve ~1 row of vertical space per label, scaled by how
            # crowded the figure is overall.
            min_spacing = max(1.0, n_rows / max(40.0, float(n_labels)))
            orig_y = [float(r) for r in keep_idx_sorted]
            display_y = list(orig_y)
            # Forward sweep: ensure top-down spacing.
            for i in range(1, len(display_y)):
                if display_y[i] - display_y[i - 1] < min_spacing:
                    display_y[i] = display_y[i - 1] + min_spacing
            # Backward sweep: if forward overshot the bottom, push earlier
            # labels upward to fit.
            for i in range(len(display_y) - 2, -1, -1):
                if display_y[i + 1] - display_y[i] < min_spacing:
                    display_y[i] = display_y[i + 1] - min_spacing
            # Final clip into axes bounds.
            top = -0.4
            bot = n_rows - 0.6
            for i, y in enumerate(display_y):
                display_y[i] = float(np.clip(y, top, bot))

            for row_idx, dy in zip(keep_idx_sorted, display_y):
                g = rows_meta_ord[row_idx].get(gene_col)
                if not g:
                    continue
                ax_lbl.annotate(
                    str(g),
                    xy=(0.0, row_idx),       # row in the heatmap
                    xytext=(0.22, dy),       # spread-out label position
                    fontsize=7, va="center", ha="left",
                    fontstyle="italic",
                    arrowprops=dict(
                        arrowstyle="-",
                        color="grey", lw=0.4, shrinkA=0, shrinkB=2,
                        connectionstyle="arc3,rad=0.0",
                    ),
                )
    else:
        ax_lbl.axis("off")

    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("Methylation level", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    if save:
        _save_fig(md, fig, save)
    return fig, {
        "heatmap": ax_heat, "dendrogram": ax_dend,
        "group_bracket": ax_group, "labels": ax_lbl, "colorbar": ax_cbar,
    }


__all__ = ["dmr_violin", "dmr_heatmap"]
