"""Interactive Plotly counterparts of the matplotlib plots used by the HTML
report.

These are thin render adapters around :mod:`epykit.pl._compute` -- the
same compute outputs feed both these and the matplotlib plots in
``pl/*.py``, so a fix in a compute function flows to both backends.

Plotly is an optional dependency. Functions import it lazily and raise a
friendly ImportError when the user lacks the report extras.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from .._style import PALETTE
from ..methyldata import MethylData
from ._compute import (
    compute_annotation_counts,
    compute_categorical_proportions,
    compute_coverage_distribution,
    compute_dmr_size_distribution,
    compute_global_methylation,
    compute_ma_data,
    compute_manhattan_data,
    compute_pca,
    compute_pvalue_histogram,
    compute_sample_correlation_matrix,
    compute_scree,
    compute_tss_metaplot,
    compute_volcano_data,
)

# Refined dashboard accents. These affect the interactive *report* figures
# only; the shared matplotlib ``_style.PALETTE`` is intentionally left
# untouched so the static plots in ``pl/*.py`` are unchanged.
_ACCENT = "#2563eb"
_VIOLET = "#7c3aed"


def _dash_layout(go, **over):
    """Shared white-background layout for dashboard report figures.

    Figures live on white plot-cards in both light and dark report themes,
    so they never need re-theming when the user toggles the theme.
    """
    base = dict(
        paper_bgcolor="white",
        plot_bgcolor="white",
        font=dict(
            family="system-ui,-apple-system,Segoe UI,Roboto,sans-serif",
            size=12,
            color="#1a2230",
        ),
        margin=dict(l=54, r=18, t=18, b=44),
        template="simple_white",
    )
    base.update(over)
    return base


def _require_plotly():
    try:
        import plotly.graph_objects as go

        return go
    except ImportError as exc:
        raise ImportError(
            "plotly is required for interactive report figures. "
            "Install with: pip install 'epykit[report]'"
        ) from exc


def volcano_plotly(
    md: MethylData, *, alpha: float = 0.05, min_abs_diff: float = 0.1, dmc=None, max_points=None
):
    go = _require_plotly()
    data = compute_volcano_data(
        md, alpha=alpha, min_abs_diff=min_abs_diff, dmc=dmc, max_points=max_points
    )
    ns = ~data.sig
    fig = go.Figure()
    fig.add_trace(
        go.Scattergl(
            x=data.meth_diff[ns],
            y=data.neg_log_p[ns],
            mode="markers",
            marker=dict(size=4, color=PALETTE["neutral"], opacity=0.4),
            name="ns",
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scattergl(
            x=data.meth_diff[data.hypo],
            y=data.neg_log_p[data.hypo],
            mode="markers",
            marker=dict(size=5, color=PALETTE["hypo"], opacity=0.7),
            name=f"hypo ({int(data.hypo.sum())})",
        )
    )
    fig.add_trace(
        go.Scattergl(
            x=data.meth_diff[data.hyper],
            y=data.neg_log_p[data.hyper],
            mode="markers",
            marker=dict(size=5, color=PALETTE["hyper"], opacity=0.7),
            name=f"hyper ({int(data.hyper.sum())})",
        )
    )
    fig.add_hline(y=-np.log10(alpha), line=dict(color="grey", dash="dash", width=1))
    fig.add_vline(x=min_abs_diff, line=dict(color="grey", dash="dash", width=1))
    fig.add_vline(x=-min_abs_diff, line=dict(color="grey", dash="dash", width=1))
    fig.update_layout(
        title="DMC volcano",
        xaxis_title="Methylation difference (treatment - control)",
        yaxis_title=f"-log_1_0({data.p_col})",
        template="simple_white",
        height=420,
    )
    return fig


def ma_plot_plotly(
    md: MethylData, *, alpha: float = 0.05, min_abs_diff: float = 0.1, dmc=None, max_points=None
):
    go = _require_plotly()
    data = compute_ma_data(
        md, alpha=alpha, min_abs_diff=min_abs_diff, dmc=dmc, max_points=max_points
    )
    ns = ~data.sig
    fig = go.Figure()
    fig.add_trace(
        go.Scattergl(
            x=data.mean_beta[ns],
            y=data.meth_diff[ns],
            mode="markers",
            marker=dict(size=4, color=PALETTE["neutral"], opacity=0.4),
            name="ns",
            hoverinfo="skip",
        )
    )
    fig.add_trace(
        go.Scattergl(
            x=data.mean_beta[data.hypo],
            y=data.meth_diff[data.hypo],
            mode="markers",
            marker=dict(size=5, color=PALETTE["hypo"], opacity=0.7),
            name=f"hypo ({int(data.hypo.sum())})",
        )
    )
    fig.add_trace(
        go.Scattergl(
            x=data.mean_beta[data.hyper],
            y=data.meth_diff[data.hyper],
            mode="markers",
            marker=dict(size=5, color=PALETTE["hyper"], opacity=0.7),
            name=f"hyper ({int(data.hyper.sum())})",
        )
    )
    fig.add_hline(y=0, line=dict(color="black", width=1))
    fig.add_hline(y=min_abs_diff, line=dict(color="grey", dash="dash", width=1))
    fig.add_hline(y=-min_abs_diff, line=dict(color="grey", dash="dash", width=1))
    fig.update_layout(
        title="MA plot",
        xaxis_title="Mean methylation",
        yaxis_title="Methylation difference (treatment - control)",
        template="simple_white",
        height=420,
    )
    return fig


def manhattan_plotly(md: MethylData, *, alpha: float = 0.05, dmc=None, max_points=None):
    go = _require_plotly()
    data = compute_manhattan_data(md, alpha=alpha, dmc=dmc, max_points=max_points)
    fig = go.Figure()
    colors = [PALETTE["hypo"], PALETTE["hyper"]]
    for i, block in enumerate(data.chrom_blocks):
        fig.add_trace(
            go.Scattergl(
                x=block["x"],
                y=block["y"],
                mode="markers",
                marker=dict(size=3, color=colors[i % 2], opacity=0.7),
                name=block["chrom"],
                hoverinfo="skip",
                showlegend=False,
            )
        )
    fig.add_hline(y=data.alpha_line_y, line=dict(color="red", dash="dash"))
    fig.update_layout(
        title="Manhattan plot",
        xaxis=dict(tickvals=data.tick_pos, ticktext=data.tick_label, title="Chromosome"),
        yaxis_title=f"-log_1_0({data.p_col})",
        template="simple_white",
        height=320,
    )
    return fig


def coverage_histogram_plotly(md: MethylData, *, bins: int = 100, max_points: int = 200_000):
    go = _require_plotly()
    cov = compute_coverage_distribution(md, max_points=max_points)
    fig = go.Figure()
    fig.add_trace(go.Histogram(x=cov, nbinsx=bins, marker_color=PALETTE["neutral"]))
    fig.update_layout(
        title="Coverage histogram",
        xaxis_title="Coverage",
        yaxis_title="CpG count",
        template="simple_white",
        height=320,
    )
    return fig


def pca_plotly(md: MethylData, *, n_sites: int = 10_000):
    go = _require_plotly()
    try:
        res = compute_pca(md, n_sites=n_sites)
    except (ImportError, ValueError):
        return None

    fig = go.Figure()
    unique = list(dict.fromkeys(res.groups))
    palette_cycle = [
        PALETTE["control"],
        PALETTE["treatment"],
        PALETTE["island"],
        PALETTE["shelf"],
        PALETTE["neutral"],
    ]
    for i, g in enumerate(unique):
        mask = np.array([gg == g for gg in res.groups])
        fig.add_trace(
            go.Scatter(
                x=res.coords[mask, 0],
                y=res.coords[mask, 1],
                mode="markers+text",
                text=[s for s, m in zip(res.samples, mask) if m],
                textposition="top center",
                marker=dict(size=12, color=palette_cycle[i % len(palette_cycle)]),
                name=str(g),
            )
        )
    fig.update_layout(
        title=f"PCA  |  n_sites={res.n_sites_used:,}",
        xaxis_title=f"PC1 ({res.explained_var[0]:.1%})",
        yaxis_title=f"PC2 ({res.explained_var[1]:.1%})",
        template="simple_white",
        height=420,
        legend_title=res.group_col,
    )
    return fig


# Human-readable suffix per annotation level, so a standalone figure stays
# self-describing even outside the report's surrounding captions.
_LEVEL_LABEL = {"dmc": "per-CpG (DMC)", "dmr": "per-region (DMR)"}


def _annot_level_table(md: MethylData, level: str):
    """Resolve the annotated table for an annotation pie, or ``None``.

    ``level="dmc"`` -> the (annotated) per-CpG table on ``md.dmc``;
    ``level="dmr"`` -> the per-region table on ``md.uns['dmr']``. Returns
    ``None`` (rather than raising) when the table is absent, so the report
    renders its graceful skip-state instead of erroring.
    """
    lvl = level.lower()
    if lvl == "dmc":
        return md.dmc
    if lvl == "dmr":
        dmr = md.uns.get("dmr")
        if isinstance(dmr, pl.DataFrame) and not dmr.is_empty():
            return dmr
        return None
    raise ValueError(f"level must be 'dmc' or 'dmr', got {level!r}")


def feature_pie_plotly(md: MethylData, *, level: str = "dmc"):
    """Pie of gene-feature distribution on the annotated DMC or DMR table.

    ``level="dmc"`` (default) weights by differential-CpG density -- "where do
    differential cytosines fall?". ``level="dmr"`` gives the region-level
    fraction -- "what fraction of DMRs hit promoters?" (the field-standard
    ChIPseeker / DSS / dmrseq view). Returns ``None`` when the chosen table or
    its ``feature_type`` column is absent.
    """
    go = _require_plotly()
    tbl = _annot_level_table(md, level)
    if tbl is None or "feature_type" not in tbl.columns:
        return None
    counts = compute_annotation_counts(tbl, annot_col="feature_type")
    fig = go.Figure(
        data=[
            go.Pie(
                labels=counts["feature_type"].to_list(),
                values=counts["count"].to_list(),
                hole=0.35,
            )
        ]
    )
    fig.update_layout(
        title=f"Gene features · {_LEVEL_LABEL.get(level.lower(), level)}",
        template="simple_white",
        height=360,
    )
    return fig


def cpg_island_pie_plotly(md: MethylData, *, level: str = "dmc"):
    """Pie of CpG-island context on the annotated DMC or DMR table.

    See :func:`feature_pie_plotly` for the ``level`` semantics. Returns
    ``None`` when the chosen table or its ``cpg_context`` column is absent.
    """
    go = _require_plotly()
    tbl = _annot_level_table(md, level)
    if tbl is None or "cpg_context" not in tbl.columns:
        return None
    counts = compute_annotation_counts(tbl, annot_col="cpg_context")
    fig = go.Figure(
        data=[
            go.Pie(
                labels=counts["cpg_context"].to_list(),
                values=counts["count"].to_list(),
                hole=0.35,
            )
        ]
    )
    fig.update_layout(
        title=f"CpG-island context · {_LEVEL_LABEL.get(level.lower(), level)}",
        template="simple_white",
        height=360,
    )
    return fig


def tss_metaplot_plotly(
    md: MethylData,
    gtf_path: str,
    *,
    window_bp: int = 2000,
    n_bins: int = 100,
    group_by: str | None = "group",
    max_genes: int | None = None,
):
    """Plotly-rendered version of :func:`pl.tss_metaplot`."""
    go = _require_plotly()
    res = compute_tss_metaplot(
        md,
        gtf_path,
        window_bp=window_bp,
        n_bins=n_bins,
        group_by=group_by,
        max_genes=max_genes,
    )
    samples = res.samples
    mean_beta = res.mean_beta
    x = res.x
    # No CpGs fell in any TSS window (e.g. GTF/store chromosome mismatch, or the
    # store path could not be resolved) -> nothing to draw. Signal "skip" so the
    # report renders its notice instead of an empty axes.
    if mean_beta.size == 0 or not np.isfinite(mean_beta).any():
        return None

    fig = go.Figure()
    palette_cycle = [
        PALETTE["control"],
        PALETTE["treatment"],
        PALETTE["island"],
        PALETTE["shelf"],
    ]
    if res.group_col is not None:
        unique = sorted(set(res.groups))
        for i, samp in enumerate(samples):
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=mean_beta[i],
                    mode="lines",
                    line=dict(
                        color=palette_cycle[unique.index(res.groups[i]) % len(palette_cycle)],
                        width=1,
                    ),
                    opacity=0.25,
                    name=str(samp),
                    showlegend=False,
                )
            )
        for j, g in enumerate(unique):
            mask = np.array([gg == g for gg in res.groups])
            if not mask.any():
                continue
            mean = np.nanmean(mean_beta[mask], axis=0)
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=mean,
                    mode="lines",
                    line=dict(color=palette_cycle[j % len(palette_cycle)], width=2.5),
                    name=str(g),
                )
            )
    else:
        for i, samp in enumerate(samples):
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=mean_beta[i],
                    mode="lines",
                    name=str(samp),
                )
            )
    fig.add_vline(x=0, line=dict(color="black", dash="dash", width=1))
    fig.update_layout(
        title=f"TSS metaplot (+/-{window_bp} bp)",
        xaxis_title="Distance from TSS (bp)",
        yaxis_title="Mean beta",
        template="simple_white",
        height=380,
    )
    return fig


def pvalue_histogram_plotly(md: MethylData, *, dmc=None):
    """Histogram of raw per-CpG p-values -- a calibration diagnostic."""
    go = _require_plotly()
    counts, edges = compute_pvalue_histogram(md, bins=30, dmc=dmc)
    centers = (edges[:-1] + edges[1:]) / 2.0
    width = float(edges[1] - edges[0]) * 0.95
    fig = go.Figure(
        [
            go.Bar(
                x=centers,
                y=counts,
                width=width,
                marker_color=_ACCENT,
                opacity=0.85,
            )
        ]
    )
    fig.update_layout(
        **_dash_layout(
            go,
            title="p-value histogram",
            xaxis_title="p-value",
            yaxis_title="CpG count",
            height=320,
            bargap=0.02,
        )
    )
    return fig


def dmr_size_hist_plotly(md: MethylData):
    """Distribution of CpGs-per-DMR."""
    go = _require_plotly()
    sizes = compute_dmr_size_distribution(md)
    fig = go.Figure(
        [
            go.Histogram(
                x=sizes,
                marker_color=_VIOLET,
                opacity=0.85,
            )
        ]
    )
    fig.update_layout(
        **_dash_layout(
            go,
            title="DMR size distribution",
            xaxis_title="CpGs per DMR",
            yaxis_title="DMR count",
            height=300,
            bargap=0.05,
        )
    )
    return fig


def global_methylation_bar_plotly(md: MethylData):
    """Per-sample global methylation, coloured by group (batch eyeball)."""
    go = _require_plotly()
    samples, values, groups = compute_global_methylation(md)
    uniq = list(dict.fromkeys(groups))
    palette_cycle = [
        PALETTE["treatment"],
        PALETTE["control"],
        PALETTE["island"],
        PALETTE["shelf"],
        PALETTE["neutral"],
    ]
    cmap = {g: palette_cycle[i % len(palette_cycle)] for i, g in enumerate(uniq)}
    colors = [cmap.get(g, PALETTE["neutral"]) for g in groups]
    fig = go.Figure([go.Bar(x=samples, y=values, marker_color=colors)])
    fig.update_layout(
        **_dash_layout(
            go,
            title="Global methylation per sample",
            yaxis_title="mean beta",
            height=300,
        )
    )
    return fig


def sample_correlation_plotly(md: MethylData):
    """Clustered all-vs-all sample correlation heatmap."""
    go = _require_plotly()
    mat, labels = compute_sample_correlation_matrix(md)
    zmin = float(min(0.7, np.nanmin(mat)))
    fig = go.Figure(
        [
            go.Heatmap(
                z=mat,
                x=labels,
                y=labels,
                colorscale=[[0, "#ffffff"], [1, _ACCENT]],
                zmin=zmin,
                zmax=1.0,
                colorbar=dict(thickness=12, len=0.75),
            )
        ]
    )
    fig.update_layout(
        **_dash_layout(
            go,
            title="Sample correlation",
            height=360,
            margin=dict(l=64, r=20, t=32, b=54),
        )
    )
    return fig


def scree_plotly(md: MethylData):
    """Variance explained by the leading principal components."""
    go = _require_plotly()
    ev = compute_scree(md, max_components=6) * 100.0
    labels = [f"PC{i + 1}" for i in range(len(ev))]
    fig = go.Figure([go.Bar(x=labels, y=ev, marker_color=_ACCENT)])
    fig.update_layout(
        **_dash_layout(
            go,
            title="Scree",
            yaxis_title="% variance explained",
            height=340,
        )
    )
    return fig


def feature_direction_stacked_plotly(md: MethylData):
    """Hyper- vs hypo-methylation proportion within each gene feature.

    Returns ``None`` when the DMC table is unannotated.
    """
    go = _require_plotly()
    dmc = md.dmc
    if dmc is None or "feature_type" not in dmc.columns or "meth_diff" not in dmc.columns:
        return None
    import polars as _pl

    work = dmc.with_columns(
        _pl.when(_pl.col("meth_diff") > 0)
        .then(_pl.lit("hyper"))
        .otherwise(_pl.lit("hypo"))
        .alias("dmr_type")
    )
    prop = compute_categorical_proportions(
        work,
        group_col="dmr_type",
        annot_col="feature_type",
        include_all_group=False,
        normalize=True,
    )
    feats = prop["feature_type"].unique(maintain_order=True).to_list()
    fig = go.Figure()
    for direction, color in (("hyper", PALETTE["hyper"]), ("hypo", PALETTE["hypo"])):
        sub = prop.filter(_pl.col("dmr_type") == direction)
        ymap = dict(zip(sub["feature_type"].to_list(), sub["proportion"].to_list()))
        fig.add_trace(
            go.Bar(
                x=feats,
                y=[ymap.get(f, 0.0) for f in feats],
                name=direction,
                marker_color=color,
            )
        )
    fig.update_layout(
        **_dash_layout(
            go,
            title="Hyper vs hypo by feature",
            barmode="stack",
            yaxis_title="proportion",
            height=300,
            showlegend=True,
            legend=dict(orientation="h"),
        )
    )
    return fig


__all__ = [
    "coverage_histogram_plotly",
    "cpg_island_pie_plotly",
    "dmr_size_hist_plotly",
    "feature_direction_stacked_plotly",
    "feature_pie_plotly",
    "global_methylation_bar_plotly",
    "ma_plot_plotly",
    "manhattan_plotly",
    "pca_plotly",
    "pvalue_histogram_plotly",
    "sample_correlation_plotly",
    "scree_plotly",
    "tss_metaplot_plotly",
    "volcano_plotly",
]
