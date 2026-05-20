"""Interactive Plotly counterparts of the matplotlib plots used by the HTML
report.

These are thin render adapters around :mod:`epykit.pl._compute` -- the
same compute outputs feed both these and the matplotlib plots in
``pl/*.py``, so a fix in a compute function flows to both backends.

Plotly is an optional dependency. Functions import it lazily and raise a
friendly ImportError when the user lacks the report extras.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import polars as pl

from .._style import PALETTE
from ._compute import (
    compute_pca,
    compute_tss_metaplot,
    compute_coverage_distribution,
    compute_volcano_data,
    compute_ma_data,
    compute_manhattan_data,
    compute_annotation_counts,
)
from ..methyldata import MethylData


def _require_plotly():
    try:
        import plotly.graph_objects as go
        return go
    except ImportError as exc:
        raise ImportError(
            "plotly is required for interactive report figures. "
            "Install with: pip install 'epykit[report]'"
        ) from exc


def volcano_plotly(md: MethylData, *, alpha: float = 0.05, min_abs_diff: float = 0.1):
    go = _require_plotly()
    data = compute_volcano_data(md, alpha=alpha, min_abs_diff=min_abs_diff)
    ns = ~data.sig
    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=data.meth_diff[ns], y=data.neg_log_p[ns], mode="markers",
        marker=dict(size=4, color=PALETTE["neutral"], opacity=0.4),
        name="ns", hoverinfo="skip",
    ))
    fig.add_trace(go.Scattergl(
        x=data.meth_diff[data.hypo], y=data.neg_log_p[data.hypo], mode="markers",
        marker=dict(size=5, color=PALETTE["hypo"], opacity=0.7),
        name=f"hypo ({int(data.hypo.sum())})",
    ))
    fig.add_trace(go.Scattergl(
        x=data.meth_diff[data.hyper], y=data.neg_log_p[data.hyper], mode="markers",
        marker=dict(size=5, color=PALETTE["hyper"], opacity=0.7),
        name=f"hyper ({int(data.hyper.sum())})",
    ))
    fig.add_hline(y=-np.log10(alpha), line=dict(color="grey", dash="dash", width=1))
    fig.add_vline(x=min_abs_diff, line=dict(color="grey", dash="dash", width=1))
    fig.add_vline(x=-min_abs_diff, line=dict(color="grey", dash="dash", width=1))
    fig.update_layout(
        title="DMC volcano",
        xaxis_title="Methylation difference (treatment - control)",
        yaxis_title=f"-log_1_0({data.p_col})",
        template="simple_white", height=420,
    )
    return fig


def ma_plot_plotly(md: MethylData, *, alpha: float = 0.05, min_abs_diff: float = 0.1):
    go = _require_plotly()
    data = compute_ma_data(md, alpha=alpha, min_abs_diff=min_abs_diff)
    ns = ~data.sig
    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=data.mean_beta[ns], y=data.meth_diff[ns], mode="markers",
        marker=dict(size=4, color=PALETTE["neutral"], opacity=0.4),
        name="ns", hoverinfo="skip",
    ))
    fig.add_trace(go.Scattergl(
        x=data.mean_beta[data.hypo], y=data.meth_diff[data.hypo], mode="markers",
        marker=dict(size=5, color=PALETTE["hypo"], opacity=0.7),
        name=f"hypo ({int(data.hypo.sum())})",
    ))
    fig.add_trace(go.Scattergl(
        x=data.mean_beta[data.hyper], y=data.meth_diff[data.hyper], mode="markers",
        marker=dict(size=5, color=PALETTE["hyper"], opacity=0.7),
        name=f"hyper ({int(data.hyper.sum())})",
    ))
    fig.add_hline(y=0, line=dict(color="black", width=1))
    fig.add_hline(y=min_abs_diff, line=dict(color="grey", dash="dash", width=1))
    fig.add_hline(y=-min_abs_diff, line=dict(color="grey", dash="dash", width=1))
    fig.update_layout(
        title="MA plot",
        xaxis_title="Mean methylation",
        yaxis_title="Methylation difference (treatment - control)",
        template="simple_white", height=420,
    )
    return fig


def manhattan_plotly(md: MethylData, *, alpha: float = 0.05):
    go = _require_plotly()
    data = compute_manhattan_data(md, alpha=alpha)
    fig = go.Figure()
    colors = [PALETTE["hypo"], PALETTE["hyper"]]
    for i, block in enumerate(data.chrom_blocks):
        fig.add_trace(go.Scattergl(
            x=block["x"], y=block["y"], mode="markers",
            marker=dict(size=3, color=colors[i % 2], opacity=0.7),
            name=block["chrom"], hoverinfo="skip", showlegend=False,
        ))
    fig.add_hline(y=data.alpha_line_y, line=dict(color="red", dash="dash"))
    fig.update_layout(
        title="Manhattan plot",
        xaxis=dict(tickvals=data.tick_pos, ticktext=data.tick_label, title="Chromosome"),
        yaxis_title=f"-log_1_0({data.p_col})",
        template="simple_white", height=320,
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
        template="simple_white", height=320,
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
        PALETTE["control"], PALETTE["treatment"],
        PALETTE["island"], PALETTE["shelf"], PALETTE["neutral"],
    ]
    for i, g in enumerate(unique):
        mask = np.array([gg == g for gg in res.groups])
        fig.add_trace(go.Scatter(
            x=res.coords[mask, 0], y=res.coords[mask, 1],
            mode="markers+text",
            text=[s for s, m in zip(res.samples, mask) if m],
            textposition="top center",
            marker=dict(size=12, color=palette_cycle[i % len(palette_cycle)]),
            name=str(g),
        ))
    fig.update_layout(
        title=f"PCA  |  n_sites={res.n_sites_used:,}",
        xaxis_title=f"PC1 ({res.explained_var[0]:.1%})",
        yaxis_title=f"PC2 ({res.explained_var[1]:.1%})",
        template="simple_white", height=420,
        legend_title=res.group_col,
    )
    return fig


def feature_pie_plotly(md: MethylData):
    """Pie chart of feature_type distribution on the annotated DMC table."""
    go = _require_plotly()
    dmc = md.dmc
    if dmc is None or "feature_type" not in dmc.columns:
        return None
    counts = compute_annotation_counts(dmc, annot_col="feature_type")
    fig = go.Figure(data=[go.Pie(
        labels=counts["feature_type"].to_list(),
        values=counts["count"].to_list(),
        hole=0.35,
    )])
    fig.update_layout(title="Genomic context (feature_type)", template="simple_white", height=360)
    return fig


def cpg_island_pie_plotly(md: MethylData):
    go = _require_plotly()
    dmc = md.dmc
    if dmc is None or "cpg_context" not in dmc.columns:
        return None
    counts = compute_annotation_counts(dmc, annot_col="cpg_context")
    fig = go.Figure(data=[go.Pie(
        labels=counts["cpg_context"].to_list(),
        values=counts["count"].to_list(),
        hole=0.35,
    )])
    fig.update_layout(title="CpG-island context", template="simple_white", height=360)
    return fig


def tss_metaplot_plotly(
    md: MethylData, gtf_path: str, *,
    window_bp: int = 2000, n_bins: int = 100,
    group_by: Optional[str] = "group", max_genes: Optional[int] = None,
):
    """Plotly-rendered version of :func:`pl.tss_metaplot`."""
    go = _require_plotly()
    res = compute_tss_metaplot(
        md, gtf_path,
        window_bp=window_bp, n_bins=n_bins,
        group_by=group_by, max_genes=max_genes,
    )
    samples = res.samples
    mean_beta = res.mean_beta
    x = res.x

    fig = go.Figure()
    palette_cycle = [
        PALETTE["control"], PALETTE["treatment"],
        PALETTE["island"], PALETTE["shelf"],
    ]
    if res.group_col is not None:
        unique = sorted(set(res.groups))
        for i, samp in enumerate(samples):
            fig.add_trace(go.Scatter(
                x=x, y=mean_beta[i], mode="lines",
                line=dict(
                    color=palette_cycle[unique.index(res.groups[i]) % len(palette_cycle)],
                    width=1,
                ),
                opacity=0.25, name=str(samp), showlegend=False,
            ))
        for j, g in enumerate(unique):
            mask = np.array([gg == g for gg in res.groups])
            if not mask.any():
                continue
            mean = np.nanmean(mean_beta[mask], axis=0)
            fig.add_trace(go.Scatter(
                x=x, y=mean, mode="lines",
                line=dict(color=palette_cycle[j % len(palette_cycle)], width=2.5),
                name=str(g),
            ))
    else:
        for i, samp in enumerate(samples):
            fig.add_trace(go.Scatter(
                x=x, y=mean_beta[i], mode="lines", name=str(samp),
            ))
    fig.add_vline(x=0, line=dict(color="black", dash="dash", width=1))
    fig.update_layout(
        title=f"TSS metaplot (+/-{window_bp} bp)",
        xaxis_title="Distance from TSS (bp)",
        yaxis_title="Mean beta",
        template="simple_white", height=380,
    )
    return fig


__all__ = [
    "volcano_plotly",
    "ma_plot_plotly",
    "manhattan_plotly",
    "coverage_histogram_plotly",
    "pca_plotly",
    "feature_pie_plotly",
    "cpg_island_pie_plotly",
    "tss_metaplot_plotly",
]
