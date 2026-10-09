from __future__ import annotations

from .._style import PALETTE, set_palette
from .._style import apply_theme as _apply_theme
from ._plotly import (
    coverage_histogram_plotly,
    cpg_island_pie_plotly,
    dmr_size_hist_plotly,
    feature_direction_stacked_plotly,
    feature_pie_plotly,
    global_methylation_bar_plotly,
    ma_plot_plotly,
    manhattan_plotly,
    pca_plotly,
    pvalue_histogram_plotly,
    sample_correlation_plotly,
    scree_plotly,
    tss_metaplot_plotly,
    volcano_plotly,
)

# Interactive Plotly aliases for standard plotting functions
volcano = volcano_plotly
ma_plot = ma_plot_plotly
manhattan = manhattan_plotly
coverage_histogram = coverage_histogram_plotly
pca = pca_plotly
feature_pie = feature_pie_plotly
cpg_island_pie = cpg_island_pie_plotly
tss_metaplot = tss_metaplot_plotly
pvalue_histogram = pvalue_histogram_plotly
dmr_size_hist = dmr_size_hist_plotly
global_methylation_bar = global_methylation_bar_plotly
sample_correlation = sample_correlation_plotly
scree = scree_plotly
feature_direction_stacked = feature_direction_stacked_plotly

__all__ = [
    "PALETTE",
    "apply_theme",
    "coverage_histogram",
    "cpg_island_pie",
    "dmr_size_hist",
    "feature_direction_stacked",
    "feature_pie",
    "global_methylation_bar",
    "ma_plot",
    "manhattan",
    "pca",
    "pvalue_histogram",
    "sample_correlation",
    "scree",
    "set_palette",
    "tss_metaplot",
    "volcano",
]


def apply_theme(context: str = "paper") -> None:
    """Re-apply the plot theme. See :func:`epykit._style.apply_theme`."""
    return _apply_theme(context)
