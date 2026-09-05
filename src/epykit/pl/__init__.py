from __future__ import annotations

from .._style import PALETTE, set_palette
from .._style import apply_theme as _apply_theme

_apply_theme()

from .annotation import (
    plot_annotation_counts,
    plot_categorical,
    plot_coannotations,
    plot_numerical_by_annotation,
)
from .clustering import pca
from .composer import figure_grid
from .correlation import sample_correlation
from .dashboard import qc_dashboard
from .differential import ma_plot, manhattan, volcano
from .dmr_boxplot import dmr_boxplot
from .dmr_summary import dmr_heatmap, dmr_violin
from .embedding import umap
from .genomic import cpg_island_pie, genomic_context_bar, karyogram
from .metaplot import gene_body_metaplot, tss_metaplot
from .overlap import dmr_overlap
from .qc import coverage_histogram, mbias_plot, methylation_heatmap

__all__ = [
    "PALETTE",
    "apply_theme",
    "coverage_histogram",
    "cpg_island_pie",
    "dmr_boxplot",
    "dmr_heatmap",
    "dmr_overlap",
    "dmr_violin",
    "figure_grid",
    "gene_body_metaplot",
    "genomic_context_bar",
    "karyogram",
    "ma_plot",
    "manhattan",
    "mbias_plot",
    "methylation_heatmap",
    "pca",
    "plot_annotation_counts",
    "plot_categorical",
    "plot_coannotations",
    "plot_numerical_by_annotation",
    "qc_dashboard",
    "sample_correlation",
    "set_palette",
    "tss_metaplot",
    "umap",
    "volcano",
]


def apply_theme(context: str = "paper") -> None:
    """Re-apply the matplotlib theme. See :func:`epykit._style.apply_theme`."""
    return _apply_theme(context)
