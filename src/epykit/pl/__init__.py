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
    "coverage_histogram",
    "methylation_heatmap",
    "mbias_plot",
    "volcano",
    "ma_plot",
    "manhattan",
    "genomic_context_bar",
    "cpg_island_pie",
    "karyogram",
    "pca",
    "tss_metaplot",
    "gene_body_metaplot",
    "umap",
    "sample_correlation",
    "qc_dashboard",
    "dmr_boxplot",
    "dmr_violin",
    "dmr_heatmap",
    "dmr_overlap",
    # annotatr-style plots
    "plot_annotation_counts",
    "plot_numerical_by_annotation",
    "plot_coannotations",
    "plot_categorical",
    # composer
    "figure_grid",
    # styling
    "apply_theme",
    "PALETTE",
    "set_palette",
]


def apply_theme(context: str = "paper") -> None:
    """Re-apply the matplotlib theme. See :func:`epykit._style.apply_theme`."""
    return _apply_theme(context)
