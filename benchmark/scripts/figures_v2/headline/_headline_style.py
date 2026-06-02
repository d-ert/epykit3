"""Shared style for the v1 headline figure set (benchmark/figures/summary/F0*).

Standalone — does not import the broken figures_v2/_style.py.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_STUDY1 = REPO_ROOT / "data" / "study1"
DATA_STUDY2 = REPO_ROOT / "data" / "study2"
DATA_STUDY3 = REPO_ROOT / "data" / "study3"
DATA_NULL = REPO_ROOT / "data" / "null_calibration"
FIG_OUT = REPO_ROOT / "figures" / "summary"

METHYLKIT_REAL = Path(
    "D:/Coding/Projeler/methyl_lib/methylkıt_realResults/scripts_and_results/methylkit_results"
)

PALETTE = {
    "epykit_lr": "#0F4C81",
    "epykit_lrplus": "#1F77B4",
    "epykit_lr+": "#1F77B4",
    "epykit_welch_t": "#2E86AB",
    "epykit_fisher": "#9BCDD2",
    "epykit_dmr_chain_merge": "#0F4C81",
    "epykit_dmr_tile": "#7FB3D5",
    "epykit_dmr_sliding_window": "#5DADE2",
    "epykit_dmr_segment": "#85C1E2",
    "methylkit": "#E07B39",
    "methylkit_tuned": "#D35400",
    "dss": "#2E8B57",
    "radmeth": "#9B59B6",
    "biseq": "#B59B6F",
    "methylsig": "#C0392B",
    "fisher": "#7F8C8D",
}

DISPLAY = {
    "epykit_lr": "epykit lr",
    "epykit_lrplus": "epykit lr+",
    "epykit_welch_t": "epykit welch_t",
    "epykit_fisher": "epykit fisher",
    "epykit_dmr_chain_merge": "epykit chain_merge",
    "epykit_dmr_tile": "epykit tile",
    "epykit_dmr_sliding_window": "epykit sliding_window",
    "epykit_dmr_segment": "epykit segment (HMM)",
    "methylkit": "methylKit",
    "methylkit_tuned": "methylKit (tuned)",
    "dss": "DSS",
    "radmeth": "RADMeth",
    "biseq": "BiSeq",
    "methylsig": "methylSig",
    "fisher": "Fisher (Piao)",
}

EPYKIT_TOOLS = {
    "epykit_lr", "epykit_lrplus", "epykit_welch_t", "epykit_fisher",
    "epykit_dmr_chain_merge", "epykit_dmr_tile",
    "epykit_dmr_sliding_window", "epykit_dmr_segment",
}


def color_for(tool: str) -> str:
    return PALETTE.get(tool, "#808080")


def label_for(tool: str) -> str:
    return DISPLAY.get(tool, tool)


def setup() -> None:
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "axes.titlesize":    11,
        "axes.titleweight":  "bold",
        "axes.labelsize":    10,
        "legend.fontsize":   9,
        "xtick.labelsize":   9,
        "ytick.labelsize":   9,
        "font.family":       "sans-serif",
        "font.sans-serif":   ["DejaVu Sans", "Arial", "Helvetica"],
        "pdf.fonttype":      42,
        "ps.fonttype":       42,
    })


def save_dual(fig, stem: str | Path) -> None:
    """Save both PNG and SVG. `stem` is a path without extension."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(stem) + ".png")
    fig.savefig(str(stem) + ".svg")
    plt.close(fig)
    print(f"  wrote {stem.name}.{{png,svg}}")
