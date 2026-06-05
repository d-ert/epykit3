"""Shared figure style + paths for figures_v2."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Self-locating: parents[2] of scripts/figures_v2/_style.py == benchmark/.
# (Was a hardcoded external benchmarkin_merges/FINAL_REPORT path pre-2026-06;
# repointed to the current committed epykit3/benchmark tree.)
BENCH_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR   = BENCH_ROOT / "data" / "study3"
# Post-rerun dis.merge sweep + multicore timing live under this sibling tree.
SWEEP_DIR  = BENCH_ROOT / "data" / "multi_thread_and_chain_sweep" / "chain_merge_dis_merge_sweep"
FIG_DIR    = BENCH_ROOT / "figures" / "study3_real_GSE263850"
THREE_WAY  = FIG_DIR / "three_way"
CM_FIG     = FIG_DIR / "chain_merge"
DSS_FIG    = FIG_DIR / "dss"

# External reference inputs that are not part of the committed benchmark tree
# (the source paper's Supp Table 5 and the methylKit-tile real-data run).
PAPER_T5_XLSX = Path(
    r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW/Paper resources/DMR_total_list.xlsx"
)
MK_TILE_DIR = Path(
    r"D:/Coding/Projeler/methyl_lib/methylkıt_realResults/scripts_and_results/methylkit_results"
)

# Color palette per caller, consistent across all figures.
CALLER_COLOR = {
    "paper":             "#2c3e50",   # dark slate — paper truth
    "paper-DSS (Supp Table 5)": "#2c3e50",
    "methylKit-tile":    "#e74c3c",   # red — outdated tile baseline
    "epykit-tile":       "#e67e22",   # orange — outdated epykit baseline
    "epykit-chain_merge-100": "#3498db",  # blue — paper-faithful chain_merge
    "epykit-chain_merge-250": "#1f77b4",  # darker blue — morphology-matched
    "DSS-from-scratch":  "#27ae60",   # green — DSS ceiling
}

# Feature-type color palette (HOMER-ish)
FEAT_COLOR = {
    "promoter-TSS": "#e74c3c",
    "5UTR":         "#f39c12",
    "exon":         "#f1c40f",
    "intron":       "#3498db",
    "3UTR":         "#9b59b6",
    "TTS":          "#e67e22",
    "non-coding":   "#95a5a6",
    "intergenic":   "#2c3e50",
}


def setup() -> None:
    """Apply a uniform matplotlib style across all v2 figures."""
    sns.set_theme(style="white", context="paper", font_scale=1.1)
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 180,
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
    })


def save_dual(fig, out_path: Path) -> None:
    """Save both PNG and SVG of a figure. out_path should be the stem
    (no extension)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path) + ".png")
    fig.savefig(str(out_path) + ".svg")
    print(f"  wrote {out_path.name}.{{png,svg}}")
