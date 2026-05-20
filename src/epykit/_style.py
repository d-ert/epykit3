from __future__ import annotations

# Default palette. Bind dynamically (plots read ``_style.PALETTE`` at call
# time, not import time) so ``set_palette`` can swap it without forcing
# a re-import.
PALETTE: dict = {
    "hyper": "#e05263",
    "hypo": "#4a90d9",
    "neutral": "#aaaaaa",
    "island": "#2ca02c",
    "shore": "#98df8a",
    "shelf": "#dbf9db",
    "open_sea": "#d3d3d3",
    "treatment": "#e05263",
    "control": "#4a90d9",

    # Okabe-Ito derived: colour-blind safe across the eight common types
    # of CVD; common choice for genomics figures.
    "promoter": "#0072B2",
    "5utr": "#56B4E9",
    "exon": "#009E73",
    "intron": "#F0E442",
    "3utr": "#E69F00",
    "intergenic": "#999999",
    "enhancer": "#CC79A7",
    "lncRNA": "#7B3294",
    "insulator": "#404040",
}

# Named palettes; switch with ``set_palette("publication")`` etc.
_PALETTES: dict[str, dict] = {
    "default": dict(PALETTE),
    "publication": {
        **PALETTE,
        # Higher-contrast hyper/hypo for print.
        "hyper": "#c0392b",
        "hypo": "#2c5aa0",
        "neutral": "#7f7f7f",
        "treatment": "#c0392b",
        "control": "#2c5aa0",
    },
    "colorblind": {
        **PALETTE,
        "hyper": "#D55E00",
        "hypo": "#0072B2",
        "treatment": "#D55E00",
        "control": "#0072B2",
    },
}


def set_palette(name: str) -> None:
    """Switch the active palette in-place.

    Plots read ``_style.PALETTE`` dynamically, so this affects every plot
    drawn after the call. Available presets: ``"default"``,
    ``"publication"``, ``"colorblind"``. Use a fresh ``dict.update``
    call for per-session custom palettes:

        >>> from epykit._style import PALETTE
        >>> PALETTE.update(hyper="#440154", hypo="#fde725")
    """
    if name not in _PALETTES:
        raise ValueError(
            f"unknown palette {name!r}; choose from {sorted(_PALETTES)}"
        )
    PALETTE.clear()
    PALETTE.update(_PALETTES[name])


def apply_theme(context: str = "paper") -> None:
    """Apply a matplotlib rcParams preset.

    Contexts
    --------
    paper
        Defaults: 150 dpi, sans-serif, small labels. Matches notebook
        rendering and the historic epykit baseline.
    publication
        300 dpi raster + vector-friendly settings (TrueType fonts embedded
        in PDF/PS, ``svg.fonttype="none"`` so text stays selectable).
        Slightly larger labels. Use once at the top of a notebook before
        drawing figures destined for a manuscript.
    talk
        Larger labels (~14pt) and lines, suitable for slide projection.
    poster
        Biggest fonts, thickest lines.
    """
    try:
        import matplotlib as mpl
    except Exception:
        return

    if context == "paper":
        rc = {
            "figure.dpi": 150,
            "figure.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "font.family": "sans-serif",
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    elif context == "publication":
        rc = {
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "figure.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
            "font.family": "sans-serif",
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "legend.fontsize": 10,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    elif context == "talk":
        rc = {
            "figure.dpi": 120,
            "figure.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.0,
            "axes.labelsize": 14,
            "axes.titlesize": 15,
            "font.family": "sans-serif",
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "lines.linewidth": 2.0,
            "legend.fontsize": 12,
        }
    elif context == "poster":
        rc = {
            "figure.dpi": 200,
            "figure.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.4,
            "axes.labelsize": 18,
            "axes.titlesize": 20,
            "font.family": "sans-serif",
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "lines.linewidth": 2.5,
            "legend.fontsize": 14,
        }
    else:
        raise ValueError(
            f"context must be one of paper / publication / talk / poster; got {context!r}"
        )

    mpl.rcParams.update(rc)


__all__ = ["PALETTE", "apply_theme", "set_palette"]
