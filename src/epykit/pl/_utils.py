from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np
from matplotlib.figure import Figure


def _get_ax(ax=None, figsize=(6, 4)) -> Tuple[Figure, object]:
    if ax is None:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)
        return fig, ax
    else:
        return ax.figure, ax


_VALID_FORMATS = ("png", "pdf", "svg", "both")


def _save_fig(
    md,
    fig: Figure,
    name: str,
    out_dir: str | None = None,
    *,
    format: str = "png",
    dpi: int = 150,
) -> str | list[str]:
    """Save ``fig`` under ``figures/<name>.<format>``.

    ``format`` accepts ``"png"``, ``"pdf"``, ``"svg"``, or ``"both"`` for
    PNG + PDF side by side (useful when you want a thumbnail and a
    vector deliverable from one call). Records the path on
    ``md.uns["figures"]`` so the report can find it later.
    """
    if format not in _VALID_FORMATS:
        raise ValueError(
            f"format must be one of {_VALID_FORMATS}; got {format!r}"
        )
    out = Path(out_dir or "figures")
    if format == "both":
        formats = ("png", "pdf")
    else:
        formats = (format,)

    written: list[str] = []
    for fmt in formats:
        path = out / f"{name}.{fmt}"
        path.parent.mkdir(parents=True, exist_ok=True)
        save_dpi = dpi if fmt == "png" else None
        fig.savefig(path, dpi=save_dpi, bbox_inches="tight")
        if hasattr(md, "uns"):
            if "figures" not in md.uns:
                md.uns["figures"] = {}
            md.uns["figures"][f"{name}.{fmt}"] = str(path)
        written.append(str(path))

    import matplotlib.pyplot as plt
    plt.close(fig)
    return written[0] if len(written) == 1 else written


def build_sample_site_matrix(md, n_sites: int = 10_000) -> tuple[np.ndarray, list[str]]:
    """Back-compat shim around :func:`epykit.pl._compute.compute_sample_site_matrix`.

    Old call sites returned ``(matrix, samples)``; the compute version
    additionally returns the number of sites kept. This wrapper drops the
    third value so existing code keeps working.
    """
    from ._compute import compute_sample_site_matrix
    matrix, samples, _n = compute_sample_site_matrix(md, n_sites=n_sites)
    return matrix, samples


__all__ = ["_get_ax", "_save_fig", "build_sample_site_matrix"]
