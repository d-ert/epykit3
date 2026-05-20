from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import polars as pl

from ._compute import compute_annotation_counts
from ._utils import _get_ax, _save_fig
from ..methyldata import MethylData


# Natural chrom sort: chr1 .. chr22, chrX, chrY, chrM. Falls back to
# lexicographic on anything outside the standard human set so the
# karyogram doesn't barf on non-human assemblies or random contigs.

def _chrom_sort_key(name: str) -> tuple:
    raw = name[3:] if name.lower().startswith("chr") else name
    if raw.isdigit():
        return (0, int(raw))
    # X / Y / MT / M get pinned after autosomes; everything else after them.
    pin = {"X": 100, "Y": 101, "M": 102, "MT": 102}
    if raw.upper() in pin:
        return (0, pin[raw.upper()])
    return (1, raw)


def genomic_context_bar(md: MethylData, ax=None, figsize=(7, 4), save: str | None = None):
    dmc = md.dmc
    if dmc is None or "feature_type" not in dmc.columns:
        raise ValueError("No feature annotations found. Run ep.tl.annotate(md, gtf=...) first.")

    counts = compute_annotation_counts(dmc, annot_col="feature_type")
    fig, ax = _get_ax(ax, figsize)
    ax.bar(counts["feature_type"].to_list(), counts["count"].to_list())
    ax.set_xlabel("Feature type")
    ax.set_ylabel("Count")
    ax.set_title("Genomic context")
    ax.tick_params(axis="x", rotation=45)

    if save:
        _save_fig(md, fig, save)
    return fig, ax


def cpg_island_pie(md: MethylData, ax=None, figsize=(5, 5), save: str | None = None):
    dmc = md.dmc
    if dmc is None or "cpg_context" not in dmc.columns:
        raise ValueError("No CpG-island annotations found. Run ep.tl.annotate(md, cpg_islands=...) first.")

    counts = compute_annotation_counts(dmc, annot_col="cpg_context")
    fig, ax = _get_ax(ax, figsize)
    ax.pie(counts["count"].to_list(), labels=counts["cpg_context"].to_list(), autopct="%1.1f%%")
    ax.set_title("CpG island context")

    if save:
        _save_fig(md, fig, save)
    return fig, ax


def karyogram(
    md: MethylData,
    value: str = "meth_diff",
    *,
    bin_size_bp: int = 1_000_000,
    chromosomes: Optional[Sequence[str]] = None,
    cmap: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    only_significant: bool = False,
    alpha: float = 0.05,
    ax=None,
    figsize=(10, 6),
    save: str | None = None,
):
    """Genome-wide chromosome-painter view of a per-CpG value.

    One horizontal bar per chromosome along the y-axis, with the
    chromosome divided into ``bin_size_bp`` tiles along the x-axis. Each
    tile is coloured by the *mean* of ``value`` (default ``meth_diff``)
    over CpGs in that tile. Empty tiles render as white. The natural use
    is to spot megabase-scale methylation differences (cancer/normal
    contrasts often show whole-arm hypomethylation, which a per-CpG
    Manhattan plot averages away).

    Parameters
    ----------
    md : MethylData
        Must have a DMC table on ``md.dmc``.
    value : str
        Column in the DMC table to summarise. Common choices:
        ``"meth_diff"`` (signed Deltabeta; RdBu_r colormap reads as
        hypo/hyper), ``"mean_beta_case"`` (absolute beta in treatment),
        ``"-log10_qvalue"`` (significance density -- autocomputed when
        the column is absent but ``qvalue`` is present).
    bin_size_bp : int
        Tile size along the chromosome axis. Default 1 Mb.
    chromosomes : sequence of str, optional
        Restrict and re-order the y-axis. Defaults to all chromosomes
        present in the DMC table, sorted by natural human order
        (chr1..chr22, chrX, chrY, chrM, others).
    cmap : str
        Matplotlib colormap. Default ``"RdBu_r"`` (good for signed
        meth_diff). Use ``"viridis"`` or ``"magma"`` for unsigned
        magnitudes.
    vmin, vmax : float, optional
        Manual colour-scale limits. Defaults pick a symmetric range
        around 0 for ``meth_diff``; otherwise the 1st / 99th percentiles
        of finite values.
    only_significant : bool
        If True, drop sites with ``qvalue >= alpha`` before binning so
        only DMC density (or signed Deltabeta of DMCs) drives the colour.
    alpha : float
        q-value threshold when ``only_significant=True``. Default 0.05.
    """
    dmc = md.dmc
    if dmc is None:
        raise ValueError(
            "md.dmc is None -- run ep.tl.dmc(md) before plotting a karyogram."
        )
    if only_significant:
        if "qvalue" not in dmc.columns:
            raise ValueError(
                "only_significant=True requires a 'qvalue' column on md.dmc."
            )
        dmc = dmc.filter(pl.col("qvalue") < alpha)
        if dmc.is_empty():
            raise ValueError(
                f"No DMCs pass qvalue<{alpha}; nothing to plot."
            )

    # Auto-compute -log10_qvalue if the user asked for it but it's absent.
    if value == "-log10_qvalue" and "-log10_qvalue" not in dmc.columns:
        if "qvalue" not in dmc.columns:
            raise ValueError(
                "value='-log10_qvalue' needs either that column or 'qvalue' "
                "on md.dmc to derive it from."
            )
        dmc = dmc.with_columns(
            (-pl.col("qvalue").log10()).alias("-log10_qvalue")
        )

    if value not in dmc.columns:
        raise ValueError(
            f"value={value!r} is not a column on md.dmc. "
            f"Available: {sorted(dmc.columns)[:15]}..."
        )
    if "chrom" not in dmc.columns or "pos" not in dmc.columns:
        raise ValueError("md.dmc must carry 'chrom' and 'pos' columns.")

    if chromosomes is None:
        chromosomes = sorted(
            set(dmc.get_column("chrom").to_list()), key=_chrom_sort_key,
        )
    else:
        chromosomes = list(chromosomes)
    if not chromosomes:
        raise ValueError("no chromosomes to plot")

    # Bin per (chrom, bin). Use lazy polars to avoid materialising the
    # full DMC table twice for large genomes.
    binned = (
        dmc.lazy()
        .filter(pl.col("chrom").is_in(chromosomes))
        .with_columns((pl.col("pos") // bin_size_bp).alias("_bin"))
        .group_by(["chrom", "_bin"])
        .agg(pl.col(value).mean().alias("_val"))
        .collect()
    )
    if binned.is_empty():
        raise ValueError("No data to bin; check chromosomes / value args.")

    # Per-chrom max bin for the heatmap width.
    chrom_max_bin = {
        row["chrom"]: int(row["_max"]) for row in
        binned.group_by("chrom").agg(pl.col("_bin").max().alias("_max"))
        .to_dicts()
    }
    n_bins = max(chrom_max_bin.values()) + 1

    grid = np.full((len(chromosomes), n_bins), np.nan, dtype=np.float64)
    chrom_to_row = {c: i for i, c in enumerate(chromosomes)}
    for row in binned.to_dicts():
        i = chrom_to_row.get(row["chrom"])
        if i is None:
            continue
        b = int(row["_bin"])
        grid[i, b] = float(row["_val"]) if row["_val"] is not None else np.nan

    # Colour limits. For signed Deltabeta-style metrics, symmetric around 0.
    if vmin is None or vmax is None:
        finite = grid[np.isfinite(grid)]
        if finite.size == 0:
            vmin, vmax = -1.0, 1.0
        elif value == "meth_diff" or (finite.min() < 0 < finite.max()):
            m = max(abs(finite.min()), abs(finite.max()))
            vmin = vmin if vmin is not None else -m
            vmax = vmax if vmax is not None else m
        else:
            lo, hi = np.percentile(finite, [1.0, 99.0])
            vmin = vmin if vmin is not None else float(lo)
            vmax = vmax if vmax is not None else float(hi)

    fig, ax = _get_ax(ax, figsize)
    im = ax.imshow(
        grid, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
        interpolation="nearest", origin="upper",
    )
    ax.set_yticks(np.arange(len(chromosomes)))
    ax.set_yticklabels(chromosomes)
    ax.set_xlabel(f"Position ({bin_size_bp / 1e6:g} Mb bins)")
    ax.set_ylabel("Chromosome")
    ax.set_title(f"Karyogram: mean {value} per {bin_size_bp / 1e6:g} Mb")
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label(value)

    if save:
        _save_fig(md, fig, save)
    return fig, ax


__all__ = ["genomic_context_bar", "cpg_island_pie", "karyogram"]
