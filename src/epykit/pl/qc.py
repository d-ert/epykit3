from __future__ import annotations

from ._compute import compute_coverage_distribution
from ._utils import _get_ax, _save_fig
from ..methyldata import MethylData
import polars as pl


def coverage_histogram(
    md: MethylData,
    bins: int = 100,
    *,
    max_points: int = 1_000_000,
    ax=None,
    figsize=(6, 4),
    save: str | None = None,
):
    """Plot histogram of coverage across all sites.

    For large stores, ``compute_coverage_distribution`` subsamples to
    ``max_points`` rows deterministically before materialising.
    """
    cov = compute_coverage_distribution(md, max_points=max_points)

    fig, ax = _get_ax(ax, figsize)
    ax.hist(cov, bins=bins, edgecolor="black")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Count")
    ax.set_title("Coverage histogram")

    if save:
        _save_fig(md, fig, save)
    return fig, ax


def methylation_heatmap(md: MethylData, n_top: int = 1000, ax=None, figsize=(8, 6), save: str | None = None):
    try:
        import seaborn as sns
    except ImportError as exc:
        raise ImportError("seaborn is required for heatmaps. Install with: pip install seaborn") from exc

    dmc = md.dmc
    if dmc is None:
        raise ValueError("No DMC results available. Run ep.tl.dmc(md) first.")

    top = (
        dmc
        .filter(pl.col("meth_diff").is_not_null())
        .with_columns(pl.col("meth_diff").abs().alias("abs_diff"))
        .sort("abs_diff", descending=True)
        .head(n_top)
        .select(["chrom", "pos"])
    )
    if len(top) == 0:
        raise ValueError("No DMC rows available to build heatmap")

    samples = md.obs.get_column("sample_id").to_list()
    
    # top is already a collected DataFrame from md.dmc, no need to call .collect()
    if len(top) == 0:
        raise ValueError("No top DMCs to build heatmap")
    
    # Process each sample separately to reduce memory footprint
    site_dfs = []
    for sample in samples:
        sample_df = (
            pl.scan_parquet(f"{md.store}/sample={sample}/chrom=*/part-*.parquet")
            .select(["chrom", "pos", "N_meth", "coverage"])
            .join(top.lazy(), on=["chrom", "pos"], how="inner")
            .collect()
        )
        if len(sample_df) > 0:
            sample_df = sample_df.with_columns(
                pl.when(pl.col("coverage") > 0)
                .then(pl.col("N_meth") / pl.col("coverage"))
                .otherwise(None)
                .alias("beta"),
                pl.lit(sample).alias("sample")
            )
            site_dfs.append(sample_df)
    
    if not site_dfs:
        raise ValueError("No sites found in store matching top DMCs")
    
    site_df = pl.concat(site_dfs)

    pivot = site_df.pivot(values="beta", index=["chrom", "pos"], on="sample", aggregate_function="mean")
    for sample in samples:
        if sample not in pivot.columns:
            pivot = pivot.with_columns(pl.lit(None).alias(sample))

    matrix = pivot.select(samples).to_numpy()
    fig, ax = _get_ax(ax, figsize)
    sns.heatmap(matrix, cmap="viridis", ax=ax)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Top DMC sites")
    ax.set_title(f"Methylation heatmap (top {n_top})")

    if save:
        _save_fig(md, fig, save)
    return fig, ax


def mbias_plot(
    mbias_data,
    *,
    context: str = "CpG",
    ax=None,
    figsize=(7, 4),
    save: str | None = None,
    md: MethylData | None = None,
):
    """Plot per-read-position methylation bias for one or more samples.

    Input shape is ``{sample_id: pl.DataFrame}`` where each frame is the
    output of :func:`epykit.nfcore_qc.parse_bismark_mbias`, or
    ``{sample_id: str | Path}`` pointing at the raw M-bias text files
    (they'll be parsed inline).

    The plot shows percent methylation against read position, one line
    per sample x read (``R1`` / ``R2``). Standard interpretation: a flat
    plateau in the middle of the read with deflections at either end
    indicates the safe trim region. A consistently sloped line is a sign
    of incomplete bisulfite conversion or fill-in artefacts at PBAT /
    EM-seq library tails.

    Parameters
    ----------
    mbias_data : dict[str, pl.DataFrame | str | Path]
        Per-sample M-bias tables or file paths.
    context : {"CpG", "CHG", "CHH"}
        Which methylation context to plot. Default ``"CpG"``.
    ax : matplotlib.axes.Axes, optional
        Existing axes to draw onto. Pass for composition into a dashboard.
    figsize : (w, h)
        Figure size when ``ax`` isn't provided. Default ``(7, 4)``.
    save : str, optional
        If a ``MethylData`` is supplied via ``md=``, the figure is saved
        under its plot directory using :func:`_save_fig`.
    md : MethylData, optional
        Used only for ``save`` path resolution; not read otherwise.

    Returns
    -------
    (matplotlib.figure.Figure, matplotlib.axes.Axes)
    """
    if not mbias_data:
        raise ValueError("mbias_data is empty; provide at least one sample.")

    from .._style import PALETTE

    fig, ax = _get_ax(ax, figsize)
    samples = sorted(mbias_data.keys())
    palette = (
        list(PALETTE.values()) if isinstance(PALETTE, dict) else list(PALETTE)
    )

    drawn = 0
    for i, sample in enumerate(samples):
        entry = mbias_data[sample]
        if not isinstance(entry, pl.DataFrame):
            from ..nfcore_qc import parse_bismark_mbias
            entry = parse_bismark_mbias(entry)
        sub = entry.filter(pl.col("context") == context)
        if sub.is_empty():
            continue
        colour = palette[i % len(palette)] if palette else None
        for read in sorted(sub.get_column("read").unique().to_list()):
            ssub = sub.filter(pl.col("read") == read).sort("position")
            if ssub.is_empty():
                continue
            ls = "-" if read == "R1" else "--"
            ax.plot(
                ssub.get_column("position").to_numpy(),
                ssub.get_column("percent").to_numpy(),
                label=f"{sample} {read}",
                color=colour, linestyle=ls, linewidth=1.2,
            )
            drawn += 1

    if drawn == 0:
        raise ValueError(
            f"No rows for context={context!r} in any provided sample. "
            "Did the M-bias report cover this context?"
        )

    ax.set_xlabel("Read position (bp)")
    ax.set_ylabel("% methylated")
    ax.set_title(f"M-bias ({context})")
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=7, ncol=2)

    if save and md is not None:
        _save_fig(md, fig, save)
    return fig, ax


__all__ = ["coverage_histogram", "methylation_heatmap", "mbias_plot"]