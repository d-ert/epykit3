"""TSS / gene-body metaplots.

Reuses the cached GTF parser in :mod:`epykit.annotate` to enumerate gene
TSS coordinates, then bins per-CpG beta values into ``n_bins`` slots from
``-window_bp`` to ``+window_bp`` around each TSS. One line per sample,
optionally grouped by ``md.obs[group_by]``.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import polars as pl

from .._style import PALETTE
from ._compute import compute_tss_metaplot
from ._utils import _get_ax, _save_fig
from ..methyldata import MethylData
from ..annotate import _parse_gtf_streaming


def _tss_table_from_gtf(gtf_path: str) -> pl.DataFrame:
    """Build a per-gene TSS DataFrame from a (cached) parsed GTF.

    Columns: chrom (Utf8), tss (Int64), strand (Utf8), gene_id (Utf8),
    gene_name (Utf8). TSS = Start for + strand, End-1 for - strand
    (GTF is 0-based half-open inside epykit's parser).
    """
    genes_pd, _ = _parse_gtf_streaming(gtf_path)
    if genes_pd is None or len(genes_pd) == 0:
        return pl.DataFrame(schema={
            "chrom": pl.Utf8, "tss": pl.Int64, "strand": pl.Utf8,
            "gene_id": pl.Utf8, "gene_name": pl.Utf8,
        })

    genes = pl.from_pandas(
        genes_pd[["Chromosome", "Start", "End", "Strand", "gene_id", "gene_name"]]
    ).rename({"Chromosome": "chrom", "Strand": "strand"})
    return (
        genes.with_columns(
            pl.when(pl.col("strand") == "-")
              .then(pl.col("End") - 1)
              .otherwise(pl.col("Start"))
              .cast(pl.Int64).alias("tss")
        )
        .select(["chrom", "tss", "strand", "gene_id", "gene_name"])
    )


def tss_metaplot(
    md: MethylData,
    gtf_path: str,
    *,
    window_bp: int = 2000,
    n_bins: int = 100,
    group_by: Optional[str] = "group",
    max_genes: Optional[int] = None,
    ax=None,
    figsize=(7, 4),
    save: str | None = None,
):
    """Plot mean beta around the TSS, averaged across genes.

    For each sample, beta values are pooled across all gene TSS in a
    ``+/-window_bp`` window, binned into ``n_bins`` slots based on relative
    position (sign-flipped for - strand so 5'->3' is left->right), and
    averaged. Each sample is drawn as a faint line, with one bold line
    per group (if ``group_by`` matches a column on ``md.obs``).

    Parameters
    ----------
    md : MethylData
        Must have a methylstore at ``md.store`` (i.e. one of
        ``read_bismark`` / ``read_nfcore_methylseq`` / ``load`` was run).
    gtf_path : str
        GTF / GFF3 (gz-allowed) for TSS coordinates. Uses epykit's
        bounded LRU GTF cache, so repeated calls within one process
        skip the streaming parse.
    window_bp : int
        Half-window width in base pairs (default 2000 -> +/-2 kb).
    n_bins : int
        Number of bins to summarise beta within the window. Default 100.
    group_by : str or None
        Optional ``md.obs`` column for per-group means. Pass ``None`` to
        draw only the per-sample lines.
    max_genes : int, optional
        Cap the number of TSS used. Useful for fast smoke tests; None
        uses every gene in the GTF.

    Returns
    -------
    (Figure, Axes)
    """
    res = compute_tss_metaplot(
        md, gtf_path,
        window_bp=window_bp, n_bins=n_bins,
        group_by=group_by, max_genes=max_genes,
    )
    samples = res.samples
    mean_beta = res.mean_beta
    x = res.x

    fig, ax = _get_ax(ax, figsize)
    if res.group_col is not None:
        groups = res.groups
        unique_groups = sorted(set(groups))
        group_palette = {
            g: PALETTE.get("treatment" if i else "control", PALETTE["neutral"])
            for i, g in enumerate(unique_groups)
        }
        for i, samp in enumerate(samples):
            ax.plot(
                x, mean_beta[i],
                color=group_palette.get(groups[i], PALETTE["neutral"]),
                alpha=0.25, linewidth=1,
            )
        for g in unique_groups:
            mask = np.array([gg == g for gg in groups])
            if not mask.any():
                continue
            grp_mean = np.nanmean(mean_beta[mask], axis=0)
            ax.plot(
                x, grp_mean,
                color=group_palette.get(g, PALETTE["neutral"]),
                linewidth=2.2, label=str(g),
            )
        ax.legend(title=res.group_col, frameon=False)
    else:
        for i, samp in enumerate(samples):
            ax.plot(x, mean_beta[i], alpha=0.7, linewidth=1.2, label=samp)
        if len(samples) <= 12:
            ax.legend(frameon=False, fontsize=8)

    ax.axvline(0, color="black", lw=0.7, ls="--", alpha=0.5)
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel("Mean beta")
    ax.set_title(f"TSS metaplot (+/-{window_bp} bp, n_bins={n_bins})")

    if save:
        _save_fig(md, fig, save)
    return fig, ax


def _gene_table_from_gtf(gtf_path: str) -> pl.DataFrame:
    """Per-gene (chrom, tx_start, tx_end, strand, gene_id, gene_name)."""
    genes_pd, _ = _parse_gtf_streaming(gtf_path)
    if genes_pd is None or len(genes_pd) == 0:
        return pl.DataFrame(schema={
            "chrom": pl.Utf8, "tx_start": pl.Int64, "tx_end": pl.Int64,
            "strand": pl.Utf8, "gene_id": pl.Utf8, "gene_name": pl.Utf8,
        })
    genes = pl.from_pandas(
        genes_pd[["Chromosome", "Start", "End", "Strand", "gene_id", "gene_name"]]
    ).rename({
        "Chromosome": "chrom", "Strand": "strand",
        "Start": "tx_start", "End": "tx_end",
    })
    return genes.select([
        "chrom", "tx_start", "tx_end", "strand", "gene_id", "gene_name",
    ])


def gene_body_metaplot(
    md: MethylData,
    gtf_path: str,
    *,
    flank_bp: int = 2000,
    n_bins_flank: int = 30,
    n_bins_body: int = 40,
    min_gene_bp: int = 500,
    group_by: Optional[str] = "group",
    max_genes: Optional[int] = None,
    ax=None,
    figsize=(8, 4),
    save: str | None = None,
):
    """Plot mean beta across a normalised gene body with TSS / TES flanks.

    The x-axis has three zones:

    * Left flank: ``-flank_bp`` to ``+0`` relative to TSS, binned into
      ``n_bins_flank`` slots.
    * Body: TSS -> TES on every gene, length-normalised into
      ``n_bins_body`` slots so a 50 kb gene and a 500 kb gene contribute
      the same number of points per bin.
    * Right flank: TES to ``+flank_bp``, binned into ``n_bins_flank`` slots.

    Genes shorter than ``min_gene_bp`` are dropped to avoid the body
    binning collapsing onto a handful of CpGs. Strand is honoured so
    5'->3' is left->right for every gene.

    Parameters mirror :func:`tss_metaplot` where they overlap.
    """
    samples = md.obs.get_column("sample_id").to_list()
    if not samples:
        raise ValueError("md.obs has no samples")

    genes = _gene_table_from_gtf(gtf_path)
    if genes.is_empty():
        raise ValueError(f"No gene records found in GTF {gtf_path!r}")
    genes = genes.filter(
        (pl.col("tx_end") - pl.col("tx_start")) >= min_gene_bp
    )
    if genes.is_empty():
        raise ValueError(
            f"No genes >={min_gene_bp} bp in the GTF -- every record is too short "
            "for body-binning."
        )
    if max_genes is not None and len(genes) > max_genes:
        genes = genes.head(max_genes)

    total_bins = n_bins_flank + n_bins_body + n_bins_flank
    sum_beta = np.zeros((len(samples), total_bins), dtype=np.float64)
    count = np.zeros((len(samples), total_bins), dtype=np.int64)
    sample_idx = {s: i for i, s in enumerate(samples)}

    body_start = n_bins_flank
    body_end = n_bins_flank + n_bins_body

    chroms = sorted(set(genes.get_column("chrom").to_list()))
    for chrom in chroms:
        g_chrom = genes.filter(pl.col("chrom") == chrom)
        if g_chrom.is_empty():
            continue
        starts = g_chrom["tx_start"].to_numpy().astype(np.int64)
        ends = g_chrom["tx_end"].to_numpy().astype(np.int64)
        strands = np.array(
            [1 if s != "-" else -1 for s in g_chrom["strand"].to_list()],
            dtype=np.int8,
        )

        pattern = f"{md.store}/sample=*/chrom={chrom}/part-*.parquet"
        try:
            chrom_df = (
                pl.scan_parquet(pattern)
                .select(["pos", "sample", "N_meth", "coverage"])
                .filter(pl.col("coverage") > 0)
                .collect()
            )
        except Exception:
            continue
        if chrom_df.is_empty():
            continue
        positions = chrom_df["pos"].to_numpy().astype(np.int64)
        samples_arr = chrom_df["sample"].to_list()
        betas = (
            chrom_df["N_meth"].to_numpy().astype(np.float64)
            / chrom_df["coverage"].to_numpy().astype(np.float64)
        )
        order = np.argsort(positions, kind="mergesort")
        positions = positions[order]
        betas = betas[order]
        samples_arr = [samples_arr[i] for i in order]
        samples_idx_arr = np.fromiter(
            (sample_idx.get(s, -1) for s in samples_arr),
            count=len(samples_arr), dtype=np.int32,
        )

        for tx_start, tx_end, strand in zip(starts, ends, strands):
            # Choose window depending on strand orientation.
            # For + strand: left flank = [tx_start - flank, tx_start),
            #               body         = [tx_start, tx_end],
            #               right flank  = (tx_end, tx_end + flank].
            # For - strand: flip so 5'->3' is left->right.
            if strand == 1:
                lo, hi = tx_start - flank_bp, tx_end + flank_bp
            else:
                lo, hi = tx_start - flank_bp, tx_end + flank_bp
            left = np.searchsorted(positions, lo, side="left")
            right = np.searchsorted(positions, hi, side="right")
            if right <= left:
                continue
            chunk_pos = positions[left:right]
            chunk_samples = samples_idx_arr[left:right]
            chunk_betas = betas[left:right]
            mask = chunk_samples >= 0
            if not mask.any():
                continue
            chunk_pos = chunk_pos[mask]
            chunk_samples = chunk_samples[mask]
            chunk_betas = chunk_betas[mask]

            # Zone assignment.
            #   pos < tx_start  -> left flank
            #   tx_start <= pos <= tx_end -> body
            #   pos > tx_end    -> right flank
            in_left = chunk_pos < tx_start
            in_right = chunk_pos > tx_end
            in_body = ~(in_left | in_right)

            bin_idx = np.zeros(chunk_pos.shape, dtype=np.int64)
            # Left flank: rel  in  [-flank, 0) -> bins [0, n_bins_flank).
            rel_left = (chunk_pos[in_left] - tx_start)  # negative
            bf = np.floor((rel_left + flank_bp) / (flank_bp / n_bins_flank))
            bin_idx[in_left] = np.clip(bf, 0, n_bins_flank - 1)
            # Body: fraction along gene -> [body_start, body_end).
            body_len = max(1, tx_end - tx_start)
            frac = (chunk_pos[in_body] - tx_start) / body_len
            bin_idx[in_body] = body_start + np.clip(
                np.floor(frac * n_bins_body).astype(np.int64),
                0, n_bins_body - 1,
            )
            # Right flank: rel  in  (0, flank] -> bins [body_end, total).
            rel_right = chunk_pos[in_right] - tx_end
            br = np.floor(rel_right / (flank_bp / n_bins_flank))
            bin_idx[in_right] = body_end + np.clip(
                br.astype(np.int64), 0, n_bins_flank - 1,
            )

            if strand == -1:
                # Flip the bin axis so 5'->3' reads left->right for - strand.
                bin_idx = (total_bins - 1) - bin_idx

            np.add.at(sum_beta, (chunk_samples, bin_idx), chunk_betas)
            np.add.at(count, (chunk_samples, bin_idx), 1)

    with np.errstate(invalid="ignore"):
        mean_beta = np.where(count > 0, sum_beta / count, np.nan)

    # X coordinate: integer bin index; we'll relabel the ticks at the
    # zone boundaries.
    x = np.arange(total_bins)

    fig, ax = _get_ax(ax, figsize)
    if group_by and group_by in md.obs.columns:
        groups = md.obs.get_column(group_by).to_list()
        unique_groups = sorted(set(groups))
        group_palette = {
            g: PALETTE.get(
                "treatment" if i else "control", PALETTE["neutral"],
            )
            for i, g in enumerate(unique_groups)
        }
        for i, samp in enumerate(samples):
            ax.plot(
                x, mean_beta[i],
                color=group_palette.get(groups[i], PALETTE["neutral"]),
                alpha=0.25, linewidth=1,
            )
        for g in unique_groups:
            mask = np.array([gg == g for gg in groups])
            if not mask.any():
                continue
            grp_mean = np.nanmean(mean_beta[mask], axis=0)
            ax.plot(
                x, grp_mean,
                color=group_palette.get(g, PALETTE["neutral"]),
                linewidth=2.2, label=str(g),
            )
        ax.legend(title=group_by, frameon=False)
    else:
        for i, samp in enumerate(samples):
            ax.plot(x, mean_beta[i], alpha=0.7, linewidth=1.2, label=samp)
        if len(samples) <= 12:
            ax.legend(frameon=False, fontsize=8)

    ax.axvline(body_start - 0.5, color="black", lw=0.7, ls="--", alpha=0.5)
    ax.axvline(body_end - 0.5, color="black", lw=0.7, ls="--", alpha=0.5)
    tick_pos = [
        0, body_start - 0.5, (body_start + body_end - 1) / 2,
        body_end - 0.5, total_bins - 1,
    ]
    tick_lbl = [
        f"-{flank_bp}", "TSS", "gene body", "TES", f"+{flank_bp}",
    ]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lbl)
    ax.set_xlabel("Position along gene")
    ax.set_ylabel("Mean beta")
    ax.set_title(
        f"Gene-body metaplot (+/-{flank_bp} bp flanks, "
        f"n_genes={len(genes):,})"
    )

    if save:
        _save_fig(md, fig, save)
    return fig, ax


__all__ = ["tss_metaplot", "gene_body_metaplot"]
