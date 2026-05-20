"""Composite QC dashboard .

A single figure summarising the most-used QC signals so users don't have
to assemble seven PNGs themselves. Reads only material that ``tl.qc`` has
already written onto ``md.obs`` / ``md.uns`` -- no extra methylstore reads.
Missing panels degrade gracefully to a placeholder caption.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from .._style import PALETTE
from ._utils import _save_fig


def _empty_panel(ax, label: str) -> None:
    ax.text(
        0.5, 0.5, label, ha="center", va="center",
        transform=ax.transAxes, fontsize=9, color="gray",
    )
    ax.set_xticks([])
    ax.set_yticks([])


def qc_dashboard(
    md,
    *,
    figsize=(15, 9),
    save: str | None = None,
):
    """Composite QC figure: conversion rate, coverage, methylation, correlation.

    Reads everything from ``md.obs`` / ``md.uns``. Missing metrics are
    rendered as captioned placeholders rather than crashing -- so
    ``pl.qc_dashboard`` works on a partially-run pipeline.
    """
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.35)

    obs = md.obs
    samples = obs.get_column("sample_id").to_list()
    groups = (
        obs.get_column("group").to_list()
        if "group" in obs.columns
        else ["sample"] * len(samples)
    )

    # Panel 1: bisulfite conversion rate (if available)
    ax1 = fig.add_subplot(gs[0, 0])
    if "bisulfite_conversion_rate" in obs.columns:
        vals = obs.get_column("bisulfite_conversion_rate").to_numpy()
        ax1.bar(range(len(samples)), vals, color=PALETTE.get("neutral", "#999"))
        ax1.axhline(0.995, color="red", linestyle="--", lw=0.8)
        ax1.set_xticks(range(len(samples)))
        ax1.set_xticklabels(samples, rotation=90, fontsize=7)
        ax1.set_ylim(0.98, 1.001)
        ax1.set_ylabel("Conversion rate")
        ax1.set_title("Bisulfite conversion")
    else:
        _empty_panel(ax1, "conversion rate\n(qc not run with CHH store)")

    # Panel 2: mean coverage
    ax2 = fig.add_subplot(gs[0, 1])
    if "mean_coverage" in obs.columns:
        vals = obs.get_column("mean_coverage").to_numpy()
        ax2.bar(range(len(samples)), vals, color=PALETTE.get("neutral", "#999"))
        ax2.set_xticks(range(len(samples)))
        ax2.set_xticklabels(samples, rotation=90, fontsize=7)
        ax2.set_ylabel("Mean coverage")
        ax2.set_title("Coverage")
    else:
        _empty_panel(ax2, "mean coverage\n(qc not run)")

    # Panel 3: global methylation grouped by treatment/group
    ax3 = fig.add_subplot(gs[0, 2])
    if "global_methylation" in obs.columns:
        vals = obs.get_column("global_methylation").to_numpy()
        unique_groups = sorted(set(groups))
        colors = [
            (PALETTE.get("treatment", "#c45") if g == "treatment"
             else PALETTE.get("control", "#456"))
            for g in groups
        ]
        ax3.bar(range(len(samples)), vals, color=colors)
        ax3.set_xticks(range(len(samples)))
        ax3.set_xticklabels(samples, rotation=90, fontsize=7)
        ax3.set_ylabel("Global beta")
        ax3.set_title("Global methylation")
    else:
        _empty_panel(ax3, "global methylation\n(qc not run)")

    # Panel 4 (row 2 span): sample correlation heatmap
    ax4 = fig.add_subplot(gs[1, :2])
    corr_df = md.uns.get("qc_sample_correlation")
    if isinstance(corr_df, pl.DataFrame) and len(corr_df) > 0:
        s_all = sorted(set(
            corr_df.get_column("sample_a").to_list()
            + corr_df.get_column("sample_b").to_list()
        ))
        n = len(s_all)
        idx = {s: i for i, s in enumerate(s_all)}
        mat = np.full((n, n), np.nan, dtype=np.float64)
        for row in corr_df.iter_rows(named=True):
            mat[idx[row["sample_a"]], idx[row["sample_b"]]] = row["correlation"]
        im = ax4.imshow(mat, vmin=-1, vmax=1, cmap="RdBu_r", aspect="equal")
        ax4.set_xticks(range(n))
        ax4.set_yticks(range(n))
        ax4.set_xticklabels(s_all, rotation=90, fontsize=7)
        ax4.set_yticklabels(s_all, fontsize=7)
        fig.colorbar(im, ax=ax4, fraction=0.046, pad=0.04)
        ax4.set_title("Sample correlation")
    else:
        _empty_panel(ax4, "sample correlation\n(run qc(..., run_sample_correlation=True))")

    # Panel 5: sex check (if available)
    ax5 = fig.add_subplot(gs[1, 2])
    sex_df = md.uns.get("qc_sex_check")
    if isinstance(sex_df, pl.DataFrame) and len(sex_df) > 0:
        ax5.scatter(
            range(len(sex_df)),
            sex_df.get_column("mean_chrx_beta").to_numpy(),
            c=[
                "tab:blue" if s == "male" else "tab:red"
                if s == "female" else "gray"
                for s in sex_df.get_column("inferred_sex").to_list()
            ],
        )
        ax5.set_xticks(range(len(sex_df)))
        ax5.set_xticklabels(
            sex_df.get_column("sample_id").to_list(),
            rotation=90, fontsize=7,
        )
        ax5.set_ylabel("Mean beta (chrX)")
        ax5.set_title("Sex check")
    else:
        _empty_panel(ax5, "sex check\n(run qc(..., run_sex_check=True))")

    # Panel 6: low-coverage flags
    ax6 = fig.add_subplot(gs[2, 0])
    if "low_coverage_flag" in obs.columns:
        flags = obs.get_column("low_coverage_flag").to_list()
        ax6.bar(
            range(len(samples)),
            [1 if f else 0 for f in flags],
            color=[("red" if f else "green") for f in flags],
        )
        ax6.set_xticks(range(len(samples)))
        ax6.set_xticklabels(samples, rotation=90, fontsize=7)
        ax6.set_ylabel("Low-cov flag")
        ax6.set_ylim(0, 1.2)
        ax6.set_title("Coverage flags")
    else:
        _empty_panel(ax6, "coverage flags\n(qc not run)")

    # Panel 7-9: free text annotations of remaining cohort-level numbers.
    ax7 = fig.add_subplot(gs[2, 1:])
    notes = []
    notes.append(f"n_samples = {len(samples)}")
    if "n_sites_filtered" in md.uns:
        notes.append(f"n_sites (filtered) = {md.uns['n_sites_filtered']:,}")
    elif "n_sites_raw" in md.uns:
        notes.append(f"n_sites (raw) = {md.uns['n_sites_raw']:,}")
    if "dmc" in md.uns and isinstance(md.uns["dmc"], dict):
        notes.append(f"DMC test = {md.uns['dmc'].get('test_used')}")
        notes.append(f"DMC n_sites tested = {md.uns['dmc'].get('n_sites')}")
    if "dmr" in md.uns:
        dmr_df = md.uns["dmr"]
        if isinstance(dmr_df, pl.DataFrame):
            notes.append(f"n_DMRs = {len(dmr_df)}")
    if "qc_sex_check" in md.uns and isinstance(md.uns["qc_sex_check"], pl.DataFrame):
        mismatches = int(md.uns["qc_sex_check"].get_column("mismatch").sum())
        notes.append(f"sex mismatches = {mismatches}")
    if "contamination_score" in obs.columns:
        top = obs.sort("contamination_score", descending=True).head(3)
        notes.append("highest contamination_score:")
        for row in top.iter_rows(named=True):
            notes.append(
                f"  {row['sample_id']}: {row['contamination_score']:.2f}"
            )
    ax7.axis("off")
    ax7.text(
        0.0, 1.0, "\n".join(notes), ha="left", va="top",
        family="monospace", fontsize=9, transform=ax7.transAxes,
    )

    fig.suptitle("epykit QC dashboard", fontsize=12, y=0.995)
    if save:
        _save_fig(md, fig, save)
    return fig, [ax1, ax2, ax3, ax4, ax5, ax6, ax7]


__all__ = ["qc_dashboard"]
