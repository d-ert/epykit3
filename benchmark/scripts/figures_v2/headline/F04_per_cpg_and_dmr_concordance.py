"""F04 — Per-CpG and DMR concordance on GSE263850.

Top left: per-CpG meth_diff scatter for sites called significant (q ≤ 0.05) by both
epykit lr and methylKit. Subsamples to ≤ 30 k points for clarity.

Top right & bottom: DMR-level meth_diff scatter, three pairs from per_dmr_stat_concordance:
(epykit_chain_merge × methylkit), (dss × epykit_chain_merge),
(dss × methylkit). Pearson r and pair-count annotated.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import polars as pl
import matplotlib.pyplot as plt

from _headline_style import (
    DATA_STUDY3, METHYLKIT_REAL, FIG_OUT, save_dual, setup,
)


def _scatter_helper(ax, x, y, *, label_x: str, label_y: str, color: str,
                    equal_aspect: bool = True) -> None:
    n = len(x)
    if n == 0:
        ax.set_visible(False)
        return
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n_finite = len(x)
    r = float(np.corrcoef(x, y)[0, 1]) if n_finite > 2 else float("nan")

    ax.scatter(x, y, s=6, c=color, alpha=0.25, edgecolors="none", rasterized=True)
    lo = float(min(x.min(), y.min()))
    hi = float(max(x.max(), y.max()))
    pad = (hi - lo) * 0.05
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
            color="#C0392B", linestyle="--", linewidth=1, label="y = x")
    ax.axhline(0, color="#cccccc", linewidth=0.5)
    ax.axvline(0, color="#cccccc", linewidth=0.5)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)
    ax.text(
        0.05, 0.95,
        f"Pearson r = {r:.3f}\nn = {n_finite:,}",
        transform=ax.transAxes,
        va="top", ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="#999999", linewidth=0.6),
    )
    if equal_aspect:
        ax.set_aspect("equal", adjustable="box")


def main() -> None:
    setup()

    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.35)

    # ---- TOP LEFT: per-CpG ----
    ax_pcpg = fig.add_subplot(gs[0, 0])

    ek = pd.read_csv(DATA_STUDY3 / "dmc_significant_qval05.csv",
                     usecols=["chrom", "pos", "meth_diff"])
    mk = pd.read_csv(METHYLKIT_REAL / "dmc_significant_qval05.csv",
                     usecols=["chrom", "pos", "meth_diff"])

    # epykit positions are 0-based (BED); methylKit positions are 1-based.
    ek = ek.assign(pos=ek["pos"] + 1)

    # Both meth_diff columns are on the same signed 0..100 scale.
    merged = pd.merge(
        ek, mk, on=["chrom", "pos"], how="inner",
        suffixes=("_ek", "_mk"),
    )
    n_inter = len(merged)
    if n_inter > 30000:
        merged = merged.sample(n=30000, random_state=42)

    _scatter_helper(
        ax_pcpg, merged["meth_diff_ek"].to_numpy(), merged["meth_diff_mk"].to_numpy(),
        label_x="epykit lr Δβ (%, sig sites)",
        label_y="methylKit Δβ (%, sig sites)",
        color="#0F4C81",
        equal_aspect=True,
    )
    ax_pcpg.set_title(
        f"Per-CpG concordance\n(sites significant in both tools, n = {n_inter:,})",
        loc="left",
        fontsize=11,
    )

    # ---- DMR pairs: top-right and bottom ----
    df = pl.read_parquet(
        DATA_STUDY3 / "comparisons_post_phase3" / "per_dmr_stat_concordance.parquet"
    ).to_pandas()

    pair_specs = [
        ("epykit_chain_merge", "methylkit", "#0F4C81",
         "epykit chain_merge", "methylKit"),
        ("dss",                "epykit_chain_merge", "#2E8B57",
         "DSS",               "epykit chain_merge"),
        ("dss",                "methylkit", "#E07B39",
         "DSS",               "methylKit"),
    ]

    subplot_positions = [(0, 1), (1, 0), (1, 1)]  # top-right, bottom-left, bottom-right

    for (row, col), (ta, tb, color, lbl_a, lbl_b) in zip(subplot_positions, pair_specs):
        ax = fig.add_subplot(gs[row, col])
        pair = df[(df["tool_a"] == ta) & (df["tool_b"] == tb)]
        x = pair["meth_diff_a"].to_numpy() / 100.0
        y = pair["meth_diff_b"].to_numpy() / 100.0
        _scatter_helper(
            ax, x, y,
            label_x=f"{lbl_a} Δβ",
            label_y=f"{lbl_b} Δβ",
            color=color,
        )
        ax.set_title(f"{lbl_a} × {lbl_b}", loc="left", fontsize=11)

    fig.suptitle(
        "F04 — Per-CpG and DMR-level Δβ concordance on GSE263850",
        fontsize=13, fontweight="bold", x=0.02, y=0.98, ha="left",
    )

    fig.text(
        0.02, -0.01,
        "Top left: union of sites called significant by both epykit lr and methylKit "
        "1.36.0, downsampled to 30 k for plotting. Other panels: DMR pairs matched by "
        "coordinate overlap.",
        fontsize=7, color="#666666", wrap=True,
    )

    plt.subplots_adjust(top=0.88, hspace=0.55)
    save_dual(fig, FIG_OUT / "F04_per_cpg_and_dmr_concordance")


if __name__ == "__main__":
    main()