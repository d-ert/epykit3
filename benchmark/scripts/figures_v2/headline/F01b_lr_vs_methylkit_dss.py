"""F01b — Three-tool TPR/FPR/F1 panel: epykit lr vs methylKit vs DSS.

Cleaner head-to-head version of F01. Same data slice (Piao 2021 simulator,
effect-size bin 0.2-0.4, n=3 vs 3, q ≤ 0.05), but only the three tools
that anchor the comparison.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import matplotlib.pyplot as plt

from _headline_style import DATA_STUDY1, FIG_OUT, setup, save_dual


PALETTE = {
    "epykit_lr": "#0F4C81",
    "methylkit": "#E07B39",
    "dss":       "#2E8B57",
}
DISPLAY = {
    "epykit_lr": "epykit lr",
    "methylkit": "methylKit",
    "dss":       "DSS",
}
TOOLS = ["epykit_lr", "methylkit", "dss"]


def main() -> None:
    setup()

    df = pl.read_parquet(DATA_STUDY1 / "eval_summary_post_phase3.parquet")
    sub = df.filter(
        (pl.col("scenario") == "dmc_coverage")
        & (pl.col("threshold_kind") == "qvalue")
        & (pl.col("meth_diff_bin") == "0.2-0.4")
    ).to_pandas()

    coverages = sorted(sub["parameter_value"].unique())

    fig, axes = plt.subplots(2, 1, figsize=(9.5, 7.5), sharex=True)
    metrics = [
        ("tpr", "True positive rate (TPR)", (0, 1.04)),
        ("fpr", "False positive rate (FPR)", (-0.002, 0.045)),
    ]

    for ax, (metric, ylabel, ylim) in zip(axes, metrics):
        for tool in TOOLS:
            rows = sub[sub["tool"] == tool].sort_values("parameter_value")
            if rows.empty or rows[metric].isna().all():
                continue
            ax.plot(
                rows["parameter_value"],
                rows[metric],
                marker="o",
                linewidth=2.4,
                markersize=9,
                color=PALETTE[tool],
                label=DISPLAY[tool],
            )
            # Value labels at each coverage
            for x, y in zip(rows["parameter_value"], rows[metric]):
                if np.isfinite(y):
                    ax.annotate(
                        f"{y:.3f}",
                        xy=(x, y),
                        xytext=(0, 8 if tool == "epykit_lr" else -14),
                        textcoords="offset points",
                        ha="center", fontsize=8,
                        color=PALETTE[tool],
                    )

        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        if metric == "fpr":
            ax.axhline(0.05, color="#C0392B", linestyle="--", linewidth=1,
                       label="nominal q = 0.05")

    axes[-1].set_xlabel("Sequencing coverage (× per CpG)")
    axes[-1].set_xticks(coverages)
    axes[0].set_title(
        "F01b — epykit lr  vs  methylKit  vs  DSS  on the Piao 2021 simulator\n"
        "(effect-size bin 0.2–0.4, n = 3 vs 3, q ≤ 0.05)",
        loc="left",
    )

    handles_t, labels_t = axes[0].get_legend_handles_labels()
    handles_f, labels_f = axes[1].get_legend_handles_labels()
    axes[0].legend(handles_t, labels_t, loc="lower right", frameon=False, fontsize=10)
    axes[1].legend(handles_f, labels_f, loc="upper right", frameon=False, fontsize=9)

    fig.text(
        0.02, -0.02,
        "epykit lr matches methylKit on TPR from 10× upward and beats it at 5× "
        "(0.835 vs 0.266). DSS underperforms in this low-effect / low-coverage "
        "slice but catches up by 15×. All three keep FPR well below the nominal "
        "0.05 cutoff at every coverage.\n"
        "Source: eval_summary_post_phase3.parquet  •  dmc_coverage  •  qvalue  •  bin 0.2–0.4.",
        fontsize=8, color="#666666", wrap=True,
    )

    plt.subplots_adjust(top=0.90, hspace=0.12)
    save_dual(fig, FIG_OUT / "F01b_lr_vs_methylkit_dss")


if __name__ == "__main__":
    main()
