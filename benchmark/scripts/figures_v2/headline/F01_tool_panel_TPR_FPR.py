"""F01 — Cross-tool TPR/FPR/F1 panel across coverage on the Piao 2021 simulator.

Three rows (TPR, FPR, F1) x five columns (coverage 5x..25x). One line per tool,
stratified by effect-size bin 0.2-0.4 (the hardest, most discriminating bin —
matches Piao 2021's reported figure).

Tools without F1 (transcribed baselines) are omitted from the F1 row.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import matplotlib.pyplot as plt

from _headline_style import (
    DATA_STUDY1, FIG_OUT, setup, save_dual, color_for, label_for,
)


def main() -> None:
    setup()

    df = pl.read_parquet(DATA_STUDY1 / "eval_summary_post_phase3.parquet")

    sub = df.filter(
        (pl.col("scenario") == "dmc_coverage")
        & (pl.col("threshold_kind") == "qvalue")
        & (pl.col("meth_diff_bin") == "0.2-0.4")
    ).to_pandas()

    sub["display_tool"] = sub["tool"].map(label_for)

    tool_order = [
        "epykit_lrplus", "epykit_lr", "epykit_welch_t", "epykit_fisher",
        "methylkit", "methylkit_tuned", "dss", "radmeth", "methylsig",
        "biseq", "fisher",
    ]
    tools_present = [t for t in tool_order if t in sub["tool"].unique()]

    coverages = sorted(sub["parameter_value"].unique())

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    metrics = [
        ("tpr", "True positive rate (TPR)", (0, 1.02)),
        ("fpr", "False positive rate (FPR)", (-0.002, 0.045)),
    ]

    for ax, (metric, ylabel, ylim) in zip(axes, metrics):
        for tool in tools_present:
            rows = sub[sub["tool"] == tool].sort_values("parameter_value")
            if rows.empty or rows[metric].isna().all():
                continue
            ax.plot(
                rows["parameter_value"],
                rows[metric],
                marker="o",
                linewidth=2,
                markersize=6,
                color=color_for(tool),
                label=label_for(tool),
            )
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        if metric == "fpr":
            ax.axhline(0.05, color="#bbbbbb", linestyle="--", linewidth=1,
                       label="nominal q=0.05")

    axes[-1].set_xlabel("Sequencing coverage (× per CpG)")
    axes[-1].set_xticks(coverages)
    axes[0].set_title(
        "F01 — Cross-tool DMC accuracy on the Piao 2021 simulator "
        "(effect-size bin 0.2–0.4, n=3 vs 3)",
        loc="left",
    )

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        loc="center right",
        bbox_to_anchor=(1.0, 0.5),
        frameon=False,
        title="tool",
    )

    fig.text(
        0.02, -0.03,
        "Source: eval_summary_post_phase3.parquet (scenario=dmc_coverage, threshold_kind=qvalue, "
        "meth_diff_bin=0.2–0.4). Baselines BiSeq/methylSig/RADMeth/DSS/Fisher transcribed from "
        "Piao et al. 2021 Table S1.",
        fontsize=7, color="#666666",
    )

    plt.subplots_adjust(right=0.78)
    save_dual(fig, FIG_OUT / "F01_tool_panel_TPR_FPR")


if __name__ == "__main__":
    main()
