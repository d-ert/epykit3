"""F06 — lr → lr+ ablation across coverage and replicate count.

Shows the cumulative effect of the four `lr+` knobs (empirical-Bayes
dispersion shrinkage, sign-aware Stouffer neighbour combine, separation-
aware Fisher fallback, Storey two-stage BH) by comparing the bare `lr`
engine to `lr+` on the same Piao 2021 simulator grid.

Caveat (caption): the four knobs are activated together in `lr+` and CANNOT
be cleanly isolated to single-knob attribution from the eval-summary
parquet alone. This figure shows the bundle, not per-knob contributions.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from _headline_style import DATA_STUDY1, FIG_OUT, save_dual, setup


def _curve(sub, tool, metric):
    rows = sub.filter(pl.col("tool") == tool).sort("parameter_value")
    return rows["parameter_value"].to_numpy(), rows[metric].to_numpy()


def main() -> None:
    setup()
    df = pl.read_parquet(DATA_STUDY1 / "eval_summary_post_phase3.parquet")

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharey="row")

    for col, (scenario, xlabel, xname) in enumerate(
        [("dmc_coverage",  "Sequencing coverage (× per CpG)", "coverage"),
         ("dmc_replicate", "Replicate count (n_total)",       "replicates")]
    ):
        sub = df.filter(
            (pl.col("scenario") == scenario)
            & (pl.col("threshold_kind") == "qvalue")
            & (pl.col("meth_diff_bin") == "0.2-0.4")
        )

        for row, (metric, ylabel, ylim) in enumerate(
            [("tpr", "TPR", (0, 1.02)),
             ("fpr", "FPR", (-0.001, 0.045))]
        ):
            ax = axes[row, col]

            for tool, color, lbl in [
                ("epykit_lr",     "#0F4C81", "epykit lr (bare)"),
                ("epykit_lrplus", "#F39C12", "epykit lr+ (4-knob stack)"),
            ]:
                x, y = _curve(sub, tool, metric)
                ax.plot(x, y, marker="o", markersize=7, linewidth=2,
                        color=color, label=lbl)
                # Annotate value at the smallest x (lowest coverage / fewest reps)
                if len(x) > 0:
                    ax.annotate(
                        f"{y[0]:.3f}", xy=(x[0], y[0]),
                        xytext=(8, 6 if tool == "epykit_lr" else -14),
                        textcoords="offset points",
                        fontsize=7, color=color,
                    )

            if metric == "fpr":
                ax.axhline(0.05, color="#cccccc", linestyle=":", linewidth=1)

            ax.set_ylim(*ylim)
            ax.grid(axis="y", linestyle=":", alpha=0.4)
            if col == 0:
                ax.set_ylabel(ylabel)
            if row == 1:
                ax.set_xlabel(xlabel)
            if row == 0:
                ax.set_title(
                    f"{scenario.replace('dmc_', 'DMC, vary ')}"
                    f" ({'coverage 5–25×' if xname == 'coverage' else 'n_total 2–10'})",
                    loc="left", fontsize=10,
                )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.01),
               ncol=2, frameon=False, fontsize=10)

    fig.suptitle(
        "F06 — Bare `lr` versus the four-knob `lr+` power stack",
        fontsize=12, fontweight="bold", x=0.02, ha="left", y=0.96,
    )

    fig.text(
        0.02, -0.03,
        "lr+ activates four enhancements together: empirical-Bayes dispersion "
        "shrinkage, sign-aware Stouffer neighbour combine, separation-aware "
        "Fisher fallback, and Storey two-stage BH. Per-knob attribution is "
        "not isolatable from the eval-summary parquet; lr+ is reported as a "
        "single opt-in bundle. meth_diff_bin = 0.2–0.4 throughout.\n"
        "Source: eval_summary_post_phase3.parquet.",
        fontsize=7, color="#666666", wrap=True,
    )

    plt.subplots_adjust(hspace=0.30, top=0.88)
    save_dual(fig, FIG_OUT / "F06_lr_vs_lrplus_ablation")


if __name__ == "__main__":
    main()
