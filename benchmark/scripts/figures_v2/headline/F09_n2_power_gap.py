"""F09 — Replicate scaling and the n = 2 regime.

At n_total = 2 (one case vs one control), per-CpG engines need to fall back
to a 2x2 contingency test because the group-level overdispersion cannot be
estimated. The honest picture on Piao 2021 simulator (effect bin 0.2–0.4):

- plain methylKit and Fisher: TPR ~ 0 (overdispersion degenerate)
- methylKit_TUNED (parameter sweep on overdispersion): TPR = 0.96
- epykit Fisher fallback: TPR = 0.29
- epykit `lr` and `lr+`: no n = 2 rows (require >= 2 per group)

epykit's strength is not the n = 2 raw-power crown — methylKit_tuned wins
that. epykit's strength is scaling, speed, and a clean Python API. This
figure shows the data as it is.

Source: eval_summary_post_phase3.parquet, scenario = dmc_replicate.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from _headline_style import DATA_STUDY1, FIG_OUT, save_dual, setup


def main() -> None:
    setup()

    df = pl.read_parquet(DATA_STUDY1 / "eval_summary_post_phase3.parquet").to_pandas()

    sub = df[
        (df["scenario"] == "dmc_replicate") &
        (df["threshold_kind"] == "qvalue") &
        (df["meth_diff_bin"] == "0.2-0.4")
    ].copy()

    tools = ["epykit_lrplus", "epykit_lr", "epykit_fisher",
             "methylkit", "methylkit_tuned", "dss",
             "radmeth", "biseq", "methylsig", "fisher"]
    tools = [t for t in tools if t in sub["tool"].unique()]
    nvals = [2, 4, 6, 8, 10]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5), sharex=True)

    palette = {
        "epykit_lrplus":   "#F39C12",
        "epykit_lr":       "#0F4C81",
        "epykit_fisher":   "#9BCDD2",
        "methylkit":       "#E07B39",
        "methylkit_tuned": "#D35400",
        "dss":             "#2E8B57",
        "radmeth":         "#9B59B6",
        "biseq":           "#B59B6F",
        "methylsig":       "#C0392B",
        "fisher":          "#7F8C8D",
    }
    display = {
        "epykit_lrplus":   "epykit lr+",
        "epykit_lr":       "epykit lr",
        "epykit_fisher":   "epykit fisher",
        "methylkit":       "methylKit (plain)",
        "methylkit_tuned": "methylKit (tuned)",
        "dss":             "DSS",
        "radmeth":         "RADMeth",
        "biseq":           "BiSeq",
        "methylsig":       "methylSig",
        "fisher":          "Fisher",
    }

    width = 0.10
    x = np.arange(len(nvals))

    for ax, metric, ylabel, ylim in [
        (ax1, "tpr", "True positive rate (TPR)", (0, 1.02)),
        (ax2, "fpr", "False positive rate (FPR)", (-0.001, 0.04)),
    ]:
        for i, tool in enumerate(tools):
            rows = sub[sub["tool"] == tool].set_index("parameter_value")
            vals = [rows.loc[n, metric] if n in rows.index else np.nan for n in nvals]
            offset = (i - (len(tools) - 1) / 2) * width
            ax.bar(x + offset, vals, width=width,
                   color=palette.get(tool, "#888888"),
                   label=display.get(tool, tool),
                   edgecolor="white", linewidth=0.5)

        ax.set_xticks(x)
        ax.set_xticklabels(nvals)
        ax.set_xlabel("n_total (samples per group × 2)")
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        ax.grid(axis="y", linestyle=":", alpha=0.4)

        # Highlight the n=2 column
        ax.axvspan(x[0] - 0.5, x[0] + 0.5, color="#FDF6E3", zorder=0)
        if metric == "fpr":
            ax.axhline(0.05, color="#cccccc", linestyle=":", linewidth=1)

    # Annotate the n=2 regime: methylKit_tuned wins; plain methylKit fails;
    # epykit's only option at n=2 is the Fisher fallback.
    ax1.annotate(
        "At n=2, plain methylKit\noverdispersion is degenerate (TPR≈0).\n"
        "methylKit_tuned wins; epykit's\nfallback is the Fisher engine.\n"
        "lr / lr+ require n ≥ 4.",
        xy=(x[0], 0.30), xytext=(x[0] + 0.5, 0.55),
        fontsize=8.5, color="#0F4C81",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#EAF1F8",
                  edgecolor="#0F4C81", linewidth=0.8),
    )

    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(tools),
               bbox_to_anchor=(0.5, -0.03), frameon=False, fontsize=8.5)

    fig.suptitle(
        "F09 — Replicate scaling on the Piao 2021 simulator (cov 10×, bin 0.2–0.4)",
        fontsize=12, fontweight="bold", x=0.02, ha="left",
    )

    fig.text(
        0.02, -0.10,
        "Highlighted n=2 column: only methylKit_tuned and epykit's Fisher fallback "
        "produce non-zero recall. All tools converge at n ≥ 4. epykit's headline "
        "value is competitive scaling (lr/lr+ match methylKit/DSS at n ≥ 4 with "
        "much lower runtime — see F02, F07, F08) and a Python-native API; the n=2 "
        "regime is the one place plain methylKit is clearly outperformed.\n"
        "Source: eval_summary_post_phase3.parquet, scenario = dmc_replicate.",
        fontsize=7, color="#666666", wrap=True,
    )

    plt.subplots_adjust(bottom=0.20, top=0.90)
    save_dual(fig, FIG_OUT / "F09_n2_power_gap")


if __name__ == "__main__":
    main()
