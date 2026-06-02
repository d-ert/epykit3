"""F08 — Pareto frontier: F1 versus runtime on the Piao 2021 simulator.

Replaces the originally proposed 13-tool heatmap (which reads as marketing).
This scatter shows wall-clock vs F1 at coverage 5×, the hardest scenario.
Only epykit engines have direct F1 measurements in eval_summary_post_phase3
(non-epykit baselines have TPR/FPR transcribed from Piao but no F1); we
approximate F1 from TPR/FPR for baselines using the canonical formula with
a fixed prior on positive-class prevalence.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl

from _headline_style import (
    DATA_STUDY1, DATA_STUDY2, FIG_OUT, color_for, label_for, save_dual, setup,
)


def main() -> None:
    setup()
    df = pl.read_parquet(DATA_STUDY1 / "eval_summary_post_phase3.parquet").to_pandas()

    # --- F1 at cov=5x for all tools (use 0.2-0.4 bin for apples-to-apples) ---
    snap = df[
        (df["scenario"] == "dmc_coverage") &
        (df["parameter_value"] == 5) &
        (df["threshold_kind"] == "qvalue") &
        (df["meth_diff_bin"] == "0.2-0.4")
    ].copy()

    # Compute F1 from TPR & FPR for baseline tools (assumes fixed positive prevalence
    # of 0.2 = 20 k true DMCs out of 100 k total in this slice).
    PREV = 0.2
    for i, row in snap.iterrows():
        if pd.isna(row.get("f1")) or row["f1"] is None:
            tpr = row["tpr"]
            fpr = row["fpr"]
            if pd.isna(tpr) or pd.isna(fpr):
                continue
            tp_rate = tpr * PREV
            fp_rate = fpr * (1 - PREV)
            denom = tp_rate + fp_rate + (1 - tpr) * PREV
            snap.at[i, "f1"] = 2 * tp_rate / denom if denom else np.nan

    # --- runtimes ---
    # epykit engines: timings_table.csv (per-engine median wall-clock on simulator)
    ek_times = pd.read_csv(DATA_STUDY1 / "timings_table.csv")
    ek_times["wall_min"] = ek_times["median_wall_s"] / 60.0
    ek_lookup = dict(zip(ek_times["tool"], ek_times["wall_min"]))

    # methylKit: study2/timings.tsv at dmc_coverage, cov=5
    mk_times = pd.read_csv(DATA_STUDY2 / "methylkit_results" / "timings.tsv", sep="\t")
    mk_row = mk_times[(mk_times["scenario"] == "dmc_coverage") & (mk_times["value"] == 5)]
    mk_wall_min = float(mk_row["wall_s"].iloc[0]) / 60.0 if not mk_row.empty else np.nan

    runtime_min = {
        **ek_lookup,
        "methylkit":       mk_wall_min,
        "methylkit_tuned": mk_wall_min,   # tuned uses same wall-clock
        # Approximate baselines from Piao 2021 (CLI/R tools, simulator-scale)
        "dss":       12.0,
        "radmeth":   20.0,
        "biseq":     90.0,
        "methylsig": 18.0,
        "fisher":     3.0,
    }
    snap["runtime_min"] = snap["tool"].map(runtime_min)

    snap = snap.dropna(subset=["runtime_min", "f1"])

    fig, ax = plt.subplots(figsize=(11, 6.5))
    for _, row in snap.iterrows():
        x, y = row["runtime_min"], row["f1"]
        ax.scatter(x, y, s=110, c=color_for(row["tool"]),
                   edgecolor="white", linewidth=1.2, zorder=3)
        ax.annotate(
            label_for(row["tool"]),
            xy=(x, y), xytext=(7, 4), textcoords="offset points",
            fontsize=9, color=color_for(row["tool"]), fontweight="bold",
        )
    ax.set_xscale("log")
    ax.set_xlabel("wall-clock per scenario (minutes, log scale)")
    ax.set_ylabel("F1 at q ≤ 0.05 (effect-size bin 0.2–0.4, cov 5×)")
    ax.set_xlim(0.005, 200)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(linestyle=":", alpha=0.4)
    ax.set_title(
        "F08 — Speed vs accuracy at the hardest scenario (5× coverage, n = 3 vs 3)",
        loc="left",
    )

    # mark the upper-left Pareto-frontier region
    ax.fill_betweenx([0.85, 1.05], 0.005, 1.5, alpha=0.07, color="#0F4C81")
    ax.text(
        0.013, 0.96, "Pareto frontier\n(fast + accurate)",
        fontsize=9, color="#0F4C81", style="italic",
    )

    fig.text(
        0.02, -0.03,
        "Baseline runtimes for DSS/RADMeth/BiSeq/methylSig are order-of-magnitude "
        "estimates from Piao 2021 reproduction notes; methylKit and epykit are "
        "measured directly. F1 is computed from TPR and FPR with prevalence = "
        "0.20 for baselines that don't report F1 directly.",
        fontsize=7, color="#666666", wrap=True,
    )

    save_dual(fig, FIG_OUT / "F08_pareto_runtime_vs_f1")


if __name__ == "__main__":
    main()
