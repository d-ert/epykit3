"""F02b — Apples-to-apples DMC runtime on the Piao 2021 simulator.

Single-panel head-to-head: epykit engines vs methylKit vs DSS at coverage 10×,
100 k CpGs / scenario, n = 3 vs 3. Bars show seed-median wall-clock; whiskers
span Q1–Q3 over 20 independent simulation seeds. Same input, same metric,
same number of seeds — the cleanest cross-tool runtime comparison the
benchmark supports.

Sources:
  data/study1b_simulator/eval_external_timings_iqr.parquet
      (methylKit + DSS IQR over 20 seeds)
  data/study1b_simulator/timings_simulator.parquet
      (per-seed epykit engines)
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl

from _headline_style import DATA_STUDY1, FIG_OUT, save_dual, setup


PALETTE = {
    "epykit_lr":      "#0F4C81",
    "epykit_lrplus":  "#F39C12",
    "epykit_welch_t": "#5DADE2",
    "epykit_fisher":  "#9BCDD2",
    "methylkit":      "#E07B39",
    "dss":            "#2E8B57",
}
DISPLAY = {
    "epykit_lr":      "epykit lr",
    "epykit_lrplus":  "epykit lr+",
    "epykit_welch_t": "epykit welch_t",
    "epykit_fisher":  "epykit fisher",
    "methylkit":      "methylKit 1.36.0",
    "dss":            "DSS (smoothed)",
}


def _iqr_at_cov10() -> pd.DataFrame:
    """Median + Q1 + Q3 + n_seeds per tool at coverage 10×."""
    base = DATA_STUDY1.parent / "study1b_simulator"

    # epykit — aggregate from per-seed file
    ek = pl.read_parquet(base / "timings_simulator.parquet").to_pandas()
    ek = ek[ek["coverage"] == 10]
    ek_agg = (
        ek.groupby("tool")["wall_s"]
        .agg(
            median="median",
            q1=lambda s: s.quantile(0.25),
            q3=lambda s: s.quantile(0.75),
            n_seeds="count",
        )
        .reset_index()
    )

    # methylKit + DSS — already aggregated
    ext = pl.read_parquet(base / "eval_external_timings_iqr.parquet").to_pandas()
    ext = ext[ext["tool"].isin(["methylkit", "dss"])]
    ext = ext[["tool", "median", "q1", "q3", "n_seeds"]]

    return pd.concat([ek_agg, ext], ignore_index=True)


def main() -> None:
    setup()
    iqr = _iqr_at_cov10()

    order = ["epykit_lr", "epykit_welch_t", "epykit_lrplus",
             "epykit_fisher", "dss", "methylkit"]
    lookup = {row["tool"]: row for _, row in iqr.iterrows()}
    rows = [lookup[t] for t in order if t in lookup]

    fig, ax = plt.subplots(figsize=(11, 6.5))

    y = np.arange(len(rows))
    medians = np.array([r["median"] for r in rows])
    q1s     = np.array([r["q1"]     for r in rows])
    q3s     = np.array([r["q3"]     for r in rows])
    err_lo  = np.maximum(medians - q1s, 0)
    err_hi  = np.maximum(q3s - medians, 0)
    colors  = [PALETTE[r["tool"]] for r in rows]
    labels  = [DISPLAY[r["tool"]] for r in rows]

    ax.barh(
        y, medians,
        xerr=[err_lo, err_hi],
        color=colors,
        edgecolor="white", linewidth=1.0, height=0.65,
        error_kw=dict(ecolor="#333", lw=1.3, capsize=4),
    )

    # Annotate value at end of each bar
    for i, r in enumerate(rows):
        ax.text(
            r["median"] * 1.08, y[i],
            f"{r['median']:.2f} s    (n = {int(r['n_seeds'])} seeds)",
            va="center", fontsize=10,
            color=PALETTE[r["tool"]], fontweight="bold",
        )

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=11.5)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("DMC step wall-clock (seconds, log scale)",
                  fontsize=11)
    ax.set_xlim(0.1, 800)
    ax.grid(axis="x", which="both", linestyle=":", alpha=0.4)

    ax.set_title(
        "F02b — DMC-step runtime on the Piao 2021 simulator  "
        "(coverage 10×, 100 k CpGs, n = 3 vs 3)",
        loc="left", fontsize=12, fontweight="bold",
    )



    plt.subplots_adjust(left=0.20, top=0.92, bottom=0.18)
    save_dual(fig, FIG_OUT / "F02b_simulator_runtime")


if __name__ == "__main__":
    main()
