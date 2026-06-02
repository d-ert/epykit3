"""F03 — Observed FDR under permutation null vs nominal q=0.05.

Forest plot: median observed FDR with Q1/Q3 across shuffles, for each
(engine × dataset) combination. Reference line at the nominal 0.05.

Substitutes the originally proposed p-value Q-Q plot because the null-
calibration data does not retain raw p-value vectors — only per-shuffle
called-site counts, which the summary parquet aggregates to FDR statistics.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from _headline_style import DATA_NULL, FIG_OUT, save_dual, setup


ENGINE_LABEL = {
    "lr": "lr",
    "lr_plus": "lr+",
    "glm": "glm",
    "welch_t": "welch_t",
    "fisher": "fisher",
}

DATASET_LABEL = {
    "gse263850":        "GSE263850 (real WGBS, 22 M sites)",
    "piao_distributed": "Piao 2021 (simulator, 100 k sites)",
    "simulator":        "epykit internal simulator (100 k sites)",
}

ENGINE_COLOR = {
    "lr":      "#0F4C81",
    "lr_plus": "#1F77B4",
    "glm":     "#2E86AB",
    "welch_t": "#5DADE2",
    "fisher":  "#9BCDD2",
}


def main() -> None:
    setup()

    df = pl.read_parquet(DATA_NULL / "summary.parquet").to_pandas()

    df["engine_lbl"] = df["engine"].map(ENGINE_LABEL).fillna(df["engine"])
    df["dataset_lbl"] = df["dataset"].map(DATASET_LABEL).fillna(df["dataset"])
    df = df.sort_values(["dataset_lbl", "engine_lbl"])
    df["row_label"] = df["dataset_lbl"] + " — " + df["engine_lbl"]

    # newest at top
    df = df.iloc[::-1].reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10.5, 5))
    y = np.arange(len(df))
    medians = df["observed_fdr_median"].values
    q1 = df["observed_fdr_q1"].values
    q3 = df["observed_fdr_q3"].values
    colors = df["engine"].map(ENGINE_COLOR).fillna("#888888").values

    err_lower = np.maximum(medians - q1, 0)
    err_upper = np.maximum(q3 - medians, 0)

    ax.errorbar(
        medians, y,
        xerr=[err_lower, err_upper],
        fmt="o",
        markersize=7,
        ecolor="#999999",
        elinewidth=1.2,
        capsize=3,
        markerfacecolor="white",
        markeredgewidth=2,
        zorder=3,
    )
    for i, c in enumerate(colors):
        ax.plot(medians[i], y[i], "o", color=c, markersize=7, zorder=4)

    ax.axvline(0.05, color="#C0392B", linestyle="--", linewidth=1.5,
               label="nominal q = 0.05")
    ax.axvline(0.0, color="#cccccc", linewidth=0.5, zorder=1)

    ax.set_yticks(y)
    ax.set_yticklabels(df["row_label"].values, fontsize=9)

    ax.set_xlabel("Observed false-discovery rate under permutation null")
    ax.set_xscale("symlog", linthresh=1e-6)
    ax.set_xlim(-1e-7, 0.1)
    ax.set_xticks([0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.05])
    ax.set_xticklabels(["0", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "0.05"])

    ax.set_title(
        "F03 — Null-shuffle FDR calibration across engines and datasets",
        loc="left",
    )

    ax.legend(loc="lower right", frameon=False)
    ax.grid(axis="x", linestyle=":", alpha=0.5)

    # Append summary count column
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(
            0.105, y[i],
            f"k={int(row['k_shuffles'])} shuffles",
            fontsize=7, color="#666666", va="center", ha="left",
            transform=ax.get_yaxis_transform(),
        )

    fig.text(
        0.02, -0.06,
        "Points are median observed FDR at q ≤ 0.05 across k label-shuffle "
        "replicates; horizontal bars span Q1–Q3. All non-trivial points sit "
        "well below 0.05 — no engine produces inflated false-discovery rates "
        "under the null hypothesis on any of the three datasets.\n"
        "Source: data/null_calibration/summary.parquet.",
        fontsize=7, color="#666666", wrap=True,
    )

    plt.subplots_adjust(left=0.34, right=0.85)
    save_dual(fig, FIG_OUT / "F03_fdr_calibration_forest")


if __name__ == "__main__":
    main()
