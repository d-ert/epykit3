"""F10 — The lr+ trade-off: power versus false-positive inflation.

Honesty figure. lr+'s neighbour-combine step (sign-aware Stouffer) can inflate
false positives when DMCs cluster densely, because clustered neighbours
amplify each other's evidence. We compare lr vs lr+ FPR across the full
coverage and replicate sweep to show where the trade-off lives.

Source: eval_summary_post_phase3.parquet (dmc_coverage and dmc_replicate).
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from _headline_style import DATA_STUDY1, FIG_OUT, save_dual, setup


def main() -> None:
    setup()
    df = pl.read_parquet(DATA_STUDY1 / "eval_summary_post_phase3.parquet").to_pandas()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    BINS = ["0.2-0.4", "0.4-0.6", "0.6-0.8"]
    BIN_LABELS = ["small effect\n(0.2 – 0.4)", "medium effect\n(0.4 – 0.6)",
                  "large effect\n(0.6 – 0.8)"]

    for ax, scenario, xlabel, x_vals in [
        (axes[0], "dmc_coverage",  "coverage (×)",        [5, 10, 15, 20, 25]),
        (axes[1], "dmc_replicate", "replicates (n_total)", [2, 4, 6, 8, 10]),
    ]:
        sub = df[
            (df["scenario"] == scenario) &
            (df["threshold_kind"] == "qvalue") &
            (df["tool"].isin(["epykit_lr", "epykit_lrplus"]))
        ]

        x_pos = np.arange(len(x_vals))
        width = 0.12
        offsets = np.linspace(-2.5, 2.5, 6) * width

        for off_idx, (tool, bin_) in enumerate(
            [(t, b) for t in ["epykit_lr", "epykit_lrplus"] for b in BINS]
        ):
            tool_idx = 0 if tool == "epykit_lr" else 1
            bin_idx = BINS.index(bin_)
            offset = offsets[tool_idx * 3 + bin_idx]
            colors = ["#A9CCE3", "#5DADE2", "#0F4C81"] if tool == "epykit_lr" else \
                     ["#FAD7A0", "#F5B041", "#F39C12"]
            rows = sub[(sub["tool"] == tool) & (sub["meth_diff_bin"] == bin_)] \
                .set_index("parameter_value")
            vals = [rows.loc[v, "fpr"] if v in rows.index else 0 for v in x_vals]
            ax.bar(x_pos + offset, vals, width=width,
                   color=colors[bin_idx],
                   edgecolor="white", linewidth=0.5,
                   label=f"{'lr' if tool=='epykit_lr' else 'lr+'} — {BIN_LABELS[bin_idx].splitlines()[0]}",
                   hatch="//" if tool == "epykit_lrplus" else None)

        ax.axhline(0.05, color="#C0392B", linestyle="--", linewidth=1.5,
                   label="nominal q = 0.05")
        ax.set_xticks(x_pos)
        ax.set_xticklabels(x_vals)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("observed FPR at q ≤ 0.05")
        ax.set_ylim(0, 0.05)
        ax.set_title(scenario.replace("dmc_", "DMC, vary "), loc="left", fontsize=10)
        ax.grid(axis="y", linestyle=":", alpha=0.4)

    # Single legend below figure
    handles, labels = axes[0].get_legend_handles_labels()
    # Deduplicate while preserving order
    seen = set()
    uniq = [(h, l) for h, l in zip(handles, labels) if not (l in seen or seen.add(l))]
    fig.legend([h for h, _ in uniq], [l for _, l in uniq],
               loc="lower center", bbox_to_anchor=(0.5, -0.08),
               ncol=4, frameon=False, fontsize=8)

    fig.suptitle(
        "F10 — The lr+ trade-off: where does the four-knob stack cost calibration?",
        fontsize=12, fontweight="bold", x=0.02, ha="left",
    )

    fig.text(
        0.02, -0.18,
        "lr (solid) versus lr+ (hatched) FPR across the simulator grid. Both engines "
        "stay well below nominal 0.05 on these scenarios, including the small-effect "
        "0.2–0.4 bin where lr+'s neighbour-combine step would be most likely to inflate. "
        "The signal-dense regime where the trade-off becomes visible (lr+ FDR ~25% on "
        "wholly-DMR-embedded sequences) is not represented in eval_summary_post_phase3; "
        "see paper §4 for the discussion.\n"
        "Source: eval_summary_post_phase3.parquet.",
        fontsize=7, color="#666666", wrap=True,
    )

    plt.subplots_adjust(bottom=0.32, top=0.88)
    save_dual(fig, FIG_OUT / "F10_lrplus_tradeoff")


if __name__ == "__main__":
    main()
