"""Cross-benchmark summary figures for the consolidated epykit report.

Loads:
* study1: eval_summary.parquet + timings.parquet
* study2: methylkit_results/timings.tsv + (study1 timings shared, since the
  Study 2 simulator and Study 1 simulator are the same input dataset)
* study3: benchmark/step_benchmarks.csv + benchmark/run_summary.csv +
  Pearson r / Jaccard hardcoded from the comparison report (no source file)

Writes to ../figures/summary/:
* S1_runtime_across_studies.png  - grouped bar, log-y
* S2_accuracy_summary.png        - 2-panel TPR vs coverage
* S3_agreement_matrix.png        - effect-size scatter + heatmap

Hardcoded numbers (methylKit GSE263850 run-summary fields) are sourced from
data/study3/benchmark/run_summary.csv for the epykit row, and from the
comparison report (review of methylkit_realResults/run_summary.csv) for the
methylKit row, since only the epykit run_summary was copied into
FINAL_REPORT/data.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
OUT = ROOT / "figures" / "summary"
OUT.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------------
# Hardcoded numbers from the GSE263850 comparison report. The methylKit row
# is not in FINAL_REPORT/data because we only copied epykit_results/. The
# numbers below are from the report's Section 1 headline table.
# ----------------------------------------------------------------------------
STUDY3_METHYLKIT = {
    "pipeline_wall_sec": 13033.0,
    "peak_rss_mb": 48001.8,
    "n_dmcs_sig_q05": 51792,
    "n_dmrs_lenient": 2661,
}
STUDY3_AGREEMENT = {
    "pearson_r_dmc": 0.9936,
    "spearman_rho_dmc": 0.9831,
    "direction_agreement_dmc": 0.9405,
    "pearson_r_dmr": 0.9970,
    "direction_agreement_dmr": 0.9182,
    "jaccard_dmc_sig": 0.234,
    "jaccard_dmr_lenient": 0.473,
    "jaccard_dmr_strict": 0.530,
    "recall_methylkit_by_epykit_dmc": 0.303,
}


def load_study1_timings() -> pd.DataFrame:
    df = pd.read_parquet(DATA / "study1" / "timings.parquet")
    df = df[df["ok"]].copy()
    return df


def load_study1_eval() -> pd.DataFrame:
    df = pd.read_parquet(DATA / "study1" / "eval_summary.parquet")
    return df


def load_study2_timings() -> pd.DataFrame:
    df = pd.read_csv(
        DATA / "study2" / "methylkit_results" / "timings.tsv", sep="\t"
    )
    return df


def load_study3_steps() -> pd.DataFrame:
    return pd.read_csv(DATA / "study3" / "benchmark" / "step_benchmarks.csv")


def load_study3_summary() -> pd.DataFrame:
    return pd.read_csv(DATA / "study3" / "benchmark" / "run_summary.csv")


# ----------------------------------------------------------------------------
# Figure S1: runtime across the three studies
# ----------------------------------------------------------------------------
def fig_s1_runtime():
    """Wall-clock totals: epykit vs methylKit on the two studies that ran
    both pipelines on the same machine with measured timings.

    Study 1 is excluded because its baseline methylKit numbers come from the
    Piao 2021 publication (transcribed from supplementary tables), without
    matching wall-clock measurements on our hardware. Study 2 is the proper
    head-to-head on the same simulated grid.

    For Study 2 we sum epykit's lr-engine timings plus all DMR timings from
    study1/timings.parquet (same simulator, same machine) and compare with
    the locally run methylKit timings from study2/methylkit_results/.

    For Study 3 we use the actual GSE263850 pipeline_wall_sec values.
    """
    t1 = load_study1_timings()
    t2 = load_study2_timings()
    s3 = load_study3_summary()

    # Study 2 epykit grid: lr DMC scenarios + all DMR scenarios.
    epykit_total_grid = (
        t1[t1["test"] == "lr"]["elapsed_s"].sum()
        + t1[t1["scenario"] == "dmr_coverage"]["elapsed_s"].sum()
    )

    # methylKit times from Study 2's run (DMC + DMR rows, summed).
    methylkit_total_grid = t2["wall_s"].sum()

    # Study 3: from the comparison report (epykit summary CSV + hardcoded mk)
    epykit_s3 = float(s3["pipeline_wall_sec"].iloc[0])
    methylkit_s3 = STUDY3_METHYLKIT["pipeline_wall_sec"]

    studies = [
        "Study 2\nsimulated grid\n(epykit lr engine + DMR)",
        "Study 3\nreal WGBS (GSE263850)\n(full pipeline)",
    ]
    epykit_vals = [epykit_total_grid, epykit_s3]
    methylkit_vals = [methylkit_total_grid, methylkit_s3]
    speedups = [
        methylkit_vals[i] / epykit_vals[i] for i in range(len(studies))
    ]

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(studies))
    w = 0.35
    bars_e = ax.bar(
        x - w / 2, epykit_vals, w, label="epykit", color="#1f77b4"
    )
    bars_m = ax.bar(
        x + w / 2,
        methylkit_vals,
        w,
        label="methylKit",
        color="#d62728",
    )
    ax.set_yscale("log")
    ax.set_ylabel("Total wall-clock time (s, log scale)")
    ax.set_xticks(x)
    ax.set_xticklabels(studies)
    ax.set_title(
        "Figure S1. Wall-clock runtime: epykit vs methylKit, per study"
    )
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))

    for bars in (bars_e, bars_m):
        for b in bars:
            h = b.get_height()
            ax.annotate(
                f"{h:.0f}s" if h > 60 else f"{h:.1f}s",
                xy=(b.get_x() + b.get_width() / 2, h),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    # Boost y axis upper bound so speedup labels fit above bars.
    top = max(methylkit_vals) * 5
    ax.set_ylim(top=top)
    for i, sp in enumerate(speedups):
        ax.annotate(
            f"{sp:.1f}× speedup",
            xy=(i, max(epykit_vals[i], methylkit_vals[i]) * 1.8),
            ha="center",
            fontsize=10,
            fontweight="bold",
            color="#2ca02c",
        )

    ax.grid(True, axis="y", which="both", alpha=0.3)
    fig.tight_layout()
    out = OUT / "S1_runtime_across_studies.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


# ----------------------------------------------------------------------------
# Figure S2: accuracy summary across the two simulated studies
# ----------------------------------------------------------------------------
def fig_s2_accuracy():
    """TPR vs coverage for epykit `lr` and methylKit, all effect sizes
    combined, q < 0.05.

    Both studies use the same Piao 2021 simulator, so the lines should
    overlap. Study 1's methylKit numbers come from Piao's published table;
    Study 2's are from a fresh local run on the same machine as epykit.
    """
    eval1 = load_study1_eval()

    # Filter to the q<0.05 cut, all effects combined, coverage grid.
    cut = (
        (eval1["scenario"] == "dmc_coverage")
        & (eval1["threshold_kind"] == "qvalue")
        & (eval1["threshold"] == 0.05)
        & (eval1["meth_diff_bin"] == "all")
    )
    e = eval1[cut].copy()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)

    # Panel A: Study 1 — epykit lr / lr+ vs methylKit (transcribed)
    axA = axes[0]
    for tool, color, marker, label in [
        ("epykit_lr", "#1f77b4", "o", "epykit lr"),
        ("epykit_lrplus", "#9467bd", "s", "epykit lr+"),
        ("methylkit", "#d62728", "^", "methylKit (Piao 2021 tables)"),
        ("dss", "#ff7f0e", "v", "DSS"),
        ("radmeth", "#2ca02c", "x", "RADMeth"),
    ]:
        sub = e[e["tool"] == tool].sort_values("parameter_value")
        if not len(sub):
            continue
        axA.plot(
            sub["parameter_value"],
            sub["tpr"],
            color=color,
            marker=marker,
            label=label,
            linewidth=1.6,
            markersize=7,
        )
    axA.set_xlabel("Coverage (×)")
    axA.set_ylabel("TPR at q < 0.05 (all effects)")
    axA.set_title("(a) Study 1 — Panel comparison")
    axA.set_ylim(0, 1.05)
    axA.legend(loc="lower right", fontsize=9)
    axA.grid(True, alpha=0.3)

    # Panel B: Study 2 — local epykit vs local methylKit
    # The HEAD_TO_HEAD numbers from the manuscript:
    cov = [5, 10, 15, 20, 25]
    tpr_epykit = [0.849, 0.944, 0.984, 0.991, 0.993]
    tpr_methylkit = [0.849, 0.944, 0.984, 0.991, 0.993]
    axB = axes[1]
    axB.plot(
        cov,
        tpr_epykit,
        color="#1f77b4",
        marker="o",
        label="epykit lr (local)",
        linewidth=1.6,
        markersize=7,
    )
    axB.plot(
        cov,
        tpr_methylkit,
        color="#d62728",
        marker="^",
        label="methylKit (local)",
        linewidth=1.6,
        markersize=7,
        linestyle="--",
    )
    axB.set_xlabel("Coverage (×)")
    axB.set_title("(b) Study 2 — Head-to-head, same machine")
    axB.set_ylim(0, 1.05)
    axB.legend(loc="lower right", fontsize=9)
    axB.grid(True, alpha=0.3)
    axB.annotate(
        "Lines overlap — TPR identical to 3 dp",
        xy=(15, 0.92),
        ha="center",
        fontsize=10,
        color="#555",
    )

    fig.suptitle(
        "Figure S2. DMC accuracy across the two simulated studies (3 vs 3 design)",
        y=1.02,
        fontsize=12,
    )
    fig.tight_layout()
    out = OUT / "S2_accuracy_summary.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


# ----------------------------------------------------------------------------
# Figure S3: agreement matrix on real data
# ----------------------------------------------------------------------------
def fig_s3_agreement():
    """For Study 3 (real GSE263850), show:
    - Left: a schematic of the effect-size scatter result (we don't load
      the 15 M raw points; we draw a stylised version with the actual r).
    - Right: a heatmap of agreement statistics (Pearson, Spearman,
      direction, Jaccard) for DMCs and DMRs.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left: stylised scatter with reported r
    ax = axes[0]
    # Generate stylised correlated points for illustration
    rng = np.random.default_rng(42)
    n = 3000
    x = rng.normal(0, 0.3, n)
    noise = rng.normal(0, 0.04, n)
    y = x + noise  # Pearson ~ 0.994
    ax.scatter(x, y, s=2, alpha=0.3, color="#1f77b4")
    ax.plot([-1, 1], [-1, 1], "k--", linewidth=0.8, alpha=0.6, label="y = x")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_xlabel("meth_diff — methylKit")
    ax.set_ylabel("meth_diff — epykit")
    ax.set_title("(a) DMC effect-size agreement, GSE263850")
    ax.annotate(
        f"Pearson r = {STUDY3_AGREEMENT['pearson_r_dmc']:.4f}\n"
        f"Spearman ρ = {STUDY3_AGREEMENT['spearman_rho_dmc']:.4f}\n"
        f"Direction agree = {STUDY3_AGREEMENT['direction_agreement_dmc']:.2%}\n"
        f"(15,597,046 shared CpGs)",
        xy=(0.05, 0.95),
        xycoords="axes fraction",
        va="top",
        ha="left",
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="#888"),
        fontsize=9,
    )
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.text(
        0.5,
        -0.18,
        "Stylised — full scatter has 15.6 M points; see study3 figure 08B for the real plot",
        transform=ax.transAxes,
        ha="center",
        fontsize=7,
        style="italic",
        color="#666",
    )

    # Right: agreement heatmap
    ax = axes[1]
    metrics = [
        "Pearson r (effect)",
        "Direction agree",
        "Jaccard (sig calls)",
        "Recall of methylKit",
    ]
    dmc_vals = [
        STUDY3_AGREEMENT["pearson_r_dmc"],
        STUDY3_AGREEMENT["direction_agreement_dmc"],
        STUDY3_AGREEMENT["jaccard_dmc_sig"],
        STUDY3_AGREEMENT["recall_methylkit_by_epykit_dmc"],
    ]
    dmr_vals = [
        STUDY3_AGREEMENT["pearson_r_dmr"],
        STUDY3_AGREEMENT["direction_agreement_dmr"],
        STUDY3_AGREEMENT["jaccard_dmr_lenient"],
        np.nan,  # recall not in summary for DMR
    ]
    data = np.array([dmc_vals, dmr_vals]).T  # rows = metrics, cols = DMC/DMR

    im = ax.imshow(
        data, cmap="RdYlGn", vmin=0.0, vmax=1.0, aspect="auto"
    )
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["DMCs (CpG-level)", "DMRs (500 bp tiles)"])
    ax.set_yticks(range(len(metrics)))
    ax.set_yticklabels(metrics)
    for i in range(len(metrics)):
        for j in range(2):
            v = data[i, j]
            label = "n/a" if np.isnan(v) else f"{v:.3f}"
            color = "white" if (not np.isnan(v) and (v < 0.3 or v > 0.85)) else "black"
            ax.text(j, i, label, ha="center", va="center", fontsize=10, color=color)
    ax.set_title("(b) Agreement statistics on real data")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="value")

    fig.suptitle(
        "Figure S3. epykit vs methylKit on real WGBS (GSE263850): "
        "same biology, different operating point",
        y=1.02,
        fontsize=12,
    )
    fig.tight_layout()
    out = OUT / "S3_agreement_matrix.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def main():
    fig_s1_runtime()
    fig_s2_accuracy()
    fig_s3_agreement()


if __name__ == "__main__":
    main()
