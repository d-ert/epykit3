"""plot_phi_sweep.py -- figures for the dispersion (phi) sweep (review M1).

Reads eval_phi_sweep_{per_cell,iqr}.parquet (written by run_phi_sweep.py) and
draws:
  * phi-vs-TPR  (per-CpG DMC tools; median line + IQR band across seeds)
  * phi-vs-FPR  (same)
  * phi-vs-peak-RSS (memory; backs the paper's peak-memory claim across phi)
  * (DMR callers dmrseq/bsmooth overlaid on TPR/FPR, dashed, with caveat)

The x-axis is the simulator ICC ``rho``; a secondary top axis shows the implied
Pearson overdispersion ``phi = 1 + (coverage-1)*rho`` at the test coverage so
the realistic-WGBS regime (phi ~ 1.5-5) is legible. Real-WGBS dispersion band
is shaded for context.

Usage:
    python plot_phi_sweep.py [--coverage 10] [--out-dir <dir>]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import polars as pl  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SIM_ROOT = ROOT / "benchmark" / "data" / "study1b_simulator"
DEFAULT_OUT = ROOT / "benchmark" / "figures" / "study1_simulated_allPackages"

# Stable colour per tool.
COLORS = {
    "epykit_lr": "#1f77b4", "epykit_lrplus": "#17becf",
    "methylkit": "#d62728", "dss": "#2ca02c", "dss_nosmooth": "#98df8a",
    "dmrseq": "#9467bd", "bsmooth": "#8c564b",
}
LABEL = {
    "epykit_lr": "epykit lr", "epykit_lrplus": "epykit lr+",
    "methylkit": "methylKit", "dss": "DSS (smooth)",
    "dss_nosmooth": "DSS (no smooth)", "dmrseq": "dmrseq (DMR→CpG)",
    "bsmooth": "BSmooth (DMR→CpG)",
}


def _phi(rho: float, coverage: int) -> float:
    return 1.0 + (coverage - 1) * rho


def _panel(ax, iqr: pl.DataFrame, metric: str, coverage: int, ylabel: str):
    rhos_all = sorted(iqr["rho"].unique().to_list())
    for tool in [t for t in COLORS if t in iqr["tool"].unique().to_list()]:
        sub = iqr.filter(pl.col("tool") == tool).sort("rho")
        if sub.height == 0 or sub[f"{metric}_median"].null_count() == sub.height:
            continue
        x = sub["rho"].to_numpy()
        y = sub[f"{metric}_median"].to_numpy()
        q1 = sub[f"{metric}_q1"].to_numpy()
        q3 = sub[f"{metric}_q3"].to_numpy()
        is_dmr = tool in ("dmrseq", "bsmooth")
        ax.plot(x, y, marker="o", ms=4, lw=1.8, color=COLORS[tool],
                ls="--" if is_dmr else "-", label=LABEL[tool], alpha=0.95)
        ax.fill_between(x, q1, q3, color=COLORS[tool], alpha=0.12)
    # realistic real-WGBS dispersion band phi ~ 1.5-5 -> rho range
    lo = (1.5 - 1) / (coverage - 1)
    hi = (5.0 - 1) / (coverage - 1)
    ax.axvspan(lo, hi, color="grey", alpha=0.08, zorder=0)
    ax.set_xlabel(r"simulator ICC  $\rho$")
    ax.set_ylabel(ylabel)
    # secondary top axis: Pearson phi
    sec = ax.secondary_xaxis("top", functions=(lambda r: _phi(r, coverage),
                                               lambda p: (p - 1) / (coverage - 1)))
    sec.set_xlabel(r"implied Pearson $\varphi = 1+(n-1)\rho$  at $n$=%d" % coverage)
    ax.grid(True, alpha=0.25)
    return rhos_all


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--coverage", type=int, default=10)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--iqr", type=Path, default=SIM_ROOT / "eval_phi_sweep_iqr.parquet")
    args = ap.parse_args(argv)

    if not args.iqr.exists():
        print(f"missing {args.iqr}; run run_phi_sweep.py first", file=sys.stderr)
        return 1
    iqr = pl.read_parquet(args.iqr)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # 3-panel: TPR, FPR, peak RSS
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    _panel(axes[0], iqr, "tpr", args.coverage, "TPR (recall) @ q<0.05")
    _panel(axes[1], iqr, "fpr", args.coverage, "FPR @ q<0.05")
    _panel(axes[2], iqr, "rss_peak_mb", args.coverage, "peak RSS (MB)")
    axes[0].legend(fontsize=7, loc="best", framealpha=0.9)
    axes[0].set_title("Sensitivity vs dispersion")
    axes[1].set_title("False-positive rate vs dispersion")
    axes[2].set_title("Peak memory vs dispersion")
    fig.suptitle("Dispersion (φ) sweep on the intrinsic-truth simulator "
                 "(coverage=%d, grey band = realistic WGBS φ≈1.5–5)" % args.coverage,
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ext in ("png", "svg"):
        p = args.out_dir / f"F9_phi_sweep.{ext}"
        fig.savefig(p, dpi=150, bbox_inches="tight")
        print(f"wrote {p}")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
