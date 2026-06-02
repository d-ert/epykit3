"""Validate the empirical-Bayes inverse-Gamma dispersion prior used by
``epykit.dmc`` (``dispersion="eb"``, ``dmc.py:775-820``).

The EB engine assumes per-site dispersion ``phi_i`` is drawn from an
inverse-Gamma prior whose parameters are estimated by method of
moments from the chromosome-wide distribution of per-site Pearson
dispersion estimates. The reviewer comment in the GB resubmission
plan flagged that this assumption is unvalidated -- there is no Q-Q
or goodness-of-fit check shipped with the package.

This script closes the gap. Given a parquet with one ``phi_site``
column (the per-site Pearson dispersion before shrinkage), it:

  1. Fits the inverse-Gamma prior by method of moments
     (mirroring ``dmc.py:814-818``).
  2. Builds a Q-Q plot of empirical phi_site quantiles against the
     theoretical inverse-Gamma quantiles.
  3. Runs a one-sample Kolmogorov-Smirnov test against the fitted
     inverse-Gamma.
  4. Optionally cross-compares against an external dispersion array
     (e.g. DSS's local-regression dispersion on the same cohort) via
     a two-sample KS test plus an overlay scatter on the Q-Q plot.

Outputs are written to ``--qq-out`` (PNG) and ``--summary-out``
(JSON with the MoM parameters, K-S statistics, n, and the
quality-of-fit interpretation).

The parquet input is produced by a small future addition to
``epykit.dmc.process_chromosomes_dmc`` (Linux Layer B). Until then
this script is callable end-to-end on any parquet with a numeric
``phi_site`` column -- the format is intentionally simple so a
researcher can also run it on an ad-hoc dispersion estimate
exported from another tool.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import polars as pl
from scipy import stats as sp_stats

logger = logging.getLogger(__name__)


# Replication of the MoM fit at dmc.py:814-818 for offline use.
def fit_inv_gamma_mom(phi_obs: np.ndarray) -> Optional[dict]:
    """Method-of-moments inverse-Gamma fit, mirroring ``dmc.py:_score_finalize``.

    Inverse-Gamma(a, b) has:
        mean = b / (a - 1)
        var  = b^2 / ((a - 1)^2 (a - 2))
    so given empirical (m, v) we recover a = m^2 / v + 2, b = m (a - 1).
    Returns ``None`` if the data is too degenerate to fit (zero
    variance, non-positive mean, or insufficient sites).
    """
    finite = phi_obs[np.isfinite(phi_obs) & (phi_obs > 0)]
    if finite.size < 100:
        return None
    m = float(np.mean(finite))
    v = float(np.var(finite))
    if v <= 1e-9 or m <= 0:
        return None
    a = m * m / v + 2.0
    b = m * (a - 1.0)
    if a <= 2.0:  # variance is undefined; fit is degenerate
        return None
    return {"a": a, "b": b, "mean": m, "variance": v, "n": int(finite.size)}


def ks_one_sample_inv_gamma(phi_obs: np.ndarray, a: float, b: float) -> dict:
    """One-sample K-S test of empirical phi_site vs InvGamma(a, b).

    scipy parameterises inverse-Gamma as ``invgamma.cdf(x, a, scale=b)``,
    matching the conventional shape-scale form. ``kstest`` accepts a
    callable CDF; we wrap the frozen distribution so we can pass
    ``scale`` without using kwargs (kstest does not forward kwargs).
    """
    finite = phi_obs[np.isfinite(phi_obs) & (phi_obs > 0)]
    if finite.size == 0:
        return {"statistic": float("nan"), "pvalue": float("nan"), "n": 0}
    frozen = sp_stats.invgamma(a, scale=b)
    res = sp_stats.kstest(finite, frozen.cdf)
    return {
        "statistic": float(res.statistic),
        "pvalue": float(res.pvalue),
        "n": int(finite.size),
    }


def ks_two_sample(a: np.ndarray, b: np.ndarray) -> dict:
    """Two-sample K-S; used to compare epykit's empirical phi against
    an external reference (e.g., DSS's local-regression dispersion)."""
    a_f = a[np.isfinite(a) & (a > 0)]
    b_f = b[np.isfinite(b) & (b > 0)]
    if a_f.size == 0 or b_f.size == 0:
        return {"statistic": float("nan"), "pvalue": float("nan"),
                "n_a": int(a_f.size), "n_b": int(b_f.size)}
    res = sp_stats.kstest(a_f, b_f)
    return {
        "statistic": float(res.statistic),
        "pvalue": float(res.pvalue),
        "n_a": int(a_f.size),
        "n_b": int(b_f.size),
    }


def _build_qq_plot(
    phi_obs: np.ndarray,
    fit: dict,
    out_path: Path,
    external: Optional[np.ndarray] = None,
    external_label: str = "external",
) -> None:
    """Save a Q-Q plot of empirical phi_site vs the fitted InvGamma(a, b).

    matplotlib is imported lazily so the script remains importable on
    runners without a display backend.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    finite = phi_obs[np.isfinite(phi_obs) & (phi_obs > 0)]
    n = finite.size
    a, b = fit["a"], fit["b"]
    # Theoretical quantiles at the empirical CDF positions.
    probs = (np.arange(1, n + 1) - 0.5) / n
    theoretical = sp_stats.invgamma.ppf(probs, a, scale=b)
    empirical = np.sort(finite)

    fig, ax = plt.subplots(figsize=(5.0, 5.0))
    ax.scatter(theoretical, empirical, s=4, alpha=0.4,
               label="epykit phi_site")
    if external is not None:
        ext_finite = external[np.isfinite(external) & (external > 0)]
        if ext_finite.size > 0:
            # Use the same theoretical quantiles for visual overlay; if
            # the external sample has different n we re-derive its
            # quantile positions.
            ext_probs = (np.arange(1, ext_finite.size + 1) - 0.5) / ext_finite.size
            ext_theoretical = sp_stats.invgamma.ppf(ext_probs, a, scale=b)
            ax.scatter(ext_theoretical, np.sort(ext_finite),
                       s=4, alpha=0.4, label=external_label)
    lim = max(theoretical.max(), empirical.max()) * 1.05
    ax.plot([0, lim], [0, lim], "k--", lw=0.5, label="y = x")
    ax.set_xlabel("Inverse-Gamma(a={:.2f}, b={:.2f}) theoretical".format(a, b))
    ax.set_ylabel("Empirical phi_site")
    ax.set_title("EB prior Q-Q (n = {:,})".format(n))
    ax.legend(loc="upper left", fontsize=8)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main(argv: Optional[list[str]] = None) -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phi-parquet", required=True, type=Path,
        help="Parquet with a 'phi_site' column (per-site Pearson dispersion)",
    )
    parser.add_argument(
        "--external-parquet", type=Path, default=None,
        help="Optional parquet with an external dispersion column for "
             "comparison (e.g. DSS local-regression). Column name "
             "configured via --external-column.",
    )
    parser.add_argument(
        "--external-column", type=str, default="phi_site",
        help="Column name in --external-parquet (default: phi_site)",
    )
    parser.add_argument(
        "--external-label", type=str, default="DSS",
        help="Label for the external sample in the Q-Q plot (default: DSS)",
    )
    parser.add_argument(
        "--qq-out", required=True, type=Path,
        help="Output PNG for the Q-Q plot",
    )
    parser.add_argument(
        "--summary-out", required=True, type=Path,
        help="Output JSON for the fit + KS + diagnostics",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    df = pl.read_parquet(str(args.phi_parquet))
    if "phi_site" not in df.columns:
        raise SystemExit(
            f"--phi-parquet must have a 'phi_site' column; "
            f"found {df.columns}"
        )
    phi = df["phi_site"].to_numpy().astype(np.float64)

    fit = fit_inv_gamma_mom(phi)
    if fit is None:
        raise SystemExit(
            "Inverse-Gamma fit is degenerate (insufficient sites or "
            "zero variance). Cannot validate the EB prior on this data."
        )

    ks_self = ks_one_sample_inv_gamma(phi, a=fit["a"], b=fit["b"])
    external_arr = None
    ks_two = None
    if args.external_parquet is not None:
        ext_df = pl.read_parquet(str(args.external_parquet))
        if args.external_column not in ext_df.columns:
            raise SystemExit(
                f"--external-parquet missing column "
                f"{args.external_column!r}; found {ext_df.columns}"
            )
        external_arr = ext_df[args.external_column].to_numpy().astype(np.float64)
        ks_two = ks_two_sample(phi, external_arr)

    _build_qq_plot(
        phi, fit, args.qq_out,
        external=external_arr,
        external_label=args.external_label,
    )

    summary = {
        "input": str(args.phi_parquet),
        "fit_method": "method_of_moments_inverse_gamma",
        "fit": fit,
        "ks_one_sample_vs_invgamma": ks_self,
        "external": (
            None if args.external_parquet is None
            else {
                "path": str(args.external_parquet),
                "column": args.external_column,
                "label": args.external_label,
                "ks_two_sample": ks_two,
            }
        ),
        # The KS one-sample test is sensitive at large n; a small p
        # does not by itself invalidate the prior, it just signals
        # a detectable systematic deviation. The Q-Q plot is the
        # primary diagnostic; KS quantifies the deviation magnitude.
        "interpretation": (
            "Use the Q-Q plot as the primary diagnostic. The KS "
            "one-sample p-value is reported for completeness but is "
            "expected to be small at WGBS scale (n > 1e6) -- look at "
            "the shape of the deviation in the Q-Q tail."
        ),
    }
    args.summary_out.parent.mkdir(parents=True, exist_ok=True)
    args.summary_out.write_text(json.dumps(summary, indent=2) + "\n")
    logger.info(
        "EB prior fit: a=%.4f b=%.4f n=%d  KS D=%.4f p=%.2g",
        fit["a"], fit["b"], fit["n"],
        ks_self["statistic"], ks_self["pvalue"],
    )

    # CLI tool, final stdout line is intentional (mirrors the
    # epykit.cli convention).
    print(
        f"EB prior fit: a={fit['a']:.4f} b={fit['b']:.4f} n={fit['n']} "
        f"(KS D={ks_self['statistic']:.4f} p={ks_self['pvalue']:.2g}); "
        f"plot={args.qq_out}, summary={args.summary_out}"
    )


if __name__ == "__main__":
    main()
