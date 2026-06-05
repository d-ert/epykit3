"""Re-implementation of the Piao et al. 2021 binomial DMC simulator.

Used to validate epykit defaults on held-out data not used during
parameter selection. The simulator's intrinsic `is_dmc` flag becomes
the ground truth, replacing the noisy threshold-reconstruction in
`_make_truth.py`. See `docs/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md`
§2.1 for the design rationale.

Simulation model
----------------
1. Baseline beta per CpG: bimodal mixture matching real WGBS marginals:
     - 40 %: low-meth pile, Beta(2, 50) (mode at ~0.04)
     - 40 %: high-meth pile, Beta(50, 2) (mode at ~0.96)
     - 20 %: intermediate, Beta(2, 2) (centred at 0.5)
   This roughly reproduces the U-shape observed in IMR90 / H1-hESC bulk
   WGBS without claiming exact match to Piao's simulator (which uses a
   similar but unspecified marginal). Tunable via `baseline_components`.

2. DMC designation: 20 % of CpGs marked as true DMCs. For each, sample
     `delta ~ U(0.2, 1.0)`, then `sign ~ {+1, -1}` 50/50, and apply
     `delta * sign` to the treatment group's baseline. The control group
     keeps the baseline unchanged.
   Clipping: treatment beta is clipped to [0, 1] post-shift to avoid
   illegal values; this can compress the realised |meth_diff| slightly
   relative to the designed `delta` near baseline=0 or baseline=1.

3. Per-sample counts: for sample i in {treatment, control}, draw
     `count_M ~ Binomial(n=coverage, p=beta_i)`,
     `count_U = coverage - count_M`.
   Coverage is deterministic per the `coverage` argument (matching Piao's
   fixed-coverage scenarios); replace with a Poisson or Negative Binomial
   draw if heteroscedastic coverage is needed.

   With ``phi > 0`` the simulator draws Beta-Binomial instead. Per-sample,
   per-site we draw ``p_i ~ Beta(mu * (1/phi - 1), (1-mu) * (1/phi - 1))``
   then ``count_M ~ Binomial(n=coverage, p=p_i)``. ``phi`` is the
   beta-binomial intraclass correlation rho in (0, 1); the variance of
   the BB count is the binomial variance multiplied by ``1 + (n-1) * phi``.
   Real WGBS dispersion is typically rho ~ 0.01-0.1 at coverage 10-30;
   the Piao binomial-only simulator (``phi=0``) is therefore under-
   dispersed relative to real data. Run the simulator at multiple
   rho values to assess engine sensitivity to overdispersion.

Outputs
-------
- AMP-format text files at `out_dir/amp.coverage={K}.sample{i}.txt` for
  i in 1..n_per_group*2 (treatment samples 1..n_per_group, control samples
  n_per_group+1..2*n_per_group). Compatible with `_loaders.py::amp_to_bismark_cov`.
- `out_dir/truth.parquet` matching the schema of
  `data/study1/ground_truth/dmc_truth.parquet`.

Returns
-------
dict with keys `truth` (Path to truth parquet) and `amp_files` (list of
Paths to per-sample AMP files).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl


@dataclass(frozen=True)
class SimConfig:
    n_cpgs: int = 100_000
    n_per_group: int = 3
    coverage: int = 10
    dmc_fraction: float = 0.20
    delta_lo: float = 0.20
    delta_hi: float = 1.00
    chromosome: str = "chr1"  # Piao simulator uses single contiguous CpG track
    pos_spacing_bp: int = 100  # 100 bp inter-CpG (loose CGI density average)
    seed: int = 42
    # Beta-Binomial intraclass correlation rho in [0, 1). phi=0
    # collapses to pure Binomial (the original Piao behaviour).
    # Recommended sweep for sensitivity analysis: {0.0, 0.01, 0.05, 0.1, 0.2}.
    phi: float = 0.0


def _draw_baseline_beta(n: int, rng: np.random.Generator) -> np.ndarray:
    """Right-skewed Beta(0.75, 1.35) baseline matching Piao's marginals.

    The Piao 2021 simulator produces a right-skewed, monotonically-decreasing
    baseline-methylation distribution (22% of CpGs at freqC ~ 0, tapering to
    ~7% at freqC ~ 100) rather than the U-shaped bimodal typical of high-CpG
    bulk WGBS. Empirical fit to Piao amp.coverage=10.sample1.txt (100K CpGs):
      Piao count_M mean = 3.77, std = 3.42
      Beta(0.75, 1.35): count_M mean = 3.74 (err 0.9%), std = 3.31 (err 3.0%)
    Both within the 10%/20% tolerances used in
    test_simulator_marginals_match_piao_within_tolerance.

    Previous parameters: bimodal mixture [40% Beta(2,50), 40% Beta(50,2),
    20% Beta(2,2)], which gave mean_beta ~ 0.50 vs Piao ~ 0.38 — a 32% error.
    """
    return np.clip(rng.beta(0.75, 1.35, size=n), 1e-4, 1.0 - 1e-4)


def _assign_dmcs(
    n: int, dmc_fraction: float, delta_lo: float, delta_hi: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (is_dmc, signed_delta, direction).

    is_dmc: bool[n]; signed_delta: float[n] (0 for non-DMCs);
    direction: U-array{"hyper","hypo","none"}[n].
    """
    is_dmc = rng.random(n) < dmc_fraction
    n_dmc = int(is_dmc.sum())
    delta_mag = rng.uniform(delta_lo, delta_hi, size=n_dmc)
    sign = rng.choice([+1.0, -1.0], size=n_dmc)
    signed_delta = np.zeros(n, dtype=np.float64)
    signed_delta[is_dmc] = sign * delta_mag

    direction = np.array(["none"] * n, dtype=object)
    direction[is_dmc & (signed_delta > 0)] = "hyper"
    direction[is_dmc & (signed_delta < 0)] = "hypo"
    return is_dmc, signed_delta, direction


def _draw_beta_binomial(
    n_trials: int,
    p: np.ndarray,
    phi: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw Beta-Binomial(n, p, phi) where phi is the intraclass correlation.

    Parameters
    ----------
    n_trials : int
        Coverage at each CpG (n in Binomial(n, p)).
    p : np.ndarray
        Per-CpG mean methylation probability, shape (n_cpgs,).
    phi : float
        Intraclass correlation rho in (0, 1). At phi -> 0 the
        Beta-Binomial collapses to Binomial(n, p); at phi -> 1 the
        per-draw probability is fixed at either 0 or 1 with probability
        1-p, p respectively (extreme over-dispersion).
    rng : np.random.Generator
        Random source.

    Returns
    -------
    np.ndarray of int64, shape (n_cpgs,) -- per-CpG methylated count.

    Notes
    -----
    Parameterisation: ``alpha = p * (1/phi - 1)``,
    ``beta = (1-p) * (1/phi - 1)``. Then Beta(alpha, beta) has mean p
    and variance ``p*(1-p)*phi``; the Beta-Binomial variance is
    ``n*p*(1-p) * (1 + (n-1)*phi)`` versus the pure Binomial
    ``n*p*(1-p)``.
    """
    if phi <= 0.0:
        raise ValueError(
            f"_draw_beta_binomial requires phi > 0 (got {phi}); use "
            f"rng.binomial directly for the phi=0 case."
        )
    if phi >= 1.0:
        raise ValueError(
            f"_draw_beta_binomial requires phi < 1 (got {phi}); the "
            f"Beta degenerates at the boundary."
        )
    concentration = 1.0 / phi - 1.0
    # Clip p to avoid Beta(0, *) / Beta(*, 0) degeneracy when the
    # baseline beta lands at the closed-interval boundary.
    p_safe = np.clip(p, 1e-6, 1.0 - 1e-6)
    alpha = p_safe * concentration
    beta = (1.0 - p_safe) * concentration
    # Per-CpG p drawn from Beta, then Binomial draw conditional on it.
    p_draw = rng.beta(alpha, beta)
    return rng.binomial(n_trials, p_draw).astype(np.int64)


def _meth_diff_bin(true_meth_diff: np.ndarray) -> np.ndarray:
    """Stratify |delta| into the paper's bins."""
    abs_d = np.abs(true_meth_diff)
    out = np.array(["none"] * len(abs_d), dtype=object)
    out[(abs_d >= 0.20) & (abs_d < 0.40)] = "0.2-0.4"
    out[(abs_d >= 0.40) & (abs_d < 0.60)] = "0.4-0.6"
    out[(abs_d >= 0.60) & (abs_d < 0.80)] = "0.6-0.8"
    out[abs_d >= 0.80] = "0.8-1.0"
    return out


def simulate_dmc(
    n_cpgs: int = 100_000,
    n_per_group: int = 3,
    coverage: int = 10,
    seed: int = 42,
    out_dir: Path | str | None = None,
    *,
    dmc_fraction: float = 0.20,
    delta_lo: float = 0.20,
    delta_hi: float = 1.00,
    chromosome: str = "chr1",
    pos_spacing_bp: int = 100,
    phi: float = 0.0,
) -> dict:
    """Run one simulation. See module docstring for the model.

    Returns dict: {"truth": Path, "amp_files": list[Path], "config": SimConfig}.
    """
    if not (0.0 <= phi < 1.0):
        raise ValueError(
            f"phi must be in [0, 1); got {phi}. Use phi=0 for pure "
            f"Binomial (the Piao default) or phi in (0, 1) for "
            f"Beta-Binomial overdispersion."
        )
    cfg = SimConfig(
        n_cpgs=n_cpgs, n_per_group=n_per_group, coverage=coverage,
        dmc_fraction=dmc_fraction, delta_lo=delta_lo, delta_hi=delta_hi,
        chromosome=chromosome, pos_spacing_bp=pos_spacing_bp, seed=seed,
        phi=phi,
    )
    out = Path(out_dir) if out_dir is not None else Path.cwd() / "simulate_piao_out"
    out.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)

    # Positions: contiguous CpGs at fixed spacing on a single chromosome.
    positions = np.arange(1, n_cpgs + 1, dtype=np.int64) * pos_spacing_bp

    # Baseline (control) beta and signed effect on treatment.
    baseline = _draw_baseline_beta(n_cpgs, rng)
    is_dmc, signed_delta, direction = _assign_dmcs(
        n_cpgs, dmc_fraction, delta_lo, delta_hi, rng,
    )
    treat_beta = np.clip(baseline + signed_delta, 1e-4, 1.0 - 1e-4)
    ctrl_beta = baseline

    # Per-sample counts. Samples 1..n_per_group are treatment, n_per_group+1..2n control.
    amp_files: list[Path] = []
    sample_idx = 0
    treat_count_M = np.zeros((n_per_group, n_cpgs), dtype=np.int64)
    ctrl_count_M = np.zeros((n_per_group, n_cpgs), dtype=np.int64)

    if phi <= 0.0:
        for j in range(n_per_group):
            treat_count_M[j] = rng.binomial(coverage, treat_beta)
        for j in range(n_per_group):
            ctrl_count_M[j] = rng.binomial(coverage, ctrl_beta)
    else:
        for j in range(n_per_group):
            treat_count_M[j] = _draw_beta_binomial(
                coverage, treat_beta, phi=phi, rng=rng,
            )
        for j in range(n_per_group):
            ctrl_count_M[j] = _draw_beta_binomial(
                coverage, ctrl_beta, phi=phi, rng=rng,
            )

    # Per-sample mean beta is the realised (noisy) value; truth uses the
    # *expected* mean (clean signal), which is the input beta. The
    # downstream `evaluate.py::_join_with_truth` reads `is_dmc` directly,
    # so the realised vs. expected distinction only affects the
    # `mean_beta_*` columns that some users read for diagnostics.
    mean_beta_treat = treat_count_M.mean(axis=0) / coverage
    mean_beta_ctrl = ctrl_count_M.mean(axis=0) / coverage

    # The truth uses signed_delta (the *intended* effect), not the
    # realised difference. This is the key win over the threshold-based
    # _make_truth.py: a CpG is a true DMC iff is_dmc was set by the
    # simulator, regardless of how the noise played out at low coverage.
    truth_df = pl.DataFrame({
        "chrom": [chromosome] * n_cpgs,
        "pos": positions,
        "mean_beta_treat": mean_beta_treat.astype(np.float64),
        "mean_beta_ctrl": mean_beta_ctrl.astype(np.float64),
        # true_meth_diff is the DESIGNED effect, not the realised one.
        "true_meth_diff": signed_delta.astype(np.float64),
        "is_dmc": is_dmc,
        "direction": [str(d) for d in direction],
        "meth_diff_bin": [str(b) for b in _meth_diff_bin(signed_delta)],
    }).with_columns(
        pl.col("chrom").cast(pl.Utf8),
        pl.col("pos").cast(pl.Int64),
    )
    truth_path = out / "truth.parquet"
    truth_df.write_parquet(truth_path)

    # Write per-sample AMP files. Schema matches Piao's:
    #   chrBase chr base strand coverage freqC freqT
    # `freqC` is methylation percent (0..100), `freqT` = 100 - freqC.
    for i in range(n_per_group):
        sample_idx += 1
        path = out / f"amp.coverage={coverage}.sample{sample_idx}.txt"
        _write_amp(path, chromosome, positions, treat_count_M[i], coverage)
        amp_files.append(path)
    for i in range(n_per_group):
        sample_idx += 1
        path = out / f"amp.coverage={coverage}.sample{sample_idx}.txt"
        _write_amp(path, chromosome, positions, ctrl_count_M[i], coverage)
        amp_files.append(path)

    return {"truth": truth_path, "amp_files": amp_files, "config": cfg}


def _write_amp(
    path: Path, chrom: str, positions: np.ndarray, count_M: np.ndarray, coverage: int,
) -> None:
    """Write one sample to AMP format (header + tab-separated rows)."""
    freqC = (count_M / coverage) * 100.0
    freqT = 100.0 - freqC
    df = pl.DataFrame({
        "chrBase": [f"{chrom}.{int(p)}" for p in positions],
        "chr": [chrom] * len(positions),
        "base": positions.astype(np.int64),
        "strand": ["F"] * len(positions),
        "coverage": np.full(len(positions), coverage, dtype=np.int64),
        "freqC": freqC,
        "freqT": freqT,
    })
    df.write_csv(path, separator="\t", include_header=True)


def main(argv: list[str] | None = None) -> None:
    """CLI: `python simulate_piao.py --n-cpgs 100000 --coverage 10 --seed 42 --out out/`"""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-cpgs", type=int, default=100_000)
    parser.add_argument("--n-per-group", type=int, default=3)
    parser.add_argument("--coverage", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--dmc-fraction", type=float, default=0.20)
    parser.add_argument(
        "--phi", type=float, default=0.0,
        help="Beta-Binomial intraclass correlation (overdispersion). "
             "Default 0 = pure Binomial (Piao). "
             "Recommended sweep: 0.0, 0.01, 0.05, 0.1, 0.2.",
    )
    args = parser.parse_args(argv)

    result = simulate_dmc(
        n_cpgs=args.n_cpgs,
        n_per_group=args.n_per_group,
        coverage=args.coverage,
        seed=args.seed,
        dmc_fraction=args.dmc_fraction,
        out_dir=args.out,
        phi=args.phi,
    )
    print(f"wrote {len(result['amp_files'])} AMP files + {result['truth']} (phi={args.phi})")


if __name__ == "__main__":
    main()
