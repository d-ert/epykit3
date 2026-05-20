"""Synthetic Bismark .cov dataset with known DMC/DMR truth.

Used by the test suite to measure power, FDR, and effect-size bias of every
DMC and DMR backend in epykit. The data-generating model is intentionally
simple so the *truth* is unambiguous:

    pi_ij  = baseline + effect[i] * is_treatment(j) + noise[i, j]
    cov_ij ~ NegativeBinomial(mean=coverage_mean, dispersion=coverage_disp)
    meth_ij ~ Binomial(cov_ij, clip(pi_ij, 0.01, 0.99))

Effect placement:

* ``n_dmrs`` contiguous regions of ``dmr_size_cpgs`` CpGs each receive the
  same signed ``dmr_effect`` (a real biological DMR).
* ``n_scattered_dmcs`` isolated CpGs outside any DMR receive +/-``dmc_effect``
  with random signs (scattered DMCs, not part of a DMR).

The remainder are null (effect == 0).

Outputs in ``out_dir/``:

  cov/<sample_id>.bismark.cov.gz   -- per-sample Bismark .cov files
  samplesheet.csv                  -- sample_id, group, path
  truth.parquet                    -- chrom, pos, is_dmc, true_meth_diff,
                                     dmr_id, in_dmr
  config.json                      -- dataclass dump for reproducibility
"""

from __future__ import annotations

import gzip
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl


@dataclass
class SimConfig:
    """Knobs for the synthetic methylation generator.

    The defaults are chosen so that a "standard" 4-vs-4 WGBS comparison
    at ~20x coverage produces detectable signal under BH-corrected
    multiple testing across ~75 k post-filter CpGs. Specifically:

    * dmc_effect = 0.40: bigger than the typical 0.20-0.30 promoter Deltabeta,
      so n=4 replicates can resolve it with a per-site count test (LR /
      score) at BH q<0.05. With Deltabeta=0.30 and n=4 the BH cutoff after
      75k tests demands an unrealistic z-statistic.
    * coverage_mean = 20 x NB(disp=5): puts most sites at >=10x after
      filter, well within the "well-powered" regime for binomial GLMs.
    * replicate_sd = 0.03: modest between-replicate variation that still
      makes the count-vs-Welch-t comparison interesting.

    Tightening any of these gives the engines an easier time and risks
    masking real bugs; loosening them makes the power thresholds in
    ``test_accuracy.py`` too optimistic to clear at honest noise levels.
    """

    n_per_group: int = 4
    chromosomes: tuple[str, ...] = ("chr1", "chr2", "chr3", "chr4", "chr5")
    cpgs_per_chrom: int = 2_000
    baseline_meth: float = 0.30
    n_scattered_dmcs: int = 500
    dmc_effect: float = 0.40
    n_dmrs: int = 10
    dmr_size_cpgs: int = 10
    dmr_effect: float = 0.40
    coverage_mean: float = 20.0
    coverage_disp: float = 5.0  # NB shape (k); larger -> less overdispersion
    replicate_sd: float = 0.03
    seed: int = 42

    # --- Multi-group / continuous-covariate extensions --------------------
    # When n_groups >= 3, samples are drawn from group_labels[:n_groups] in
    # equal-size blocks. The first half of scattered DMCs become "multi-
    # group" DMCs whose effect ramps with the group index (multiplied by
    # `dmc_effect_multigroup_step`). The other half remain binary "treat-
    # ment vs others" DMCs.
    n_groups: int = 2
    group_labels: tuple[str, ...] = ("control", "treatment", "extra1", "extra2")
    dmc_effect_multigroup_step: float = 0.20
    # When set, adds an `age` column to the samplesheet drawn U(age_low,
    # age_high) and injects a linear age x beta effect at `n_age_dmcs`
    # sites with slope `age_effect_per_year` (per-year Deltabeta).
    continuous_covariate: bool = False
    age_low: float = 20.0
    age_high: float = 80.0
    n_age_dmcs: int = 200
    age_effect_per_year: float = 0.005  # 0.5pp Deltabeta per year of age, 60yr span = 0.30

    @property
    def total_samples(self) -> int:
        return self.n_per_group * self.n_groups

    @property
    def n_total_sites(self) -> int:
        return len(self.chromosomes) * self.cpgs_per_chrom


def _nb_coverage(rng: np.random.Generator, mean: float, k: float, n: int) -> np.ndarray:
    """Negative-binomial coverage with mean ``mean`` and dispersion ``k``.

    numpy's parametrisation: NB(k, p) where mean = k * (1 - p) / p. We solve
    for p given the requested mean and k. Returns int64, minimum 1 to avoid
    accidental zero-coverage sites at low means.
    """
    p = k / (k + mean)
    out = rng.negative_binomial(k, p, size=n).astype(np.int64)
    return np.maximum(out, 1)


def _place_effects(cfg: SimConfig, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Decide per-site effect size + DMR membership.

    Returns (effects, dmr_id, chrom_of_site) where:
      effects[i]    = true Deltabeta at site i (0 for null sites)
      dmr_id[i]     = DMR index (0..n_dmrs-1) or -1 if not in a DMR
      chrom_of_site = chromosome name per site (for joining)
    """
    n = cfg.n_total_sites
    effects = np.zeros(n, dtype=np.float64)
    dmr_id = np.full(n, -1, dtype=np.int32)

    chroms_arr = np.repeat(np.array(cfg.chromosomes, dtype=object), cfg.cpgs_per_chrom)
    n_chroms = len(cfg.chromosomes)

    # 1. Place DMRs first: each DMR sits entirely on one chromosome.
    # Sign assignment is *deterministically* balanced (half hyper, half hypo)
    # so the truth table has zero mean effect over DMRs. Without this the
    # random +/-1 draw produces a small fixture-level imbalance (e.g. 7/3 at
    # n_dmrs=10) that shows up as "bias" in test assertions.
    dmr_signs = np.concatenate([
        np.full(cfg.n_dmrs // 2, 1.0),
        np.full(cfg.n_dmrs - cfg.n_dmrs // 2, -1.0),
    ])
    rng.shuffle(dmr_signs)
    for d in range(cfg.n_dmrs):
        chrom_idx = rng.integers(0, n_chroms)
        within = rng.integers(0, cfg.cpgs_per_chrom - cfg.dmr_size_cpgs)
        start = chrom_idx * cfg.cpgs_per_chrom + int(within)
        end = start + cfg.dmr_size_cpgs
        effects[start:end] = float(dmr_signs[d]) * cfg.dmr_effect
        dmr_id[start:end] = d

    # 2. Scatter remaining DMCs into not-in-DMR positions. Signs balanced
    # deterministically (same reasoning as DMRs).
    pool = np.where(dmr_id == -1)[0]
    n_chosen = min(cfg.n_scattered_dmcs, len(pool))
    chosen = rng.choice(pool, size=n_chosen, replace=False)
    scatter_signs = np.concatenate([
        np.full(n_chosen // 2, 1.0),
        np.full(n_chosen - n_chosen // 2, -1.0),
    ])
    rng.shuffle(scatter_signs)
    effects[chosen] = scatter_signs * cfg.dmc_effect

    return effects, dmr_id, chroms_arr


def _positions(cfg: SimConfig, rng: np.random.Generator) -> np.ndarray:
    """Generate sorted positions per chromosome, spaced 50-400 bp apart."""
    parts = []
    base = 1_000
    for _ in cfg.chromosomes:
        gaps = rng.integers(50, 400, cfg.cpgs_per_chrom)
        pos = base + np.cumsum(gaps).astype(np.int64)
        parts.append(pos)
    return np.concatenate(parts)


def _write_cov_gz(path: Path, chroms: np.ndarray, positions: np.ndarray,
                  N_meth: np.ndarray, N_unmeth: np.ndarray) -> None:
    """Write a Bismark .cov.gz file.

    Format (tab-separated, no header):
        chrom  start  end  methyl_percent  N_meth  N_unmeth

    epykit's converter treats ``start`` as 0-based (BED-format), matching
    nf-core/methylseq's bismark2bedGraph output. We write start = pos,
    end = pos + 1 (single-CpG interval).
    """
    coverage = (N_meth + N_unmeth).astype(np.int64)
    methyl_pct = 100.0 * N_meth / np.maximum(coverage, 1)
    df = pd.DataFrame({
        "chrom": chroms,
        "start": positions.astype(np.int64),
        "end": (positions + 1).astype(np.int64),
        "methyl_pct": methyl_pct.astype(np.float64),
        "N_meth": N_meth.astype(np.int64),
        "N_unmeth": N_unmeth.astype(np.int64),
    })
    with gzip.open(path, "wt", newline="") as fh:
        df.to_csv(fh, sep="\t", header=False, index=False, float_format="%.4f",
                  lineterminator="\n")


def generate(cfg: SimConfig, out_dir: str | Path) -> dict:
    """Generate Bismark .cov.gz files + samplesheet + truth table.

    Returns a dict with paths and summary stats so callers can hand
    ``samplesheet`` straight into ``ep.read_bismark`` and join ``truth``
    onto DMC results.
    """
    out_dir = Path(out_dir)
    cov_dir = out_dir / "cov"
    cov_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(cfg.seed)

    # Per-site truth (same across all samples).
    effects, dmr_id, chroms_arr = _place_effects(cfg, rng)
    positions = _positions(cfg, rng)
    n = cfg.n_total_sites

    # ---- Optional multi-group / continuous-covariate truths --------------
    multigroup_sites = np.zeros(n, dtype=bool)
    age_sites = np.zeros(n, dtype=bool)
    age_slope = np.zeros(n, dtype=np.float64)
    if cfg.n_groups >= 3:
        # Carve a fresh band of sites (outside DMRs and scattered DMCs) for
        # multi-group effects: the per-site effect grows linearly with the
        # group index, so a joint F-test on group should detect them and a
        # binary test (group_0 vs everyone else) should miss most of them.
        pool = np.where((dmr_id == -1) & (effects == 0))[0]
        n_mg = min(200, len(pool))
        if n_mg > 0:
            chosen = rng.choice(pool, size=n_mg, replace=False)
            multigroup_sites[chosen] = True

    if cfg.continuous_covariate:
        # Place age-effect sites on the remaining null positions.
        pool = np.where(
            (dmr_id == -1) & (effects == 0) & (~multigroup_sites)
        )[0]
        n_age = min(cfg.n_age_dmcs, len(pool))
        if n_age > 0:
            chosen = rng.choice(pool, size=n_age, replace=False)
            age_sites[chosen] = True
            # Half positive, half negative slope.
            slopes = np.concatenate([
                np.full(n_age // 2, cfg.age_effect_per_year),
                np.full(n_age - n_age // 2, -cfg.age_effect_per_year),
            ])
            rng.shuffle(slopes)
            age_slope[chosen] = slopes

    # Resolve per-sample group assignment + optional continuous covariate.
    group_labels = list(cfg.group_labels[:cfg.n_groups])
    # When n_groups=2 the legacy layout is "treatment first, control second";
    # we keep that ordering by reversing the group_labels list. Beyond that
    # we just take the first n_groups labels.
    if cfg.n_groups == 2:
        ordered_labels = ["treatment", "control"]
    else:
        ordered_labels = group_labels[:cfg.n_groups]

    # Generate per-sample read counts.
    sample_records: list[dict] = []
    for sample_idx in range(cfg.total_samples):
        # Determine which group this sample belongs to (block layout).
        group_idx = sample_idx // cfg.n_per_group
        group = ordered_labels[group_idx]
        sid = f"{group}_{(sample_idx % cfg.n_per_group) + 1}"
        age = (
            float(rng.uniform(cfg.age_low, cfg.age_high))
            if cfg.continuous_covariate else None
        )

        # True per-site beta for this sample. Decompose into:
        #   - binary effect (existing): only "treatment" group gets it
        #   - multi-group effect: scales with group_idx (centered)
        #   - age effect: linear with age at age_sites
        if group == "treatment":
            per_site_effect = effects.copy()
        else:
            per_site_effect = np.zeros_like(effects)
        if cfg.n_groups >= 3 and multigroup_sites.any():
            # Center group index around (n_groups - 1) / 2 so the
            # baseline-shifted F-test sees the full range.
            centered = (group_idx - (cfg.n_groups - 1) / 2.0)
            per_site_effect = per_site_effect + (
                multigroup_sites.astype(np.float64)
                * centered
                * cfg.dmc_effect_multigroup_step
            )
        if cfg.continuous_covariate and age_sites.any():
            per_site_effect = per_site_effect + age_slope * (age or 0.0)

        # Independent replicate noise around the per-sample mean.
        rep_noise = rng.normal(0.0, cfg.replicate_sd, size=n)
        beta = np.clip(cfg.baseline_meth + per_site_effect + rep_noise, 0.01, 0.99)

        cov = _nb_coverage(rng, cfg.coverage_mean, cfg.coverage_disp, n)
        meth = rng.binomial(cov, beta).astype(np.int64)
        unmeth = cov - meth

        cov_path = cov_dir / f"{sid}.bismark.cov.gz"
        _write_cov_gz(cov_path, chroms_arr, positions, meth, unmeth)
        rec = {
            "sample_id": sid,
            "group": group,
            "path": str(cov_path),
        }
        if cfg.continuous_covariate:
            rec["age"] = age
        sample_records.append(rec)

    # Samplesheet.
    samplesheet_path = out_dir / "samplesheet.csv"
    pd.DataFrame(sample_records).to_csv(samplesheet_path, index=False)

    # Truth table.
    #
    # ``true_meth_diff`` stores the **post-clip effective Deltabeta**, not the
    # raw ``effects[i]``. With baseline=0.30 + effect=+/-0.40, the intended
    # hypo beta_treatment of -0.10 gets clipped to 0.01 by the sampling
    # model, so the actual realised Deltabeta is 0.01 - 0.30 = -0.29 (not -0.40).
    # If we stored the unclipped intended effect, the estimator would
    # appear ~5pp positively biased on hypo sites even though it's
    # correctly recovering what was sampled.
    pi_treat = np.clip(cfg.baseline_meth + effects, 0.01, 0.99)
    pi_ctrl  = np.full_like(effects, cfg.baseline_meth)
    true_meth_diff_postclip = pi_treat - pi_ctrl
    truth_df = pl.DataFrame({
        "chrom": [str(c) for c in chroms_arr.tolist()],
        "pos": positions.astype(np.int64).tolist(),
        "is_dmc": (effects != 0).tolist(),
        "true_meth_diff": true_meth_diff_postclip.tolist(),
        "dmr_id": dmr_id.tolist(),
        "in_dmr": (dmr_id >= 0).tolist(),
        "is_multigroup_dmc": multigroup_sites.tolist(),
        "is_age_dmc": age_sites.tolist(),
        "age_slope": age_slope.tolist(),
    })
    truth_path = out_dir / "truth.parquet"
    truth_df.write_parquet(truth_path)

    # Config dump for reproducibility.
    (out_dir / "config.json").write_text(json.dumps(asdict(cfg), indent=2, default=str))

    return {
        "samplesheet": str(samplesheet_path),
        "truth": truth_path,
        "out_dir": str(out_dir),
        "cov_dir": str(cov_dir),
        "n_total_sites": n,
        "n_dmcs_true": int((effects != 0).sum()),
        "n_dmrs": cfg.n_dmrs,
        "n_multigroup_dmcs": int(multigroup_sites.sum()),
        "n_age_dmcs": int(age_sites.sum()),
        "chromosomes": list(cfg.chromosomes),
        "sample_ids": [r["sample_id"] for r in sample_records],
        "treatment_ids": [r["sample_id"] for r in sample_records if r["group"] == "treatment"],
        "control_ids": [r["sample_id"] for r in sample_records if r["group"] == "control"],
        "group_ids": {
            label: [r["sample_id"] for r in sample_records if r["group"] == label]
            for label in ordered_labels
        },
        "config": cfg,
    }
