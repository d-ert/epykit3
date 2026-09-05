# Phase 2: Benchmark Scripts (simulator + null calibration + CIs) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the three benchmark-side gaps from §2 of the design spec — method-shaped truth, no null calibration, no CIs on headline numbers — by adding a Piao re-implementation simulator, a label-shuffle null-calibration runner, and a Wilson/bootstrap CI helper. Also consolidate the `benchmarkin_merges/FINAL_REPORT/` work into `epykit3/benchmark/` for single-repo reproducibility.

**Architecture:** All new scripts live under `epykit3/benchmark/scripts/` once consolidation completes. The simulator is a small standalone module (~300 lines) that writes Piao-compatible AMP-format files and an intrinsic-truth parquet matching the existing `dmc_truth.parquet` schema (`chrom, pos, mean_beta_treat, mean_beta_ctrl, true_meth_diff, is_dmc, direction, meth_diff_bin`). Null calibration reuses the existing per-engine runners by re-labelling treatment / control IDs before dispatch. CIs operate on the existing `eval_summary.parquet` schema without touching the engines.

**Tech Stack:** Python 3.12, polars, numpy, scipy.stats (binomial, hypergeom, Wilson via `scipy.stats.binomtest`), statsmodels (for `proportion_confint` if more convenient), pytest. R is not touched.

**Companion spec:** [docs/history/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md](../specs/2026-05-27-paper-defendable-benchmark-design.md) §2 (gaps), §5 (script catalogue).

---

## Scope decisions (anchored to spec)

- **In scope (this plan):** consolidation, `simulate_piao.py`, `run_null_calibration.py`, `wilson_bootstrap_ci.py`.
- **Out of scope (Phase 3+):** `methylkit_stouffer_combine.R`, `regen_all.py`, `bug_fix_audit.py`. These pair with the locked benchmark re-run, which depends on this Phase 2 producing the simulator and CIs first.
- **The simulator is for held-out validation, not to replace Piao-as-distributed.** §2.1 of the spec: tune defaults on Piao, freeze, then run once on simulator output. This plan produces the simulator and its self-tests; the held-out re-run is Phase 3.
- **Null calibration runs on Piao-as-distributed first**, simulator output second. §2.2 of the spec.

---

## File structure

| File | Why touched | Phase 2 task |
|---|---|---|
| `epykit3/_legacy_benchmark/` (new, was `epykit3/benchmark/`) | Archive of pre-consolidation benchmark content | Task 1 |
| `epykit3/benchmark/` (replaced) | New canonical location, mirrors `benchmarkin_merges/FINAL_REPORT/` | Task 1 |
| `epykit3/.gitignore` | Exclude regenerable artefacts under the new `benchmark/` path | Task 1 |
| `epykit3/benchmark/scripts/simulate_piao.py` (new) | The Python re-implementation of Piao 2021's binomial DMC simulator | Tasks 2, 3 |
| `epykit3/benchmark/scripts/tests/test_simulate_piao.py` (new) | Unit tests for the simulator (determinism, marginals, file schema) | Tasks 2, 3 |
| `epykit3/benchmark/scripts/wilson_bootstrap_ci.py` (new) | Wilson CIs on TPR/FPR; bootstrap CIs on AUROC/F1 | Tasks 4, 5 |
| `epykit3/benchmark/scripts/tests/test_wilson_bootstrap_ci.py` (new) | CI math correctness on known reference values | Tasks 4, 5 |
| `epykit3/benchmark/scripts/run_null_calibration.py` (new) | Label-shuffle empirical FDR runner per (engine, scenario) | Task 6 |
| `epykit3/benchmark/scripts/tests/test_run_null_calibration.py` (new) | Mocked-engine test of the shuffle loop + Wilson-CI integration | Task 6 |
| `epykit3/benchmark/scripts/conftest.py` (new) | pytest config for the benchmark test suite | Task 2 |
| `CHANGELOG.md` | New `### Added` section for Phase 2 scripts | Tasks per-commit |

**Note on test layout:** the existing epykit package tests live at `epykit3/tests/`. The new benchmark-script tests live at `epykit3/benchmark/scripts/tests/` because the scripts are deliberately decoupled from the package import path. The benchmark test suite gets its own pytest invocation: `uv run pytest benchmark/scripts/tests/`. The main `tests/` suite stays untouched.

---

## Pre-flight (once, before Task 1)

- [ ] **Step 1: Verify on `p0-fixes` branch, clean tree, tag intact**

Run:
```
git status --short
git log --oneline -3
git tag --list "v0.7.3-p0-complete"
```

Expected: working tree shows only untracked dirs (`.github/`, `CLAUDE.md`, `benchmark/`, `docs/history/superpowers/plans/`). HEAD is `8f447a9` or later. Tag exists.

- [ ] **Step 2: Confirm benchmark tests baseline**

Run: `uv run pytest -m "not slow" --strict-markers -x -q`
Expected: 229 passed, 5 skipped, 0 failed (matches Phase 1 wrap-up).

- [ ] **Step 3: Read the existing `dmc_truth.parquet` schema you'll be matching**

The simulator must write a truth parquet matching:
```
chrom: String
pos: Int64
mean_beta_treat: Float64
mean_beta_ctrl: Float64
true_meth_diff: Float64
is_dmc: Boolean
direction: String  # "hyper" | "hypo" | "none"
meth_diff_bin: String  # "none" | "0.2-0.4" | "0.4-0.6" | "0.6-0.8" | "0.8-1.0"
```

This is the schema `evaluate.py::_join_with_truth` reads. The simulator's truth must drop into evaluate.py without changes to evaluate.py.

---

## Task 1: Consolidation — move `benchmarkin_merges/FINAL_REPORT/` into `epykit3/benchmark/`

**Risk:** large directory move with mixed git-tracked and untracked content. Stage in three explicit chunks rather than a `mv` that destroys traceability.

**Files:**
- Rename: `epykit3/benchmark/` → `epykit3/_legacy_benchmark/`
- Create: `epykit3/benchmark/` (populated from `D:\Coding\Projeler\methyl_lib\benchmarkin_merges\FINAL_REPORT\`)
- Modify: `epykit3/.gitignore`
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Sanity-check the source before any move**

Run:
```
uv run python -c "
from pathlib import Path
src = Path(r'D:/Coding/Projeler/methyl_lib/benchmarkin_merges/FINAL_REPORT')
total = sum(f.stat().st_size for f in src.rglob('*') if f.is_file())
n_files = sum(1 for _ in src.rglob('*') if _.is_file())
print(f'{n_files:,} files, {total / 1e9:.2f} GB total')
for top in sorted(src.iterdir()):
    if top.is_dir():
        size = sum(f.stat().st_size for f in top.rglob('*') if f.is_file()) / 1e6
        print(f'  {top.name}/  ({size:,.1f} MB)')
"
```

Expected: prints file count + MB per top-level subdir. Confirm sizes are reasonable (FINAL_REPORT is mostly text/parquet, total well under 5 GB). If the total exceeds 5 GB or one subdir is >2 GB, STOP and report — large binaries should not enter git.

- [ ] **Step 2: Archive the legacy benchmark directory**

```
git mv benchmark _legacy_benchmark
git status
```

Expected: `_legacy_benchmark/` shows as renamed. The benchmark dir is currently untracked (per Phase 1 wrap-up `git status --short`), so `git mv` may fall back to a regular `mv` and report nothing tracked. In that case:

```
mv benchmark _legacy_benchmark
```

(plain shell move). Confirm `epykit3/benchmark/` no longer exists.

- [ ] **Step 3: Copy FINAL_REPORT into `epykit3/benchmark/`**

```
uv run python -c "
import shutil
from pathlib import Path
src = Path(r'D:/Coding/Projeler/methyl_lib/benchmarkin_merges/FINAL_REPORT')
dst = Path(r'D:/Coding/Projeler/methyl_lib/epykit3/benchmark')
shutil.copytree(src, dst)
print(f'copied tree to {dst}')
"
```

Expected: prints "copied tree". Verify with `ls benchmark/` — should show `README.md`, `EXECUTIVE_SUMMARY.md`, `PROTOCOL.md`, `paper/`, `report/`, `figures/`, `data/`, `scripts/`.

- [ ] **Step 4: Update `.gitignore` for the new benchmark path**

Open `epykit3/.gitignore`. The existing block at the top reads:

```
# benchmark regenerables (not bundled in PyPI wheel; rebuilt by benchmark/scripts/*)
benchmark/raw_sim_data/
benchmark/ground_truth/
benchmark/_converted/
benchmark/_runs/
benchmark/results/
```

That's still correct for the consolidated layout. **No change needed** unless `benchmark/data/study{1,2,3}/_runs/` or similar regenerable subdirs exist — confirm with `ls benchmark/data/study*/` and add patterns if needed. For example, if `benchmark/data/study3/chain_merge/` or similar holds heavy regenerables not bundled, append:

```
benchmark/data/study*/_runs/
benchmark/data/study*/_converted/
benchmark/data/study*/chain_merge/
benchmark/data/study*/dss/
```

Document the addition in a comment line: `# Phase 2: consolidated benchmark regenerables, added 2026-05-27`.

- [ ] **Step 5: Check overall size of the staged content**

Run:
```
uv run python -c "
from pathlib import Path
import subprocess
out = subprocess.check_output(['git', 'add', '--dry-run', 'benchmark', '_legacy_benchmark', '.gitignore'], text=True)
n = len(out.splitlines())
print(f'would add {n:,} files')
"
```

Expected: prints file count. If >10,000 or single files >50 MB, STOP and investigate — likely a regenerable directory slipped past .gitignore.

- [ ] **Step 6: Stage and commit consolidation**

```
git add benchmark _legacy_benchmark .gitignore
git status --short | head -20
```

Then verify with `git diff --cached --stat | tail -5` that the byte totals look like text+parquet (not GB-scale binaries).

If clean, commit:
```
git commit -m "$(cat <<'EOF'
chore(benchmark): consolidate benchmarkin_merges/FINAL_REPORT into epykit3/benchmark

Per Open Q1 of the design spec, the FINAL_REPORT directory at
benchmarkin_merges/ (not under version control) is now the canonical
benchmark home and moves into the epykit3 repo for single-repo
reproducibility (one Zenodo DOI, one tag).

- Archives the previous outdated content of epykit3/benchmark/ as
  epykit3/_legacy_benchmark/ (Piao-only, pre-2026-05-22 work).
- Copies benchmarkin_merges/FINAL_REPORT/ to epykit3/benchmark/.
- Extends .gitignore to cover the consolidated regenerable paths.

Future work in the benchmark suite happens under epykit3/benchmark/
exclusively. benchmarkin_merges/ on disk remains as the source-of-truth
copy until Phase 3's locked re-run, after which it can be deleted.

Refs docs/history/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 7: Append CHANGELOG entry**

Append under the existing `## Unreleased` heading, in a NEW `### Added` subsection (create the subsection if it doesn't exist yet):

```markdown
### Added

- **Benchmark consolidation**: the canonical benchmark suite is now
  `epykit3/benchmark/` (was previously split between `epykit3/benchmark/`
  and `benchmarkin_merges/FINAL_REPORT/`). One repo, one tag, reproducible
  via `uv run pytest benchmark/scripts/tests/` for the new script tests.
```

Commit this small CHANGELOG addition with the consolidation commit if you remember; otherwise as a follow-up:

```
git add CHANGELOG.md
git commit -m "docs(changelog): note benchmark consolidation under Unreleased

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: `simulate_piao.py` core — binomial simulator with intrinsic truth

**What it does:** Generates a Piao-style DMC simulation entirely in Python. Per-CpG baseline methylation is drawn from a bimodal distribution matching real WGBS (heavy at β=0 and β=1, light intermediate). ~20% of CpGs are designated true DMCs with `|meth_diff| ~ U(0.2, 1.0)`, sign 50/50. Per-sample counts are drawn binomially at the requested coverage. Outputs (a) AMP-format `amp.coverage=K.sampleN.txt` files compatible with `_loaders.py::amp_to_bismark_cov` and (b) a `truth.parquet` matching `dmc_truth.parquet`'s schema with **intrinsic** `is_dmc` flag.

**Files:**
- Create: `epykit3/benchmark/scripts/simulate_piao.py`
- Create: `epykit3/benchmark/scripts/tests/__init__.py` (empty)
- Create: `epykit3/benchmark/scripts/tests/conftest.py`
- Create: `epykit3/benchmark/scripts/tests/test_simulate_piao.py`

- [ ] **Step 1: Create `conftest.py` for the benchmark-script test suite**

`epykit3/benchmark/scripts/tests/conftest.py`:

```python
"""Pytest configuration for the benchmark-script test suite.

These tests are independent of the main epykit test suite — run them via
`uv run pytest benchmark/scripts/tests/` from the repo root. They exercise
the simulator, null-calibration runner, and CI helpers without touching
the epykit package internals.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make the scripts directory importable as a flat package so tests can do
# `from simulate_piao import simulate_dmc` rather than messing with PYTHONPATH.
SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
```

And the empty `__init__.py`:
```bash
touch benchmark/scripts/tests/__init__.py
```

- [ ] **Step 2: Write the failing test for the simulator's truth schema and determinism**

`epykit3/benchmark/scripts/tests/test_simulate_piao.py`:

```python
"""Tests for the Piao 2021 binomial simulator re-implementation."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest


def test_simulate_dmc_truth_schema_matches_dmc_truth_parquet(tmp_path):
    """The simulator's truth parquet must match the columns + dtypes of
    evaluate.py's expected `dmc_truth.parquet`."""
    from simulate_piao import simulate_dmc

    result = simulate_dmc(
        n_cpgs=2000,
        n_per_group=3,
        coverage=10,
        seed=42,
        out_dir=tmp_path,
    )
    truth = pl.read_parquet(result["truth"])
    assert set(truth.columns) == {
        "chrom", "pos", "mean_beta_treat", "mean_beta_ctrl",
        "true_meth_diff", "is_dmc", "direction", "meth_diff_bin",
    }, f"Schema mismatch: {truth.columns}"
    # Dtype check on the columns that downstream code casts.
    assert truth["chrom"].dtype == pl.Utf8
    assert truth["pos"].dtype == pl.Int64
    assert truth["true_meth_diff"].dtype == pl.Float64
    assert truth["is_dmc"].dtype == pl.Boolean


def test_simulate_dmc_is_deterministic_with_seed(tmp_path):
    """Two calls with the same seed must produce bit-identical truth + reads."""
    from simulate_piao import simulate_dmc

    a = simulate_dmc(n_cpgs=1000, n_per_group=3, coverage=10, seed=7, out_dir=tmp_path / "a")
    b = simulate_dmc(n_cpgs=1000, n_per_group=3, coverage=10, seed=7, out_dir=tmp_path / "b")
    truth_a = pl.read_parquet(a["truth"])
    truth_b = pl.read_parquet(b["truth"])
    # Same is_dmc vector and same true_meth_diff vector.
    np.testing.assert_array_equal(truth_a["is_dmc"].to_numpy(), truth_b["is_dmc"].to_numpy())
    np.testing.assert_array_equal(truth_a["true_meth_diff"].to_numpy(), truth_b["true_meth_diff"].to_numpy())
    # Same per-sample AMP files.
    for i in range(1, 7):
        f_a = tmp_path / "a" / f"amp.coverage=10.sample{i}.txt"
        f_b = tmp_path / "b" / f"amp.coverage=10.sample{i}.txt"
        assert f_a.read_bytes() == f_b.read_bytes(), f"sample{i} differs"


def test_simulate_dmc_marginals_match_design(tmp_path):
    """~20% true DMCs with |meth_diff| in [0.2, 1.0] and 50/50 direction split."""
    from simulate_piao import simulate_dmc

    result = simulate_dmc(
        n_cpgs=10000,
        n_per_group=3,
        coverage=10,
        seed=1,
        out_dir=tmp_path,
    )
    truth = pl.read_parquet(result["truth"])

    n_dmc = int(truth["is_dmc"].sum())
    frac = n_dmc / len(truth)
    assert 0.18 <= frac <= 0.22, f"DMC fraction {frac:.3f} outside design 0.20 ± 0.02"

    # Among true DMCs, |true_meth_diff| ~ U(0.2, 1.0): mean ≈ 0.6, min >= 0.2, max <= 1.0.
    dmc_only = truth.filter(pl.col("is_dmc"))
    abs_diff = dmc_only["true_meth_diff"].abs().to_numpy()
    assert abs_diff.min() >= 0.20, f"min |meth_diff| {abs_diff.min():.3f} below 0.2"
    assert abs_diff.max() <= 1.00, f"max |meth_diff| {abs_diff.max():.3f} above 1.0"
    assert 0.55 <= abs_diff.mean() <= 0.65, f"mean |meth_diff| {abs_diff.mean():.3f} outside U(0.2,1.0) expectation 0.6"

    # Direction split is ~50/50 among true DMCs.
    n_hyper = int((dmc_only["direction"] == "hyper").sum())
    n_hypo = int((dmc_only["direction"] == "hypo").sum())
    assert abs(n_hyper - n_hypo) / n_dmc < 0.05, (
        f"direction split unbalanced: {n_hyper} hyper vs {n_hypo} hypo (ratio "
        f"{n_hyper / n_dmc:.3f})"
    )
```

- [ ] **Step 3: Run the tests to confirm they fail (no module yet)**

Run: `uv run pytest benchmark/scripts/tests/test_simulate_piao.py -v`
Expected: all 3 tests FAIL with `ModuleNotFoundError: No module named 'simulate_piao'`.

- [ ] **Step 4: Write the simulator module**

`epykit3/benchmark/scripts/simulate_piao.py`:

```python
"""Re-implementation of the Piao et al. 2021 binomial DMC simulator.

Used to validate epykit defaults on held-out data not used during
parameter selection. The simulator's intrinsic `is_dmc` flag becomes
the ground truth, replacing the noisy threshold-reconstruction in
`_make_truth.py`. See `docs/history/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md`
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


def _draw_baseline_beta(n: int, rng: np.random.Generator) -> np.ndarray:
    """Bimodal mixture: 40 % low, 40 % high, 20 % intermediate."""
    component = rng.choice([0, 1, 2], size=n, p=[0.4, 0.4, 0.2])
    out = np.empty(n, dtype=np.float64)
    out[component == 0] = rng.beta(2.0, 50.0, size=int((component == 0).sum()))
    out[component == 1] = rng.beta(50.0, 2.0, size=int((component == 1).sum()))
    out[component == 2] = rng.beta(2.0, 2.0, size=int((component == 2).sum()))
    return np.clip(out, 1e-4, 1.0 - 1e-4)


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
) -> dict:
    """Run one simulation. See module docstring for the model.

    Returns dict: {"truth": Path, "amp_files": list[Path], "config": SimConfig}.
    """
    cfg = SimConfig(
        n_cpgs=n_cpgs, n_per_group=n_per_group, coverage=coverage,
        dmc_fraction=dmc_fraction, delta_lo=delta_lo, delta_hi=delta_hi,
        chromosome=chromosome, pos_spacing_bp=pos_spacing_bp, seed=seed,
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

    for j in range(n_per_group):
        treat_count_M[j] = rng.binomial(coverage, treat_beta)
    for j in range(n_per_group):
        ctrl_count_M[j] = rng.binomial(coverage, ctrl_beta)

    # Per-sample mean beta is the realised (noisy) value; truth uses the
    # *expected* mean (clean signal), which is the input beta. The
    # downstream `evaluate.py::_join_with_truth` reads `is_dmc` directly,
    # so the realised vs. expected distinction only affects the
    # `mean_beta_*` columns that some users read for diagnostics.
    mean_beta_treat = treat_count_M.mean(axis=0) / coverage
    mean_beta_ctrl = ctrl_count_M.mean(axis=0) / coverage
    realised_diff = mean_beta_treat - mean_beta_ctrl

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
    args = parser.parse_args(argv)

    result = simulate_dmc(
        n_cpgs=args.n_cpgs,
        n_per_group=args.n_per_group,
        coverage=args.coverage,
        seed=args.seed,
        dmc_fraction=args.dmc_fraction,
        out_dir=args.out,
    )
    print(f"wrote {len(result['amp_files'])} AMP files + {result['truth']}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Run the three tests to verify pass**

Run: `uv run pytest benchmark/scripts/tests/test_simulate_piao.py -v`
Expected: all 3 PASS.

- [ ] **Step 6: Smoke-run the simulator from the CLI**

Run:
```
uv run python benchmark/scripts/simulate_piao.py --n-cpgs 5000 --coverage 10 --seed 99 --out /tmp/simpiao_smoke
```

Expected: prints `wrote 6 AMP files + /tmp/simpiao_smoke/truth.parquet`. Inspect:
```
uv run python -c "
import polars as pl
t = pl.read_parquet(r'/tmp/simpiao_smoke/truth.parquet')
print(t.describe())
print('DMC fraction:', t['is_dmc'].mean())
"
```

Expected: DMC fraction in [0.18, 0.22]; `true_meth_diff` describe shows non-trivial spread.

- [ ] **Step 7: Commit**

```
git add benchmark/scripts/simulate_piao.py benchmark/scripts/tests/__init__.py benchmark/scripts/tests/conftest.py benchmark/scripts/tests/test_simulate_piao.py
git commit -m "$(cat <<'EOF'
feat(benchmark): simulate_piao.py -- Python re-implementation of Piao 2021 simulator

Closes the AUROC tautology in the existing benchmark (truth was a
methylation-difference threshold on the same simulator output the
methods read). This simulator has an *intrinsic* is_dmc flag — a CpG
is a true DMC iff designated as such by the simulator, regardless of
how the binomial noise played out at low coverage.

Outputs:
- AMP-format text files (Piao schema, drop-in to _loaders.py).
- truth.parquet matching evaluate.py's expected schema
  (chrom, pos, mean_beta_treat, mean_beta_ctrl, true_meth_diff,
  is_dmc, direction, meth_diff_bin).

Tests verify: deterministic with seed, truth-schema match, ~20% DMCs
with |delta| ~ U(0.2, 1.0), 50/50 direction split.

Refs docs/history/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md §2.1

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: `simulate_piao.py` validation — marginals match real Piao data

**Goal:** verify the re-implementation produces marginal statistics within tolerance of Piao's distributed simulator output. This is the "did I rebuild the right simulator?" check.

**Files:**
- Modify: `epykit3/benchmark/scripts/tests/test_simulate_piao.py` (append marginal-match test)

- [ ] **Step 1: Identify a Piao-as-distributed sample file for comparison**

The Piao raw data lives at `benchmark/data/study1/...` if it was bundled, or `benchmarkin_merges/epykit_vs_allPackages(simulated_approxData)/raw_sim_data/simulated_datasets/dmc_simulation/coverage/amp.coverage=10.sample1.txt` if not. Check:

```
ls benchmark/data/study1/ 2>/dev/null
ls "D:/Coding/Projeler/methyl_lib/benchmarkin_merges/epykit_vs_allPackages(simulated_approxData)/raw_sim_data/simulated_datasets/dmc_simulation/coverage/" 2>/dev/null | head -5
```

Use whichever exists. The validation test will skip gracefully if neither is available (so the test suite still runs without raw Piao data).

- [ ] **Step 2: Add the marginal-match test**

Append to `tests/test_simulate_piao.py`:

```python
def _piao_sample_path() -> Path | None:
    """Return the path to a Piao coverage=10 sample if it exists locally."""
    candidates = [
        Path("benchmark/data/study1/raw_sim_data/simulated_datasets/dmc_simulation/coverage/amp.coverage=10.sample1.txt"),
        Path("D:/Coding/Projeler/methyl_lib/benchmarkin_merges/epykit_vs_allPackages(simulated_approxData)/raw_sim_data/simulated_datasets/dmc_simulation/coverage/amp.coverage=10.sample1.txt"),
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def test_simulator_marginals_match_piao_within_tolerance(tmp_path):
    """The simulator's per-CpG count_M distribution at coverage=10 should
    have the same first two moments as Piao's distributed sample within
    Monte Carlo noise on ~100K CpGs. This is the 'did I rebuild the right
    simulator?' check.

    Tolerance is loose because we don't know Piao's exact baseline model:
    mean within ±10%, std within ±20%. Failure here means investigate the
    `_draw_baseline_beta` mixture parameters.
    """
    from simulate_piao import simulate_dmc

    piao = _piao_sample_path()
    if piao is None:
        pytest.skip("Piao raw data not available locally; skipping marginal match")

    # Read Piao's count_M distribution.
    piao_df = pl.read_csv(piao, separator="\t")
    piao_count_M = (piao_df["coverage"].cast(pl.Float64) * piao_df["freqC"] / 100.0).round().cast(pl.Int64)
    piao_mean = float(piao_count_M.mean())
    piao_std = float(piao_count_M.std())

    # Match Piao's CpG count (100K for DMC sim) and coverage (10).
    res = simulate_dmc(n_cpgs=len(piao_count_M), n_per_group=3, coverage=10,
                       seed=12345, out_dir=tmp_path)
    sim_amp = res["amp_files"][0]
    sim_df = pl.read_csv(sim_amp, separator="\t")
    sim_count_M = (sim_df["coverage"].cast(pl.Float64) * sim_df["freqC"] / 100.0).round().cast(pl.Int64)
    sim_mean = float(sim_count_M.mean())
    sim_std = float(sim_count_M.std())

    # Tolerances are loose: we don't claim Piao's exact baseline model.
    rel_mean_err = abs(sim_mean - piao_mean) / piao_mean
    rel_std_err = abs(sim_std - piao_std) / piao_std
    assert rel_mean_err < 0.10, (
        f"simulator count_M mean {sim_mean:.2f} vs Piao {piao_mean:.2f}: "
        f"rel error {rel_mean_err:.3f} > 10%"
    )
    assert rel_std_err < 0.20, (
        f"simulator count_M std {sim_std:.2f} vs Piao {piao_std:.2f}: "
        f"rel error {rel_std_err:.3f} > 20%"
    )


def test_simulator_truth_dmc_count_close_to_piao_design(tmp_path):
    """Piao's design has exactly 20,000 / 100,000 = 20% true DMCs.
    The simulator with default dmc_fraction=0.2 should land at 20% ± 0.5%
    on 100K CpGs (~50 std error)."""
    from simulate_piao import simulate_dmc

    res = simulate_dmc(n_cpgs=100_000, n_per_group=3, coverage=10,
                       seed=2026, out_dir=tmp_path)
    truth = pl.read_parquet(res["truth"])
    n_dmc = int(truth["is_dmc"].sum())
    assert 19_500 <= n_dmc <= 20_500, (
        f"n_dmc = {n_dmc:,}; design is 20,000 (20% of 100,000). "
        f"Outside ±0.5% tolerance suggests a bug in _assign_dmcs."
    )
```

- [ ] **Step 3: Run the new tests**

Run: `uv run pytest benchmark/scripts/tests/test_simulate_piao.py -v`
Expected: 5 tests total (3 from Task 2 + 2 new). The marginal-match test SKIPS if Piao raw data isn't local; the design-count test passes.

- [ ] **Step 4: If marginal-match fails, tune the baseline mixture**

If `test_simulator_marginals_match_piao_within_tolerance` fails (and Piao data was available), inspect:

```
uv run python -c "
import polars as pl
from pathlib import Path
# Adjust the path to whichever exists on your machine.
p = Path(r'D:/Coding/Projeler/methyl_lib/benchmarkin_merges/epykit_vs_allPackages(simulated_approxData)/raw_sim_data/simulated_datasets/dmc_simulation/coverage/amp.coverage=10.sample1.txt')
df = pl.read_csv(p, separator='\t')
print('coverage describe:', df['coverage'].describe())
print('freqC describe:', df['freqC'].describe())
print('freqC histogram:')
import numpy as np
hist, edges = np.histogram(df['freqC'].to_numpy(), bins=10)
for h, e in zip(hist, edges[:-1]):
    print(f'  [{e:5.1f}, {e+10:5.1f}): {h:>6,}')
"
```

The histogram tells you the baseline-beta marginal you're trying to match. Adjust the mixture probabilities (`[0.4, 0.4, 0.2]`) and/or the Beta parameters in `_draw_baseline_beta` to match. Then re-run.

- [ ] **Step 5: Commit (only if test passes or if you adjusted mixture)**

```
git add benchmark/scripts/tests/test_simulate_piao.py benchmark/scripts/simulate_piao.py
git commit -m "$(cat <<'EOF'
test(benchmark): simulate_piao validates marginals match Piao sample

Adds two marginal-match tests:
- count_M first two moments within 10% / 20% of Piao's coverage=10
  sample1.txt (skips gracefully if Piao raw data not local).
- truth DMC count = 20,000 ± 500 on 100K CpGs (design check).

If the marginal-match test fails locally, tune the
_draw_baseline_beta mixture parameters and document the change in
the commit message; record the comparison stats (Piao mean/std vs
simulator mean/std) for the bug-fix audit table.

Refs docs/history/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md §2.1

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: `wilson_bootstrap_ci.py` — Wilson CI on TPR/FPR

**What it does:** Adds Wilson 95 % CI lo/hi columns to any DataFrame with `tp, fp, tn, fn` counts. Operates on the existing `eval_summary.parquet` schema. The Wilson interval is closed-form via `scipy.stats.binomtest(k, n).proportion_ci(method="wilson")` (or equivalently `statsmodels.stats.proportion.proportion_confint`).

**Files:**
- Create: `epykit3/benchmark/scripts/wilson_bootstrap_ci.py`
- Create: `epykit3/benchmark/scripts/tests/test_wilson_bootstrap_ci.py`

- [ ] **Step 1: Write failing tests**

`epykit3/benchmark/scripts/tests/test_wilson_bootstrap_ci.py`:

```python
"""Tests for Wilson CI / bootstrap CI helpers."""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest


def test_add_wilson_ci_for_tpr_matches_scipy_reference():
    """Wilson 95% CI on a known proportion. Reference: scipy directly."""
    from wilson_bootstrap_ci import add_wilson_ci

    # 90 successes out of 100 -> p_hat = 0.9
    df = pl.DataFrame({
        "tp": [90], "fp": [10], "fn": [10], "tn": [90],
        "tpr": [0.9], "fpr": [10 / 100],
    })
    out = add_wilson_ci(df, rate="tpr", k_col="tp", n_col_expr=lambda d: d["tp"] + d["fn"])
    assert "tpr_ci_lo" in out.columns
    assert "tpr_ci_hi" in out.columns

    # Scipy reference for p_hat=0.9, n=100:
    from scipy.stats import binomtest
    ref = binomtest(90, 100).proportion_ci(method="wilson", confidence_level=0.95)
    assert abs(out["tpr_ci_lo"][0] - ref.low) < 1e-10
    assert abs(out["tpr_ci_hi"][0] - ref.high) < 1e-10


def test_add_wilson_ci_handles_zero_count_edges():
    """k=0 and k=n must not crash and must produce sensible intervals."""
    from wilson_bootstrap_ci import add_wilson_ci

    df = pl.DataFrame({
        "tp": [0, 100], "fp": [0, 0], "fn": [100, 0], "tn": [100, 100],
        "tpr": [0.0, 1.0], "fpr": [0.0, 0.0],
    })
    out = add_wilson_ci(df, rate="tpr", k_col="tp", n_col_expr=lambda d: d["tp"] + d["fn"])
    # k=0/n=100 Wilson lo = 0, hi ~ 0.037
    assert out["tpr_ci_lo"][0] == 0.0
    assert 0.02 < out["tpr_ci_hi"][0] < 0.05
    # k=100/n=100 Wilson lo ~ 0.963, hi = 1.0
    assert 0.95 < out["tpr_ci_lo"][1] < 0.98
    assert out["tpr_ci_hi"][1] == 1.0


def test_add_wilson_ci_zero_denominator_returns_nan():
    """If tp+fn = 0 (no positives in the truth set), CI is NaN."""
    from wilson_bootstrap_ci import add_wilson_ci

    df = pl.DataFrame({
        "tp": [0], "fp": [5], "fn": [0], "tn": [100],
        "tpr": [0.0], "fpr": [5 / 105],
    })
    out = add_wilson_ci(df, rate="tpr", k_col="tp", n_col_expr=lambda d: d["tp"] + d["fn"])
    assert np.isnan(out["tpr_ci_lo"][0])
    assert np.isnan(out["tpr_ci_hi"][0])
```

- [ ] **Step 2: Run failing tests**

Run: `uv run pytest benchmark/scripts/tests/test_wilson_bootstrap_ci.py -v`
Expected: 3 tests FAIL with ModuleNotFoundError.

- [ ] **Step 3: Write the Wilson-CI module (bootstrap added in Task 5)**

`epykit3/benchmark/scripts/wilson_bootstrap_ci.py`:

```python
"""Wilson 95% CIs on per-cell TPR/FPR and bootstrap CIs on AUROC/F1.

Reads the existing eval_summary.parquet schema:
  tool, scenario, parameter, parameter_value, test, meth_diff_bin,
  threshold_kind, threshold, tp, fp, tn, fn, tpr, fpr, precision, f1, auroc

Adds:
  tpr_ci_lo, tpr_ci_hi (Wilson, 95%)
  fpr_ci_lo, fpr_ci_hi (Wilson, 95%)
  auroc_ci_lo, auroc_ci_hi (bootstrap, 95%, B=1000) -- Task 5
  f1_ci_lo, f1_ci_hi (bootstrap, 95%, B=1000) -- Task 5

This module is pure-Python; it does not call any epykit engine. Re-runs
of the engines are not needed -- CIs operate on the counts already in
eval_summary.parquet.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import polars as pl
from scipy.stats import binomtest


def _wilson_single(k: int, n: int, confidence: float = 0.95) -> tuple[float, float]:
    """Return Wilson lo/hi for k successes / n trials. NaN/NaN if n == 0."""
    if n == 0:
        return (float("nan"), float("nan"))
    res = binomtest(k, n).proportion_ci(method="wilson", confidence_level=confidence)
    return (float(res.low), float(res.high))


def add_wilson_ci(
    df: pl.DataFrame,
    rate: str,
    k_col: str,
    n_col_expr: Callable[[pl.DataFrame], pl.Series],
    confidence: float = 0.95,
) -> pl.DataFrame:
    """Add `<rate>_ci_lo` and `<rate>_ci_hi` columns via Wilson interval.

    Parameters
    ----------
    df : input DataFrame with `<rate>`, `<k_col>` columns.
    rate : the rate column name (e.g. "tpr", "fpr").
    k_col : the integer-count column for successes (e.g. "tp" for TPR).
    n_col_expr : callable mapping the input DataFrame to a polars Series
        of trial counts (denominators). For TPR: lambda d: d["tp"] + d["fn"].
        For FPR: lambda d: d["fp"] + d["tn"].
    confidence : confidence level, default 0.95.
    """
    k = df[k_col].to_numpy().astype(np.int64)
    n = n_col_expr(df).to_numpy().astype(np.int64)
    lo = np.empty(len(k), dtype=np.float64)
    hi = np.empty(len(k), dtype=np.float64)
    for i in range(len(k)):
        lo[i], hi[i] = _wilson_single(int(k[i]), int(n[i]), confidence)
    return df.with_columns([
        pl.Series(f"{rate}_ci_lo", lo),
        pl.Series(f"{rate}_ci_hi", hi),
    ])


def add_wilson_ci_for_tpr_fpr(
    df: pl.DataFrame, confidence: float = 0.95,
) -> pl.DataFrame:
    """Convenience: add both tpr and fpr Wilson CIs to an eval_summary frame."""
    df = add_wilson_ci(
        df, rate="tpr", k_col="tp",
        n_col_expr=lambda d: d["tp"] + d["fn"],
        confidence=confidence,
    )
    df = add_wilson_ci(
        df, rate="fpr", k_col="fp",
        n_col_expr=lambda d: d["fp"] + d["tn"],
        confidence=confidence,
    )
    return df


def main(argv: list[str] | None = None) -> None:
    """CLI: `python wilson_bootstrap_ci.py --eval eval_summary.parquet --out eval_summary_with_ci.parquet`"""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eval", required=True, type=str, help="Path to eval_summary.parquet")
    parser.add_argument("--out", required=True, type=str, help="Output parquet path")
    args = parser.parse_args(argv)

    df = pl.read_parquet(args.eval)
    df = add_wilson_ci_for_tpr_fpr(df)
    # Bootstrap CIs for AUROC/F1 land in Task 5; this CLI stub leaves them.
    df.write_parquet(args.out)
    print(f"wrote {args.out} with TPR/FPR Wilson CIs added")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests, expect pass**

Run: `uv run pytest benchmark/scripts/tests/test_wilson_bootstrap_ci.py -v`
Expected: 3 PASS.

- [ ] **Step 5: Commit**

```
git add benchmark/scripts/wilson_bootstrap_ci.py benchmark/scripts/tests/test_wilson_bootstrap_ci.py
git commit -m "$(cat <<'EOF'
feat(benchmark): wilson_bootstrap_ci.py -- Wilson CIs on TPR/FPR

Closes the 'no CIs on headline numbers' gap from §2.3 of the design
spec for the proportion-based metrics. Adds <rate>_ci_lo /
<rate>_ci_hi columns to any DataFrame with tp/fp/tn/fn counts. CLI
reads/writes eval_summary.parquet schemas.

Bootstrap CIs for AUROC / F1 land in the next commit.

Refs docs/history/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md §2.3

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: `wilson_bootstrap_ci.py` — Bootstrap CI on AUROC/F1

**What it does:** Bootstraps CpGs with replacement (B=1000) and recomputes AUROC / F1 per bootstrap draw. Returns the 2.5 / 97.5 percentiles as the CI. **Requires the raw per-CpG joined DataFrame** (truth + predicted p-values), not just the per-cell counts. The CLI invocation therefore takes a different signature for this metric — it operates per-tool-per-scenario rather than per-row.

**Files:**
- Modify: `epykit3/benchmark/scripts/wilson_bootstrap_ci.py` (append bootstrap functions)
- Modify: `epykit3/benchmark/scripts/tests/test_wilson_bootstrap_ci.py` (append bootstrap tests)

- [ ] **Step 1: Write failing tests for bootstrap CI**

Append to `test_wilson_bootstrap_ci.py`:

```python
def test_bootstrap_auroc_ci_contains_point_estimate():
    """Bootstrap CI for AUROC must contain the original point estimate."""
    from wilson_bootstrap_ci import bootstrap_auroc_ci

    rng = np.random.default_rng(42)
    n = 1000
    is_dmc = rng.random(n) < 0.2
    # Pvalues correlated with is_dmc but noisy.
    pvalues = np.where(is_dmc, rng.beta(0.5, 5.0, n), rng.beta(5.0, 0.5, n))
    point = _auroc_reference(is_dmc, pvalues)

    lo, hi = bootstrap_auroc_ci(
        is_dmc=is_dmc, pvalues=pvalues, B=200, seed=42, confidence=0.95,
    )
    assert lo < point < hi, f"CI [{lo:.4f}, {hi:.4f}] does not contain point {point:.4f}"
    assert (hi - lo) < 0.10, f"CI [{lo:.4f}, {hi:.4f}] too wide (>0.1) for n=1000"


def test_bootstrap_auroc_ci_is_deterministic_with_seed():
    """Same seed -> same CI bounds."""
    from wilson_bootstrap_ci import bootstrap_auroc_ci

    rng = np.random.default_rng(7)
    n = 500
    is_dmc = rng.random(n) < 0.2
    pvalues = np.where(is_dmc, rng.beta(0.5, 5.0, n), rng.beta(5.0, 0.5, n))

    lo_a, hi_a = bootstrap_auroc_ci(is_dmc=is_dmc, pvalues=pvalues, B=100, seed=99)
    lo_b, hi_b = bootstrap_auroc_ci(is_dmc=is_dmc, pvalues=pvalues, B=100, seed=99)
    assert lo_a == lo_b
    assert hi_a == hi_b


def test_bootstrap_f1_ci_contains_point_estimate():
    """Same shape for F1 at a fixed q-threshold."""
    from wilson_bootstrap_ci import bootstrap_f1_ci

    rng = np.random.default_rng(1)
    n = 1000
    is_dmc = rng.random(n) < 0.2
    qvalues = np.where(is_dmc, rng.beta(0.5, 5.0, n), rng.beta(5.0, 0.5, n))

    # Compute point F1 at q < 0.05.
    pred = qvalues < 0.05
    tp = int((pred & is_dmc).sum())
    fp = int((pred & ~is_dmc).sum())
    fn = int((~pred & is_dmc).sum())
    point_f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else 0.0

    lo, hi = bootstrap_f1_ci(
        is_dmc=is_dmc, qvalues=qvalues, threshold=0.05, B=200, seed=1,
    )
    assert lo < point_f1 < hi, f"CI [{lo:.4f}, {hi:.4f}] does not contain point {point_f1:.4f}"


# --- helper for tests -------------------------------------------------------


def _auroc_reference(is_dmc: np.ndarray, pvalues: np.ndarray) -> float:
    """Reference AUROC via Mann-Whitney U."""
    score = 1.0 - pvalues
    n_pos = int(is_dmc.sum())
    n_neg = len(is_dmc) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    # Average ranks for ties.
    order = np.argsort(score)
    ranks = np.empty_like(order, dtype=np.float64)
    ranks[order] = np.arange(1, len(order) + 1)
    # Handle ties (rare here; use scipy if needed).
    sum_ranks_pos = ranks[is_dmc].sum()
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2
    return u / (n_pos * n_neg)
```

- [ ] **Step 2: Run failing tests**

Run: `uv run pytest benchmark/scripts/tests/test_wilson_bootstrap_ci.py -v`
Expected: 3 new tests FAIL with AttributeError (`bootstrap_auroc_ci` not defined).

- [ ] **Step 3: Append bootstrap implementations**

Append to `wilson_bootstrap_ci.py`:

```python
# --- Bootstrap CIs for AUROC and F1 ----------------------------------------


def _auroc_mwu(is_dmc: np.ndarray, score: np.ndarray) -> float:
    """AUROC via Mann-Whitney U with average ranks for ties."""
    n_pos = int(is_dmc.sum())
    n_neg = len(is_dmc) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(score, kind="mergesort")
    ranks = np.empty_like(order, dtype=np.float64)
    ranks[order] = np.arange(1, len(order) + 1, dtype=np.float64)
    sum_ranks_pos = ranks[is_dmc].sum()
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


def bootstrap_auroc_ci(
    is_dmc: np.ndarray,
    pvalues: np.ndarray,
    B: int = 1000,
    seed: int = 0,
    confidence: float = 0.95,
) -> tuple[float, float]:
    """Bootstrap AUROC CI by resampling CpGs with replacement.

    Returns (lo, hi) at the given confidence level (two-sided percentile).
    """
    rng = np.random.default_rng(seed)
    n = len(is_dmc)
    score = 1.0 - np.asarray(pvalues, dtype=np.float64)
    is_dmc = np.asarray(is_dmc, dtype=bool)

    boot = np.empty(B, dtype=np.float64)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        boot[b] = _auroc_mwu(is_dmc[idx], score[idx])
    alpha = (1.0 - confidence) / 2.0
    lo = float(np.nanpercentile(boot, 100 * alpha))
    hi = float(np.nanpercentile(boot, 100 * (1.0 - alpha)))
    return (lo, hi)


def bootstrap_f1_ci(
    is_dmc: np.ndarray,
    qvalues: np.ndarray,
    threshold: float = 0.05,
    B: int = 1000,
    seed: int = 0,
    confidence: float = 0.95,
) -> tuple[float, float]:
    """Bootstrap F1 CI at a fixed q-threshold, resampling CpGs with replacement."""
    rng = np.random.default_rng(seed)
    n = len(is_dmc)
    is_dmc = np.asarray(is_dmc, dtype=bool)
    pred = np.asarray(qvalues, dtype=np.float64) < threshold

    boot = np.empty(B, dtype=np.float64)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        tp = int((pred[idx] & is_dmc[idx]).sum())
        fp = int((pred[idx] & ~is_dmc[idx]).sum())
        fn = int((~pred[idx] & is_dmc[idx]).sum())
        denom = 2 * tp + fp + fn
        boot[b] = (2 * tp / denom) if denom > 0 else 0.0
    alpha = (1.0 - confidence) / 2.0
    lo = float(np.nanpercentile(boot, 100 * alpha))
    hi = float(np.nanpercentile(boot, 100 * (1.0 - alpha)))
    return (lo, hi)
```

- [ ] **Step 4: Run all tests in the file**

Run: `uv run pytest benchmark/scripts/tests/test_wilson_bootstrap_ci.py -v`
Expected: 6 tests PASS.

- [ ] **Step 5: Commit**

```
git add benchmark/scripts/wilson_bootstrap_ci.py benchmark/scripts/tests/test_wilson_bootstrap_ci.py
git commit -m "$(cat <<'EOF'
feat(benchmark): wilson_bootstrap_ci -- bootstrap CIs for AUROC and F1

Closes the threshold-free metric CI gap. Bootstraps CpGs with
replacement (B=1000 default) and reports 2.5/97.5 percentile bounds
for AUROC and F1. Operates on per-CpG joined frames (is_dmc +
pvalue/qvalue arrays) rather than the per-cell counts that Wilson
uses, so the calling convention differs.

Determinism: seeded RNG; two calls with the same seed yield
identical bounds.

Refs docs/history/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md §2.3

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: `run_null_calibration.py` — engine label-shuffle calibration

**What it does:** For each (engine, scenario), runs `k` label-shuffles where treatment / control labels are randomly assigned to the same samples (so no true DMCs exist), runs the engine on the shuffled labels, and records the observed FDR at nominal `q < 0.05`. Outputs a table with `engine, scenario, k_shuffle, observed_fdr, n_called, n_total` plus Wilson CIs on observed FDR via the helper from Task 4.

**Files:**
- Create: `epykit3/benchmark/scripts/run_null_calibration.py`
- Create: `epykit3/benchmark/scripts/tests/test_run_null_calibration.py`

- [ ] **Step 1: Write failing tests with a mocked engine**

`epykit3/benchmark/scripts/tests/test_run_null_calibration.py`:

```python
"""Tests for run_null_calibration.py.

The real null-calibration runner dispatches to epykit engines via
`ep.tl.dmc`, which is heavy and integration-test-only. These tests use
a mock engine that returns prebuilt q-value arrays so we can verify the
shuffle loop, the FDR computation, and the Wilson CI integration without
running real DMC calls.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest


def test_compute_observed_fdr_matches_definition():
    """Observed FDR at nominal q<0.05 = (# q<0.05) / total. On a fixed
    array of qvalues, the function must compute this exactly."""
    from run_null_calibration import compute_observed_fdr

    qvals = np.array([0.001, 0.01, 0.04, 0.05, 0.1, 0.5, 0.9])
    # 3 sites with q<0.05 out of 7.
    out = compute_observed_fdr(qvals, q_thresh=0.05)
    assert out["n_called"] == 3
    assert out["n_total"] == 7
    assert abs(out["observed_fdr"] - 3 / 7) < 1e-12


def test_compute_observed_fdr_ignores_nan():
    """NaN q-values are excluded from the denominator."""
    from run_null_calibration import compute_observed_fdr

    qvals = np.array([0.001, 0.01, np.nan, np.nan, 0.5])
    out = compute_observed_fdr(qvals, q_thresh=0.05)
    assert out["n_called"] == 2
    assert out["n_total"] == 3
    assert abs(out["observed_fdr"] - 2 / 3) < 1e-12


def test_run_shuffles_returns_one_row_per_shuffle_with_ci(tmp_path):
    """The runner returns a frame with k rows per (engine, scenario),
    each carrying observed_fdr and tpr_ci_lo/tpr_ci_hi (here repurposed
    as Wilson CI bounds on the observed FDR proportion)."""
    from run_null_calibration import run_null_calibration

    # Mock engine: returns deterministic q-values per shuffle seed.
    def mock_engine(samples_treatment, samples_control, seed=0, n_sites=200):
        rng = np.random.default_rng(seed)
        # Under-null engine: q-values uniform-ish on [0, 1].
        return rng.uniform(0, 1, size=n_sites)

    samples = [f"s{i}" for i in range(1, 7)]
    out = run_null_calibration(
        engine_fn=mock_engine,
        engine_name="mock_lr",
        scenario_name="cov10_3v3",
        samples=samples,
        n_per_group=3,
        k_shuffles=10,
        q_thresh=0.05,
        seed=42,
    )
    assert isinstance(out, pl.DataFrame)
    assert set(out.columns) >= {
        "engine", "scenario", "k_shuffle", "observed_fdr",
        "n_called", "n_total",
        "observed_fdr_ci_lo", "observed_fdr_ci_hi",
    }
    assert len(out) == 10
    # All entries reference the same engine + scenario.
    assert out["engine"].unique().to_list() == ["mock_lr"]
    assert out["scenario"].unique().to_list() == ["cov10_3v3"]
    # Wilson CI bounds bracket observed_fdr.
    assert (out["observed_fdr_ci_lo"] <= out["observed_fdr"]).all()
    assert (out["observed_fdr_ci_hi"] >= out["observed_fdr"]).all()


def test_run_shuffles_is_deterministic_with_seed(tmp_path):
    """Two runs with same seed give identical observed_fdr column."""
    from run_null_calibration import run_null_calibration

    def mock_engine(samples_treatment, samples_control, seed=0, n_sites=200):
        rng = np.random.default_rng(seed)
        return rng.uniform(0, 1, size=n_sites)

    samples = [f"s{i}" for i in range(1, 7)]
    a = run_null_calibration(
        engine_fn=mock_engine, engine_name="mock", scenario_name="s1",
        samples=samples, n_per_group=3, k_shuffles=5, q_thresh=0.05, seed=7,
    )
    b = run_null_calibration(
        engine_fn=mock_engine, engine_name="mock", scenario_name="s1",
        samples=samples, n_per_group=3, k_shuffles=5, q_thresh=0.05, seed=7,
    )
    np.testing.assert_array_equal(
        a["observed_fdr"].to_numpy(), b["observed_fdr"].to_numpy(),
    )
```

- [ ] **Step 2: Run failing tests**

Run: `uv run pytest benchmark/scripts/tests/test_run_null_calibration.py -v`
Expected: 4 tests FAIL with ModuleNotFoundError.

- [ ] **Step 3: Write the module**

`epykit3/benchmark/scripts/run_null_calibration.py`:

```python
"""Null-calibration runner: empirical FDR under label-shuffled data.

For each (engine, scenario), randomly shuffles treatment/control labels
over the same samples, re-runs the engine, and records the observed
proportion of sites called significant at the nominal threshold. With
no true DMCs in the shuffled design, observed FDR at nominal q < 0.05
should be ~ 0.05 if the test is well-calibrated, OR much lower if the
test is conservative on the input data's noise regime.

See `docs/history/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md`
§2.2.

This module exposes:

- `compute_observed_fdr(qvalues, q_thresh)`: the pure-arithmetic kernel.
- `run_null_calibration(...)`: orchestrates k shuffles, calls the engine
  closure once per shuffle, returns a polars DataFrame with one row per
  shuffle and Wilson CI bounds on observed FDR.

The engine closure has signature
    engine_fn(samples_treatment, samples_control, seed=int) -> np.ndarray
returning per-site q-values for that shuffle. Real callers wrap epykit's
``ep.tl.dmc`` (or equivalent) in such a closure; tests use a fake.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import polars as pl

from wilson_bootstrap_ci import _wilson_single


def compute_observed_fdr(qvalues: np.ndarray, q_thresh: float = 0.05) -> dict:
    """Per-shuffle stats: n_called, n_total, observed_fdr."""
    q = np.asarray(qvalues, dtype=np.float64)
    finite = np.isfinite(q)
    q_clean = q[finite]
    n_called = int((q_clean < q_thresh).sum())
    n_total = int(len(q_clean))
    if n_total == 0:
        return {"n_called": 0, "n_total": 0, "observed_fdr": float("nan")}
    return {
        "n_called": n_called,
        "n_total": n_total,
        "observed_fdr": float(n_called / n_total),
    }


def run_null_calibration(
    engine_fn: Callable,
    engine_name: str,
    scenario_name: str,
    samples: list[str],
    n_per_group: int,
    k_shuffles: int = 20,
    q_thresh: float = 0.05,
    seed: int = 0,
) -> pl.DataFrame:
    """Run k label-shuffles and aggregate observed FDR per shuffle.

    Returns a DataFrame with columns:
      engine, scenario, k_shuffle, n_called, n_total, observed_fdr,
      observed_fdr_ci_lo, observed_fdr_ci_hi.

    Wilson CIs treat each shuffle's observed FDR as a binomial proportion
    of `n_called` / `n_total`. They quantify the within-shuffle estimation
    uncertainty; for across-shuffle variability, compute the median + IQR
    over rows externally.
    """
    rng = np.random.default_rng(seed)
    n = len(samples)
    if n_per_group * 2 != n:
        raise ValueError(
            f"need {n_per_group * 2} samples for {n_per_group}v{n_per_group}, got {n}"
        )

    rows = []
    for k in range(k_shuffles):
        # Local RNG per shuffle so a single global seed reproduces the run.
        shuffled = rng.permutation(samples).tolist()
        treat = shuffled[:n_per_group]
        ctrl = shuffled[n_per_group:]

        qvals = engine_fn(samples_treatment=treat, samples_control=ctrl, seed=seed + k + 1)
        stats = compute_observed_fdr(qvals, q_thresh=q_thresh)
        lo, hi = _wilson_single(stats["n_called"], stats["n_total"])
        rows.append({
            "engine": engine_name,
            "scenario": scenario_name,
            "k_shuffle": k,
            "n_called": stats["n_called"],
            "n_total": stats["n_total"],
            "observed_fdr": stats["observed_fdr"],
            "observed_fdr_ci_lo": lo,
            "observed_fdr_ci_hi": hi,
        })
    return pl.DataFrame(rows)


def main(argv: list[str] | None = None) -> None:
    """CLI stub. Real callers wire ``ep.tl.dmc`` here; this stub demonstrates
    the integration surface with a deterministic noise engine."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", default="noise", help="Engine label for reporting")
    parser.add_argument("--scenario", default="demo", help="Scenario label")
    parser.add_argument("--k-shuffles", type=int, default=20)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True, type=str, help="Output parquet path")
    args = parser.parse_args(argv)

    def noise_engine(samples_treatment, samples_control, seed=0):
        rng = np.random.default_rng(seed)
        return rng.uniform(0, 1, size=1000)

    df = run_null_calibration(
        engine_fn=noise_engine,
        engine_name=args.engine,
        scenario_name=args.scenario,
        samples=["s1", "s2", "s3", "s4", "s5", "s6"],
        n_per_group=3,
        k_shuffles=args.k_shuffles,
        seed=args.seed,
    )
    df.write_parquet(args.out)
    print(f"wrote {args.out} with {len(df)} shuffle rows")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests, expect pass**

Run: `uv run pytest benchmark/scripts/tests/test_run_null_calibration.py -v`
Expected: 4 PASS.

- [ ] **Step 5: Smoke-run the CLI**

Run:
```
uv run python benchmark/scripts/run_null_calibration.py --engine demo --scenario test --k-shuffles 5 --seed 1 --out /tmp/null_calib_smoke.parquet
uv run python -c "
import polars as pl
df = pl.read_parquet('/tmp/null_calib_smoke.parquet')
print(df)
"
```

Expected: prints a 5-row DataFrame with `engine=demo`, `scenario=test`, observed FDR around 0.05 (since the noise engine returns uniform q-values).

- [ ] **Step 6: Commit**

```
git add benchmark/scripts/run_null_calibration.py benchmark/scripts/tests/test_run_null_calibration.py
git commit -m "$(cat <<'EOF'
feat(benchmark): run_null_calibration.py -- engine label-shuffle calibration

Closes the 'no null calibration' gap from §2.2 of the design spec.
For each (engine, scenario), runs k label-shuffles of treatment /
control over the same samples (no true DMCs by construction) and
reports observed FDR at nominal q<0.05 with Wilson CI per shuffle.

If observed FDR ~ 0.05 +/- CI on a well-powered shuffle, the test is
well-calibrated. If observed FDR is much smaller (e.g. 1e-5), the
test is conservative on the input data's noise regime (e.g. Piao's
underdispersed simulator + lr default).

Module decoupled from epykit -- takes an `engine_fn` closure so it
tests with a mock noise engine. Real callers wire ep.tl.dmc.

Refs docs/history/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md §2.2

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Wrap-up (after Task 6)

- [ ] **Step 1: Full benchmark test pass**

Run: `uv run pytest benchmark/scripts/tests/ -v`
Expected: all tests PASS (3 + 2 + 3 + 3 + 4 = 15 minimum, may be more if marginal-match test ran).

- [ ] **Step 2: Main epykit test pass (no regressions from consolidation)**

Run: `uv run pytest -m "not slow" --strict-markers -x -q`
Expected: 229 passed, 5 skipped, 0 failed (matches Phase 1 wrap-up baseline).

- [ ] **Step 3: ruff on the new code**

Run: `uv run ruff check benchmark/scripts/`
Expected: no F-level issues in the new files. Pre-existing issues elsewhere are unchanged.

- [ ] **Step 4: CHANGELOG sanity check**

Open `CHANGELOG.md`. The Unreleased section should now contain:

- `### Added` with bullets for: benchmark consolidation, simulate_piao.py, wilson_bootstrap_ci.py, run_null_calibration.py.
- `### Fixed (P0 manifest, paper preparation)` (from Phase 1) — unchanged.
- `### Changed (breaking on the lr+ schema)` (from Phase 1) — unchanged.

If any bullets are missing, add them now.

- [ ] **Step 5: Tag**

```
git tag -a v0.7.4-phase2-scripts -m "Phase 2 benchmark scripts complete

simulate_piao.py + wilson_bootstrap_ci.py + run_null_calibration.py
landed; benchmarkin_merges/FINAL_REPORT consolidated into
epykit3/benchmark/. Next: Phase 3 -- locked benchmark re-run with
post-P0 epykit, multi-seed simulator runs, CI columns on every
headline cell, null-calibration table."
```

Hold the tag until the user confirms (same rule as Phase 1).

- [ ] **Step 6: Brief out**

Summarise (one paragraph): which scripts landed, what each is for, what's still needed before the locked re-run. Add the summary as a comment on `c53178d`'s tag, or append to the design spec's §7 as a Phase 2 closeout note.

---

## Self-review (writer's pass, completed inline)

- **Spec coverage:** §2.1 truth tautology → simulator (Tasks 2-3); §2.2 null calibration → Task 6; §2.3 CIs → Tasks 4-5; Open Q1 consolidation → Task 1. All four §2 gaps covered.
- **Out-of-scope confirmed:** `methylkit_stouffer_combine.R`, `regen_all.py`, `bug_fix_audit.py` are not present in any task.
- **Placeholder scan:** no TBDs, no "implement later", no "add error handling" without showing the handler, no "similar to Task N" without inline code.
- **Type consistency:** truth parquet schema (`chrom: Utf8, pos: Int64, ...`) referenced consistently across Tasks 2-3 and matches the existing `dmc_truth.parquet` it's drop-in compatible with. `engine_fn` callable signature (`samples_treatment, samples_control, seed`) consistent across Task 6 tests and implementation. `_wilson_single` helper imported from `wilson_bootstrap_ci` into `run_null_calibration` — name matches.
- **Test layout decision:** new tests at `benchmark/scripts/tests/` independent of `tests/`. `conftest.py` injects the scripts dir into `sys.path` so tests do bare `from simulate_piao import ...`. This avoids needing to install the scripts as a package.
- **Risk acknowledged in Task 1 Step 5:** large directory move size check before commit — prevents accidental GB-scale commits.
- **Risk acknowledged in Task 3 Step 4:** marginal-match tolerance is wide because we don't know Piao's exact baseline model. The test skips when Piao raw data isn't local, so the suite is portable.
