# DMC DMR Power Calibration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Improve epykit DMC/DMR recall at low coverage and small effect sizes while preserving explicit FPR calibration gates.

**Architecture:** Treat this as a calibration change, not a pure parameter edit. Each default change gets a focused unit/regression test, then a small synthetic power/FPR benchmark gate, then documentation/benchmark metadata updates so claims cannot drift from behavior.

**Tech Stack:** Python 3.9+, `polars`, `numpy`, `scipy`, `statsmodels`, `pytest`, `uv`; no new runtime dependencies.

## Global Constraints

- Preserve the methylstore streaming contract: no full-genome materialization inside DMC/DMR engines except the existing global FDR vector step.
- Preserve Windows compatibility; all new tests must pass under the CI matrix `{ubuntu-latest, windows-latest} x {py3.9, py3.12}`.
- Library code under `epykit.*` must use stdlib `logging`, not `print()`.
- Do not make `power_stack="lr+"` the default; benchmark notes show it can inflate FPR on real WGBS.
- Any `DF_PHI_FLOOR` change must update `tests/test_principled_df.py` and include a null-calibration benchmark artifact.
- Use `uv run pytest -m "not slow" --strict-markers -ra` for CI-equivalent verification.
- Do not update benchmark headline claims until the benchmark scripts have been re-run and the report tables regenerated.

---

## File Structure

- Modify `src/epykit/tl.py`: high-level `ep.tl.dmc` and `ep.tl.dmr` defaults, parameter metadata recorded in `md.uns`.
- Modify `src/epykit/cli.py`: CLI parity for `dmc --fdr-method`, `dmr --min-cpgs`, and any DMC smoothing default if accepted.
- Modify `src/epykit/dmc.py`: lower-level DMC defaults only if API/CLI parity requires it; `DF_PHI_FLOOR` only after calibration passes.
- Modify `src/epykit/dmr.py`: chain-merge preset/default resolution and tuning guidance.
- Modify `tests/test_cli_api_parity.py`: API/CLI default parity for DMC dispersion/FDR behavior.
- Modify `tests/test_lr_improvements.py`: FDR method default and `power_stack` non-default guards.
- Modify `tests/test_smoothed_dmc.py`: DSS-style count smoothing behavior and cache-key coverage.
- Modify `tests/test_principled_df.py`: `DF_PHI_FLOOR` calibration expectation, only if the floor changes.
- Modify `tests/test_dmr_min_cpgs_parity.py`: chain-merge high-level default now resolves to 3.
- Modify `tests/test_dmr_presets_and_diagnose.py`: preset tuning expectations.
- Add `tests/test_power_calibration_defaults.py`: lightweight synthetic TPR/FPR guards for selected default profiles.
- Add `benchmark/scripts/run_default_calibration_grid.py`: reproducible benchmark grid over low coverage, small effect, and null scenarios.
- Add `benchmark/data/default_calibration/README.md`: schema and regeneration command for calibration outputs.

---

### Task 1: Add Calibration Tests Before Changing Defaults

**Files:**
- Create: `tests/test_power_calibration_defaults.py`
- Modify: none

**Interfaces:**
- Consumes: `tests.fixtures.synth`, `ep.tl.dmc`, `ep.tl.dmr`
- Produces: regression tests named `test_lr_default_small_effect_power_floor`, `test_lr_default_null_fpr_ceiling`, `test_chain_merge_default_recovers_three_cpg_region`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_power_calibration_defaults.py`:

```python
from __future__ import annotations

import copy

import polars as pl

import epykit as ep
from epykit.dmr import call_dmr_chain_merge


def _score_dmc(df: pl.DataFrame, truth: pl.DataFrame, *, alpha: float = 0.05) -> tuple[float, float]:
    joined = df.select("chrom", "pos", "qvalue").join(
        truth.select("chrom", "pos", "is_dmc"),
        on=["chrom", "pos"],
        how="inner",
    )
    called = joined["qvalue"].fill_null(1.0) < alpha
    truth_mask = joined["is_dmc"]
    tp = int((called & truth_mask).sum())
    fp = int((called & ~truth_mask).sum())
    fn = int((~called & truth_mask).sum())
    tn = int((~called & ~truth_mask).sum())
    tpr = tp / max(tp + fn, 1)
    fpr = fp / max(fp + tn, 1)
    return tpr, fpr


def test_lr_default_small_effect_power_floor(synth_md_filtered, synth_bundle):
    """Default lr should not collapse on the low-coverage small-effect slice."""
    md = copy.deepcopy(synth_md_filtered)
    ep.tl.dmc(md, test="lr")
    tpr, _fpr = _score_dmc(md.get_dmc(test="lr"), synth_bundle.truth)
    assert tpr >= 0.08


def test_lr_default_null_fpr_ceiling(synth_md_filtered, synth_bundle):
    """Default lr tuning must preserve a strict null-side ceiling."""
    md = copy.deepcopy(synth_md_filtered)
    ep.tl.dmc(md, test="lr")
    _tpr, fpr = _score_dmc(md.get_dmc(test="lr"), synth_bundle.truth)
    assert fpr <= 0.02


def test_chain_merge_default_recovers_three_cpg_region():
    dmc = pl.DataFrame(
        {
            "chrom": ["chr1"] * 3,
            "pos": [100, 180, 260],
            "pvalue": [1e-7, 1e-7, 1e-7],
            "qvalue": [1e-7, 1e-7, 1e-7],
            "meth_diff": [0.22, 0.24, 0.20],
        }
    )
    out = call_dmr_chain_merge(dmc, alpha=0.05, min_abs_meth_diff=0.1)
    assert len(out) == 1
    assert out.item(0, "n_cpgs") == 3
```

- [ ] **Step 2: Run tests to verify current behavior**

Run:

```bash
uv run pytest tests/test_power_calibration_defaults.py -ra
```

Expected: the two DMC tests pass or establish the current baseline; the DMR direct-engine test should pass because `call_dmr_chain_merge()` already falls back to `min_cpgs=3`.

- [ ] **Step 3: Commit**

```bash
git add tests/test_power_calibration_defaults.py
git commit -m "test: add default power calibration guards"
```

---

### Task 2: Tune Chain-Merge High-Level Defaults

**Files:**
- Modify: `src/epykit/dmr.py`
- Modify: `src/epykit/tl.py`
- Modify: `src/epykit/cli.py`
- Modify: `tests/test_dmr_min_cpgs_parity.py`
- Modify: `tests/test_dmr_presets_and_diagnose.py`

**Interfaces:**
- Consumes: `resolve_layer_min_cpgs(min_cpgs: int | None, preset: str | None) -> int`
- Produces: high-level bare `ep.tl.dmr(method="chain_merge")` and CLI bare `epykit dmr --method chain_merge` resolve `min_cpgs=3`

- [ ] **Step 1: Update tests for the new high-level default**

In `tests/test_dmr_min_cpgs_parity.py`, change the bare API/CLI assertions from `5` to `3`:

```python
assert md.uns["dmr_params"]["min_cpgs"] == 3
```

and:

```python
assert got == 3
```

- [ ] **Step 2: Run the affected tests and verify they fail**

Run:

```bash
uv run pytest tests/test_dmr_min_cpgs_parity.py tests/test_dmr_presets_and_diagnose.py -ra
```

Expected: failures show the current high-level bare default is still `5`.

- [ ] **Step 3: Change the high-level default resolver**

In `src/epykit/dmr.py`, change:

```python
_DMR_DEFAULT_MIN_CPGS = 5
```

to:

```python
_DMR_DEFAULT_MIN_CPGS = 3
```

Update the adjacent docstring to say this now matches the direct DSS-style engine default.

- [ ] **Step 4: Keep `pct_sig=0.5` as the balanced default and make permissive explicitly recall-oriented**

In `src/epykit/dmr.py`, keep `DMR_PRESETS["default"]["pct_sig"]` at `0.5`. Change only `DMR_PRESETS["permissive"]["pct_sig"]`:

```python
"permissive": dict(
    alpha=1e-4, min_abs_meth_diff=0.05, dis_merge_bp=1000,
    min_cpgs=3, pct_sig=0.4, minlen_bp=50,
),
```

Update the preset summary in `src/epykit/tl.py` to mention `pct_sig=0.4` for `permissive`.

- [ ] **Step 5: Run chain-merge parity tests**

Run:

```bash
uv run pytest tests/test_dmr_chain_merge.py tests/test_dmr_min_cpgs_parity.py tests/test_dmr_presets_and_diagnose.py tests/test_cli_dmr_qfilter.py -ra
```

Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add src/epykit/dmr.py src/epykit/tl.py src/epykit/cli.py tests/test_dmr_min_cpgs_parity.py tests/test_dmr_presets_and_diagnose.py
git commit -m "feat: tune chain-merge defaults for recall"
```

---

### Task 3: Switch Default DMC FDR to Two-Stage BH With FPR Guardrails

**Files:**
- Modify: `src/epykit/tl.py`
- Modify: `src/epykit/cli.py`
- Modify: `tests/test_lr_improvements.py`
- Modify: `tests/test_cli_dmc_contrast_forwarding.py`
- Modify: `tests/test_cli_api_parity.py`
- Modify: `tests/test_power_stack_equivalence.py`

**Interfaces:**
- Consumes: `apply_multiple_testing_correction(dmc_results, method="fdr_tsbh")`
- Produces: default `ep.tl.dmc(..., fdr_method="fdr_tsbh")` and CLI `epykit dmc` parity

- [ ] **Step 1: Update default tests**

In `tests/test_lr_improvements.py`, add:

```python
def test_tl_dmc_default_fdr_method_is_tsbh():
    import inspect
    import epykit as ep

    assert inspect.signature(ep.tl.dmc).parameters["fdr_method"].default == "fdr_tsbh"
```

In `tests/test_cli_api_parity.py`, add a CLI parser assertion:

```python
def test_cli_dmc_default_fdr_method_is_tsbh():
    from epykit.cli import build_parser

    args = build_parser().parse_args(["dmc", "--methylstore", "store", "--output", "out.parquet"])
    assert args.fdr_method == "fdr_tsbh"
```

- [ ] **Step 2: Run tests and verify they fail**

Run:

```bash
uv run pytest tests/test_lr_improvements.py::test_tl_dmc_default_fdr_method_is_tsbh tests/test_cli_api_parity.py::test_cli_dmc_default_fdr_method_is_tsbh -ra
```

Expected: both fail with current `fdr_bh` defaults.

- [ ] **Step 3: Change API and CLI defaults**

In `src/epykit/tl.py`, change:

```python
fdr_method: str = "fdr_bh",
```

to:

```python
fdr_method: str = "fdr_tsbh",
```

In `src/epykit/cli.py`, change the `--fdr-method` default:

```python
default="fdr_tsbh",
help="Multiple-testing correction method (default: fdr_tsbh).",
```

- [ ] **Step 4: Preserve power-stack semantics**

Update tests that assumed `power_stack="off"` means `fdr_bh` by default. The invariant should become:

```python
ep.tl.dmc(md, test="lr", power_stack="off", neighbour_combine=False)
assert md.uns["dmc"]["fdr_method"] == "fdr_tsbh"
```

Keep `power_stack="lr+"` equivalence pinned to explicit `fdr_method="fdr_tsbh"`.

- [ ] **Step 5: Run FDR and parity tests**

Run:

```bash
uv run pytest tests/test_fdr_nan_exclusion.py tests/test_lr_improvements.py tests/test_power_stack_equivalence.py tests/test_cli_api_parity.py tests/test_cli_dmc_contrast_forwarding.py -ra
```

Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add src/epykit/tl.py src/epykit/cli.py tests/test_lr_improvements.py tests/test_cli_api_parity.py tests/test_cli_dmc_contrast_forwarding.py tests/test_power_stack_equivalence.py
git commit -m "feat: default dmc fdr to two-stage bh"
```

---

### Task 4: Decide Whether DSS-Style Count Smoothing Can Be a Default

**Files:**
- Modify: `src/epykit/tl.py`
- Modify: `src/epykit/cli.py`
- Modify: `src/epykit/dmc.py`
- Modify: `tests/test_smoothed_dmc.py`
- Modify: `tests/test_empirical_fdr_sep_fallback_parity.py`

**Interfaces:**
- Consumes: `process_chromosomes_dmc(..., smoothing: bool, smoothing_span_bp: int)`
- Produces: either default smoothing for `test="lr"` or a named opt-in profile if FPR/TPR gates fail

- [ ] **Step 1: Add explicit tests for the desired default**

In `tests/test_smoothed_dmc.py`, add:

```python
def test_tl_dmc_lr_default_records_smoothing_enabled(synth_md_filtered):
    ep.tl.dmc(synth_md_filtered, test="lr")
    assert synth_md_filtered.uns["dmc"]["smoothing"] is True
    assert synth_md_filtered.uns["dmc"]["smoothing_span_bp"] == 500
```

- [ ] **Step 2: Run the new test and verify it fails**

Run:

```bash
uv run pytest tests/test_smoothed_dmc.py::test_tl_dmc_lr_default_records_smoothing_enabled -ra
```

Expected: fails because `smoothing=False` today.

- [ ] **Step 3: Change only the high-level API/CLI default first**

In `src/epykit/tl.py`, change:

```python
smoothing: bool = False,
```

to:

```python
smoothing: bool = True,
```

In `src/epykit/cli.py`, expose the same behavior by making the CLI DMC smoothing flag default true and adding a negative flag:

```python
p_dmc.set_defaults(smoothing=True)
p_dmc.add_argument("--no-smoothing", dest="smoothing", action="store_false")
```

Do not change lower-level `process_chromosomes_dmc(..., smoothing=False)` until API/CLI parity and calibration are proven; direct engine callers should remain explicit.

- [ ] **Step 4: Run calibration tests**

Run:

```bash
uv run pytest tests/test_smoothed_dmc.py tests/test_empirical_fdr_sep_fallback_parity.py tests/test_power_calibration_defaults.py -ra
```

Expected: pass only if smoothing does not breach the FPR guard. If FPR exceeds the guard, revert this task and implement a named `dmc_profile="dss_smooth"` instead of changing the default.

- [ ] **Step 5: Commit only if FPR guard passes**

```bash
git add src/epykit/tl.py src/epykit/cli.py tests/test_smoothed_dmc.py tests/test_empirical_fdr_sep_fallback_parity.py
git commit -m "feat: enable lr count smoothing by default"
```

---

### Task 5: Revisit `DF_PHI_FLOOR` With Calibration Evidence

**Files:**
- Modify: `src/epykit/dmc.py`
- Modify: `tests/test_principled_df.py`
- Modify: `tests/test_p0_eb_adaptive_f_floor.py`
- Modify: `tests/test_overdispersed_calibration.py`

**Interfaces:**
- Consumes: `DF_PHI_FLOOR` and `_score_finalize(..., reference="adaptive")`
- Produces: either `DF_PHI_FLOOR=20.0` or a rejected-change note in benchmark outputs

- [ ] **Step 1: Update the pinning test to express the proposed floor**

In `tests/test_principled_df.py`, change:

```python
assert DF_PHI_FLOOR == 50.0
```

to:

```python
assert DF_PHI_FLOOR == 20.0
```

Change `_DF_PHI_FLOOR_F_VS_CHI2_TOL_AT_P05` expectation to the actual relative excess for F(1, 20) versus chi-square at p=0.05:

```python
assert 0.0 <= relative_excess <= 0.35
```

- [ ] **Step 2: Run the pinning test and verify it fails**

Run:

```bash
uv run pytest tests/test_principled_df.py -ra
```

Expected: fails because `DF_PHI_FLOOR` remains `50.0`.

- [ ] **Step 3: Change the constant**

In `src/epykit/dmc.py`, change:

```python
DF_PHI_FLOOR: float = 50.0
_DF_PHI_FLOOR_F_VS_CHI2_TOL_AT_P05: float = 0.12
```

to:

```python
DF_PHI_FLOOR: float = 20.0
_DF_PHI_FLOOR_F_VS_CHI2_TOL_AT_P05: float = 0.35
```

- [ ] **Step 4: Run overdispersed calibration**

Run:

```bash
uv run pytest tests/test_principled_df.py tests/test_p0_eb_adaptive_f_floor.py tests/test_overdispersed_calibration.py -ra
```

Expected: all pass. If `test_overdispersed_calibration.py` fails, revert the floor to `50.0` and do not ship this task.

- [ ] **Step 5: Commit only if calibration passes**

```bash
git add src/epykit/dmc.py tests/test_principled_df.py tests/test_p0_eb_adaptive_f_floor.py tests/test_overdispersed_calibration.py
git commit -m "feat: lower adaptive dispersion df floor"
```

---

### Task 6: Repair or Demote `welch_t` for Small-N DMC

**Files:**
- Modify: `src/epykit/dmc.py`
- Modify: `src/epykit/tl.py`
- Modify: `tests/test_accuracy.py`
- Modify: `tests/test_engine_oracles.py`
- Modify: `tests/test_overdispersed_calibration.py`
- Modify: `tests/test_cli_allow_n1_fisher.py`

**Interfaces:**
- Consumes: existing `test="welch_t"` path
- Produces: explicit behavior for `welch_t` at `n<4`: either hard error or automatic documented fallback

- [ ] **Step 1: Add a small-N behavior test**

In `tests/test_cli_allow_n1_fisher.py`, add:

```python
def test_welch_t_small_n_is_rejected_without_force():
    from epykit.dmc import _warn_if_low_power_test

    with pytest.raises(ValueError, match="welch_t requires at least 4 samples per group"):
        _warn_if_low_power_test("welch_t", min_n=3, force=False)
```

- [ ] **Step 2: Run the test and verify it fails**

Run:

```bash
uv run pytest tests/test_cli_allow_n1_fisher.py::test_welch_t_small_n_is_rejected_without_force -ra
```

Expected: fails because current behavior warns but allows the engine.

- [ ] **Step 3: Implement explicit small-N rejection or force flag**

In `src/epykit/dmc.py`, update the low-power guard so `welch_t` with `min_n < 4` raises `ValueError` unless an existing or new force flag is passed:

```python
raise ValueError(
    "welch_t requires at least 4 samples per group for usable variance; "
    "use test='lr' for small-n count-aware DMC calling."
)
```

- [ ] **Step 4: Run Welch tests**

Run:

```bash
uv run pytest tests/test_accuracy.py::test_dmc_welch_t_power_and_fdr tests/test_engine_oracles.py tests/test_overdispersed_calibration.py -ra
```

Expected: existing oracle behavior remains valid for supported sample sizes.

- [ ] **Step 5: Commit**

```bash
git add src/epykit/dmc.py src/epykit/tl.py tests/test_accuracy.py tests/test_engine_oracles.py tests/test_overdispersed_calibration.py tests/test_cli_allow_n1_fisher.py
git commit -m "fix: guard welch t at unsupported small sample sizes"
```

---

### Task 7: Add a Reproducible Default Calibration Grid

**Files:**
- Create: `benchmark/scripts/run_default_calibration_grid.py`
- Create: `benchmark/data/default_calibration/README.md`
- Modify: `benchmark/README.md`

**Interfaces:**
- Consumes: epykit public API and existing benchmark data layout
- Produces: parquet/TSV outputs with `tool`, `profile`, `coverage`, `effect`, `n_per_group`, `tpr`, `fpr`, `f1`, `runtime_s`

- [ ] **Step 1: Create the benchmark script**

Create `benchmark/scripts/run_default_calibration_grid.py`:

```python
from __future__ import annotations

import argparse
import time
from pathlib import Path

import polars as pl


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, default=Path("benchmark/data/default_calibration/results.tsv"))
    p.add_argument("--quick", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()
    rows: list[dict] = []
    coverages = [10] if args.quick else [5, 10, 15, 20]
    effects = [0.2, 0.3, 0.4] if args.quick else [0.1, 0.2, 0.3, 0.4, 0.6]
    n_groups = [3] if args.quick else [2, 3, 5]
    for coverage in coverages:
        for effect in effects:
            for n_per_group in n_groups:
                start = time.perf_counter()
                rows.append(
                    {
                        "tool": "epykit",
                        "profile": "default",
                        "coverage": coverage,
                        "effect": effect,
                        "n_per_group": n_per_group,
                        "tpr": None,
                        "fpr": None,
                        "f1": None,
                        "runtime_s": time.perf_counter() - start,
                        "status": "pending_simulator_hook",
                    }
                )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows).write_csv(args.out, separator="\\t")


if __name__ == "__main__":
    main()
```

This script starts with the output schema and CLI. In implementation, replace the `pending_simulator_hook` row body with the existing benchmark simulator runner used by the paper grid.

- [ ] **Step 2: Add calibration output README**

Create `benchmark/data/default_calibration/README.md`:

```markdown
# Default Calibration Grid

Generated by:

```bash
uv run python benchmark/scripts/run_default_calibration_grid.py
```

Required columns:

- `tool`
- `profile`
- `coverage`
- `effect`
- `n_per_group`
- `tpr`
- `fpr`
- `f1`
- `runtime_s`

The default profile may ship only when the null FPR ceiling remains within the agreed guardrail and the small-effect TPR improves over the pre-tuning baseline.
```

- [ ] **Step 3: Run the quick command**

Run:

```bash
uv run python benchmark/scripts/run_default_calibration_grid.py --quick
```

Expected: writes `benchmark/data/default_calibration/results.tsv`.

- [ ] **Step 4: Commit**

```bash
git add benchmark/scripts/run_default_calibration_grid.py benchmark/data/default_calibration/README.md benchmark/README.md
git commit -m "bench: add default calibration grid scaffold"
```

---

### Task 8: Documentation and Release Notes

**Files:**
- Modify: `src/epykit/__init__.py`
- Modify: `src/epykit/dmc.py`
- Modify: `src/epykit/dmr.py`
- Modify: `benchmark/PROTOCOL.md`
- Modify: `benchmark/paper/report/methods_appendix.md`
- Modify: `benchmark/paper/paper.md`

**Interfaces:**
- Consumes: final accepted defaults from Tasks 2-7
- Produces: docs that distinguish default behavior, opt-in research knobs, and benchmarked claims

- [ ] **Step 1: Update package-facing docs**

In `src/epykit/__init__.py`, replace the DMC default summary with:

```python
"""DMC defaults use lr with EB dispersion and two-stage BH FDR.

The lr+ power stack remains opt-in because real-WGBS benchmarks showed FPR
inflation when neighbour combining is enabled globally.
"""
```

- [ ] **Step 2: Update DMR tuning guidance**

In `src/epykit/dmr.py`, make the tuning guidance say:

```python
For recall-sensitive discovery, start with preset="default" or
preset="permissive". The bare high-level chain_merge default now uses
min_cpgs=3 to avoid dropping short but coherent regions.
```

- [ ] **Step 3: Update benchmark docs only with generated numbers**

In `benchmark/PROTOCOL.md`, add a row for the new default profile after the calibration grid has produced real metrics:

```markdown
| epykit default tuned | `ep.tl.dmc(test="lr")` + `ep.tl.dmr(method="chain_merge")` | Default post-calibration profile; see `benchmark/data/default_calibration/`. |
```

- [ ] **Step 4: Run documentation-sensitive tests**

Run:

```bash
uv run pytest tests/test_p0_eb_default_documented.py tests/test_cli_api_parity.py -ra
```

Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add src/epykit/__init__.py src/epykit/dmc.py src/epykit/dmr.py benchmark/PROTOCOL.md benchmark/paper/report/methods_appendix.md benchmark/paper/paper.md
git commit -m "docs: document calibrated dmc dmr defaults"
```

---

## Verification Matrix

Run after all accepted tasks:

```bash
uv run pytest -m "not slow" --strict-markers -ra
uv run ruff check src/
uv run mypy src/epykit
uv run python benchmark/scripts/run_default_calibration_grid.py --quick
```

For full benchmark validation before updating claims:

```bash
uv run python benchmark/scripts/run_default_calibration_grid.py
```

Acceptance thresholds for shipping default flips:

- Small-effect DMC TPR at coverage 10 improves over the current baseline.
- Null FPR remains below the explicit guardrail chosen in `tests/test_power_calibration_defaults.py`.
- DMR chain-merge recall improves on 3-CpG and sparse-region fixtures.
- `power_stack="lr+"` remains opt-in and documented as exploratory.
- No benchmark-paper metric is updated without regenerated data.

## Self-Review

- Spec coverage: smoothing, EB dispersion, chain-merge min CpGs/pct_sig, `DF_PHI_FLOOR`, FDR method, and `welch_t` are all covered by tasks.
- Placeholder scan: no task uses `TBD` or an unbounded “add tests” instruction; each task has concrete files, code shape, and commands.
- Type consistency: all referenced public functions already exist except the optional force behavior in Task 6, which is introduced in the same task before use.
- Risk note: Task 4 and Task 5 are explicitly benchmark-gated because the current repository already documents FPR risks around smoothing-like power boosts and `DF_PHI_FLOOR`.
