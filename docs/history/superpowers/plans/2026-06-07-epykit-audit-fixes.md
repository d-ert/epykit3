# epykit 2026-06-07 audit-fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remediate the defects found in the 2026-06-07 peer-review-grade audit (`docs/review/2026-06-07-epykit-codebase-audit.md`) so epykit 1.0 ships a correct, reproducible, methods-paper-defensible release.

**Architecture:** Phased remediation in four batches by topic and review-priority (Batch 1 = permutation/empirical-FDR; Batch 2 = CLI/API parity + ingest defaults; Batch 3 = aggregate_regions + GLM consistency + welch CI; Batch 4 = read-level/export/packaging + downgraded D-bundle + calibration tests). Every fix lands as an atomic commit with a regression test that would have caught the bug. Streaming/logging/Windows-CI invariants from `CLAUDE.md` are preserved.

**Tech Stack:** polars, numpy, scipy, statsmodels, pytest (with `slow` marker), `uv` lockfile in CI.

**Source of findings:** [docs/review/2026-06-07-epykit-codebase-audit.md](../../review/2026-06-07-epykit-codebase-audit.md) (18 module/dimension reviewers + adversarial verification pass; "core-engine recovery pass" findings flagged unverified — M13 is the strongest of those).

**Prior batch:** [2026-06-06-epykit-review-fixes-design.md](2026-06-06-epykit-review-fixes-design.md). The current branch (`review-fixes-batch-2`) already landed batch-2 commits for the 2026-06-06 audit. This plan is the **third batch**, covering everything new + the items the 06-07 audit confirms still open.

---

## Non-goals (explicit)

- No `lr+` knob-default changes (would invalidate the paper's ablations; `lr+` stays a research opt-in per `CLAUDE.md`).
- No re-running of benchmark headlines in Batch 1; benchmarks re-run only if a fix changes a default that the paper claims.
- No bug-fix manifest in the manuscript (per `feedback_paper_no_bugfix_history.md`).
- Packaging URL/`[all]` fixes already landed in batch-2 are not re-done; this plan only touches what's still open.

## Guiding constraints (pre-flight)

- Windows + Linux, py3.9–3.12 must stay green (CI matrix).
- Library code emits via `logging`, never `print`; CLI owns stdout.
- Preserve O(largest-chromosome) streaming contract; don't materialize whole-genome frames.
- Every fix = one atomic commit with a test that would have caught the bug. Where the existing test encodes the bug, fix the test as part of the fix and reference it explicitly in the commit message.

## Batch ordering (audit §4 sequence)

| Batch | Theme | Findings | Posture |
|---|---|---|---|
| **1** | **Permutation/empirical-FDR correctness** | M1, M2, M3, m-perm (failed-perm denominator bias, no `merge_adjacent`/`backend` propagation, untested-NaN-as-neighbour) | **Written in full below — TDD step-by-step.** Audit names these "fix before submission". |
| **2** | **CLI/API parity + ingest defaults** | M5 (merge_strands no-op), M10 (CLI `min_cpgs` drop), D10 (CLI merge default), D11 (region q-filter), D12 (CLI `--allow-n1` help), `fdr_method` ignored in contrast path, `--unite` divergence, `tl.contrast` CLI dispersion/reference forwarding | Scope + acceptance enumerated. |
| **3** | **Aggregate_regions + GLM consistency + welch CI** | M4 (GLM DF_PHI_FLOOR), M6 (aggregate_regions stale output), M7 (overlapping BED), M9 (methylKit 1-based off-by-one), M13 (welch CI uses z but test uses t), D2 (GPU IRLS NaN mask), D8 (multi-group `df2` audit), D9 (neighbour combine NaN audit), `_per_site_variance_test` cleanup | Scope + acceptance enumerated. |
| **4** | **Read-level modules + exports + downgraded D-bundle + calibration coverage** | M11 (BAM/ASM/entropy/clock bugs — experimental gate + targeted fixes), M12 (export streaming), D1 (segment unsigned Stouffer), D3-D7, D13-D17, D18 (φ-overdispersed FPR tests for glm/welch_t), assorted Minors/Nitpicks | Scope + acceptance enumerated; some items deferred to "experimental" namespace. |

Each batch is a tier gate: full `uv run pytest -m "not slow"` green on both OSes + the batch's specific acceptance criteria + at least one new regression test per finding before the next batch begins.

---

# BATCH 1 — Permutation/empirical-FDR correctness

**Why first.** Audit §4 #1 names the empirical-FDR story the single biggest reviewer-exposure surface. Three confirmed Majors compose into one coherent fix because they all live in the `tl.dmr` / `tl.dmc` permutation harness and share callers, fixtures, and tests:

- **M1** — stratified permutation confounds groups (`tl.py:1409-1420`, consumed at `dmr.py:1445-1454`); paired/batch designs get a degenerate null.
- **M2** — `empirical_fdr=True` is silently a no-op for `chain_merge` (default), `sliding_window`, `segment` (`tl.py:1402-1442`).
- **M3** — empirical DMC null doesn't reproduce `sep_fallback`/`smoothing` from the observed run (`tl.py:920-936`) → anti-conservative.

Two related Minors fold in cleanly:
- **m-perm-1** — failed permutations counted in the `n_perm+1` denominator (`dmc.py:3086-3090`, `dmr.py:1466-1471`).
- **m-perm-2** — permutation-null tile run doesn't inherit `merge_adjacent`/`backend` from the observed run (`tl.py:1422-1442`).

## Batch-1 file structure

| File | Responsibility | Verb |
|---|---|---|
| `src/epykit/dmr.py` | `empirical_fdr_for_dmr`: within-stratum k-of-n shuffle; failed-perm denominator correction; method dispatch to non-tile callers. | Modify |
| `src/epykit/tl.py` | `tl.dmr` empirical branch: lift out of `method=="tile"`; thread `sep_fallback`/`smoothing` into permutation kwargs; propagate `merge_adjacent`/`backend`. `tl.dmc` empirical branch: thread `sep_fallback`/`smoothing`. | Modify |
| `src/epykit/dmc.py` | `empirical_fdr_for_dmc`: forward `sep_fallback`/`smoothing`/`sep_threshold`/`smoothing_span_bp`; failed-perm denominator correction. | Modify |
| `tests/test_empirical_fdr_stratified.py` | Within-stratum permutation test (k-of-n). **New file.** | Create |
| `tests/test_empirical_fdr_method_coverage.py` | Asserts `empirical_pvalue`/`empirical_qvalue` columns appear for every DMR method (or that `NotImplementedError` is raised). **New file.** | Create |
| `tests/test_empirical_fdr_sep_fallback_parity.py` | Asserts observed-run + null-run use the same `sep_fallback`/`smoothing` settings. **New file.** | Create |
| `tests/test_empirical_fdr_denominator.py` | Asserts failed permutations do not count in denominator. **New file.** | Create |

## Task 1.1 — Within-stratum k-of-n permutation (M1)

**Files:**
- Modify: `src/epykit/dmr.py:1442-1471` (`empirical_fdr_for_dmr._run_one_perm`)
- Modify: `src/epykit/tl.py:1409-1420` (strata_map construction — record original treatment count per stratum)
- Test: `tests/test_empirical_fdr_stratified.py` (create)

- [ ] **Step 1.1.1 — Write the failing test (paired-design degeneracy).**

```python
# tests/test_empirical_fdr_stratified.py
"""M1: stratified empirical FDR must permute within strata while preserving
each stratum's original treatment/control split. The pre-fix implementation
shuffled in dict order then split globally, so a paired 1T+1C design sent
ALL pairs' first elements to treatment every permutation."""
import numpy as np
import polars as pl
import pytest
from epykit.dmr import _stratified_permutation_assignment


def test_stratified_permutation_preserves_per_stratum_counts():
    # 4 paired strata, each 1 treatment + 1 control. n_treat = 4, n_ctrl = 4.
    strata_map = {
        "pair_A": ["A_T", "A_C"],
        "pair_B": ["B_T", "B_C"],
        "pair_C": ["C_T", "C_C"],
        "pair_D": ["D_T", "D_C"],
    }
    original_treatment = ["A_T", "B_T", "C_T", "D_T"]
    original_control = ["A_C", "B_C", "C_C", "D_C"]
    rng = np.random.default_rng(0)

    # Run many permutations; each must produce exactly one treatment sample
    # per stratum (k_treat = 1 per pair).
    for _ in range(200):
        perm_t, perm_c = _stratified_permutation_assignment(
            strata_map=strata_map,
            samples_treatment=original_treatment,
            samples_control=original_control,
            rng=rng,
        )
        assert len(perm_t) == 4
        assert len(perm_c) == 4
        for stratum, members in strata_map.items():
            n_in_treat = sum(1 for s in perm_t if s in members)
            assert n_in_treat == 1, (
                f"Stratum {stratum} must contribute exactly 1 sample to "
                f"treatment (its original count); got {n_in_treat}."
            )


def test_stratified_permutation_unequal_strata_preserves_counts():
    # Stratum sizes 2T+1C and 1T+2C.
    strata_map = {"S1": ["s1a", "s1b", "s1c"], "S2": ["s2a", "s2b", "s2c"]}
    original_treatment = ["s1a", "s1b", "s2a"]   # 2 from S1, 1 from S2
    original_control = ["s1c", "s2b", "s2c"]
    rng = np.random.default_rng(1)

    for _ in range(200):
        perm_t, _ = _stratified_permutation_assignment(
            strata_map=strata_map,
            samples_treatment=original_treatment,
            samples_control=original_control,
            rng=rng,
        )
        s1_in_t = sum(1 for s in perm_t if s in strata_map["S1"])
        s2_in_t = sum(1 for s in perm_t if s in strata_map["S2"])
        assert s1_in_t == 2 and s2_in_t == 1
```

- [ ] **Step 1.1.2 — Run, confirm it fails on the import.**

```powershell
uv run pytest tests/test_empirical_fdr_stratified.py -v
```
Expected: `ImportError: cannot import name '_stratified_permutation_assignment' from 'epykit.dmr'`.

- [ ] **Step 1.1.3 — Extract the new helper in `dmr.py` and use it from `_run_one_perm`.**

Add at module scope in `src/epykit/dmr.py` (above `empirical_fdr_for_dmr`):

```python
def _stratified_permutation_assignment(
    *,
    strata_map: dict[str, list[str]],
    samples_treatment: list[str],
    samples_control: list[str],
    rng: np.random.Generator,
) -> tuple[list[str], list[str]]:
    """Per-stratum k-of-n permutation.

    For each stratum, randomly select k samples as treatment, where k equals
    the original number of treatment samples in that stratum. The remaining
    samples become control. Preserves per-stratum group sizes -- the
    invariant pre-fix code violated by shuffling within strata and then
    splitting globally.
    """
    treat_set = set(samples_treatment)
    perm_treat: list[str] = []
    perm_ctrl: list[str] = []
    for stratum_samples in strata_map.values():
        k_treat = sum(1 for s in stratum_samples if s in treat_set)
        shuffled = list(rng.permutation(stratum_samples))
        perm_treat.extend(shuffled[:k_treat])
        perm_ctrl.extend(shuffled[k_treat:])
    return perm_treat, perm_ctrl
```

Then replace the body of the `if empirical_strata is not None:` branch inside `_run_one_perm` (currently `dmr.py:1445-1454`):

```python
        if empirical_strata is not None:
            perm_treat, perm_ctrl = _stratified_permutation_assignment(
                strata_map=empirical_strata,
                samples_treatment=samples_treatment,
                samples_control=samples_control,
                rng=local_rng,
            )
        else:
            shuffled = pool.copy()
            local_rng.shuffle(shuffled)
            perm_treat = shuffled[:n_treat]
            perm_ctrl = shuffled[n_treat:]
```

Remove the now-dead `pool = list(samples_treatment) + list(samples_control)` line at `dmr.py:1439` only if it is no longer referenced by the unstratified branch (keep it — the `else` branch above still uses `pool` and `n_treat` from the enclosing scope; verify by reading the function before deleting).

- [ ] **Step 1.1.4 — Run the new test; confirm pass.**

```powershell
uv run pytest tests/test_empirical_fdr_stratified.py -v
```
Expected: PASS.

- [ ] **Step 1.1.5 — Mirror the strata_map construction so `tl.dmr` passes original treatment-count information through.**

`tl.py:1411-1420` already groups by `empirical_strata` column. No semantic change needed — the *helper* derives `k_treat` from `samples_treatment` ∩ stratum members, which is the post-fix contract. Verify by reading: the strata_map values are sample-ids; the helper consults `samples_treatment` for membership. No edit required here. **Document this in a comment on the helper definition** so it's not silently broken later by someone refactoring `tl.py` to drop the original treatment set.

- [ ] **Step 1.1.6 — Run the full empirical-FDR test suite.**

```powershell
uv run pytest tests/ -k "empirical" -v
```
Expected: all green (the existing smoke tests pass because the contract — column exists, treats/ctrls non-empty — still holds).

- [ ] **Step 1.1.7 — Commit.**

```powershell
git add src/epykit/dmr.py tests/test_empirical_fdr_stratified.py
git commit -m "fix(dmr): within-stratum k-of-n permutation for empirical FDR (M1)

Pre-fix: empirical_fdr_for_dmr's stratified branch shuffled within each
stratum then split globally, sending the first ceil(n/2) of the dict-
ordered pool to treatment every permutation. For paired 1T+1C designs
this produced a degenerate, maximally-confounded null and silently
invalidated empirical_pvalue/empirical_qvalue for exactly the batch/
paired designs the feature exists to serve.

Now: for each stratum, randomly select k samples as treatment where k
is that stratum's original treatment count; pool into perm_treat/
perm_ctrl; never split globally.

Audit: docs/review/2026-06-07-epykit-codebase-audit.md (M1)."
```

## Task 1.2 — Failed-permutation denominator correction (m-perm-1)

**Files:**
- Modify: `src/epykit/dmr.py:1466-1496` (count non-failed perms; use that count in the denominator)
- Modify: `src/epykit/dmc.py:3086-3090` (same fix on the DMC side)
- Test: `tests/test_empirical_fdr_denominator.py` (create)

- [ ] **Step 1.2.1 — Write the failing test.**

```python
# tests/test_empirical_fdr_denominator.py
"""m-perm-1: failed permutations must NOT count in the empirical-p
denominator. Pre-fix: denominator was always n_perm+1, biasing empirical p
downward by the failure rate (which can be non-trivial under small-n DMR
calling where some permuted splits produce zero candidate regions)."""
import numpy as np
import polars as pl
from epykit.dmr import _empirical_pvalues_from_null_pool


def test_failed_permutations_excluded_from_denominator():
    observed = np.array([1e-6, 1e-3, 0.05])
    null_pool = np.array([1e-4, 1e-2, 0.5])  # only 3 successful perms
    n_perm_attempted = 10
    n_perm_successful = 3

    emp = _empirical_pvalues_from_null_pool(
        observed_pvalues=observed,
        null_pvalues_pool=null_pool,
        n_perm_used=n_perm_successful,
    )
    # For observed[0]=1e-6: 0 nulls <= it -> emp = (0+1)/(3+1) = 0.25.
    # The pre-fix code used (10+1) in the denominator -> 0.0909, biased low.
    assert abs(emp[0] - 0.25) < 1e-9, (
        f"Denominator must be n_perm_successful+1=4, got {1.0/emp[0] - 1}"
    )
```

- [ ] **Step 1.2.2 — Run, confirm it fails on the import.**

```powershell
uv run pytest tests/test_empirical_fdr_denominator.py -v
```
Expected: ImportError on `_empirical_pvalues_from_null_pool`.

- [ ] **Step 1.2.3 — Extract the helper used by both `dmr.py` and `dmc.py`.**

If the current code inlines this computation, extract a module-private helper into `src/epykit/dmr.py`:

```python
def _empirical_pvalues_from_null_pool(
    *,
    observed_pvalues: np.ndarray,
    null_pvalues_pool: np.ndarray,
    n_perm_used: int,
) -> np.ndarray:
    """Pooled-null Westfall-Young empirical p-values.

    Parameters
    ----------
    n_perm_used : int
        Number of permutations that produced at least one usable null
        p-value. Failed permutations (zero null regions, exception in the
        engine) MUST be excluded. The denominator is ``n_perm_used + 1``,
        not ``n_perm_requested + 1``.
    """
    if null_pvalues_pool.size == 0 or n_perm_used <= 0:
        return np.ones_like(observed_pvalues, dtype=np.float64)
    sorted_null = np.sort(null_pvalues_pool)
    counts = np.searchsorted(sorted_null, observed_pvalues, side="right")
    # Standard Westfall-Young adjustment: (count + 1) / (n + 1).
    return (counts.astype(np.float64) + 1.0) / (n_perm_used + 1.0)
```

Find the current denominator at `dmr.py:1488-1496` (or wherever the empirical p-value is computed; grep `n_perm + 1`) and replace with a call to the helper. Compute `n_perm_used = sum(1 for arr in null_pvals_list if arr.size > 0)`.

Repeat in `dmc.py:3086-3090` (the analogous block in `empirical_fdr_for_dmc`).

- [ ] **Step 1.2.4 — Run the test and the full empirical suite.**

```powershell
uv run pytest tests/test_empirical_fdr_denominator.py tests/ -k "empirical" -v
```
Expected: all green.

- [ ] **Step 1.2.5 — Commit.**

```powershell
git add src/epykit/dmr.py src/epykit/dmc.py tests/test_empirical_fdr_denominator.py
git commit -m "fix(empirical): exclude failed permutations from denominator (m-perm-1)

Pre-fix: failed permutations (zero null regions, engine exception) were
counted in the n_perm+1 denominator, biasing empirical p downward by the
failure rate. Now: empirical p = (n_null<=obs + 1) / (n_perm_used + 1),
where n_perm_used counts only permutations that produced >=1 null value.

Audit: docs/review/2026-06-07-epykit-codebase-audit.md (Minors / perm-1)."
```

## Task 1.3 — Forward `sep_fallback`/`smoothing` into the DMC permutation null (M3)

**Files:**
- Modify: `src/epykit/tl.py:920-936` (forward `sep_fallback`, `sep_threshold`, `smoothing`, `smoothing_span_bp` into `empirical_fdr_for_dmc`)
- Modify: `src/epykit/dmc.py` (signature of `empirical_fdr_for_dmc`; thread through to its inner permutation runner)
- Test: `tests/test_empirical_fdr_sep_fallback_parity.py` (create)

- [ ] **Step 1.3.1 — Locate `empirical_fdr_for_dmc` and its inner permutation runner.**

```powershell
uv run python -c "import epykit.dmc as d, inspect; print(inspect.getsourcefile(d.empirical_fdr_for_dmc)); print(inspect.signature(d.empirical_fdr_for_dmc))"
```
Record the current signature; the new kwargs need to be added without breaking existing callers.

- [ ] **Step 1.3.2 — Write the failing test.**

```python
# tests/test_empirical_fdr_sep_fallback_parity.py
"""M3: empirical DMC FDR must run permutations with the SAME sep_fallback
and smoothing settings as the observed run; otherwise the Westfall-Young
statistic compares deflated observed p-values against an un-deflated null
pool, producing anti-conservative empirical p."""
import polars as pl
from unittest.mock import patch
from epykit.dmc import empirical_fdr_for_dmc


def test_sep_fallback_propagates_to_permutations(tmp_path):
    # Create a tiny synthetic store with one chromosome.
    from tests.fixtures.synth import build_two_group_store
    store_path = build_two_group_store(
        tmp_path=tmp_path, n_treat=3, n_ctrl=3, n_sites=200,
        chrom="chr1", seed=0,
    )
    observed = pl.DataFrame({
        "chrom": ["chr1"] * 5,
        "pos": [100, 200, 300, 400, 500],
        "pvalue": [1e-3, 1e-2, 0.1, 0.3, 0.5],
        "meth_diff": [0.4, 0.3, 0.1, 0.05, 0.0],
    })

    captured_kwargs: list[dict] = []
    original = empirical_fdr_for_dmc.__wrapped__ if hasattr(
        empirical_fdr_for_dmc, "__wrapped__"
    ) else None

    # Patch the per-permutation worker to capture its kwargs.
    from epykit.dmc import _run_one_dmc_perm
    with patch("epykit.dmc._run_one_dmc_perm", autospec=True) as mock_runner:
        mock_runner.return_value = pl.DataFrame({"pvalue": [0.5]})
        empirical_fdr_for_dmc(
            methylstore_path=store_path,
            samples_treatment=["sample_t0", "sample_t1", "sample_t2"],
            samples_control=["sample_c0", "sample_c1", "sample_c2"],
            observed_dmc=observed,
            n_perm=2, seed=0, n_jobs=1, test="lr",
            sep_fallback=True, sep_threshold=0.05,
            smoothing=True, smoothing_span_bp=500,
        )
        # Every permutation must have been invoked WITH the same sep/smoothing.
        assert mock_runner.call_count == 2
        for call in mock_runner.call_args_list:
            kwargs = call.kwargs
            assert kwargs.get("sep_fallback") is True, (
                f"sep_fallback not forwarded: {kwargs}"
            )
            assert kwargs.get("smoothing") is True
            assert kwargs.get("smoothing_span_bp") == 500
```

(If `_run_one_dmc_perm` doesn't exist as a separate symbol, the test should patch the inner `process_chromosomes_dmc` call inside the perm worker and assert on its kwargs.)

- [ ] **Step 1.3.3 — Run, confirm failure.**

```powershell
uv run pytest tests/test_empirical_fdr_sep_fallback_parity.py -v
```
Expected: FAIL — `sep_fallback` not in `kwargs` (the function signature doesn't accept it yet).

- [ ] **Step 1.3.4 — Add `sep_fallback`, `sep_threshold`, `smoothing`, `smoothing_span_bp` to `empirical_fdr_for_dmc`.**

In `src/epykit/dmc.py`, add to the signature (after `reference`):

```python
def empirical_fdr_for_dmc(
    *,
    methylstore_path: str,
    samples_treatment: list[str],
    samples_control: list[str],
    observed_dmc: pl.DataFrame,
    n_perm: int = 100,
    seed: int = 0,
    n_jobs: int = 1,
    test: str = "lr",
    chromosomes: Optional[Sequence[str]] = None,
    unite: str = "intersect",
    min_samples_treatment: int = 2,
    min_samples_control: int = 2,
    dispersion: str = "eb",
    reference: str = "adaptive",
    sep_fallback: bool = False,
    sep_threshold: float = 0.05,
    smoothing: bool = False,
    smoothing_span_bp: int = 1000,
) -> pl.DataFrame:
```

Forward all four into the per-permutation `process_chromosomes_dmc` call (find the call inside the perm loop / worker — grep `process_chromosomes_dmc` inside `empirical_fdr_for_dmc`).

- [ ] **Step 1.3.5 — Forward from `tl.dmc` at `tl.py:920-936`.**

Replace the existing call:

```python
        if empirical_fdr and len(result) > 0:
            result = empirical_fdr_for_dmc(
                methylstore_path=_dmc_store,
                samples_treatment=md.treatment_ids,
                samples_control=md.control_ids,
                observed_dmc=result,
                n_perm=n_perm,
                seed=perm_seed,
                n_jobs=perm_n_jobs,
                test=selected_test,
                chromosomes=chromosomes,
                unite=unite,
                min_samples_treatment=min_samples_treatment,
                min_samples_control=min_samples_control,
                dispersion=dispersion,
                reference=reference,
                sep_fallback=sep_fallback,
                sep_threshold=sep_threshold,
                smoothing=use_smoothed,
                smoothing_span_bp=smoothing_span_bp,
            )
```

(Verify `use_smoothed`, `smoothing_span_bp`, `sep_fallback`, `sep_threshold` are in scope at this point in `tl.dmc`; if `smoothing_span_bp` lives under a different local name, read its surrounding block to find the canonical local.)

- [ ] **Step 1.3.6 — Run the failing test → pass.**

```powershell
uv run pytest tests/test_empirical_fdr_sep_fallback_parity.py -v
```
Expected: PASS.

- [ ] **Step 1.3.7 — Run the full DMC suite.**

```powershell
uv run pytest tests/test_dmc.py tests/test_dmc_multigroup.py tests/ -k "empirical or dmc" -v
```
Expected: all green.

- [ ] **Step 1.3.8 — Commit.**

```powershell
git add src/epykit/dmc.py src/epykit/tl.py tests/test_empirical_fdr_sep_fallback_parity.py
git commit -m "fix(empirical): propagate sep_fallback/smoothing into DMC permutations (M3)

Pre-fix: tl.dmc's empirical_fdr path forwarded only test/unite/min_samples/
dispersion/reference to the permutation null, while the observed run also
applied sep_fallback and smoothing. The Westfall-Young statistic then
compared deflated observed p-values against an un-deflated null pool ->
anti-conservative empirical_pvalue (silent, reachable via power_stack='lr+'
or smoothing=True).

Now: sep_fallback/sep_threshold/smoothing/smoothing_span_bp are forwarded
into empirical_fdr_for_dmc and through to each permutation's
process_chromosomes_dmc call. Observed and null runs share the engine
configuration.

Audit: docs/review/2026-06-07-epykit-codebase-audit.md (M3)."
```

## Task 1.4 — Lift empirical-FDR out of the `method=="tile"` branch (M2)

**Decision:** the audit gives a binary — implement empirical FDR for `chain_merge`/`sliding_window`/`segment`, or raise `NotImplementedError`. **We raise** in Batch 1 (one-line correctness fix; full implementation deferred to a follow-up because each non-tile caller has a different region-definition story and the permutation harness for them is non-trivial). This eliminates the silent no-op + suppressed calibration warning, which is the actual reviewer-exposure surface. Add a `# TODO(batch-4-followup)` referencing this plan for the implement-for-real path.

**Files:**
- Modify: `src/epykit/tl.py:1402-1466` (move the `empirical_fdr` check before the `method=="tile"` branch; raise for non-tile; do not suppress the calibration warning)
- Test: `tests/test_empirical_fdr_method_coverage.py` (create)

- [ ] **Step 1.4.1 — Write the failing test.**

```python
# tests/test_empirical_fdr_method_coverage.py
"""M2: empirical_fdr=True was silently a no-op for chain_merge (the
default), sliding_window, and segment. The calibration-warning note was
also suppressed in those branches, leaving users with no signal that
combined_qvalue is anti-conservative. Batch-1 contract: empirical_fdr=True
must produce columns (tile) OR raise NotImplementedError (others); never
silently no-op."""
import polars as pl
import pytest
from epykit.methyldata import MethylData


def test_chain_merge_empirical_fdr_raises_notimplemented(small_md_with_dmc):
    md = small_md_with_dmc
    from epykit import tl
    with pytest.raises(NotImplementedError, match="empirical_fdr.*tile"):
        tl.dmr(md, method="chain_merge", empirical_fdr=True, n_perm=10)


def test_sliding_window_empirical_fdr_raises_notimplemented(small_md_with_dmc):
    md = small_md_with_dmc
    from epykit import tl
    with pytest.raises(NotImplementedError, match="empirical_fdr.*tile"):
        tl.dmr(md, method="sliding_window", empirical_fdr=True, n_perm=10)


def test_segment_empirical_fdr_raises_notimplemented(small_md_with_dmc):
    md = small_md_with_dmc
    from epykit import tl
    with pytest.raises(NotImplementedError, match="empirical_fdr.*tile"):
        tl.dmr(md, method="segment", empirical_fdr=True, n_perm=10)


def test_tile_empirical_fdr_still_works(small_md_with_dmc):
    md = small_md_with_dmc
    from epykit import tl
    tl.dmr(md, method="tile", empirical_fdr=True, n_perm=5, perm_seed=0)
    dmr = md.uns["dmr"]
    assert "empirical_pvalue" in dmr.columns
    assert "empirical_qvalue" in dmr.columns


def test_calibration_warning_not_suppressed_for_empirical_fdr_request(
    small_md_with_dmc, caplog,
):
    """Pre-fix bug: the 'your q-value may be anti-conservative' note fired
    only when `empirical_fdr` was False, so a user passing empirical_fdr=True
    on chain_merge got (a) no empirical columns, (b) no error, AND (c) the
    warning suppressed. Now an explicit NotImplementedError surfaces the
    issue; the warning is no longer the user's only signal so its
    suppression-by-flag bug becomes moot."""
    # If the NotImplementedError is raised before reaching the warning,
    # there's nothing to suppress. This test pins the contract.
    md = small_md_with_dmc
    from epykit import tl
    with pytest.raises(NotImplementedError):
        tl.dmr(md, method="chain_merge", empirical_fdr=True)


@pytest.fixture
def small_md_with_dmc(tmp_path):
    """Minimal MethylData with a DMC result populated -- enough for tl.dmr to
    dispatch each method. Uses the existing synth helpers."""
    from tests.fixtures.synth import build_two_group_store
    from epykit.io import read_bismark  # noqa: F401  (just ensures import path)
    store = build_two_group_store(
        tmp_path=tmp_path, n_treat=3, n_ctrl=3, n_sites=200, chrom="chr1", seed=0,
    )
    md = MethylData.from_methylstore(store)
    from epykit import tl as _tl
    _tl.dmc(md, test="lr")
    return md
```

(If `build_two_group_store` doesn't exist in `tests/fixtures/synth.py`, locate the existing two-group fixture — `tests/test_empirical_fdr.py` already uses one — and reuse it. The fixture's exact name is implementation detail; the test contract is what matters.)

- [ ] **Step 1.4.2 — Run, confirm failure (silent no-op rather than NotImplementedError).**

```powershell
uv run pytest tests/test_empirical_fdr_method_coverage.py -v
```
Expected: most tests FAIL — `chain_merge`/`sliding_window`/`segment` silently return with no `empirical_pvalue` column instead of raising.

- [ ] **Step 1.4.3 — Move the empirical-FDR check outside the tile branch.**

In `src/epykit/tl.py`, restructure `tl.dmr` so the `empirical_fdr` gate is **before** the method dispatch:

```python
    if empirical_fdr and method != "tile":
        raise NotImplementedError(
            f"empirical_fdr=True is currently implemented only for "
            f"method='tile'. Got method={method!r}. Use method='tile' or "
            f"omit empirical_fdr=True. (Follow-up: implement permutation "
            f"FDR for chain_merge/sliding_window/segment -- tracked in "
            f"docs/history/superpowers/plans/2026-06-07-epykit-audit-fixes.md "
            f"Batch-4-followup.)"
        )
```

Place this right after the method validation and before the `if method == "tile":` block. Verify the existing tile path still threads `empirical_fdr` through unchanged (it does — that block remains).

Then audit any `if not empirical_fdr:` calibration-warning suppression and **remove the `not empirical_fdr` condition**. The warning should fire on `method in {"chain_merge","sliding_window","segment"}` regardless of the `empirical_fdr` flag (since empirical FDR doesn't actually apply there).

Grep for the suppression: `grep -n "not empirical_fdr" src/epykit/tl.py` — verify and fix each occurrence.

- [ ] **Step 1.4.4 — Update `tl.dmr` docstring** so it no longer claims `sliding_window` supports empirical FDR.

Find the docstring's `empirical_fdr` parameter description (likely near the function signature in `tl.py`). Replace any "supported for chain_merge / sliding_window" wording with: `"empirical_fdr : currently supported only for method='tile'. Raises NotImplementedError otherwise."`

Also update `CLAUDE.md` line 50 (`"All [DMR callers] support optional permutation empirical FDR"`) → `"Permutation empirical FDR is implemented for method='tile'; other callers raise NotImplementedError pending Batch-4 follow-up."`

- [ ] **Step 1.4.5 — Run failing test → pass; run full suite.**

```powershell
uv run pytest tests/test_empirical_fdr_method_coverage.py -v
uv run pytest tests/ -k "dmr or empirical" -v
```
Expected: all green.

- [ ] **Step 1.4.6 — Commit.**

```powershell
git add src/epykit/tl.py CLAUDE.md tests/test_empirical_fdr_method_coverage.py
git commit -m "fix(tl): raise NotImplementedError for empirical_fdr on non-tile DMR (M2)

Pre-fix: empirical_fdr=True was silently a no-op for chain_merge (default),
sliding_window, and segment -- no error, no warning, no
empirical_pvalue/empirical_qvalue columns. The calibration-warning note
that says 'combined_qvalue is anti-conservative' was suppressed when
empirical_fdr=True, so a user calling tl.dmr(md, empirical_fdr=True) on
the DEFAULT caller got nothing AND lost the only signal of trouble.

Batch-1 posture: raise NotImplementedError for non-tile callers (silent
no-op is a worse failure mode than an explicit raise). Implementing real
per-method permutation FDR for chain_merge/sliding_window/segment is
deferred to Batch-4 (each caller's region-definition needs its own
permutation harness).

The calibration-warning suppression is also removed: the note now fires
on chain_merge/sliding_window/segment regardless of empirical_fdr.

Updated tl.dmr docstring and CLAUDE.md to match.

Audit: docs/review/2026-06-07-epykit-codebase-audit.md (M2)."
```

## Task 1.5 — Propagate `merge_adjacent`/`backend` into tile permutations (m-perm-2)

**Files:**
- Modify: `src/epykit/tl.py:1422-1442` (forward `merge_adjacent` and `backend` to `empirical_fdr_for_dmr`)
- Modify: `src/epykit/dmr.py` (`empirical_fdr_for_dmr` signature + per-perm `call_dmr_tile_based` call)
- Test: extend `tests/test_empirical_fdr_method_coverage.py` with one assertion

- [ ] **Step 1.5.1 — Extend the existing test.**

Append to `tests/test_empirical_fdr_method_coverage.py`:

```python
def test_tile_empirical_fdr_propagates_merge_adjacent(small_md_with_dmc):
    """m-perm-2: permutation null must use the same merge_adjacent and
    backend as the observed run; otherwise the observed vs null
    distributions are computed under different region definitions."""
    from unittest.mock import patch
    md = small_md_with_dmc
    with patch("epykit.dmr.call_dmr_tile_based") as mock:
        # Make the mock return a non-empty DMR frame so the empirical step runs.
        mock.return_value = pl.DataFrame({
            "chrom": ["chr1"], "start": [100], "end": [200],
            "pvalue": [0.01], "meth_diff": [0.3],
        })
        from epykit import tl
        tl.dmr(
            md, method="tile", empirical_fdr=True, n_perm=2, perm_seed=0,
            merge_adjacent=False, backend="sequential",
        )
        # The observed call + at least 2 null calls should all carry
        # merge_adjacent=False, backend='sequential'.
        for call in mock.call_args_list:
            assert call.kwargs.get("merge_adjacent") is False
            assert call.kwargs.get("backend") == "sequential"
```

- [ ] **Step 1.5.2 — Add `merge_adjacent`/`backend` to `empirical_fdr_for_dmr`'s signature.**

In `src/epykit/dmr.py`:

```python
def empirical_fdr_for_dmr(
    *,
    # ... existing kwargs ...
    merge_adjacent: bool = True,
    backend: str = "sequential",
) -> pl.DataFrame:
```

Thread both into the `call_dmr_tile_based(...)` invocation inside `_run_one_perm`.

In `src/epykit/tl.py:1422-1442`, forward `merge_adjacent=merge_adjacent, backend=backend` in the `empirical_fdr_for_dmr(...)` call.

- [ ] **Step 1.5.3 — Run, confirm pass; commit.**

```powershell
uv run pytest tests/test_empirical_fdr_method_coverage.py -v
```

```powershell
git add src/epykit/dmr.py src/epykit/tl.py tests/test_empirical_fdr_method_coverage.py
git commit -m "fix(empirical): propagate merge_adjacent/backend into tile permutations (m-perm-2)

Observed tile DMR run used merge_adjacent/backend; permutation nulls
silently fell back to defaults, so observed and null distributions were
computed under different region definitions. Now both share configuration.

Audit: docs/review/2026-06-07-epykit-codebase-audit.md (Minors / perm-2)."
```

## Task 1.6 — Batch-1 gate (full suite + Windows-equivalence check)

- [ ] **Step 1.6.1 — Full `not slow` suite.**

```powershell
uv run pytest -m "not slow" --strict-markers -ra
```
Expected: all green.

- [ ] **Step 1.6.2 — Slow tier locally (gates the calibration test landed in batch 2).**

```powershell
uv run pytest -m slow -ra
```
Expected: all green.

- [ ] **Step 1.6.3 — Update batch tracker.**

Append to the bottom of this plan file under a new "Batch 1 — DONE" section a one-liner per merged commit with its hash. The next batch starts only after this gate is green on CI for both OSes.

---

# BATCH 2 — CLI/API parity + ingest defaults (enumerated)

**Why second.** Audit §4 #3 + #4. Several CLI commands silently use different defaults than `tl.*`. Reproducing the paper via CLI gives a different DMR set than via API. `merge_strands` default behaviour on the namesake input is the single most likely "disqualifying" finding a methylation reviewer will catch.

## Findings in this batch

| ID | Location | Fix |
|---|---|---|
| **M5** | `convert.py:542-552`; `io.py` `read_bismark` signature; `read_bismark` docstring | Merge CpG dyads from the +/- position offset (pos N with N+1) without needing a reference — implement a strand-free pair-merge helper; route to it when `reference_fasta is None` and `merge_strands=True`. Document `reference_fasta` as the *better* (strand-aware) path in the docstring; keep the warning for fully-unmerge-able cases. Add regression test: two-strand `.cov` → `read_bismark(..., merge_strands=True)` (no ref) → tested-site count halved, coverages summed pairwise. |
| **M10** | `cli.py:308-317, 846` (chain_merge branch); engine default `dmr.py:684`; `tl.py:1210, 1582` | Add `min_cpgs` forwarding to `_cmd_dmr`'s chain_merge branch using a *sentinel* (`None` means "let the engine use its preset-aware default"); explicit `--min-cpgs N` propagates. Also fix `tl.dmr`'s hardcoded `5` to pass `None` when the user didn't set it so engine presets keep working. Regression test: `epykit dmr ... --min-cpgs 10 --method chain_merge` vs unset → DMR-count difference observable. |
| **D10** | `cli.py:638, 142` | Resolve the `--merge-strands` `None` sentinel to `True` (match API default). |
| **D11** | `cli.py:293-414` | Add `--min-mean-qvalue` to `epykit dmr`; default `0.05` for chain_merge/sliding_window paths so CLI matches `tl.dmr`'s post-filter. Tile path already filters at α. |
| **D12** | `cli.py:169-174, 709-716` | If `--allow-n1` and `min(n)<2` and `--test` is unset or `lr`, resolve to `fisher` (the engine that actually has the n=1 fallback) with a `logger.warning`. Update `--allow-n1` help text accordingly. |
| **Minor — contrast path** | `tl.py:1163-1165`; `cli.py:213-221` | Forward `fdr_method`, `dispersion`, `reference` through contrast/GLM DMC path; stop hardcoding `fdr_bh`. |
| **Minor — `--unite`** | `cli.py` | Resolve `--unite` to match `tl.dmc`'s default (currently the CLI defaults to `intersect` while `tl.dmc` derives union from unset `md.uns["unite"]`). |

## Batch-2 acceptance

- New `tests/test_cli_api_parity.py` (extension of the existing batch-2 file if present) runs both `epykit dmr --method chain_merge` and `tl.dmr(method='chain_merge')` on the same fixture with the same args; asserts identical DMR sets row-for-row.
- New `tests/test_merge_strands_no_reference.py` ingests a synthetic two-strand `.cov` with `merge_strands=True, reference_fasta=None`; asserts `pos` count is halved, per-pair coverages summed.
- `epykit dmr --help` lists `--min-mean-qvalue`, `--min-cpgs` (already there) forwards to chain_merge.
- README/CLAUDE.md updated: `read_bismark` defaults documented; `merge_strands` is genuinely effective without a reference.
- Windows + Linux `not slow` suites green.

---

# BATCH 3 — Aggregate_regions + GLM/CI consistency (enumerated)

**Why third.** Audit §4 #2 (calibration) + correctness clean-ups that surface when reviewers compute by hand.

## Findings in this batch

| ID | Location | Fix |
|---|---|---|
| **M4** | `_glm.py:739-744, 928-940`; callers `dmc.py:1615-1617, 1704-1707` | Apply `np.maximum(df_resid, DF_PHI_FLOOR)` inside the GLM `reference_pvalues`/`wald_test`, matching the `lr` path. Document the floor's rationale in one central comment (and fix the contradictory `dmc.py:211` docstring claiming "FPR inflation" — the math says the F tail is *more* conservative). Test: GLM p-value vs `lr` p-value at the same statistic across a range of `df_phi ∈ {2, 4, 8, 50, 200}` — assert agreement to within a documented tolerance for `df_phi ≥ floor`. |
| **M6** | `pp.py:313-314` | `rmtree` the regions store before re-writing (mirror `filter_sites` and `normalize_coverage_store`). Test: run `aggregate_regions` twice with different BEDs; assert downstream `tl.dmc` doesn't see stale chrom partitions. |
| **M7** | `pp.py:404-439` (`_assign_cpgs_to_regions`) | Replace `searchsorted+filter` with `bioframe.overlap` (already a dep); each CpG contributes to every overlapping region. Test: nested BED (outer region containing an inner region) — assert CpGs inside the inner are aggregated into both regions; assert overlapping non-nested regions both receive their shared CpGs. |
| **M9** | `methylkit_io.py:70-80` | Use `pos + 1` for both `base` and `chrBase`; add a coordinate-correctness regression test that writes a known 0-based store row and asserts the exported `base` is 1-based. |
| **M13** | `_glm.py:1017`; caller plumbing in `dmc.py:1139-area` | Add a `dof` kwarg to `welch_meth_diff_ci`; use `sp_stats.t.isf(alpha/2, dof)` instead of `norm.isf`. Plumb the Satterthwaite `dof` from the welch_t test code into the CI call. Test: at n=3 per group, assert CI width ≈ `t_{0.025,4}/1.96` × current-z-CI-width (≈1.418×). |
| **D2** | `_glm_gpu.py:220` vs `_glm.py:491-492` | Add NaN-masking of non-converged sites to the GPU IRLS path. New CPU/GPU parity test that includes deliberately non-converging sites. (Gated on `[gpu]` extra; CI may xfail if no GPU.) |
| **D8** | `dmc.py:1704-1740` | Emit the `df_phi` actually used by the multi-group F-test in a new audit column; deprecate `df2` for that engine. Document the adaptive F/χ² switch in the engine's docstring. |
| **D9** | `dmc.py:2540-2547, 2570-2584` | In neighbour Stouffer combine, exclude NaN-p neighbours from the mask (drop the `z=0` line). Update the `pvalue_combined_n_neighbours` audit count to reflect contributing-not-counted neighbours. Test: a window with 1 valid + 4 NaN-p neighbours should now report `n_neighbours=1`. |
| **Minor — dead-code** | `dvc.py:58, 185` (`_per_site_variance_test`/dead `alpha`) | Remove dead path; or, if intended, wire `alpha` through. |
| **Minor — `_solve_weighted_lsq` docstring** | `_glm.py:549` | Correct docstring: ridge-regularizes (does not return NaN on singular). |

## Batch-3 acceptance

- All `tl.dmc(test="glm")` p-values within ~1e-12 of `tl.dmc(test="lr")` p-values when the test statistic is identical and `df_phi ≥ DF_PHI_FLOOR`.
- `aggregate_regions` re-run with stricter `min_cpgs_per_region` or different BED produces a clean store; integration test passes.
- methylKit cross-validation round-trip: epykit `tl.dmc` → `to_methylkit_tabix` → methylKit `methRead` → coordinates align with annotation (manual smoke or skipped on CI; coord regression test in pytest).
- `welch_meth_diff_ci` test/CI parity: a per-site Welch p > α now never coincides with a CI excluding 0 at the same α (modulo numerical tolerance).
- Windows + Linux `not slow` suites green.

---

# BATCH 4 — Read-level + exports + downgraded D-bundle + calibration tests (enumerated)

**Why fourth.** Audit §4 #5 (read-level paper-surface), §4 #2 (calibration coverage), plus the D-bundle (real but narrower bugs) and Minors that are best fixed together.

## Major: M11 — read-level modules (BAM/ASM/entropy/clocks)

**Decision posture.** These modules are paper-surface but research-grade. We take a two-step approach:

1. **Gate them behind `epykit.experimental.*`** with `DeprecationWarning` from the current public names (one release window). This matches the existing design doc's M-PKG5 stance and stops new users from depending on shaky code.
2. **Fix the silent-wrong-result bugs in place** (reverse-strand MM/ML; ASM bisulfite-SNV confound; clock 1-based manifest) and **NotImplementedError or doc-as-experimental** the won't-scale bugs (full-BAM materialization in `read_methylation_calls`, O(n_windows × n_reads) entropy scan) so they don't silently OOM on real data.

| Sub-finding | Location | Fix |
|---|---|---|
| MM/ML reverse-strand miscall | `bam_io.py:236-301` | Replace hand-rolled MM/ML parser with `read.modified_bases` (pysam ≥0.22). Regression test on a 2-read BAM (one forward, one reverse). |
| ASM bisulfite-SNV confound | `asm.py:134-176` | At anchor selection, drop C/T and A/G transition SNVs by phasing strand. Document the restriction; flag tested haplotypes. Regression test with a deliberately conversion-confounded SNV. |
| `read_methylation_calls` full-BAM materialization | `bam_io.py:118-141`; callers `asm.py:83-86`, `entropy.py:93-95` | Iterate `regions` per chromosome by default; raise if `regions=None` and the BAM is large (>X reads). Stream-per-chromosome. |
| Entropy O(n×r) scan | `entropy.py:146-166` | Index reads by covered CpGs; raise `NotImplementedError("genome-scale entropy not yet supported; pass a regions list")` if the input would otherwise OOM. |
| Clock manifest coordinate base | `clocks.py:73-77, 205-218` | New `coordinate_base: Literal["auto","one_based","zero_based"]` param mirroring the convert-side fix from batch-2 (2026-06-06 design); default `"auto"` with detection by overlap fraction. Warn when resolved-fraction < 50% even with `impute_missing=True`. |
| MethylDackel context filter ignored | `bam_io.py:133-136` | CpG-filter explicitly or raise `ValueError("MethylDackel mode currently emits all contexts; pass context='all'")`. Downgraded from Major → Minor; do the simple raise. |

## Major-leaning: D1 (segment unsigned Stouffer)

`src/epykit/dmr_segment.py:31-39`. Reuse `_stouffer_combine_signed` from `dmr.py`. Regression test: a region with all-hyper sites must produce p ≈ that of a one-sided combine, not the unsigned two-sided shrinkage that pre-fix code emits.

## Major-leaning: M12 (export streaming)

`src/epykit/export.py:69, 100-116, 179-196`. Iterate `chrom=*` partitions; write one chromosome at a time; never materialize a whole-sample list. Test: `to_bedgraph` peak RSS on a synthetic whole-genome store stays below a documented bound (use `tracemalloc` or `resource.getrusage`).

## D-bundle (downgraded; quick wins)

| ID | Fix |
|---|---|
| D3 | `qc.sex_check`: when diptest absent, fall back to fixed threshold + `logger.warning`. Test with a synthetic single-sex cohort. |
| D4 | `qc.sample_correlation`: reservoir-sample during accumulation; honor `max_sites` before building the n_samples × n_sites matrix. |
| D5 | `convert.convert_sample`: avoid the eager `df = lf.collect()` on the no-reference path. |
| D6 | `nfcore_qc.py:24-31`: replace regexes with MultiQC's two patterns; add a real-report fixture. |
| D7 | `methylkit_io.py:59, 121-126`: switch to `pysam.tabix_compress` + `line_skip=1`, or drop the "tabix-indexed" claim from the docstring. Detect failure and surface it (don't `except: pass`). |
| D13 | `pl/_compute.py:638`: subtract `pos.min()` from chromosome tick positions. |
| D14 | `pl/metaplot.py:233-240`: rewrite `gene_body_metaplot` to mirror `compute_tss_metaplot`'s per-sample + window-filter loop. |
| D15 | `_compute.py:149-162, 201-209`: submit dask/ray tasks in chunks of `n_workers`; workers write parquet, return paths. |
| D16 | `tl.py:419`; `_dmc_store.py:135-171`: document the `materialize=True` ~700 MB-1.4 GB tradeoff in the public docstring (the streaming opt-out is already wired). |
| D17 | `dmc.py:1502-1503, 2133-2136`: correct the `lr` `meth_diff`/`mean_beta` docstring (coverage-weighted, not Welford; doesn't cite the removed `score` engine). |
| D18 | `tests/test_overdispersed_calibration.py` (extend): add **two-sided** φ-overdispersed FPR tests at φ ∈ {1.5, 2, 3, 5} for `glm` and `welch_t` (the `lr` test already covers φ≈4). Tighten any loose one-sided bounds. **If a test fails, that's a real finding — surface it and rescope.** |

## Minors/Nitpicks pulled into Batch 4 (one commit each)

- `convert._infer_strand` negative-pos guard (`convert.py:329-335`)
- `normalize_coverage` rounding floor (`filter.py:489-498`)
- `intersect_sites` on `(chrom,pos)` not `(chrom,pos,strand)` (`filter.py:555-569`)
- Sliding-window DMR boundary overshoot + silent >10 kb drop (`dmr.py:333-334, 562-563`)
- `empirical_fdr_for_dmr` docstring rewording (`dmr.py:1384-1388`) — pooled-null vs Westfall-Young max-T language
- EB dispersion moments on clamped φ (`_glm.py:671-684`)
- ASM BH pools across samples (`asm.py:113-115`)
- Single-sample clock `impute_missing=True` no-op (`clocks.py:236-273`)
- Multi-mod MM tags mis-align ML probabilities (`bam_io.py:238-289`)
- Scree vs PCA differing site panels (`report.py:915-916`)
- `compute_pca` cache key (`pl/_compute.py:239-241`)
- Coverage-histogram y-axis label (`pl/_compute.py:486-495`)
- `nearest_tss_distance` int32.min sentinel in Int32 (`annotate.py:1333-1334`)
- `annotate_cpg_islands` empty-input schema (`annotate.py:1395-1397`)
- Per-chromosome `except Exception` over-broad handlers (`annotate.py:626-654, 1247-1262`)
- `pl.umap` seed not forwarded (`pl/embedding.py:51`)
- `coverage_histogram` filesystem-traversal-order subsampling (`pl/_compute.py:188-216, 488-493`)
- Parquet `part-*` concat order (`filter.py`, `pp.py`)
- `region_beta` predicate pushdown (`methyldata.py:570-573`)
- DMCStore pointer not rebased on `load()` (`methyldata.py:268-289`)
- `set_tmp_dir(None)` env-var restore (`_config.py:47-53`)
- Checkpoint/resume API doc cleanup (`_cache.py`, `methyldata.py:123-184`)
- `uns` JSON `default=str` silent stringification (`methyldata.py:392-420`)
- `power()` n=10000 silent cap (`qc.py:866-871`)
- `EPYKIT_NO_AUTO_TSV` env-var case-match (`tl.py:254`)
- `materialize=False` resumable-cache hit override (`tl.py:771-811`)
- `reference_level` patsy formula injection (`_glm.py:178`)
- Design-matrix rank check (`_glm.py:237-243`)
- Samplesheet duplicate-id / path-exists check (`io.py:60-82`)
- `to_bigwig` missing-chrom drop (`export.py:186-188`)
- BedGraph `%.6g` scientific notation on int coverage (`export.py:115`)
- DMC-BED `qvalue_combined` gate (`export.py:241-246`)
- nf-core samplesheet `sample` vs `sample_id` (`nfcore_qc.py:169-171`)
- Report "not run" vs "failed" distinction (`report.py:761-768, 902-941`)
- Manhattan banding palette collision (`pl/_plotly.py:146-152`)
- M-bias plot semantic palette (`pl/qc.py:155`)
- Methylation heatmap no sample names (`pl/qc.py:94`)
- `annot_log` ignored in proportion mode (`pl/annotation.py:265`)
- `feature_direction_stacked` mixed normalization (`pl/_plotly.py:399-410`)
- CLI ValueError/RuntimeError → traceback (`cli.py:1029-1032`)
- All nitpick docstring/dead-code items from audit §2.NITPICK (one commit batched)

## Batch-4 acceptance

- `epykit.experimental.{bam_io,asm,entropy,clocks}` exist; importing the old names emits a `DeprecationWarning` referencing the new path.
- New BAM regression tests cover MM/ML reverse-strand parity and bisulfite-SNV exclusion.
- `to_bedgraph`/`to_bigwig` peak RSS bounded by O(largest chromosome) under a measured assertion.
- `tests/test_overdispersed_calibration.py` includes `glm` and `welch_t` at φ ∈ {1.5, 2, 3, 5} — and any failure is documented and scoped, not loosened.
- All D-bundle items have at least one regression test or a documentation change citing this plan.

---

## Cross-batch deliverables

- **CHANGELOG.md** updated per batch with a "Reviewer-visible behavior changes" section calling out: M5 default behaviour change (no-reference merging now works), M9 methylKit off-by-one fix (existing exports were 1 bp left), M4 glm p-value change (now matches lr at floored df).
- **`docs/review/2026-06-07-epykit-codebase-audit.md`** annotated at the top with a "Status: remediated" header once Batch 4 completes; each finding ID gets a commit-hash backlink.
- **No re-running of benchmark headlines** unless one of M4/M5 actually changes a default the paper claims around. The audit's adversarial verification already established M4 is *conservative* and M5 only affects unmerged-strand inputs; the paper uses `combined_strand_bed` per `CLAUDE.md`. If the φ-FPR tests in Batch 4 fail, that's a paper-level finding requiring escalation.

---

## Risks & rollback

| Risk | Mitigation |
|---|---|
| M2 NotImplementedError breaks existing user code | Branch tag before merge; add migration note to CHANGELOG; documented in tl.dmr docstring; the silent no-op was strictly worse. |
| M4 glm fix changes published glm p-values | Bias is purely toward lower (more powerful, less conservative); document; users who want the old behaviour can pass a custom `df_phi`. |
| M5 strandless pair-merge has edge cases on non-standard `.cov` formats | Conservative implementation: only merge pairs where `pos[i+1] - pos[i] == 1`; otherwise leave as-is and log a count of unmerged. |
| Batch-4 calibration tests fail for `glm`/`welch_t` at high φ | Surface as a real finding; potentially rescope to a §2.3-style "recovery pass" investigation; do not loosen bounds to make tests pass. |
| Read-level deprecation breaks downstream notebooks | One-release deprecation window with explicit warning + import shims at the old names. |

---

## Self-review notes

- **Spec coverage:** all 13 Majors (M1-M13), 18 D-items (D1-D18), Minors/Nitpicks have an explicit fix location and acceptance criterion. Refuted item ("paired-end ASM double-count") is intentionally absent. §2.3 "recovery pass" items are deferred *as a research-pass investigation* rather than direct fixes — `dispersion="site"` anti-conservativeness needs a φ-stratified FPR run before action, which I've scheduled within Batch 4's D18 work.
- **Placeholders:** none — every batch-1 step has runnable code; later batches list file + finding + fix.
- **Type consistency:** `_stratified_permutation_assignment` (Task 1.1) and `_empirical_pvalues_from_null_pool` (Task 1.2) are the two new helpers; both are referenced by their full names throughout. `welch_meth_diff_ci`'s new `dof` kwarg (Task 3-M13) matches the existing `dof` local in `dmc.py:1136`.
- **Out-of-scope items that surfaced:** none — every audit finding is either in a batch or explicitly deferred with reason.

---

## Batch 1 — DONE (2026-06-07)

Landed on `review-fixes-batch-2` (10 commits). Gate green: 500 tests passed across `not slow` + `slow` tiers (10 skips, all pre-existing optional-dep gates: pysam, CuPy, pyBigWig, dask-missing-extra).

| Task | Finding | Commit |
|---|---|---|
| 1.1 | M1 within-stratum k-of-n permutation | `298f258` + `855ee3c` lint polish |
| 1.2 | m-perm-1 failed-perm denominator | `8c107a2` + `78df7fb` docstring polish |
| 1.3 | M3 sep_fallback/smoothing into DMC perms | `80ba94f` + `d8e6575` negative-path test |
| 1.4 | M2 NotImplementedError non-tile + CLI gate | `0969565` + `bd1ee86` log message + `70da835` CLI mirror |
| 1.5 | m-perm-2 merge_adjacent/backend into tile perms | `a455194` |

### Follow-ups surfaced during Batch 1 (deferred to Batch 4)

- **M3-DMR** — `empirical_fdr_for_dmr` (tile path) has the same M3 gap as DMC: `sep_fallback`/`smoothing` are not forwarded into permutation `call_dmr_tile_based` invocations. Surfaced by Task 1.3 code-quality review.
- **m-perm-2-DMC** — `empirical_fdr_for_dmc` has no explicit `backend`/`n_workers` parameters and `tl.dmc` doesn't forward them. Performance/parity issue (defaults match), not correctness. Surfaced by Task 1.5 code-quality review.

Both tracked as Tasks #7 and #8 in the session task list.
