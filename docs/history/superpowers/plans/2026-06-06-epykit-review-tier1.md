# epykit Review Remediation — Tier 1 Implementation Plan

> **For agentic workers:** Use TDD where a behavioral oracle exists. Each task = atomic commit with its test. Steps use checkbox (`- [ ]`) syntax. Full finding detail: [docs/review/2026-06-06-epykit-peer-review.md](../../review/2026-06-06-epykit-peer-review.md). Design rationale + decisions: [docs/history/superpowers/specs/2026-06-06-epykit-review-fixes-design.md](../specs/2026-06-06-epykit-review-fixes-design.md).

**Goal:** Remediate the paper-blocking + silent-wrong-science Critical findings (repro cluster + C1/C2/C3) without disturbing the benchmark paper's numbers.

**Architecture:** Surgical edits to `convert.py`/`io.py` (C1), `annotate.py` (C2), `clocks.py` (C3), `cli.py`+`pyproject.toml`+CI (repro cluster). Each lands with a test that would have caught the bug; where a test currently encodes the bug, it is corrected.

**Tech stack:** polars, numpy, scipy, pytest, uv, GitHub Actions.

**Branch:** `review-remediation` (off `main` @ `7e51dc4`). Decisions baked in: phased; fix-the-science-properly; C1 via start/end auto-detect with `coordinate_base` override; C4/C5 deferred to Tier 2.

**Tier 1 gate:** `uv run pytest -m "not slow"` green on Windows (local) + the new slow job green on Linux; CLI and API produce identical DMC q-values on the parity fixture; fresh resolve installs polars≥1.0.

---

## Task ordering (smallest/safest → largest)
T1.7 (Horvath) → T1.1 (polars floor) → T1.3 (uv.lock) → T1.2 (CLI/API parity) → T1.6 (C2 chrom guard) → T1.5 (C1 coordinate) → T1.4 (calibration test + slow CI).

---

### Task T1.7 — C3: Horvath clock anti-transform

**Files:** Modify `src/epykit/clocks.py` (line 238 + docstring ~153); Test `tests/test_stats_new.py:230-240`.

- [ ] Extract pure helper in `clocks.py` (module scope):
```python
def _horvath_anti_transform(linear, adult_age: float = 20.0):
    """Horvath (2013) anti.trafo: (1+adult_age)*exp(x)-1 for x<0 else (1+adult_age)*x+adult_age."""
    linear = np.asarray(linear, dtype=float)
    return np.where(
        linear < 0,
        (1.0 + adult_age) * np.exp(linear) - 1.0,
        (1.0 + adult_age) * linear + adult_age,
    )
```
- [ ] Replace `clocks.py:238` body with `age = _horvath_anti_transform(linear)`.
- [ ] Fix the docstring at ~153 (`age = exp(linear) - 1 if linear < 0`) → `age = (1+adult_age)*exp(linear) - 1 if linear < 0 else (1+adult_age)*linear + adult_age`.
- [ ] Rewrite `test_age_clock_horvath_transform_branches` to call the real helper:
```python
def test_age_clock_horvath_transform_branches():
    from epykit.clocks import _horvath_anti_transform
    # negative branch: 21*exp(-0.5) - 1 = 11.7371...
    assert abs(float(_horvath_anti_transform(-0.5)) - 11.7371438) < 1e-5
    # positive branch: 1.0*21 + 20 = 41.0
    assert float(_horvath_anti_transform(1.0)) == 41.0
    # vector + adult_age override
    import numpy as np
    out = _horvath_anti_transform(np.array([-0.5, 1.0]))
    assert out.shape == (2,)
```
- [ ] Run: `uv run pytest tests/test_stats_new.py -q` → PASS. Commit `fix(clocks): correct Horvath anti-transform negative branch (C3)`.

**Acceptance:** young/under-20 ages correct; the test exercises the real transform.

---

### Task T1.1 — M-PKG1: polars floor + dependency-floor audit

**Files:** Modify `pyproject.toml`; Test `tests/test_compute_backends.py` (or new pivot smoke).

- [ ] `pyproject.toml`: `polars>=0.20.0` → `polars>=1.0,<2`. In the `all` extra, the duplicated `scipy>=1.11` is fine. Reconcile `numba>=0.59`→`numba>=0.60` (numpy-2 support) since numpy is uncapped.
- [ ] Add a pivot smoke test asserting a `DataFrame.pivot(on=...)` path runs (guards the floor):
```python
def test_polars_pivot_on_keyword_available():
    import polars as pl
    df = pl.DataFrame({"i": [1, 1], "k": ["a", "b"], "v": [1.0, 2.0]})
    wide = df.pivot(values="v", index="i", on="k")
    assert "a" in wide.columns and "b" in wide.columns
```
- [ ] Run: `uv run pytest -q -k "pivot or compute_backends"` → PASS. Commit `fix(deps): bump polars floor to >=1.0 (pivot on= requires 1.0) (M-PKG1)`.

**Acceptance:** floor matches actual API usage; smoke test green.

---

### Task T1.3 — M-PKG4: commit uv.lock + document thread pinning

**Files:** Modify `.gitignore` (remove line 273 `uv.lock`); `.github/workflows/test.yml` (cache comment + glob); `benchmark/README.md` (thread pinning note); add `uv.lock` to tracking.

- [ ] Remove `uv.lock` from `.gitignore` (line 273); leave the explanatory comment trimmed to note it is now tracked.
- [ ] `test.yml`: change `cache-dependency-glob: "pyproject.toml"` → `"uv.lock"` (3 occurrences) and update the line-26 comment.
- [ ] Append a "Reproducing the published numbers" note to `benchmark/README.md` documenting `uv sync --frozen` and `POLARS_MAX_THREADS=1` / `OMP_NUM_THREADS=1` pinning.
- [ ] `git add uv.lock .gitignore .github/workflows/test.yml benchmark/README.md` → commit `chore(repro): commit uv.lock + document thread pinning (M-PKG4)`.

**Acceptance:** `uv.lock` tracked; `uv sync --frozen` path documented.

---

### Task T1.2 — M-PKG2: CLI/API q-value parity

**Files:** Modify `src/epykit/cli.py` (`p_dmc` subparser ~668; `_cmd_dmc` non-contrast call ~219-229); Test new `tests/test_cli_api_parity.py`.

- [ ] Add to `p_dmc` (after `--allow-n1`):
```python
p_dmc.add_argument("--dispersion", choices=["site", "eb", "shrink", "chrom"], default="eb",
                   help="Dispersion estimator (default: eb, matching ep.tl.dmc).")
p_dmc.add_argument("--reference", choices=["adaptive", "chi2", "f"], default="adaptive",
                   help="Reference distribution for the LR statistic (default: adaptive).")
p_dmc.add_argument("--fdr-method", dest="fdr_method", default="fdr_bh",
                   help="Multiple-testing method (default: fdr_bh).")
```
  (Confirm the exact `choices` against `process_chromosomes_dmc`'s accepted `dispersion`/`reference` values when implementing.)
- [ ] In `_cmd_dmc` non-contrast path: pass `dispersion=args.dispersion, reference=args.reference` to `process_chromosomes_dmc(...)` and `method=args.fdr_method` to `apply_multiple_testing_correction(...)`.
- [ ] Write `tests/test_cli_api_parity.py`: build a small store, run `tl.dmc(test="lr")` (eb default) and the CLI handler/`process_chromosomes_dmc` with matched args; assert `qvalue` columns equal within 1e-12.
- [ ] Run: `uv run pytest tests/test_cli_api_parity.py -q` → PASS. Commit `fix(cli): dmc defaults match tl.dmc (dispersion=eb) + expose knobs (M-PKG2)`.

**Acceptance:** CLI and API q-values identical on the fixture; `--dispersion/--reference/--fdr-method` in `--help`.

---

### Task T1.6 — C2: annotation chromosome-name mismatch guard

**Files:** Modify `src/epykit/annotate.py` (after feature/TSS/island maps are built); Test `tests/test_annotate_multi.py`.

- [ ] After `features_by_chrom`/`tss_by_chrom` (and island BED load), add a guard helper that computes `set(site_chroms) ∩ set(feature_chroms)`: raise `ValueError` with sampled names from each side when empty; `logger.warning` when coverage <50%. Apply to both gene-feature and island joins.
- [ ] Optional `harmonize_chr: bool = False` param to add/strip `chr` prefix on the feature side.
- [ ] Test: annotate a `chr1`-named site frame against a `1`-named GTF → assert `ValueError` (and harmonize path succeeds if implemented).
- [ ] Run: `uv run pytest tests/test_annotate_multi.py -q` → PASS. Commit `fix(annotate): guard chromosome-name mismatch -> no silent 100%-intergenic (C2)`.

**Acceptance:** mismatch raises/warns instead of silently returning all-intergenic.

---

### Task T1.5 — C1: Bismark coordinate auto-detect (largest)

**Files:** Modify `src/epykit/convert.py` (else branch ~395-413; manifest dict ~100-110; docstring 3-9), `src/epykit/io.py` (`read_bismark` + `read_methyldackel` + `_read_methylation_samplesheet` to thread `coordinate_base`); Tests `tests/fixtures/synth.py`, `tests/test_methyldackel.py`, new end-to-end regression.

- [ ] Add `coordinate_base: Literal["auto","one_based","zero_based"] = "auto"` to `convert`, `read_bismark`, `read_methyldackel`, `_read_methylation_samplesheet`; thread through.
- [ ] In `convert.py` else branch: read `end`; resolve offset:
  - `methyldackel` → always 0-based (`pos = start`).
  - `bismark` + `auto`: sample first N rows; if `start == end` predominates → 1-based → `pos = start - 1` (one `logger.info`); if `start == end - 1` → 0-based → `pos = start`.
  - explicit `one_based`/`zero_based` override bypasses detection.
  - `combined_strand_bed` branch unchanged (0-based).
- [ ] Bump manifest payload (`convert.py:107` area): add `"manifest_version": 2` and `"coordinate_base": <resolved>`. On load of a store with no/old version written from real Bismark, `logger.warning` that positions may be +1 (handled in `_cache`/io load path).
- [ ] Rewrite `convert.py:3-9` docstring to the truth (Bismark `.cov` 1-based `start==end`; bedGraph/combined 0-based; auto-detect + override).
- [ ] Fix `tests/fixtures/synth.py` `_write_cov_gz`: emit real 1-based Bismark `.cov` (`start = pos + 1`, `end = start`); update the comment. (Verify the session-scoped store still yields the same `pos` the truth table uses — truth `pos` stays the 0-based store coordinate, so the writer must shift +1 on disk and the converter shifts −1 back.)
- [ ] Rewrite `tests/test_methyldackel.py` equivalence: bismark rows at `start = mdk_start + 1, end = start` (1-based) vs methyldackel rows at `(mdk_start, mdk_start+1)` (0-based); assert identical `pos`.
- [ ] New end-to-end regression (`tests/test_combined_strand_bed.py` or a new `tests/test_coordinate_convention.py`): a `.cov` row at a known 1-based pos → `read_bismark` → assert `pos` is the 0-based coord → `annotate_features` against a known feature → assert expected `feature_type`/distance.
- [ ] Run: `uv run pytest tests/test_methyldackel.py tests/test_coordinate_convention.py -q` then full `uv run pytest -m "not slow" -q` → PASS. Commit `fix(convert): ingest Bismark .cov as 1-based via start/end auto-detect (C1)`.

**Acceptance:** real 1-based `.cov` → correct 0-based `pos`; pre-shifted 0-based detected (no double shift); coord→annotation chain pinned. Migration: manifest-version bump + load warning; CHANGELOG note.

**Risk:** changes outputs for real-Bismark users. The auto-detect keeps existing 0-based fixtures/pre-shifters working (their files have `start==end-1`).

---

### Task T1.4 — M-PKG3: slow-tier CI job + real null-calibration test

**Files:** Modify `.github/workflows/test.yml` (new job); new `tests/test_null_calibration.py`.

- [ ] `tests/test_null_calibration.py` (slow-marked): build a true-null store (both groups from one beta-binomial at φ≈2, or label-permuted) using the `synth` generator with `dmc_effect=0`/`n_dmrs=0`/`n_scattered_dmcs=0`; run `ep.tl.dmc(test=t)` for `t in ("lr","glm","welch_t")`; assert (i) `scipy.stats.kstest(pvalues, "uniform")` not rejected at α=0.01; (ii) empirical FPR at raw `pvalue<0.05` ∈ [0.03, 0.07].
- [ ] `test.yml`: add a `slow` job (ubuntu-latest, py3.12) running `uv run pytest -m slow --strict-markers -ra`. Bound cost with `if: github.event_name == 'push'`.
- [ ] Run locally: `uv run pytest tests/test_null_calibration.py -q` → PASS for all three engines. **If an engine fails calibration, that is a real finding — record it, do not loosen the test.**
- [ ] Commit `test(calibration): real-engine null FPR/KS test + slow CI job (M-PKG3)`.

**Acceptance:** calibration test passes for all three engines on the new CI job; or a genuine calibration failure is documented for Tier 2.

---

## Self-review
- Spec coverage: 1.1✓(T1.1) 1.2✓(T1.2) 1.3✓(T1.3) 1.4✓(T1.4) 1.5✓(T1.5) 1.6✓(T1.6) 1.7✓(T1.7).
- Type/name consistency: `_horvath_anti_transform`, `coordinate_base`, `--dispersion/--reference/--fdr-method` used consistently across tasks.
- Open confirmations (resolve at implementation): exact `choices` accepted by `process_chromosomes_dmc` for `dispersion`/`reference`; the manifest-write call site line; `_read_methylation_samplesheet`→`convert` call signature.
