# Testing Patterns

**Analysis Date:** 2026-06-06

## Test Framework

**Runner:**
- Framework: pytest >= 8.4.2
- Config: `pyproject.toml [tool.pytest.ini_options]` (lines 132–144)
- Test paths: `tests/` directory (auto-discovered)

**Assertion Library:**
- Standard pytest assertions (`assert`, `assert ... in ...`)
- Custom helpers in `tests/conftest.py` for accuracy metrics:
  - `power_at_threshold(dmc_df, truth, alpha=0.05)` — sensitivity (TPR)
  - `fdr_at_threshold(dmc_df, truth, alpha=0.05)` — empirical false-discovery rate
  - `meth_diff_bias(dmc_df, truth)` — effect-size bias and MAE
  - `dmr_recovery(dmr_df, truth, cfg, alpha=0.05)` — seeded DMR recovery rate

**Run Commands:**
```bash
# CI invocation (default; excludes slow tests)
uv run pytest -m "not slow" --strict-markers -ra

# Slow tests only (opt-in; >~5s tests)
uv run pytest -m slow

# Single test
uv run pytest tests/test_dmc_multigroup.py::test_name

# With coverage
uv run pytest -m "not slow" --cov=epykit --cov-report=term-missing
```

## Test File Organization

**Location:**
- Tests co-located in `tests/` directory (not next to source files)
- Fixtures and shared utilities in `tests/conftest.py` (session-scoped synthetic dataset)
- Test data generators: `tests/fixtures/synth.py` (synthetic Bismark .cov with known truth)

**Naming:**
- Test files: `test_*.py` prefix (auto-discovery)
- Test functions: `test_*` prefix
- Fixtures: named after what they provide (e.g., `synth_md`, `synth_md_filtered`, `synth_bundle`)

**Structure:**
```
tests/
├── conftest.py                        # pytest fixtures (session-scoped synthetic data)
├── fixtures/
│   ├── __init__.py
│   └── synth.py                       # SimConfig + generate() for test methylstores
├── test_api.py                        # MethylData contract, state derivation, save/load
├── test_dmc_lr.py                     # DMC tests: LR engine specifics (Newcombe CI)
├── test_dmc_multigroup.py             # Multi-group / F-test contrasts
├── test_dmc_smooth_dispersion.py      # Dispersion shrinkage (EB, median-unbiased)
├── test_dmc_empirical_fdr.py          # Permutation-based empirical FDR
├── test_cli.py                        # CLI argument parsing, subcommand dispatch
├── test_calibration.py                # @slow: null calibration (K-S test)
└── test_viz_new.py                    # @slow: plotting (matplotlib + plotly)
```

## Test Structure

**Suite Organization:**

```python
"""Module docstring describing test scope.

Layer X: [Purpose]
- Tests what they verify
- Assumptions about fixtures
"""

from __future__ import annotations
import warnings
import pytest
import epykit as ep

# Example from test_dmc_lr.py
def test_lr_emits_log2_odds_ratio_pooled(synth_md_filtered):
    """Test title line: one-liner describing the assertion.
    
    Optional: longer explanation of the test's purpose and any nuance.
    """
    md = synth_md_filtered
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(md, test="lr")
    
    df = md.dmc
    assert "log2_odds_ratio_pooled" in df.columns, (
        f"missing log2_odds_ratio_pooled; got {df.columns}"
    )
    # ... more assertions ...
    
    # Check for expected warnings
    fut = [w for w in caught if issubclass(w.category, FutureWarning)]
    assert fut and "log2_odds_ratio" in str(fut[0].message)
```

**Setup Pattern:**
```python
@pytest.fixture
def synth_md_filtered(synth_md):
    """Fresh MethylData after filtering; ready for DMC.
    
    Fixture docstring explains what the object has undergone.
    """
    import epykit as ep
    ep.pp.filter_coverage(synth_md, lo_count=5, hi_perc=99.9)
    ep.pp.set_unite_type(synth_md, type="intersect")
    return synth_md
```

**Teardown Pattern:**
- Fixtures with temporary directories use pytest's `tmp_path` and `tmp_path_factory`
- No explicit teardown needed; pytest cleans up temp dirs on test completion
- Example from `conftest.py` line 44: `@pytest.fixture(scope="session")` uses `tmp_path_factory` for session-level reuse

**Assertion Pattern:**
```python
# Multi-line assertions with clear context
assert value > threshold, (
    f"Expected {label} > {threshold}; got {value}"
)

# Membership checks
assert item in collection, f"Expected {item} in {collection}"

# Floating-point comparisons with tolerance
assert np.allclose(result, expected, rtol=1e-6, atol=1e-10)

# NaN checks (numpy)
assert np.isnan(value).all(), f"Expected all NaN; got {value[:5]}"
```

## Mocking

**Framework:** unittest.mock (standard library)
- Not heavily used in core tests
- When mocking external I/O: wrap with context managers in test functions
- Example: `unittest.mock.patch()` for file paths, subprocess calls

**Patterns:**
- Mock file systems using `tmp_path` (pytest) instead of full mocking
- Mock subprocess calls for CLI tests (see `test_cli.py` line 20):
  ```python
  def _epykit(*args, check: bool = True, **kwargs) -> subprocess.CompletedProcess:
      cmd = [sys.executable, "-m", "epykit.cli", *map(str, args)]
      return subprocess.run(cmd, capture_output=True, text=True, check=check, **kwargs)
  ```

**What to Mock:**
- External subprocesses (CLI invocations)
- File I/O operations if you need to test error paths (use `tmp_path` for actual data)

**What NOT to Mock:**
- DMC engines — test them against real synthetic data
- Parquet I/O — always use actual files (uses pytest temp dirs)
- Logging calls — assert on log output via `caplog` fixture if needed

## Fixtures and Factories

**Test Data:**

The synthetic dataset builder in `tests/fixtures/synth.py` generates:
```python
@dataclass
class SimConfig:
    """Knobs for the synthetic methylation generator.
    
    Defaults chosen so 4-vs-4 WGBS at ~20x coverage produces
    detectable signal under BH-correction across ~75k post-filter CpGs.
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
    coverage_disp: float = 5.0
    replicate_sd: float = 0.03
    seed: int = 42
    # Multi-group / continuous covariate extensions...
```

Output: `cov/*.bismark.cov.gz`, `samplesheet.csv`, `truth.parquet`, `config.json`

**Location:**
- Fixture generation: `tests/fixtures/synth.py` function `generate(cfg, out_dir)`
- Per-test fresh `MethylData` objects: `synth_md` fixture in `conftest.py` (line 62)
- Session-level synthetic store reuse: `synth_bundle` fixture (line 44, scope="session")

**Usage Example:**
```python
def test_dmc_on_synthetic_data(synth_md_filtered):
    """Test against the 4-vs-4 synthetic fixture."""
    import epykit as ep
    ep.tl.dmc(synth_md_filtered, test="lr")
    # Assertions on md.dmc (the DMC result frame)
```

**SynthBundle Dataclass** (lines 29–40 in `conftest.py`):
```python
@dataclass
class SynthBundle:
    """Bundle of paths + truth table + ids passed around by fixtures."""
    samplesheet: str
    truth: pl.DataFrame                    # Ground truth: is_dmc, true_meth_diff, dmr_id
    store_root: str
    treatment_ids: list[str]               # e.g., ["treatment_0", "treatment_1", ...]
    control_ids: list[str]
    n_dmcs_true: int                       # Count of seeded DMC sites
    n_dmrs: int                            # Count of seeded DMR regions
    config: SimConfig
```

## Coverage

**Requirements:**
- No hard enforcement in CI (not gated)
- Optional: run with `--cov=epykit --cov-report=term-missing` to generate coverage reports

**View Coverage:**
```bash
uv run pytest -m "not slow" --cov=epykit --cov-report=html
# Generates htmlcov/index.html
```

## Test Types

**Unit Tests:**
- Scope: Single module or function (e.g., `test_dmc_lr.py` tests DMC LR engine specifics)
- Approach: Synthetic data fixture with known ground truth; verify accuracy metrics (power, FDR, bias)
- Example: `test_lr_meth_diff_ci_is_asymmetric_near_boundary()` (lines 53–100 in `test_dmc_lr.py`) checks that Newcombe CIs are meaningfully asymmetric near β=0/1

**Integration Tests:**
- Scope: Multi-step pipelines (e.g., filter → dmc → dmr → annotate)
- Approach: Run the full epykit API flow on synthetic data; check end-to-end consistency
- Example: `test_api.py` verifies MethylData round-trip through save/load and state derivation from `_store_history`

**E2E Tests:**
- Scope: CLI argument parsing, subprocess invocation, exit codes
- Framework: subprocess.run() to invoke `python -m epykit.cli`
- Example: `test_cli.py` (lines 34–80) tests `--help` output, subcommand routing, argument validation
- Tests run the CLI as a real subprocess so they catch import-time side effects and argparse wiring

## Common Patterns

**Async Testing:**
Not used in epykit (synchronous I/O only). Pytest's `pytest-asyncio` plugin not in dependencies.

**Error Testing:**
```python
def test_invalid_input_raises_valueerror():
    """Test that invalid parameters are rejected."""
    with pytest.raises(ValueError, match="Invalid hi_perc"):
        ep.pp.filter_coverage(md, hi_perc=-1)
```

**Warning Testing:**
```python
def test_fisher_emits_user_warning(synth_md_filtered):
    """Test that using fisher test emits a warning about bias."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(synth_md_filtered, test="fisher")
    
    # Filter by warning category
    fisher_warns = [w for w in caught if issubclass(w.category, UserWarning)]
    assert fisher_warns, "Expected UserWarning about fisher bias"
```

**Deprecation Testing:**
```python
def test_removed_engine_raises_with_migration_hint():
    """Old engines raise ValueError with a helpful message."""
    with pytest.raises(ValueError, match="logit_t.*use.*welch_t"):
        ep.tl.dmc(md, test="logit_t")
```

## Filtering Warnings

**Configuration:**
```toml
# pyproject.toml [tool.pytest.ini_options]
filterwarnings = [
    # Allow epykit's own warnings through so tests can assert on them
    "default::DeprecationWarning:epykit",
    "default::UserWarning:epykit",
    # Silence third-party noise that flooded test output
    "ignore::DeprecationWarning:pyarrow",
]
```

**Rationale:**
- Load-bearing: Tests need to see and assert on epykit's own DeprecationWarnings (migration hints for removed engines)
- Pyarrow deprecation warnings are suppressed because they're frequent but not actionable by epykit users

## Slow Marker

**Registration:**
```toml
# pyproject.toml [tool.pytest.ini_options]
markers = [
    "slow: tests that take more than ~5s; deselect with -m 'not slow'",
]
```

**Usage:**
```python
@pytest.mark.slow
def test_calibration_on_large_dataset():
    """End-to-end test of null calibration scaffolding.
    
    Runs K-S test on 20k simulated p-values and permutation-based
    empirical FDR. Takes >5s; only run when explicitly requested.
    """
    # ... test implementation ...
```

**CI Behavior:**
- Default CI run: `pytest -m "not slow" --strict-markers -ra` (excludes slow tests)
- Slow tests opt-in for CI runs or local development with `pytest -m slow`
- `--strict-markers` enforces that all markers used in tests are registered in `pyproject.toml`

## CI Matrix

**Platforms and Python Versions:**
- Platforms: `ubuntu-latest`, `windows-latest` (Windows compatibility is load-bearing)
- Python: `3.9`, `3.12` (CI matrix covers these; mypy checks against 3.10)
- Defined in `.github/workflows/test.yml`

**Invocation:**
```bash
uv sync --extra dev --extra all
uv run pytest -m "not slow" --strict-markers -ra
```

**Windows-Specific Considerations:**
- Extras gated by `sys_platform != 'win32'`: `pyBigWig`, `pysam` (methylkit, bam)
- All tests must pass on Windows without those extras
- Path handling must use `pathlib.Path` (not hardcoded `/`)

## Example Test Session

```bash
# Install with dev + all optional extras
$ uv sync --extra dev --extra all

# Run fast tests (CI default)
$ uv run pytest -m "not slow" --strict-markers -ra
tests/test_api.py::test_read_bismark_produces_well_formed_methyldata PASSED
tests/test_dmc_lr.py::test_lr_emits_log2_odds_ratio_pooled PASSED
tests/test_dmc_lr.py::test_lr_meth_diff_ci_is_asymmetric_near_boundary PASSED
tests/test_cli.py::test_cli_dmc_help_shows_lr_default PASSED
... (20+ tests pass in ~30s)

# Run slow tests separately
$ uv run pytest -m slow
tests/test_calibration.py::test_uniform_engine_passes_ks_calibration PASSED [slow]
... (~5 slow tests take 2+ minutes)

# Run coverage report
$ uv run pytest -m "not slow" --cov=epykit --cov-report=term-missing
```

---

*Testing analysis: 2026-06-06*
