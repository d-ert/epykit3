# Coding Conventions

**Analysis Date:** 2026-06-06

## Naming Patterns

**Files:**
- Modules use lowercase with underscores: `methyldata.py`, `dmc.py`, `_glm.py`, `_dmc_store.py`
- Private/internal modules prefixed with underscore: `_config.py`, `_cache.py`, `_compute.py`, `_glm_gpu.py`, `_hmm.py`, `_smoothed_store.py`
- Subpackages (e.g., `pl/`, `templates/`) lowercase
- Re-export files explicitly (e.g., `src/epykit/__init__.py` F401-ignored for re-exports)

**Functions:**
- Snake case: `read_bismark()`, `filter_coverage()`, `apply_multiple_testing_correction()`, `bisulfite_conversion_rate()`
- Prefix with underscore for internal/private functions: `_append_store_history()`, `_auto_test()`, `_welford_init()`
- Module-level constants ALL_CAPS: `_EMPTY_SCHEMA`, `_SMOOTH_BOX_NJIT_FN`, `_BED_BASE_COLS`, `_GLM_BACKENDS`, `_FISHER_WARNED`

**Variables:**
- Snake case throughout: `samplesheet`, `treatment_ids`, `control_ids`, `n_samples`, `store_path`
- Loop variables, accumulators: `M2`, `S0_g`, `S1_g`, `Sigmam^2` (per Welford nomenclature in `dmc.py`)
- Dataframe columns use underscore: `mean_beta_case`, `mean_beta_control`, `meth_diff_ci_lo`, `meth_diff_ci_hi`, `log2_odds_ratio_pooled`

**Types:**
- Literal unions for test engine names: `Literal["lr", "glm", "welch_t", "fisher"]`
- Generic type hints: `dict[str, pl.DataFrame]`, `list[str]`, `Optional[str]` (from typing module)
- Use `from __future__ import annotations` in all files for forward references

## Code Style

**Formatting:**
- Line length: 100 characters (enforced by ruff configuration in `pyproject.toml`)
- Linter: ruff with F-only baseline (undefined names, unused imports only)
- Indentation: 4 spaces
- Type stub: `py.typed` marker present in `src/epykit/py.typed`

**Linting:**
- Baseline check: `uv run ruff check src/` runs F-only rules
- Config in `pyproject.toml [tool.ruff]` sets `select = ["F"]` as the default
- Per-file ignores: F401 in test files and `__init__.py` (re-exports), F811 in tests (redefinition)
- Target version: Python 3.9+ via `target-version = "py39"` (Windows compatibility required)

**Type Checking:**
- Tool: mypy (configured in `pyproject.toml [tool.mypy]`)
- Run: `uv run mypy src/epykit`
- Settings: `python_version = "3.10"`, `follow_imports = "silent"`, `ignore_missing_imports = true`, `no_strict_optional = true`, `warn_unused_ignores = true`

## Import Organization

**Order:**
1. `from __future__ import annotations` (always first)
2. Standard library (`logging`, `json`, `os`, `sys`, `pathlib.Path`, `tempfile`, `time`, `warnings`)
3. Third-party (`numpy`, `polars`, `scipy`, `statsmodels`, `pyarrow`, `patsy`, `bioframe`)
4. Local package imports (relative imports with `.` prefix)

**Path Aliases:**
- First-party: `epykit` (configured as `known-first-party` in ruff isort settings)
- Relative imports preferred within package: `from . import filter`, `from ._cache import count_store_rows`
- Re-exports in `__init__.py` use absolute imports: `from .methyldata import MethylData`

**Pattern - Import Submodules:**
- `import epykit as ep` in user code and tests
- Within library: `from . import pp, tl, pl, query` (register namespaces in `__init__.py`)
- Lazy imports for heavy dependencies inside functions (e.g., `jinja2`, `plotly` in `report.py`)

## Error Handling

**Patterns:**
- Raise `ValueError` for invalid user input: `raise ValueError(f"Invalid hi_perc/quantile value: {hi_perc}")`
- Raise `FileNotFoundError` or use `Path.resolve()` for missing files
- Check preconditions and fail early: `if quantile <= 0 or quantile > 1: raise ValueError(...)`
- Preserve deprecation shims with clear migration hints (see "Deprecation pattern" below)

**Warnings:**
- Use stdlib `warnings.warn()` with `DeprecationWarning` or `UserWarning`
- One-shot warning gates (module-level flag): `_FISHER_WARNED` in `tl.py` prevents spam across chromosomes
- Tests can assert on warnings via `warnings.catch_warnings(record=True)` (see `test_dmc_lr.py`)

## Logging

**Framework:** Standard library `logging` module

**Pattern (Library Code):**
```python
import logging
logger = logging.getLogger(__name__)
```
- Every module under `epykit.*` (except `epykit.cli`) declares a module-level logger
- Use `logger.info()`, `logger.debug()`, `logger.warning()` for progress/diagnostics
- **CRITICAL: Never call `print()` in library code** — stdout is reserved for CLI output

**Pattern (CLI Code):**
- File: `src/epykit/cli.py`
- CLI entry point reserves `print()` for final user-facing result lines on stdout
- Structured progress logs flow through the same logging hierarchy, controlled via `-v` / `-q` flags
- This split allows host applications and notebooks to consume epykit without stdout pollution

**Example from `src/epykit/dmc.py`:**
```python
import logging
logger = logging.getLogger(__name__)

def process_chromosomes_dmc(...):
    logger.info(f"Processing {len(chromosomes)} chromosomes...")
    # ... no print() calls ...
```

**Example from `src/epykit/cli.py`:**
```python
import logging
# Configure logging for CLI output
def main(args):
    # ... call library code ...
    # Final result:
    print(f"DMC results written to {output_path}")  # OK for CLI
```

## Comments

**When to Comment:**
- Explain *why* non-obvious logic exists (e.g., why Newcombe CIs vs Wald, why Welford's algorithm is used)
- Pinpoint complex algorithms or statistical choices
- Reference papers/issues: `P1-11: log2_odds_ratio renamed per backend` (see `test_dmc_lr.py`)
- Avoid restating what the code clearly shows

**JSDoc/TSDoc:**
- Use docstrings for public functions and classes
- Triple-quote format for docstring blocks
- Include Parameters, Returns, Raises sections for non-trivial APIs
- Example from `src/epykit/pp.py`:
  ```python
  def normalize_coverage(md: MethylData, method: str = "median") -> None:
      """Per-sample coverage normalisation, in-place on a MethylData.
      
      Computes a per-sample scaling factor so that each sample's central
      coverage statistic (median by default, or mean) matches a common
      target...
      
      Parameters
      ----------
      md : MethylData
          Object whose store has been ``filter_coverage``'d.
      method : {"median", "mean"}
          Central statistic to align.
      
      Raises
      ------
      ValueError
          If ``filter_coverage`` has not been called yet.
      """
  ```

## Function Design

**Size:** Aim for functions that fit on one screen (~50 lines); refactor helpers for clarity over brevity

**Parameters:**
- Use descriptive names; avoid single letters except loop vars and math conventions
- Defaults should be safe/sensible (e.g., `lo_count=10`, `hi_perc=99.9` for coverage filtering)
- Keyword-only args for clarity on complex functions: `def filter_coverage(..., blacklist_bed: str | None = None, output_store: str | None = None)`
- Type hints mandatory on all public API functions

**Return Values:**
- Mutating functions return `None` (e.g., `filter_coverage(md, ...)` mutates `md` in place, returns nothing)
- Data-returning functions return typed values (e.g., `power_at_threshold(...) -> float`)
- When multiple related values are returned, use a dataclass (`@dataclass`) not a tuple (e.g., `SynthBundle` in `tests/conftest.py`)

## Module Design

**Exports:**
- Public API in `__all__` list in `__init__.py` (e.g., lines 102–161)
- Private/internal functions prefixed with `_` and not exported
- Deprecation shims use `__getattr__` for backward compat (lines 173–192)

**Barrel Files:**
- `src/epykit/__init__.py` re-exports: main entry points (`read_bismark`, `MethylData`), namespaces (`pp`, `tl`, `pl`), and top-level tools
- Per-file ignores for F401 (unused imports in re-exports) in ruff config
- Deprecation shims documented in module docstring

## Platform Compatibility

**Windows Load-Bearing:**
- CI matrix: `{ubuntu-latest, windows-latest} × {py3.9, py3.12}`
- Extras gated by `sys_platform != 'win32'`: `pyBigWig` (export), `pysam` (methylkit, bam) in `pyproject.toml`
- Path handling: Always use `pathlib.Path`, never hardcoded `/` separators
- Temp directory: Use `ep.set_tmp_dir(path)` to redirect tempfile on Windows (where `C:\` may be too small for whole-genome staging)
- All functions must work cross-platform unless explicitly documented as platform-specific

## Canonical DMC Output Schema

Every DMC engine (`"lr"`, `"welch_t"`, `"fisher"`, `"glm"`) outputs the same base schema:
- `chrom` (Utf8) — chromosome name
- `pos` (Int32) — genomic position
- `strand` (Utf8) — strand (if applicable)
- `n_case` (Int32) — number of treatment/case samples with coverage
- `n_control` (Int32) — number of control samples with coverage
- `mean_beta_case` (Float32) — mean methylation fraction in case group
- `mean_beta_control` (Float32) — mean methylation fraction in control group
- `pvalue` (Float64) — raw p-value (per-site significance test)
- `meth_diff` (Float32) — methylation difference (treatment − control, in [−1, 1])
- `meth_diff_ci_lo` (Float32) — lower confidence interval bound
- `meth_diff_ci_hi` (Float32) — upper confidence interval bound
- `log2_odds_ratio_pooled` (Float64) — log2(odds ratio) for lr/fisher backends
- `log2_odds_ratio` (Float64) — **deprecated transitional column, NaN-filled with FutureWarning**

**Engine-Specific Extras:**
- GLM: adds `coef_treatment`, `coef_treatment_log2` (logit coefficient, not pooled log2-OR)
- Multi-group F-tests: adds `f_stat`, `df1`, `df2`
- When `neighbour_combine=True`: adds `pvalue_combined`, `qvalue_combined`, `pvalue_combined_n_neighbours`, `qvalue_combined_reject` (raw `pvalue`/`qvalue` remain unchanged)

See `dmc.py` lines 54–68 for the `_EMPTY_SCHEMA` definition.

## Store History Pattern

**Load-Bearing Convention:**
Preprocessing state (filtered, united, smoothed) is **derived** from `uns["_store_history"]` rather than stored as independent booleans. This prevents state drift.

**Pattern - Adding a Pipeline Step:**
```python
def new_preprocessing_step(md: MethylData, ...) -> None:
    """Your step description."""
    # ... compute output store ...
    md.store = output_dir
    
    # Record the step in history (load-bearing for md._filtered / md._united / etc.)
    n_sites = _count_parquet_rows(output_dir)  # optional but recommended
    _append_store_history(md, "step_name", output_dir, n_sites)
    
    # Optional: store step-specific config for audit trail
    md.uns["step_name"] = {"param1": value1, "param2": value2}
```

Helper: `_append_store_history(md, step, path, n_sites)` in `src/epykit/pp.py` appends to the history list.

**Why:**
- Allows any downstream code to check `md.state` to see which steps have run
- Properties like `md._filtered` read the history: `any(h.get("step") == "filtered" for h in md.uns.get("_store_history", []))`
- Prevents bugs where boolean flags drift from actual store state

## Deprecation & Migration Pattern

**For Removed Engines:**
When an engine is no longer supported, use a clear migration hint:

```python
def _canonicalise_test_name(test: str) -> str:
    """Resolve shorthand and deprecated engine names; raise if invalid."""
    if test == "logit_t":
        raise ValueError(
            "test='logit_t' was deprecated in 0.7.5 and removed in 1.0. "
            "Use test='welch_t' instead (both apply Welch t to logit-transformed betas)."
        )
    if test == "bb_lr":
        raise ValueError(
            "test='bb_lr' (beta-binomial LR) was removed in 1.0. "
            "Use test='lr' (quasi-binomial LR) instead."
        )
    # ... etc ...
```

Removed engines (as of 0.7.5 / 1.0):
- `"logit_t"` → use `"welch_t"`
- `"bb_lr"` → use `"lr"`
- `"score"` → use `"lr"`
- `"cmh"` → use formula-based GLM: `formula='~ group + batch'`

Documented in `CLAUDE.md` and `src/epykit/cli.py` docstring.

---

*Convention analysis: 2026-06-06*
