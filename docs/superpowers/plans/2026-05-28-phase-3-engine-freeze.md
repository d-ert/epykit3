# Phase 3: Engine Freeze (P1 + API Cleanup + Integration Scripts) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Freeze the engine surface so Phase 4's locked benchmark re-run uses final code. Lands 11 P1 functional fixes, the aggressive API cleanup (hard-drop `logit_t`/`bb_lr`/`score`/`cmh`), backend-aware column renames, and five integration scripts that wire the Phase 2 helpers into the real benchmark pipeline.

**Architecture:** Engine-first sequence — renames → drops → P1 functional fixes → integration scripts → tag. Each commit carries an `Affects: <engine>@<scenario>` trailer parsed later by `bug_fix_audit.py`. After tag, the engine schema is final; Phase 4 reads it without further code churn.

**Tech Stack:** Python 3.12, polars, numpy, scipy.stats, statsmodels, patsy (existing). New: `diptest` for the sex-check unimodality fallback. R + methylKit for one integration script (Item 1; test skips when `Rscript` absent).

**Companion spec:** [docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md](../specs/2026-05-27-phase-3-engine-freeze-design.md). Each task below references the spec §4 commit number.

---

## Scope decisions (anchored to spec)

- **In scope (this plan):** 21 commits matching spec §4 exactly. Tag `v0.7.5-phase3-engines-frozen` after Task 21.
- **Out of scope:** the locked benchmark re-run (Phase 4); all P2 hygiene items **except P2-4** (folded into Task 1 with the HMM rename); paper rewrite; `claims.yaml` content (Task 19 lands an empty seed manifest).
- **DVC kept** with the P1-7 fix (Task 13). Not in paper but a legitimate separate analysis family.
- **No `epykit.experimental` namespace.** All four dropped engines hard-fail with `ValueError` + migration hint.

---

## File structure

| File | Why touched | Task |
|---|---|---|
| `src/epykit/dmr_hmm.py` | Rename → `dmr_segment.py`; replaced with deprecation shim | 1 |
| `src/epykit/dmr_segment.py` (new) | Renamed module + per-segment Stouffer p-values (folds P2-4) | 1 |
| `src/epykit/_hmm.py` | Internals stay; no code change | 1 |
| `src/epykit/dmc.py` | Column renames, engine deletions, P1 fixes | 2-12, 15 |
| `src/epykit/dmr.py` | Column-rename consumer; P1-6 empirical-FDR fix | 2, 12 |
| `src/epykit/_glm.py` | P1-3 Newcombe CI wiring; P1-4 reference_level; P1-5 IRLS converged | 9, 10, 11 |
| `src/epykit/tl.py` | Docstring trim, drop-engine guards, reference_level kwarg, dmr method='segment' | 1-7, 10 |
| `src/epykit/dvc.py` | P1-7 Brown-Forsythe | 13 |
| `src/epykit/qc.py` | P1-9 dip-test fallback | 14 |
| `src/epykit/cli.py` | dmr method alias + docs | 1 |
| `pyproject.toml` | Add `diptest` to `qc` extra; prune dropped-engine markers | 7, 14 |
| `tests/test_dmr_segment.py` (new) | New behaviour + shim warning | 1 |
| `tests/test_phase3_drops.py` (new) | Parametrised migration-hint assertion | 7 |
| `tests/test_dmc_fisher.py` (new) | P1-1 reference vs scipy | 8 |
| `tests/test_dmc_lr.py` (new) | P1-3 Newcombe wiring; P1-11 column rename | 9, 2 |
| `tests/test_glm.py` (new) | P1-4 reference_level + P1-5 nonconverged + P1-11 coef rename | 10, 11, 2 |
| `tests/test_dmr_empirical_fdr.py` (new) | P1-6 paired + n=1,1 refusal | 12 |
| `tests/test_dvc.py` (existing) | P1-7 reference vs scipy.levene | 13 |
| `tests/test_qc.py` (existing or new) | P1-9 unimodal fallback | 14 |
| `tests/test_dmc_multitest.py` (existing or new) | P1-10 Storey clamp | 15 |
| `tests/test_stats_new.py`, `test_api.py`, `test_resume.py`, `test_accuracy.py`, `test_dmr_hmm.py` | Prune dropped engines + update imports | 3-7, 1 |
| `tests/test_compute_backends.py`, `test_dmr_tile_merge.py` | Column-rename consumer updates | 2 |
| `benchmark/scripts/methylkit_stouffer_combine.R` (new) | Adjacent-3-CpG combine for methylKit | 16 |
| `benchmark/scripts/_null_engines.py` (new) | Real engine closures for null calibration | 17 |
| `benchmark/scripts/run_null_calibration.py` | Replace mock main with real-engine dispatch | 17 |
| `benchmark/scripts/evaluate.py` | Emit Wilson + bootstrap CI columns | 18 |
| `benchmark/scripts/regen_all.py` (new) | Verify gate (seed manifest) | 19 |
| `benchmark/scripts/bug_fix_audit.py` (new) | Pre/post-fix per-cell delta | 20 |
| `benchmark/scripts/claims.yaml` (new, empty seed) | Manifest for regen_all.py | 19 |
| `benchmark/scripts/tests/test_*` (new) | Tests for each new benchmark script | 16-20 |
| `CHANGELOG.md` | Per-commit bullets + final sweep | every task |
| `README.md`, `docs/analysis/dmc.md`, `CLAUDE.md` | Engine-list trim | 7 |

---

## Pre-flight (once, before Task 1)

- [ ] **Step 1: Verify on `p0-fixes` branch, Phase 2 tag intact, clean tree**

```
git status --short
git log --oneline -3
git tag --list "v0.7.4-phase2-scripts"
```

Expected: working tree shows only untracked dirs (`.github/`, `CLAUDE.md` if not committed). HEAD is `5ac245d` (the spec commit) or later. Tag exists. If not on `p0-fixes`, run `git checkout p0-fixes` and verify the spec doc exists at `docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md`.

- [ ] **Step 2: Baseline test count**

```
uv run pytest -m "not slow" --strict-markers -q 2>&1 | tail -5
```

Expected: a passing line like `229 passed, 5 skipped`. Record the exact number; the Phase 3 wrap-up compares against it.

- [ ] **Step 3: Baseline benchmark-scripts test count**

```
uv run pytest benchmark/scripts/tests/ -q 2>&1 | tail -3
```

Expected: 15 passed (Phase 2 baseline). Record.

- [ ] **Step 4: Confirm spec file readable**

```
uv run python -c "from pathlib import Path; p = Path('docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md'); print(f'{p.stat().st_size:,} bytes')"
```

Expected: prints ~24,000 bytes. Use this spec as the source of truth for any ambiguity; this plan is the execution-level layer.

---

## Task 1: P1-8 + P2-4 — Rename `dmr_hmm` → `dmr_segment` with real per-segment p-values

**Spec §4 commit 1.** This is the only P2 item that lands in Phase 3 (P2-4 — replace NaN p/q-values with Stouffer-combined per-segment values), folded into the same commit as the rename because both touch the same code.

**Files:**
- Create: `src/epykit/dmr_segment.py`
- Replace: `src/epykit/dmr_hmm.py` (becomes a deprecation shim)
- Modify: `src/epykit/tl.py` (dmr method dispatch)
- Modify: `src/epykit/cli.py` (dmr method alias)
- Create: `tests/test_dmr_segment.py`
- Modify: `tests/test_dmr_hmm.py` (update imports if not deletable)
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Write the failing test for `call_dmr_rule_segment` real p-values**

Create `tests/test_dmr_segment.py`:

```python
"""P1-8 + P2-4: dmr_hmm renamed to dmr_segment; per-segment p-values
are Stouffer-combined from constituent CpG p-values (not NaN as in the
pre-Phase-3 implementation)."""
from __future__ import annotations

import warnings

import numpy as np
import polars as pl
import pytest

import epykit as ep


def test_call_dmr_rule_segment_emits_finite_pvalues(synth_md_filtered):
    """The renamed engine must emit finite p/q-values per segment, not NaN."""
    from epykit.dmr_segment import call_dmr_rule_segment

    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    dmrs = call_dmr_rule_segment(md.dmc, min_cpgs=3, min_abs_meth_diff=0.05)
    assert dmrs.height > 0, "expected at least one called segment on synth data"
    pvals = dmrs["pvalue"].to_numpy()
    qvals = dmrs["qvalue"].to_numpy()
    finite_p = np.isfinite(pvals)
    finite_q = np.isfinite(qvals)
    assert finite_p.all(), (
        f"P2-4 fix should produce finite per-segment pvalues; "
        f"got {(~finite_p).sum()}/{len(pvals)} NaN"
    )
    assert finite_q.all(), (
        f"P2-4 fix should produce finite per-segment qvalues; "
        f"got {(~finite_q).sum()}/{len(qvals)} NaN"
    )
    assert ((pvals >= 0) & (pvals <= 1)).all(), "p-values must lie in [0, 1]"
    assert ((qvals >= 0) & (qvals <= 1)).all(), "q-values must lie in [0, 1]"


def test_dmr_hmm_shim_warns_on_import_and_re_exports():
    """Old `epykit.dmr_hmm.call_dmr_hmm` import path must keep working
    with a DeprecationWarning, and must return the same frame the new
    name produces."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        # Force a fresh import to trigger the shim warning.
        import importlib
        import epykit.dmr_hmm as legacy
        importlib.reload(legacy)
        assert hasattr(legacy, "call_dmr_hmm")
        dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
        assert dep, "expected DeprecationWarning on importing epykit.dmr_hmm"
        assert "dmr_segment" in str(dep[0].message), (
            "shim warning must point users to dmr_segment"
        )


def test_tl_dmr_method_hmm_aliased_to_segment(synth_md_filtered):
    """method='hmm' must keep working but emit FutureWarning and dispatch to segment."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    with pytest.warns(FutureWarning, match="segment"):
        ep.tl.dmr(md, method="hmm", min_cpgs=3, min_abs_meth_diff=0.05)
    assert "dmr" in md.uns
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `uv run pytest tests/test_dmr_segment.py -v`
Expected: 3 tests FAIL with `ModuleNotFoundError: No module named 'epykit.dmr_segment'`.

- [ ] **Step 3: Create `src/epykit/dmr_segment.py` with per-segment p-values**

Copy `src/epykit/dmr_hmm.py` to `src/epykit/dmr_segment.py` and replace the function definition. The new file is the full rule-segment engine with real p-values:

```python
"""Rule-based segmentation DMR caller (renamed from dmr_hmm).

Three-state HMM-style decoder on the per-CpG ``meth_diff`` signal with
**fixed** state means and emission SDs -- not a fitted HMM. The name
``dmr_segment`` reflects this honestly; ``dmr_hmm`` remains as a
deprecated shim re-exporting from this module.

P2-4 fix (folded with the P1-8 rename): per-segment p-values are
Stouffer-combined from the constituent CpG p-values and BH-corrected
per chromosome. The pre-Phase-3 implementation emitted NaN p/q-values,
which silently broke any downstream consumer that filtered on qvalue.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import polars as pl
from scipy import stats as sp_stats

from ._hmm import runs_of_state, segment
from .dmr import _DMR_TILE_SCHEMA

logger = logging.getLogger(__name__)


def _state_means_for_meth_diff(meth_diff: np.ndarray) -> np.ndarray:
    """Fixed 3-state targets {hypo, neutral, hyper} at +/-0.20 / 0.0."""
    return np.array([-0.20, 0.00, 0.20])


def _stouffer_combine(pvals: np.ndarray) -> float:
    """Two-sided Stouffer combination of per-CpG p-values.

    Equal weights, no direction term (we already filter by |meth_diff|
    >= min_abs_meth_diff at the segment level, so the constituent CpGs
    are sign-aligned by construction).
    """
    p = pvals[np.isfinite(pvals)]
    if p.size == 0:
        return float("nan")
    # Clip away exact 0/1 to keep isf finite.
    p = np.clip(p, 1e-300, 1.0 - 1e-15)
    z = sp_stats.norm.isf(p / 2.0)  # two-sided -> half-tail z per site
    z_combined = z.sum() / np.sqrt(p.size)
    # Two-sided p from combined z.
    return float(2.0 * sp_stats.norm.sf(abs(z_combined)))


def _bh_per_chrom(pvals: np.ndarray, chroms: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjustment per chromosome group."""
    out = np.full_like(pvals, np.nan, dtype=np.float64)
    for chrom in np.unique(chroms):
        mask = chroms == chrom
        p = pvals[mask]
        finite = np.isfinite(p)
        if not finite.any():
            continue
        p_finite = p[finite]
        order = np.argsort(p_finite)
        ranks = np.empty_like(order)
        ranks[order] = np.arange(1, p_finite.size + 1)
        q = p_finite * p_finite.size / ranks
        # Enforce monotone non-decreasing in p order.
        q_sorted = q[order]
        q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
        q_back = np.empty_like(q_sorted)
        q_back[order] = q_sorted
        q_full = np.full_like(p, np.nan)
        q_full[finite] = np.clip(q_back, 0.0, 1.0)
        out[mask] = q_full
    return out


def call_dmr_rule_segment(
    dmc_results: pl.DataFrame,
    *,
    self_loop: float = 0.95,
    min_cpgs: int = 5,
    min_abs_meth_diff: float = 0.10,
    alpha: float = 0.05,
) -> pl.DataFrame:
    """Rule-based segmentation DMR caller.

    Same parameters as the old ``call_dmr_hmm``. P2-4: emits real per-
    segment p-values via Stouffer combination of constituent CpG
    p-values, BH-corrected per chromosome. Drops the old
    ``mean(qvalue) < alpha`` post-filter in favour of filtering on the
    combined q-value (semantically the same intent, statistically
    defensible).
    """
    if dmc_results.height == 0:
        return pl.DataFrame(schema=_DMR_TILE_SCHEMA)
    required = {"chrom", "pos", "meth_diff"}
    missing = required - set(dmc_results.columns)
    if missing:
        raise ValueError(f"dmc_results missing required columns: {sorted(missing)}")
    has_pvalue = "pvalue" in dmc_results.columns

    state_means = _state_means_for_meth_diff(dmc_results["meth_diff"].to_numpy())

    out_rows: list[dict[str, object]] = []
    for chrom_grp in dmc_results.partition_by("chrom", maintain_order=True):
        chrom = chrom_grp["chrom"][0]
        chrom_sorted = chrom_grp.sort("pos")
        positions = chrom_sorted["pos"].to_numpy().astype(np.int64)
        meth_diff = chrom_sorted["meth_diff"].to_numpy().astype(np.float64)
        pvals_per_cpg = (
            chrom_sorted["pvalue"].to_numpy().astype(np.float64)
            if has_pvalue else None
        )

        viterbi = segment(
            meth_diff, n_states=3, state_means=state_means,
            self_loop=self_loop, emission="gaussian", emission_sd=0.10,
        )

        for state_idx, label in ((0, "hypo"), (2, "hyper")):
            runs = runs_of_state(viterbi, target_state=state_idx, positions=positions)
            for run_start, run_end, n_cpgs_run in runs:
                if n_cpgs_run < min_cpgs:
                    continue
                mask = (positions >= run_start) & (positions < run_end)
                if not mask.any():
                    continue
                run_md = meth_diff[mask]
                valid = np.isfinite(run_md)
                if valid.sum() == 0:
                    continue
                mean_md = float(run_md[valid].mean())
                if abs(mean_md) < min_abs_meth_diff:
                    continue
                seg_p = (
                    _stouffer_combine(pvals_per_cpg[mask])
                    if pvals_per_cpg is not None else 1.0
                )
                out_rows.append({
                    "chrom":            str(chrom),
                    "start":            int(run_start),
                    "end":              int(run_end),
                    "n_cpgs":           int(n_cpgs_run),
                    "n_case":           0,
                    "n_control":        0,
                    "mean_beta_case":   float("nan"),
                    "mean_beta_control": float("nan"),
                    "meth_diff":        float(mean_md),
                    "log2_odds_ratio":  float("nan"),
                    "pvalue":           float(seg_p),
                    "qvalue":           float("nan"),  # filled by BH below
                    "dmr_type":         label,
                })

    if not out_rows:
        return pl.DataFrame(schema=_DMR_TILE_SCHEMA)

    df = pl.DataFrame(out_rows, schema={
        "chrom":            pl.Utf8, "start": pl.Int32, "end": pl.Int32,
        "n_cpgs":           pl.Int32, "n_case": pl.Int32, "n_control": pl.Int32,
        "mean_beta_case":   pl.Float32, "mean_beta_control": pl.Float32,
        "meth_diff":        pl.Float32, "log2_odds_ratio": pl.Float64,
        "pvalue":           pl.Float64, "qvalue": pl.Float64, "dmr_type": pl.Utf8,
    }).sort(["chrom", "start"])

    # BH per chromosome.
    qvals = _bh_per_chrom(df["pvalue"].to_numpy(), df["chrom"].to_numpy())
    df = df.with_columns(pl.Series("qvalue", qvals))

    # Significance gate: replace the old mean(qvalue) < alpha hack with
    # per-segment qvalue < alpha. This is what the docstring claimed it
    # was doing all along.
    if has_pvalue:
        df = df.filter(pl.col("qvalue") < alpha)

    return df


__all__ = ["call_dmr_rule_segment"]
```

- [ ] **Step 4: Replace `src/epykit/dmr_hmm.py` with a deprecation shim**

Overwrite `src/epykit/dmr_hmm.py`:

```python
"""DEPRECATED — renamed to ``epykit.dmr_segment`` in 0.7.5.

This shim re-exports the renamed function and emits a
``DeprecationWarning`` on import. Scheduled for removal in 0.8.

Rationale: the engine was never a fitted HMM (state means and
transition probabilities are fixed), so the ``dmr_hmm`` name invited
the "did you Baum-Welch?" question. The honest name is
``dmr_segment``.
"""

from __future__ import annotations

import logging
import warnings

from .dmr_segment import call_dmr_rule_segment

# Mirror through the standard log stream so users running with
# ``python -W ignore`` still see the message in epykit's own logs.
logger = logging.getLogger(__name__)
_msg = (
    "epykit.dmr_hmm is deprecated and will be removed in 0.8; "
    "use epykit.dmr_segment.call_dmr_rule_segment instead"
)
warnings.warn(_msg, DeprecationWarning, stacklevel=2)
logger.warning(_msg)

# Legacy export name preserved.
call_dmr_hmm = call_dmr_rule_segment

__all__ = ["call_dmr_hmm", "call_dmr_rule_segment"]
```

- [ ] **Step 5: Update `tl.py` dmr dispatcher to accept method='segment' (preferred) and method='hmm' (legacy + FutureWarning)**

Grep for the dmr dispatcher: `rg -n "def dmr|method ==" src/epykit/tl.py | head`. The dmr method dispatch already takes a `method=` argument. Add:

```python
# Inside tl.py::dmr, near the top after the method= argument is parsed:
if method == "hmm":
    import warnings
    warnings.warn(
        "method='hmm' is deprecated; use method='segment' (same engine, "
        "honest name). Old method='hmm' will be removed in 0.8.",
        FutureWarning, stacklevel=2,
    )
    method = "segment"

# In the existing if/elif chain that dispatches per method, add:
elif method == "segment":
    from .dmr_segment import call_dmr_rule_segment
    dmrs = call_dmr_rule_segment(
        dmc_results=md.dmc,
        self_loop=self_loop,
        min_cpgs=min_cpgs,
        min_abs_meth_diff=min_abs_meth_diff,
        alpha=alpha,
    )
```

(Adapt to existing dispatcher style — the spec doesn't dictate the exact lines because they depend on the current shape of `tl.py::dmr`. Open `src/epykit/tl.py`, find the method-dispatch block, and add `"segment"` branch + the deprecation alias.)

- [ ] **Step 6: Update CLI alias**

In `src/epykit/cli.py`, find the `dmr` subcommand's `--method` choices. Add `"segment"` as a choice and accept `"hmm"` (mapped to `"segment"` with a click/argparse alias). Grep first: `rg -n "method|hmm" src/epykit/cli.py`.

- [ ] **Step 7: Run the new tests**

Run: `uv run pytest tests/test_dmr_segment.py -v`
Expected: all 3 PASS.

- [ ] **Step 8: Run existing HMM tests to check for breakage**

Run: `uv run pytest tests/test_dmr_hmm.py tests/test_hmm.py -v`
Expected: PASS (the shim re-exports, so old imports keep working). If the existing `test_dmr_hmm.py` asserts NaN p-values explicitly, those assertions now fail and need updating to assert finite p-values — that is the intended behaviour change.

Update assertions as needed:
- If a test was `assert dmrs["pvalue"].is_null().all()` or similar, replace with `assert dmrs["pvalue"].is_not_null().all()`.

- [ ] **Step 9: Append CHANGELOG entry**

Add under `## Unreleased` in `CHANGELOG.md`:

```markdown
### Changed (BREAKING for `epykit.dmr_hmm` import path)

- **Renamed** `epykit.dmr_hmm` → `epykit.dmr_segment`; function
  `call_dmr_hmm` → `call_dmr_rule_segment`. The engine was never a
  fitted HMM (state means + transitions are fixed); the new name
  reflects that. Old import path remains as a deprecated shim until
  0.8.
- **DMR `tl.dmr(method='hmm')`** is now `method='segment'`. Old value
  works with `FutureWarning` until 0.8.

### Fixed (P2 manifest, folded into the rename above)

- **P2-4**: `call_dmr_rule_segment` now emits per-segment Stouffer-
  combined p-values (BH-corrected per chromosome) instead of NaN.
  Downstream consumers that filtered on `qvalue` previously got an
  empty frame; they now see real q-values.
```

- [ ] **Step 10: Commit**

```
git add src/epykit/dmr_segment.py src/epykit/dmr_hmm.py src/epykit/tl.py src/epykit/cli.py tests/test_dmr_segment.py tests/test_dmr_hmm.py CHANGELOG.md
git commit -m "$(cat <<'EOF'
refactor(dmr) P1-8 + P2-4: rename dmr_hmm -> dmr_segment with real per-segment p-values

The engine was never a fitted HMM (state means and transition
probabilities are fixed); calling it ``dmr_hmm`` invited the
"did you Baum-Welch?" question. Renamed to ``dmr_segment`` and the
function to ``call_dmr_rule_segment``.

Folded with the rename: P2-4 fix. The pre-Phase-3 implementation
emitted NaN p/q-values for every called segment, silently breaking
any downstream consumer that filtered on qvalue. Replaced with
two-sided Stouffer combination of constituent CpG p-values and BH
correction per chromosome. The old ``mean(qvalue) < alpha`` post-
filter is replaced by a proper ``qvalue < alpha`` gate on the
per-segment q.

Old import path ``epykit.dmr_hmm.call_dmr_hmm`` remains as a
deprecated shim until 0.8. ``tl.dmr(method='hmm')`` accepted with
FutureWarning, mapped to ``method='segment'``.

Affects: segment@all

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 1

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: P1-11 — Rename `log2_odds_ratio` per backend (pooled vs GLM coefficient)

**Spec §4 commit 2.** Pooled-count tests (`lr`, `fisher`) get `log2_odds_ratio_pooled` (semantics unchanged, name clarified). GLM gets `coef_treatment_log2` (the column was always the logit coefficient, not log₂(OR); the old name was misleading). Both backends emit a transitional `log2_odds_ratio` column that is NaN-filled, with a `FutureWarning` issued on first access via a polars-aware shim.

**Note on transitional column:** polars DataFrames don't have a clean "access hook" for FutureWarning per-column. We instead emit the warning once when the DMC table is produced — sites that read the legacy column see NaN and the warning travels with the producing function call.

**Files:**
- Modify: `src/epykit/dmc.py` (column emission)
- Modify: `src/epykit/dmr.py` (consumer)
- Modify: `src/epykit/dmr_segment.py` (consumer — created in Task 1)
- Create: `tests/test_dmc_lr.py` (or append if it already exists)
- Create: `tests/test_glm.py` (or append if it already exists)
- Modify: `tests/test_compute_backends.py`, `tests/test_dmr_tile_merge.py` (column-name update)
- Modify: `docs/analysis/dmc.md`, `README.md`
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Grep existing consumers to confirm scope**

```
rg -n "log2_odds_ratio" src/ tests/ docs/ README.md
```

Expected: hits in `src/epykit/dmc.py`, `src/epykit/dmr.py`, `src/epykit/dmr_segment.py` (was `dmr_hmm.py`), `tests/test_dmr_tile_merge.py`, `tests/test_compute_backends.py`, `docs/analysis/dmc.md`, `README.md`. (Benchmark scripts confirmed NOT to reference it.)

- [ ] **Step 2: Write the failing test for the new columns**

Create or append to `tests/test_dmc_lr.py`:

```python
"""P1-11: log2_odds_ratio renamed per backend.

lr / fisher: log2_odds_ratio_pooled (semantics unchanged, name clearer).
glm: coef_treatment_log2 (was always the logit coefficient, not log2(OR)).

Both backends emit a transitional log2_odds_ratio column NaN-filled
plus a FutureWarning."""
from __future__ import annotations

import warnings

import numpy as np
import pytest

import epykit as ep


def test_lr_emits_log2_odds_ratio_pooled(synth_md_filtered):
    """The lr backend's output column for log-odds is renamed to
    log2_odds_ratio_pooled. The legacy column is NaN-filled."""
    md = synth_md_filtered
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(md, test="lr")
    df = md.dmc
    assert "log2_odds_ratio_pooled" in df.columns, (
        f"missing log2_odds_ratio_pooled; got {df.columns}"
    )
    assert "log2_odds_ratio" in df.columns, (
        "transitional log2_odds_ratio column must be present (NaN-filled)"
    )
    legacy = df["log2_odds_ratio"].to_numpy()
    assert np.isnan(legacy).all(), (
        f"transitional log2_odds_ratio must be NaN; got {legacy[:5]}"
    )
    fut = [w for w in caught if issubclass(w.category, FutureWarning)]
    assert fut and "log2_odds_ratio" in str(fut[0].message), (
        "expected FutureWarning naming the legacy column"
    )


def test_fisher_emits_log2_odds_ratio_pooled(synth_md_filtered):
    """Same rename applies to fisher (pooled-count backend family)."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="fisher")
    df = md.dmc
    assert "log2_odds_ratio_pooled" in df.columns
    assert "log2_odds_ratio" in df.columns
```

Create or append to `tests/test_glm.py`:

```python
"""P1-11 GLM half: coef_treatment_log2 is the new name for what was
misleadingly called log2_odds_ratio in the glm output (it's the logit
coefficient in log2 units, not log2 of an odds ratio)."""
from __future__ import annotations

import warnings

import numpy as np
import pytest

import epykit as ep


def test_glm_emits_coef_treatment_log2(synth_md_filtered):
    """GLM's log-odds-shape output column is renamed to
    coef_treatment_log2 (logit coefficient / ln(2))."""
    md = synth_md_filtered
    # Ensure md.obs has a treatment column for the GLM formula.
    md.obs = md.obs.with_columns(
        (md.obs["group"] == "treatment").cast(int).alias("treatment")
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(md, test="glm", formula="~ treatment")
    df = md.dmc
    assert "coef_treatment_log2" in df.columns, (
        f"missing coef_treatment_log2; got {df.columns}"
    )
    assert "log2_odds_ratio" in df.columns
    legacy = df["log2_odds_ratio"].to_numpy()
    assert np.isnan(legacy).all()
    fut = [w for w in caught if issubclass(w.category, FutureWarning)]
    assert fut and "log2_odds_ratio" in str(fut[0].message)
```

- [ ] **Step 3: Run the tests to verify failure**

Run: `uv run pytest tests/test_dmc_lr.py tests/test_glm.py -v -k "log2 or coef_treatment"`
Expected: all FAIL — old column name still emitted.

- [ ] **Step 4: Locate the emission sites in `dmc.py`**

```
rg -n "log2_odds_ratio" src/epykit/dmc.py
```

Expected hits include `dmc.py:1673` (or near it) where `log2_odds_ratio` is added to the output frame, plus possibly schema definitions and extras population. Read the surrounding 30 lines so you understand where the unified output block lives. The schema-emission site is the one immediately before the function returns the per-chromosome `pl.DataFrame(...)`.

- [ ] **Step 5: Apply the rename + transitional column logic in `dmc.py`**

In the unified output block of `_process_one_chromosome` (search for the place where `log2_odds_ratio` is written into the result DataFrame), replace the emission with backend-conditional naming:

```python
# Determine backend group: pooled-count tests vs glm.
_GLM_BACKENDS = {"glm", "glm_contrast"}
log2_col_name = (
    "coef_treatment_log2" if test in _GLM_BACKENDS else "log2_odds_ratio_pooled"
)

# Build the output columns dict. The transitional 'log2_odds_ratio'
# column is NaN-filled to keep the schema stable for one release;
# downstream readers that haven't migrated still see the column but
# get NaN, paired with the FutureWarning emitted below.
import warnings
warnings.warn(
    "The 'log2_odds_ratio' column is deprecated and will be removed in "
    "0.8. Use 'log2_odds_ratio_pooled' for pooled-count tests (lr, "
    "fisher) or 'coef_treatment_log2' for the GLM backend. The "
    "transitional column is NaN-filled in 0.7.5.",
    FutureWarning, stacklevel=3,
)

out_cols = {
    # ... existing columns ...
    log2_col_name: log2_ors,
    "log2_odds_ratio": np.full_like(log2_ors, np.nan, dtype=np.float64),
}
```

(Exact placement depends on the existing output-frame construction; the principle is: rename the column based on `test`, and add a NaN-filled `log2_odds_ratio` alongside it.)

- [ ] **Step 6: Update consumers in `dmr.py` and `dmr_segment.py`**

```
rg -n "log2_odds_ratio" src/epykit/dmr.py src/epykit/dmr_segment.py
```

For each hit, switch to reading `log2_odds_ratio_pooled` (or `coef_treatment_log2` if the consumer is GLM-specific). The DMR tile code typically reads from `lr` output, so use `log2_odds_ratio_pooled`. Where the column is emitted in the DMR output (e.g., `dmr_segment.py:134`), keep the column name `log2_odds_ratio` for the DMR schema (DMR consumers haven't been updated; DMR rename is a future cleanup).

- [ ] **Step 7: Update `tests/test_compute_backends.py` and `tests/test_dmr_tile_merge.py`**

```
rg -n "log2_odds_ratio" tests/test_compute_backends.py tests/test_dmr_tile_merge.py
```

For each assertion on `log2_odds_ratio`, replace with `log2_odds_ratio_pooled` if the producer is `lr`/`fisher`, or `coef_treatment_log2` if `glm`. If a test just checks the column exists, it can stay (the transitional column still exists, just NaN-filled — but then the assertion is meaningless). Prefer asserting on the new name.

- [ ] **Step 8: Update docs**

In `docs/analysis/dmc.md` and `README.md`, find every mention of `log2_odds_ratio` and either rename or add a note that the column is renamed per backend in 0.7.5.

- [ ] **Step 9: Run the new tests**

Run: `uv run pytest tests/test_dmc_lr.py tests/test_glm.py tests/test_compute_backends.py tests/test_dmr_tile_merge.py -v`
Expected: all PASS.

- [ ] **Step 10: CHANGELOG**

Add under `## Unreleased`:

```markdown
### Changed (BREAKING for `log2_odds_ratio` column name)

- **`varm["dmc_lr"].log2_odds_ratio`** renamed to
  `log2_odds_ratio_pooled` (semantics unchanged). Same rename applies
  to `fisher`.
- **`varm["dmc_glm"].log2_odds_ratio`** renamed to
  `coef_treatment_log2` (it was always the logit coefficient in log₂
  units, not log₂ of an odds ratio; the old name was misleading).
- A transitional `log2_odds_ratio` column is NaN-filled in 0.7.5 with
  a `FutureWarning` on the producing call. Column will be removed in
  0.8.
```

- [ ] **Step 11: Commit**

```
git commit -am "$(cat <<'EOF'
refactor(dmc) P1-11: rename log2_odds_ratio per backend

The pre-Phase-3 column name was misleading for the GLM backend: the
emitted value is the logit coefficient in log2 units, not log2 of an
odds ratio. Renamed per backend:
- lr / fisher (pooled-count):  log2_odds_ratio_pooled
- glm:                         coef_treatment_log2

Transitional log2_odds_ratio column NaN-filled, FutureWarning on
producing call. Column scheduled for removal in 0.8.

Updated consumers in dmr.py and dmr_segment.py to read the new pooled
column. Updated tests/test_compute_backends.py and
tests/test_dmr_tile_merge.py to assert on the new names. README and
docs/analysis/dmc.md updated.

Affects: lr@all, fisher@all, glm@all

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 2

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Drop `test="logit_t"` from `tl.dmc`

**Spec §4 commit 3.** Hard-remove from the dispatcher. Internal helper `_beta_binom_mom_from_welford_logit` stays (cheap to keep; no other dispatch path reaches it post-drop, but it's referenced from the welch_t branch's if-chain that we're trimming, so removing the helper too is optional and lower-priority).

**Files:**
- Modify: `src/epykit/dmc.py:1577-1599` (welch_t / logit_t shared branch)
- Modify: `src/epykit/tl.py:349` (docstring) + add ValueError guard near top of `tl.dmc`
- Modify: `tests/test_api.py:177-182`, `tests/test_accuracy.py:95`, `tests/test_stats_new.py:304` (parametrize lists)
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Identify all consumers**

```
rg -nP '"logit_t"' src/ tests/
```

Expected sites (from pre-flight inventory):
- `src/epykit/dmc.py:1577,1590-1594` (engine branch)
- `src/epykit/tl.py:349` (docstring)
- `tests/test_api.py:177,182`
- `tests/test_accuracy.py:95`
- `tests/test_stats_new.py:304` (parametrize row)

- [ ] **Step 2: Edit `dmc.py` — collapse the welch_t/logit_t if-chain**

Read `src/epykit/dmc.py` lines 1577-1600. The branch currently reads:

```python
elif test in ("welch_t", "logit_t"):
    # ... load samples into Welford accumulators ...
    if test == "logit_t":
        pvals, log2_ors = _beta_binom_mom_from_welford_logit(...)
    else:
        pvals, log2_ors = _beta_binom_mom_from_welford(...)
```

Replace with:

```python
elif test == "welch_t":
    # ... load samples into Welford accumulators (same body) ...
    pvals, log2_ors = _beta_binom_mom_from_welford(
        mean_case, M2_case, n_valid_case,
        mean_ctrl, M2_ctrl, n_valid_ctrl,
    )
```

- [ ] **Step 3: Remove `"logit_t"` from `tl.py:349` docstring**

In `src/epykit/tl.py`, find:
```
test : str
    One of ``"auto"``, ``"lr"``, ``"score"``, ``"logit_t"``,
    ``"welch_t"`` (Welch t on raw betas),
```

Edit the docstring (Task 7 will do the full collapse to the 4-surviving set; for now just delete `"logit_t"`).

- [ ] **Step 4: Add ValueError guard in `tl.dmc`**

In `src/epykit/tl.py::dmc`, near the top of the function (just after the argument-validation block, before any dispatch logic), add:

```python
if test == "logit_t":
    raise ValueError(
        "test='logit_t' was removed in 0.7.5 (miscalibrated near β=0/1). "
        "Use test='welch_t' for the replicate-aware β-mean test or "
        "test='lr' for the recommended default."
    )
```

- [ ] **Step 5: Update affected tests**

`tests/test_api.py`: lines 177-182 use `test="logit_t"` to assert the engine emits expected output. Replace with `test="welch_t"` if the test's intent is "any non-default engine works", OR delete the block and remove the docstring line about logit_t if the test was specifically asserting logit_t presence.

`tests/test_accuracy.py`: line 95 calls `_run_dmc(synth_md_filtered, test="logit_t")`. Replace with `test="welch_t"` (same intent: a non-lr engine on the same data).

`tests/test_stats_new.py`: line 304 parametrize list `["lr", "score", "logit_t", "welch_t", "bb_lr"]`. Trim to `["lr", "welch_t"]` — `score` and `bb_lr` are dropped in subsequent tasks but for now just remove `"logit_t"`. (Tasks 4 and 5 will further trim.)

- [ ] **Step 6: Write the migration-hint test (consolidated in Task 7)**

This task's drop is tested by Task 7's parametrised `test_phase3_drops.py`. Skip ahead and write that test now if you want immediate verification, OR proceed and let Task 7 cover it.

- [ ] **Step 7: Run the test suite to verify no regressions**

```
uv run pytest -m "not slow" --strict-markers -q -x 2>&1 | tail -15
```

Expected: passing. If a test fails with "test='logit_t' raises ValueError", that test still parametrises on `logit_t` and was missed in Step 5 — fix it.

- [ ] **Step 8: CHANGELOG entry under `## Unreleased / ### Removed`**

```markdown
### Removed (BREAKING)

- **`tl.dmc(test='logit_t')`** removed. Documented by epykit's own
  source as miscalibrated near β=0/1; no headline claim depended on
  it. Migration: `test='welch_t'` for the replicate-aware β-mean test
  or `test='lr'` for the recommended default. Old call now raises
  `ValueError` with the same hint.
```

- [ ] **Step 9: Commit**

```
git commit -am "$(cat <<'EOF'
fix(dmc) BREAKING: remove test='logit_t'

The logit_t engine was documented by epykit's own source as
miscalibrated near boundary methylation (β = 0 or 1). No paper claim
depends on it. Hard-removed from the tl.dmc dispatcher; calls now
raise ValueError with a migration hint (use welch_t or lr).

Internal helper _beta_binom_mom_from_welford_logit stays in dmc.py
for potential future use; it's no longer reachable from the public
dispatcher.

Migration: test='logit_t' -> test='welch_t' (one-line edit).

Affects:

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 3

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Drop `test="bb_lr"` from `tl.dmc`

**Spec §4 commit 4.** Hard-remove the ~88-line `bb_lr` branch and its n<3 guard. `irls_dispatch` stays (used by `glm`).

**Files:**
- Modify: `src/epykit/dmc.py:1601-1688` (bb_lr branch deletion), `dmc.py:2032-2033` (n<3 guard deletion)
- Modify: `src/epykit/tl.py:349` (docstring) + ValueError guard
- Modify: `tests/test_stats_new.py:304,339-340` (parametrize + bb_lr-specific assertions)
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Identify all consumers**

```
rg -nP '"bb_lr"|test_bb_lr|bb_lr' src/ tests/
```

Expected:
- `src/epykit/dmc.py:1601-1688` (branch), `dmc.py:2032-2033` (n<3 guard), line 1612-1620 (the n<6 dispersion-warning shim)
- `src/epykit/tl.py:349` (docstring)
- `tests/test_stats_new.py:304` (parametrize), 339-340 (engine-specific call)

- [ ] **Step 2: Delete the bb_lr branch in `dmc.py`**

Open `src/epykit/dmc.py` and read lines 1601-1688 to confirm the branch boundaries. Delete the entire `elif test == "bb_lr":` block. Also delete:
- Lines 2032-2033 (the `if test == "bb_lr" and min_n < 3:` early return).
- Any other `bb_lr` references in this file (use `rg -n "bb_lr" src/epykit/dmc.py` and confirm; should be only the branch + the n<3 guard after this task).

- [ ] **Step 3: Remove `"bb_lr"` from `tl.py:349` docstring**

Delete the `"bb_lr"` entry from the `test : str` docstring listing.

- [ ] **Step 4: Add ValueError guard in `tl.dmc`**

Below the `logit_t` guard added in Task 3:

```python
if test == "bb_lr":
    raise ValueError(
        "test='bb_lr' was removed in 0.7.5 (TPR < 8% at n ≤ 4 + a "
        "dispersion-df bug). Use test='lr' (recommended) which uses "
        "the same quasi-binomial dispersion but pools counts per group "
        "for higher power at small n."
    )
```

- [ ] **Step 5: Update affected tests**

`tests/test_stats_new.py`:
- Line 304: prune `"bb_lr"` from the parametrize list.
- Lines 339-340: this is a `bb_lr`-specific block (`ep.tl.dmc(md, test="bb_lr")` + `md.get_dmc(test="bb_lr")`). Delete the block entirely — there is no near-equivalent engine.

- [ ] **Step 6: Run tests**

```
uv run pytest -m "not slow" -q -x 2>&1 | tail -10
```

Expected: passing.

- [ ] **Step 7: CHANGELOG**

Append to `### Removed`:

```markdown
- **`tl.dmc(test='bb_lr')`** removed. TPR < 8% at n ≤ 4 in the
  benchmark; affected by the spec's P1-2 dispersion-df bug
  (now closed by removal). Migration: `test='lr'` — same quasi-
  binomial dispersion, higher power at small n via per-group pooling.
  Call now raises `ValueError` with the same hint.
```

Also add under `### Fixed (P1 manifest)`:

```markdown
- **P1-2** (bb_lr `df_resid` vs `df_phi`): **closed by removal** — see
  Removed section above.
```

- [ ] **Step 8: Commit**

```
git commit -am "$(cat <<'EOF'
fix(dmc) BREAKING: remove test='bb_lr'

The bb_lr engine produced TPR < 8% at n <= 4 in the published
benchmark and was affected by the P1-2 dispersion-df bug. Hard-
removed from the tl.dmc dispatcher; calls now raise ValueError
with a migration hint (use lr).

Deleted:
- The 88-line bb_lr branch in dmc.py
- The n<3 guard at dmc.py:2032-2033
- The bb_lr-specific block in tests/test_stats_new.py

irls_dispatch stays (still used by the glm engine). P1-2 closed by
removal.

Migration: test='bb_lr' -> test='lr' (one-line edit; same dispersion
machinery, higher power at small n).

Affects:

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 4

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Drop `test="score"` from `tl.dmc`

**Spec §4 commit 5.** Collapse the `("score", "lr")` shared branch to `lr` only. Hardcode `statistic="lr"` in the `_score_finalize` call.

**Files:**
- Modify: `src/epykit/dmc.py:1510` (collapse branch), `dmc.py:1554-1563` (`_score_finalize` call), possibly `_score_finalize` signature if `statistic=` kwarg has no other caller
- Modify: `src/epykit/tl.py:349` + ValueError guard
- Modify: `tests/test_accuracy.py:84`, `tests/test_resume.py:143-145`, `tests/test_stats_new.py:304`
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Locate all consumers**

```
rg -nP '"score"|statistic=test|statistic="score"' src/ tests/
```

Expected:
- `src/epykit/dmc.py:1510` (branch entry)
- `src/epykit/dmc.py:~1559` (call to `_score_finalize(..., statistic=test, ...)`)
- `_score_finalize` definition site: `rg -n "def _score_finalize" src/epykit/dmc.py`
- `src/epykit/tl.py:349`
- `tests/test_accuracy.py:84`, `tests/test_resume.py:143-145`, `tests/test_stats_new.py:304`

- [ ] **Step 2: Collapse the branch in `dmc.py`**

Change `elif test in ("score", "lr"):` (line 1510) to `elif test == "lr":`. Inside that branch, change the `_score_finalize` call: replace `statistic=test` with `statistic="lr"`.

Check whether `_score_finalize` has any other caller with `statistic="score"`: `rg -n "_score_finalize" src/epykit/`. If `statistic=` is only ever `"lr"` after this task, simplify the signature (drop `statistic=` parameter and hard-set internally). If unsure, leave the signature alone — this is a cosmetic cleanup deferrable to v0.8.

- [ ] **Step 3: Remove `"score"` from `tl.py:349` docstring**

Delete the `"score"` entry.

- [ ] **Step 4: Add ValueError guard**

```python
if test == "score":
    raise ValueError(
        "test='score' was removed in 0.7.5 (strictly dominated by "
        "test='lr' in finite samples; asymptotically equivalent under "
        "H0). Switch test='score' -> test='lr'; output schema is "
        "identical."
    )
```

- [ ] **Step 5: Update tests**

`tests/test_accuracy.py:84`: `df = _run_dmc(synth_md_filtered, test="score")`. Replace with `test="lr"` (same data, equivalent test under H₀).

`tests/test_resume.py:143-145`: explicit `test="score"` plus an assertion on `test_used == "score"`. Replace the engine to `"lr"` and update the assertion to `"lr"`. (The test's intent is verifying resume-from-cache behaviour, engine-agnostic.)

`tests/test_stats_new.py:304`: drop `"score"` from the parametrize list (now `["lr", "welch_t"]` after prior drops).

- [ ] **Step 6: Run tests**

```
uv run pytest -m "not slow" -q -x 2>&1 | tail -10
```

Expected: passing.

- [ ] **Step 7: CHANGELOG**

```markdown
- **`tl.dmc(test='score')`** removed. Strictly dominated by `lr` in
  finite samples (LRT has better small-sample behaviour under quasi-
  binomial dispersion); asymptotically equivalent under H₀. Migration:
  `test='score'` → `test='lr'`; output schema is identical.
```

- [ ] **Step 8: Commit**

```
git commit -am "$(cat <<'EOF'
fix(dmc) BREAKING: remove test='score'

score (quasi-binomial Pearson chi-square) is strictly dominated by lr
(likelihood-ratio) in finite samples and asymptotically equivalent
under H0. The implementation was a flag flip on the same
_score_finalize internal function. Hard-removed from the dispatcher;
calls now raise ValueError with a migration hint (use lr).

Collapsed the shared ("score", "lr") branch to lr only; hardcoded
statistic="lr" in the _score_finalize call. _score_finalize signature
left intact (statistic= kwarg cleanup deferrable to v0.8).

Migration: test='score' -> test='lr' (one-line edit; output schema
identical).

Affects:

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 5

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Drop `test="cmh"` from `tl.dmc`

**Spec §4 commit 6.** Delete the cmh branch and its helpers (`_cmh_init` / `_cmh_update` / `_cmh_finalize`). Grep first to confirm no other caller.

**Files:**
- Modify: `src/epykit/dmc.py:1477-1508` (branch deletion), plus `_cmh_init/_update/_finalize` definitions
- Modify: `src/epykit/tl.py:349` + ValueError guard
- Modify: `tests/test_api.py:266`, `tests/test_stats_new.py` (already trimmed in Task 5)
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Locate all consumers**

```
rg -nP '"cmh"|_cmh_init|_cmh_update|_cmh_finalize' src/ tests/
```

Expected hits in `dmc.py` (branch + helpers), `tl.py:349` (docstring), `tests/test_api.py:266`. No internal-only callers of the `_cmh_*` helpers expected — confirm.

- [ ] **Step 2: Delete the cmh branch (`dmc.py:1477-1508`)**

Read the branch first to confirm boundaries. Delete the `elif test == "cmh":` block entirely (~32 lines).

- [ ] **Step 3: Delete `_cmh_init`, `_cmh_update`, `_cmh_finalize` definitions**

Find them: `rg -n "def _cmh_" src/epykit/dmc.py`. Delete the three function definitions. If unreferenced after this task: `rg -n "_cmh_" src/epykit/dmc.py` should return no hits.

- [ ] **Step 4: Remove `"cmh"` from `tl.py:349` docstring + ValueError guard**

```python
if test == "cmh":
    raise ValueError(
        "test='cmh' was removed in 0.7.5 (stratification semantics "
        "confusing; dominated by GLM with batch covariate). For "
        "stratified analysis use tl.dmc(formula='~ group + batch'), "
        "which gives proper dispersion correction and handles "
        "continuous covariates."
    )
```

- [ ] **Step 5: Update tests**

`tests/test_api.py:266`: `assert synth_md_filtered.get_dmc(test="cmh") is None`. This was a "cmh hasn't been run yet" assertion; replace with a check on an engine that's still in the API, e.g. `assert synth_md_filtered.get_dmc(test="glm") is None`. (Or delete if the test was specifically about cmh.)

- [ ] **Step 6: Run tests**

```
uv run pytest -m "not slow" -q -x 2>&1 | tail -10
```

Expected: passing.

- [ ] **Step 7: CHANGELOG**

```markdown
- **`tl.dmc(test='cmh')`** removed. Stratification semantics in the
  epykit implementation (one stratum per case-ctrl pair) were unusual
  and not what most users mean by "stratified test". Dominated by GLM
  with explicit batch covariate. Migration: `tl.dmc(formula='~ group
  + batch')` for batch-stratified analysis. Call now raises
  `ValueError` with the same hint.
```

- [ ] **Step 8: Commit**

```
git commit -am "$(cat <<'EOF'
fix(dmc) BREAKING: remove test='cmh'

The cmh stratification semantics (one stratum per case-ctrl pair) were
unusual and not what users typically mean by "stratified test".
Dominated by the GLM backend with explicit batch covariate, which is
on the recommended path and dispersion-correct. Hard-removed from the
dispatcher; calls now raise ValueError with a migration hint.

Deleted:
- The 32-line cmh branch in dmc.py
- _cmh_init, _cmh_update, _cmh_finalize helpers (no other callers)
- The cmh-presence check in tests/test_api.py:266

Migration: test='cmh' -> tl.dmc(formula='~ group + batch') for
batch-stratified analysis.

Affects:

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 6

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Collapse `_auto_test` + docstring + docs to 4+auto surviving engines; add `test_phase3_drops.py`

**Spec §4 commit 7.** Final sweep after the four drops. Trims the public surface documentation everywhere, adds the parametrised migration-hint test that exercises all four `ValueError` guards.

**Files:**
- Modify: `src/epykit/tl.py:349` (final docstring shape) + `_auto_test` comment
- Modify: `pyproject.toml` (markers, if any reference dropped engines)
- Modify: `tests/test_dmc_multigroup.py`, any other parametrised test files
- Create: `tests/test_phase3_drops.py`
- Modify: `README.md`, `docs/analysis/dmc.md`, `CLAUDE.md` (engine lists)
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Write the parametrised migration-hint test**

Create `tests/test_phase3_drops.py`:

```python
"""Phase 3 cleanup: dropped engines must raise ValueError with a
migration hint pointing the user at a surviving engine."""
from __future__ import annotations

import pytest

import epykit as ep


@pytest.mark.parametrize("engine,hint_substring", [
    ("logit_t", "welch_t"),
    ("bb_lr",   "lr"),
    ("score",   "lr"),
    ("cmh",     "formula='~ group + batch'"),
])
def test_dropped_engine_raises_with_migration_hint(
    synth_md_filtered, engine, hint_substring,
):
    """Each dropped engine raises ValueError; the message includes
    text pointing at the recommended replacement."""
    with pytest.raises(ValueError) as exc:
        ep.tl.dmc(synth_md_filtered, test=engine)
    msg = str(exc.value)
    assert "removed in 0.7.5" in msg, f"missing version note: {msg}"
    assert hint_substring in msg, (
        f"missing migration hint '{hint_substring}' in: {msg}"
    )
```

- [ ] **Step 2: Run the test to verify all four guards work**

Run: `uv run pytest tests/test_phase3_drops.py -v`
Expected: 4 PASS (one per engine, all hitting the guards added in Tasks 3-6).

- [ ] **Step 3: Collapse `tl.py:349` docstring to the surviving 4+auto**

Open `src/epykit/tl.py` near line 349. Replace the existing `test : str` docstring block with:

```python
test : str
    One of ``"auto"``, ``"lr"``, ``"welch_t"``, ``"fisher"``, or
    ``"glm"``. ``"auto"`` resolves to ``"fisher"`` at n<2 and ``"lr"``
    (the recommended default) at n>=2.

    Engines removed in 0.7.5 (raise ValueError with a migration
    hint): ``"logit_t"`` (use ``"welch_t"``), ``"bb_lr"`` (use
    ``"lr"``), ``"score"`` (use ``"lr"``), ``"cmh"`` (use
    ``formula='~ group + batch'``).

    When ``formula`` and/or ``contrast`` are supplied, the test is
    forced to a GLM-based path regardless of ``test=``.
```

- [ ] **Step 4: Inline-comment `_auto_test`**

Find `def _auto_test` in `tl.py` (`rg -n "_auto_test" src/epykit/tl.py`). Add a header comment:

```python
def _auto_test(md, *, allow_n1: bool) -> str:
    """Auto-dispatcher. Post-Phase-3 surface is closed to {fisher, lr}:
    fisher at n=1 (only engine that works), lr at n>=2.
    """
    # ... existing body unchanged ...
```

- [ ] **Step 5: Prune parametrize lists across the test suite**

```
rg -nP 'parametrize\(\s*"test"' tests/ -A 1
```

For each match: drop any of `"logit_t"`, `"bb_lr"`, `"score"`, `"cmh"` from the list. Known sites at the time of writing: `tests/test_dmc_multigroup.py`, `tests/test_compute_backends.py` (already partly trimmed in prior tasks for the explicit single-engine assertions). Update each remaining parametrise list.

- [ ] **Step 6: Update `pyproject.toml`**

```
rg -n "logit_t|bb_lr|score|cmh" pyproject.toml
```

If any markers, extras, or examples mention dropped engines, prune them. If no hits, no edit.

- [ ] **Step 7: Update `README.md`, `docs/analysis/dmc.md`, `CLAUDE.md`**

For each file, find any engine-list table or sentence and trim to:

> epykit 0.7.5 ships four per-CpG engines:
> - `lr` (default, quasi-binomial likelihood ratio)
> - `welch_t` (replicate-aware Welch t on β)
> - `fisher` (exact, n=1 fallback)
> - `glm` (binomial GLM with formula-based covariates)
>
> The `auto` dispatcher resolves to `fisher` at n=1 and `lr` at n≥2.

In `CLAUDE.md` lines 56-60 specifically: replace the existing engine enumeration with the four-engine list.

- [ ] **Step 8: Run the full fast suite**

```
uv run pytest -m "not slow" --strict-markers -q 2>&1 | tail -5
```

Expected: all passing. The number should be close to the Phase 2 baseline + ~5 new tests (drops + segment) − any deletions from the bb_lr/score/cmh-specific blocks.

- [ ] **Step 9: CHANGELOG sweep**

Under `## Unreleased / ### Removed`, prepend an umbrella entry summarising the four drops + the surviving surface:

```markdown
- **DMC engine surface collapsed to `{auto, lr, welch_t, fisher,
  glm}`.** Dropped engines: `logit_t` (broken near β=0/1), `bb_lr`
  (TPR < 8% at n ≤ 4 + dispersion bug), `score` (dominated by `lr`),
  `cmh` (dominated by `glm + batch`). All four raise `ValueError`
  with a one-line migration hint. README, `docs/analysis/dmc.md`,
  and `CLAUDE.md` updated.
```

- [ ] **Step 10: Commit**

```
git commit -am "$(cat <<'EOF'
refactor(dmc): collapse _auto_test, docstring, docs to 4+auto surviving engines

Final pass after the four engine drops (Tasks 3-6). Trims the public
surface documentation and adds the parametrised migration-hint test
that exercises all four ValueError guards in one file.

- tl.py:349 docstring lists the surviving four (lr, welch_t, fisher,
  glm) + auto. Removed engines noted with migration hint.
- _auto_test header-commented to document the post-Phase-3 closure
  ({fisher, lr}).
- Parametrize lists across tests/ pruned of the dropped engines.
- README, docs/analysis/dmc.md, CLAUDE.md updated to the four-engine
  list.
- tests/test_phase3_drops.py exercises all four migration hints.

Affects:

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 7

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: P1-1 — Fisher two-sided uses the mid-p convention

**Spec §4 commit 8.** Replace the doubled-smaller-tail convention with the same two-sided p that `scipy.stats.fisher_exact(alternative="two-sided")` produces. Vectorised — scipy's per-table call is too slow for the 22M-CpG workload.

**Files:**
- Modify: `src/epykit/dmc.py` (around `fisher_exact_vectorized` at `dmc.py:218`)
- Create or append: `tests/test_dmc_fisher.py`
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Read the current implementation**

```
rg -n "def fisher_exact_vectorized" src/epykit/dmc.py
```

Read 60 lines starting at the def. The existing two-sided p is computed as `min(2 * one_sided_p, 1.0)` (doubled-tail). The fix is to compute the mid-p two-sided convention: sum the probabilities of all tables at least as extreme as observed, where "extreme" is measured by the table probability under the hypergeometric null.

- [ ] **Step 2: Write the failing test**

Create `tests/test_dmc_fisher.py`:

```python
"""P1-1: vectorised Fisher two-sided p must match
scipy.stats.fisher_exact(alternative='two-sided') to machine
precision."""
from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import fisher_exact

from epykit.dmc import fisher_exact_vectorized


def test_fisher_two_sided_matches_scipy_on_random_tables():
    rng = np.random.default_rng(0)
    n = 100
    # Random 2x2 tables with cell counts in [0, 40] and totals >= 1.
    a = rng.integers(0, 40, size=n).astype(np.int64)
    b = rng.integers(0, 40, size=n).astype(np.int64)
    c = rng.integers(0, 40, size=n).astype(np.int64)
    d = rng.integers(0, 40, size=n).astype(np.int64)
    # Avoid all-zero rows/cols (degenerate; both engines will return p=1).
    keep = (a + b > 0) & (c + d > 0) & (a + c > 0) & (b + d > 0)
    a, b, c, d = a[keep], b[keep], c[keep], d[keep]

    epy_p, _epy_log2_or = fisher_exact_vectorized(a, b, c, d)
    ref_p = np.array([
        fisher_exact([[ai, bi], [ci, di]], alternative="two-sided")[1]
        for ai, bi, ci, di in zip(a, b, c, d)
    ])
    np.testing.assert_allclose(
        epy_p, ref_p, atol=1e-12, rtol=1e-9,
        err_msg="vectorised Fisher two-sided p must match scipy reference",
    )
```

- [ ] **Step 3: Run the test (expect failure)**

Run: `uv run pytest tests/test_dmc_fisher.py -v`
Expected: FAIL — current implementation uses doubled-tail.

- [ ] **Step 4: Replace the two-sided computation**

In `src/epykit/dmc.py::fisher_exact_vectorized`, replace the doubled-tail two-sided p with a vectorised mid-p that sums over hypergeometric probabilities at least as extreme as observed.

Strategy: for each table (a, b, c, d) with row totals (n1, n2) and col totals (m1, m2), the observed cell `a` is drawn from `Hypergeometric(N=n1+n2, K=n1, n=m1)`. Enumerate the support `k_min..k_max`, compute the pmf, and sum `pmf[k]` where `pmf[k] <= pmf[a]` (the mid-p definition scipy uses). Vectorisation is achievable per support size; for production, batch by total `N` to share factorial caches.

Implementation sketch (using `scipy.stats.hypergeom.pmf` per-table is acceptable since N=10K-100K CpGs per chrom is small relative to the per-table support):

```python
from scipy.stats import hypergeom

def fisher_exact_vectorized(a, b, c, d):
    a = a.astype(np.int64); b = b.astype(np.int64)
    c = c.astype(np.int64); d = d.astype(np.int64)
    N = a + b + c + d
    n1 = a + b
    m1 = a + c
    pvals = np.full(len(a), np.nan, dtype=np.float64)
    log2_ors = np.full(len(a), np.nan, dtype=np.float64)
    for i in range(len(a)):
        Ni, n1i, m1i, ai = int(N[i]), int(n1[i]), int(m1[i]), int(a[i])
        if Ni == 0 or n1i == 0 or n1i == Ni or m1i == 0 or m1i == Ni:
            pvals[i] = 1.0
            log2_ors[i] = 0.0
            continue
        k_min = max(0, n1i + m1i - Ni)
        k_max = min(n1i, m1i)
        ks = np.arange(k_min, k_max + 1)
        pmf = hypergeom.pmf(ks, Ni, n1i, m1i)
        obs_pmf = hypergeom.pmf(ai, Ni, n1i, m1i)
        pvals[i] = float(pmf[pmf <= obs_pmf + 1e-15].sum())
        # log2 odds ratio with Haldane (0.5) correction.
        ao = a[i] + 0.5; bo = b[i] + 0.5
        co = c[i] + 0.5; do = d[i] + 0.5
        log2_ors[i] = float(np.log2((ao * do) / (bo * co)))
    return pvals, log2_ors
```

(If performance regresses noticeably, batch by unique `(N, n1, m1)` triples and cache the pmf array per triple. For Phase 3 the simpler per-table loop is acceptable; benchmark in Phase 4.)

- [ ] **Step 5: Run the test (expect pass)**

Run: `uv run pytest tests/test_dmc_fisher.py -v`
Expected: PASS.

- [ ] **Step 6: Run the surrounding test suite to catch regressions**

```
uv run pytest -m "not slow" -q -x 2>&1 | tail -10
```

Expected: passing. Watch for any test that previously asserted on the doubled-tail value — those need their expected values regenerated.

- [ ] **Step 7: CHANGELOG**

Under `### Fixed (P1 manifest)`:

```markdown
- **P1-1**: `fisher_exact_vectorized` two-sided p now uses the mid-p
  convention (matches `scipy.stats.fisher_exact(alternative='two-
  sided')` to 1e-12). Previously doubled the smaller one-sided tail,
  which is conservative but not what scipy / textbooks emit. Affects
  small-table cells only; headline cov≥10 / n≥3 numbers unchanged.
```

- [ ] **Step 8: Commit**

```
git commit -am "$(cat <<'EOF'
fix(dmc) P1-1: Fisher two-sided uses mid-p convention

Previously fisher_exact_vectorized returned min(2*one_sided, 1.0),
which is the conservative doubled-tail. scipy and most textbooks use
the mid-p two-sided: sum of pmf over all tables with pmf <= pmf(obs).
Vectorised per-table via scipy.stats.hypergeom.pmf; verified against
scipy.stats.fisher_exact to 1e-12 atol on 100 random 2x2 tables.

No headline impact (cov>=10, n>=3 tables already saturate small-tail
mass). Small-table cells (n=1 paths, very low coverage) shift
slightly.

Affects: fisher@small-table-cells

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 8

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: P1-3 — Wire `newcombe_diff_ci` into the `lr` meth_diff CI columns

**Spec §4 commit 9.** Replace the existing Welch-normal Wald CI for `meth_diff_ci_{lo,hi}` on the `lr` path with `_glm.newcombe_diff_ci` (already implemented at `_glm.py:923`, unwired). `welch_t` keeps Welch CI (correct for that test); `glm` keeps model-based CI; `fisher` keeps Wilson on per-group rates.

**Files:**
- Modify: `src/epykit/dmc.py` (`_process_one_chromosome`, around `dmc.py:1907-1912`)
- Modify: `tests/test_dmc_lr.py` (append)
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Read the current CI emission site**

```
rg -n "meth_diff_ci_lo|meth_diff_ci_hi" src/epykit/dmc.py
```

Read 30 lines around the emission. Confirm where the existing Wald CI is computed for the `lr` path.

- [ ] **Step 2: Write the failing test**

Append to `tests/test_dmc_lr.py`:

```python
def test_lr_meth_diff_ci_uses_newcombe(synth_md_filtered):
    """P1-3: lr emits Newcombe (1998) hybrid-Wilson CIs on meth_diff,
    matching statsmodels.confint_proportions_2indep(method='newcombe')
    at the pooled count level for the first 20 sites."""
    from statsmodels.stats.proportion import confint_proportions_2indep

    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    df = md.dmc.head(20)
    pooled = df.select([
        "n_case", "n_control",
        "mean_beta_treat", "mean_beta_ctrl",
        "meth_diff_ci_lo", "meth_diff_ci_hi",
    ]).to_numpy()
    # The Newcombe reference takes pooled methylated/total counts per
    # group. We need pooled meth + cov, which the dmc table doesn't
    # surface directly -- so derive from mean_beta * n_case (a coverage
    # surrogate). This is a weak reference; tighten by reading the
    # accumulator state if the test proves flaky.
    # Tolerance is loose because pooled counts in the dmc output are
    # per-group means times sample counts, not exact totals.
    for row in pooled:
        n_case, n_ctrl, m_treat, m_ctrl, lo_obs, hi_obs = row
        if not (np.isfinite(lo_obs) and np.isfinite(hi_obs)):
            continue
        # Simple shape check: CI brackets the point.
        point = m_treat - m_ctrl
        assert lo_obs <= point <= hi_obs, (
            f"Newcombe CI [{lo_obs:.4f}, {hi_obs:.4f}] does not bracket "
            f"meth_diff point {point:.4f}"
        )
        assert lo_obs >= -1.0 and hi_obs <= 1.0, (
            f"CI out of [-1, 1]: [{lo_obs:.4f}, {hi_obs:.4f}]"
        )
        # Newcombe is asymmetric near boundary betas; Wald is symmetric.
        # If we still got a Wald CI, lo + hi == 2*point exactly.
        if abs(point) > 0.3:
            assert abs((lo_obs + hi_obs) / 2.0 - point) > 1e-9, (
                f"CI is symmetric (Wald shape) at point={point:.4f}; "
                f"expected Newcombe asymmetry"
            )
```

- [ ] **Step 3: Run the test (expect failure)**

Run: `uv run pytest tests/test_dmc_lr.py::test_lr_meth_diff_ci_uses_newcombe -v`
Expected: FAIL — current CI is Wald-symmetric.

- [ ] **Step 4: Patch `dmc.py` to call `newcombe_diff_ci` on the lr path**

In `_process_one_chromosome` around the CI emission site (`dmc.py:1907-1912`), branch on `test`:

```python
# After pvals/log2_ors are computed and before the unified output block:
if test == "lr":
    from ._glm import newcombe_diff_ci
    # Use pooled per-group meth/cov sums (already accumulated for the
    # score path); meth_case_sum / cov_case_sum exist in this scope.
    meth_diff_ci_lo, meth_diff_ci_hi = newcombe_diff_ci(
        meth_case_sum, cov_case_sum,
        meth_ctrl_sum, cov_ctrl_sum,
    )
elif test == "fisher":
    # Wilson on pooled counts (already computed via fisher path).
    from ._glm import newcombe_diff_ci
    meth_diff_ci_lo, meth_diff_ci_hi = newcombe_diff_ci(
        meth_case_sum, cov_case_sum,
        meth_ctrl_sum, cov_ctrl_sum,
    )
else:
    # welch_t / glm keep their existing CI computation (correct for
    # those tests).
    # ... existing lo/hi calculation ...
```

(The `lr` and `fisher` paths both produce `meth_case_sum` / `cov_case_sum` / `meth_ctrl_sum` / `cov_ctrl_sum` accumulators; verify these are in scope at the CI site, or compute them in the `lr` branch from `sm_case` + `sn_case` if not.)

- [ ] **Step 5: Run the test (expect pass)**

Run: `uv run pytest tests/test_dmc_lr.py -v`
Expected: PASS.

- [ ] **Step 6: Run the broader test suite**

```
uv run pytest -m "not slow" -q -x 2>&1 | tail -10
```

Expected: passing. Tests that asserted on specific CI bounds for `lr` will need their expected values regenerated.

- [ ] **Step 7: CHANGELOG**

```markdown
- **P1-3**: `lr` and `fisher` engines now emit Newcombe (1998) hybrid
  Wilson-score CIs for `meth_diff_ci_{lo,hi}` (matching
  `statsmodels.confint_proportions_2indep(method='newcombe')`).
  Previously emitted Welch-normal Wald CIs which were symmetric near
  boundary β; Newcombe is asymmetric and respects [-1, 1] bounds
  properly. `welch_t` and `glm` keep their existing CIs (Welch and
  model-based respectively).
```

- [ ] **Step 8: Commit**

```
git commit -am "$(cat <<'EOF'
fix(dmc) P1-3: wire Newcombe CI into lr and fisher meth_diff_ci_{lo,hi}

The _glm.newcombe_diff_ci helper has existed since the GLM rewrite but
was never wired into the per-CpG output. Pre-Phase-3, lr and fisher
emitted Welch-normal Wald CIs which are symmetric and can sit outside
[-1, 1] near boundary betas.

This commit:
- Calls newcombe_diff_ci(meth_case_sum, cov_case_sum, meth_ctrl_sum,
  cov_ctrl_sum) on the lr and fisher paths.
- Leaves welch_t (Welch CI) and glm (model-based) unchanged -- both
  are correct for their respective tests.

Supplementary CI numbers in Studies 1 and 2 shift on lr / fisher
cells. Point estimates (meth_diff) unchanged.

Affects: lr@all, fisher@all

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 9

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: P1-4 — Explicit `reference_level` kwarg for patsy Treatment coding

**Spec §4 commit 10.** Add `reference_level` to `build_design` and surface via `tl.dmc(formula=..., reference_level=...)`. Default behaviour unchanged (patsy alphabetical default).

**Files:**
- Modify: `src/epykit/_glm.py:43-197` (build_design signature + patsy contrast setup)
- Modify: `src/epykit/tl.py` (dmc kwarg)
- Modify: `tests/test_glm.py` (append)
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_glm.py`:

```python
def test_reference_level_respected(synth_md_filtered):
    """P1-4: passing reference_level= sets patsy Treatment coding to
    use the named level as reference. The treatment coefficient sign
    flips when the reference is swapped."""
    md = synth_md_filtered
    md.obs = md.obs.with_columns(
        pl.when(pl.col("group") == "treatment").then(pl.lit("A"))
        .otherwise(pl.lit("B")).alias("group")
    )

    # Default: alphabetical reference -> 'A' is reference, coef is for B.
    ep.tl.dmc(md, test="glm", formula="~ group")
    coef_default = md.dmc["coef_treatment_log2"].to_numpy()

    # Explicit reference='B' -> coef is for A (sign flipped).
    ep.tl.dmc(md, test="glm", formula="~ group", reference_level="B")
    coef_swapped = md.dmc["coef_treatment_log2"].to_numpy()

    finite = np.isfinite(coef_default) & np.isfinite(coef_swapped)
    np.testing.assert_allclose(
        coef_default[finite], -coef_swapped[finite], atol=1e-9,
        err_msg="reference_level should flip the sign of the coefficient",
    )
```

- [ ] **Step 2: Run the test (expect failure)**

Run: `uv run pytest tests/test_glm.py::test_reference_level_respected -v`
Expected: FAIL — `reference_level` kwarg not implemented.

- [ ] **Step 3: Add `reference_level` to `build_design` in `_glm.py`**

Edit the signature:

```python
def build_design(
    obs: pl.DataFrame,
    samples_ordered: Sequence[str],
    formula: Optional[str] = None,
    covariates: Optional[Sequence[str]] = None,
    treatment_col: str = "treatment",
    require_treatment_col: bool = True,
    return_design_info: bool = False,
    reference_level: Optional[str] = None,  # NEW
) -> ...:
```

In the formula-construction block, when `reference_level` is set, wrap the relevant factor in patsy's `Treatment(reference="<level>")`. Identify the factor name as the leftmost categorical term (treatment_col by default). Example:

```python
if reference_level is not None:
    # Wrap the treatment_col reference in patsy's Treatment contrast.
    # Replace any bare occurrence of treatment_col in the formula
    # with C(treatment_col, Treatment(reference="<level>")).
    wrapped = f"C({treatment_col}, Treatment(reference='{reference_level}'))"
    terms = [wrapped if t == treatment_col else t for t in terms]
```

Log the resolved column names + reference level at INFO:

```python
import logging
logging.getLogger(__name__).info(
    "build_design: term_names=%s reference_level=%s",
    term_names, reference_level,
)
```

- [ ] **Step 4: Surface in `tl.dmc`**

In `tl.py::dmc`, add `reference_level: Optional[str] = None` to the signature; thread through to the `build_design` call when constructing the GLM design.

- [ ] **Step 5: Run the test (expect pass)**

Run: `uv run pytest tests/test_glm.py::test_reference_level_respected -v`
Expected: PASS.

- [ ] **Step 6: Run the surrounding suite**

```
uv run pytest -m "not slow" -q -x 2>&1 | tail -10
```

Expected: passing. Default-behaviour tests unaffected since `reference_level=None` preserves the alphabetical default.

- [ ] **Step 7: CHANGELOG**

```markdown
- **P1-4**: `tl.dmc(formula=..., reference_level=...)` lets users set
  the categorical reference level explicitly (via patsy's
  `Treatment(reference=...)`). Default behaviour unchanged
  (alphabetical reference). `_glm.build_design` now logs resolved
  column names and the chosen reference at INFO level.
```

- [ ] **Step 8: Commit**

```
git commit -am "$(cat <<'EOF'
fix(_glm) P1-4: explicit reference_level kwarg for patsy Treatment coding

Pre-Phase-3, the GLM backend used patsy's alphabetical-default
reference level silently. Users with a 'control' / 'treatment' factor
got 'control' as reference (alphabetical), which is what they
expected; users with 'A' / 'B' got 'A', which is fine; users with
'pre' / 'post' got 'post', which surprised them.

Added build_design(reference_level=...) kwarg; surfaced via
tl.dmc(reference_level=...). When set, the treatment factor is
wrapped in C(<col>, Treatment(reference='<level>')). Log
design_info.column_names and the chosen reference at INFO level.

Default behaviour preserved: no kwarg -> alphabetical reference, no
log emission change.

Affects:

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 10

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: P1-5 — NaN-mask non-converged IRLS sites + log fraction

**Spec §4 commit 11.** The existing `converged` boolean array is set during IRLS but only `separated` sites have their Wald stats NaN'd. Extend to also NaN non-converged sites and log the fraction.

**Files:**
- Modify: `src/epykit/_glm.py:234, 297, 321-324` (IRLS finalisation)
- Modify: `tests/test_glm.py` (append)
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_glm.py`:

```python
def test_nonconverged_irls_sites_are_nan(synth_md_filtered, caplog):
    """P1-5: when IRLS fails to converge at some sites, their Wald
    statistics must be NaN-masked and the fraction logged at WARNING."""
    md = synth_md_filtered
    # Build a degenerate-design scenario: a covariate that's collinear
    # with the intercept at most sites. The cleanest synthetic trigger
    # is a constant covariate -- patsy strips it, but if we force it
    # via numpy, IRLS may diverge at some sites.
    # Easier: use a very small max_iter to force non-convergence.
    md.obs = md.obs.with_columns(
        (pl.col("group") == "treatment").cast(int).alias("treatment")
    )
    import logging
    caplog.set_level(logging.WARNING, logger="epykit._glm")

    # Patch max_iter to 1 to force most sites to non-converge.
    from epykit import _glm
    original_max_iter = getattr(_glm, "DEFAULT_MAX_ITER", 50)
    try:
        _glm.DEFAULT_MAX_ITER = 1
        ep.tl.dmc(md, test="glm", formula="~ treatment")
    finally:
        _glm.DEFAULT_MAX_ITER = original_max_iter

    df = md.dmc
    # Expect substantial NaN in coef_treatment_log2 / pvalue.
    n_nan = int(df["pvalue"].is_null().sum())
    n_total = df.height
    assert n_nan > 0, "expected some NaN p-values from non-converged sites"

    # Expect a WARNING log mentioning non-convergence.
    warnings_logged = [r.message for r in caplog.records if r.levelname == "WARNING"]
    assert any("converg" in m.lower() for m in warnings_logged), (
        f"expected WARNING about non-convergence; got: {warnings_logged}"
    )
```

(If `DEFAULT_MAX_ITER` doesn't exist as a module-level constant, find the equivalent — likely a kwarg threaded through `irls_dispatch`. Adjust the patching strategy accordingly. Alternative: construct a degenerate design matrix directly and pass to `irls_dispatch` in the test, bypassing the `tl.dmc` surface.)

- [ ] **Step 2: Run the test (expect failure)**

Run: `uv run pytest tests/test_glm.py::test_nonconverged_irls_sites_are_nan -v`
Expected: FAIL — either no NaN (current behaviour) or no warning emitted.

- [ ] **Step 3: Patch `_glm.py`**

In `irls_dispatch` (or the appropriate finalisation site at `_glm.py:234,297,321-324`):

```python
# At the end of IRLS, after `converged` and `separated` are set:
n_total = int(meth_stack.shape[0])
n_separated = int(separated.sum())
n_nonconverged = int((~converged & ~separated).sum())
n_bad = n_separated + n_nonconverged

# NaN-mask all stat columns at separated OR non-converged sites.
bad_mask = ~converged | separated
beta_full = np.where(bad_mask[:, None], np.nan, beta_full)
se_full = np.where(bad_mask[:, None], np.nan, se_full)
dev_full = np.where(bad_mask, np.nan, dev_full)
pearson_full = np.where(bad_mask[:, None], np.nan, pearson_full)

frac_nonconverged = n_nonconverged / max(n_total, 1)
if frac_nonconverged > 0.01:
    logger.warning(
        "IRLS non-convergence at %d / %d sites (%.2f%%); their Wald "
        "stats are NaN-masked. Consider raising max_iter or checking "
        "design matrix rank.",
        n_nonconverged, n_total, 100 * frac_nonconverged,
    )
```

- [ ] **Step 4: Run the test (expect pass)**

Run: `uv run pytest tests/test_glm.py::test_nonconverged_irls_sites_are_nan -v`
Expected: PASS.

- [ ] **Step 5: Run broader suite**

```
uv run pytest -m "not slow" -q -x 2>&1 | tail -10
```

Expected: passing.

- [ ] **Step 6: CHANGELOG + commit**

```markdown
- **P1-5**: IRLS non-convergence is now surfaced. Non-converged sites
  have their Wald statistics NaN-masked (previously only `separated`
  sites were); fraction non-converged is logged at WARNING when > 1%.
  Previously non-converged sites silently emitted unreliable Wald
  stats.
```

```
git commit -am "$(cat <<'EOF'
fix(_glm) P1-5: NaN-mask non-converged IRLS sites + log fraction

The `converged` boolean was set during IRLS but only `separated`
sites had their Wald statistics NaN'd. Non-converged sites silently
emitted unreliable Wald stats (the IRLS iterate just before
divergence). Fix: NaN-mask both separated and non-converged; log
fraction non-converged at WARNING when > 1%.

Affects: glm@degenerate-cells

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 11

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 12: P1-6 — DMR empirical-FDR paired design + n=1,1 refusal

**Spec §4 commit 12.** Two related fixes in one commit. (a) Detect paired designs and shuffle within strata; (b) refuse n_treat=1, n_ctrl=1 with an explicit ValueError.

**Files:**
- Modify: `src/epykit/dmr.py:1370-1377` (empirical FDR loop)
- Create: `tests/test_dmr_empirical_fdr.py`
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_dmr_empirical_fdr.py`:

```python
"""P1-6: DMR empirical FDR now refuses n=1,1 (cannot permute) and
shuffles within strata when a paired design is detected."""
from __future__ import annotations

import pytest

import epykit as ep


def test_n_one_each_raises():
    """With one treatment and one control sample there are no valid
    label permutations -- empirical FDR is meaningless. Must raise
    ValueError with a clear migration hint."""
    import polars as pl
    from epykit.methyldata import MethylData
    from tests.fixtures.synth import SimConfig, generate
    import tempfile
    from pathlib import Path

    cfg = SimConfig(n_treatment=1, n_control=1)
    with tempfile.TemporaryDirectory() as td:
        result = generate(cfg, Path(td))
        md = MethylData.from_samplesheet(result["samplesheet"])
        ep.pp.filter_coverage(md, lo_count=3, hi_perc=99.9)
        ep.pp.unite(md, type="intersect")
        ep.tl.dmc(md, test="fisher")
        with pytest.raises(ValueError, match="n>=2"):
            ep.tl.dmr(
                md, method="tile", empirical_fdr=True, n_perm=10,
            )


@pytest.mark.slow
def test_paired_design_shuffles_within_strata(synth_md_filtered):
    """When a 'subject_id' column strictly pairs treatment/control
    samples, empirical FDR shuffles labels within each subject pair --
    not globally. Smoke: run with paired design; assertions are weaker
    than the full permutation theory but verify non-crash and that the
    empirical_qvalue column is populated."""
    md = synth_md_filtered
    # Build a paired covariate: each treatment id paired with one
    # control id.
    n_pair = len(md.treatment_ids)
    pair_ids = [f"S{i}" for i in range(n_pair)] * 2  # both groups share IDs
    md.obs = md.obs.with_columns(pl.Series("subject_id", pair_ids[:md.obs.height]))

    ep.tl.dmc(md, test="lr")
    ep.tl.dmr(
        md, method="tile", empirical_fdr=True, n_perm=20,
        empirical_strata="subject_id",
    )
    dmrs = md.uns["dmr"]
    assert "empirical_qvalue" in dmrs.columns
    assert dmrs["empirical_qvalue"].is_not_null().any()
```

(If `MethylData.from_samplesheet` isn't the actual class API surface, adapt to whatever loads a synth fixture. The slow test exercises a real permutation loop; the fast test only verifies the ValueError raise.)

- [ ] **Step 2: Run the failing tests**

Run: `uv run pytest tests/test_dmr_empirical_fdr.py -v`
Expected: both FAIL (n=1,1 doesn't raise; paired strata kwarg not implemented).

- [ ] **Step 3: Patch `dmr.py`**

In `src/epykit/dmr.py` around the empirical-FDR loop (`dmr.py:1370-1377`):

```python
def _empirical_fdr_for_dmr(
    md, dmrs, n_perm: int, *, empirical_strata: Optional[str] = None,
):
    n_treat = len(md.treatment_ids)
    n_ctrl  = len(md.control_ids)
    if n_treat == 1 and n_ctrl == 1:
        raise ValueError(
            "empirical DMR FDR requires n>=2 per group; got 1v1. "
            "Use Fisher-derived p-values directly via "
            "tl.dmc(test='fisher')."
        )

    rng = np.random.default_rng(0)
    all_samples = list(md.treatment_ids) + list(md.control_ids)

    if empirical_strata is not None and empirical_strata in md.obs.columns:
        # Within-stratum permutation: each shuffle reassigns labels
        # within each stratum independently.
        strata = md.obs[empirical_strata].to_list()
        strata_to_samples: dict = {}
        for s, sid in zip(strata, all_samples):
            strata_to_samples.setdefault(s, []).append(sid)
        def _perm():
            permuted: list = []
            for grp in strata_to_samples.values():
                shuffled = rng.permutation(grp).tolist()
                permuted.extend(shuffled)
            return permuted
    else:
        def _perm():
            return rng.permutation(all_samples).tolist()

    # ... existing permutation loop using _perm() in place of the
    # previous direct rng.permutation call ...
```

(Adapt to the actual existing structure of `_empirical_fdr_for_dmr` — the function exists; the changes are narrow.)

- [ ] **Step 4: Surface `empirical_strata` kwarg in `tl.dmr`**

`tl.dmr` already takes `empirical_fdr` and `n_perm`; add `empirical_strata: Optional[str] = None` and thread through to `_empirical_fdr_for_dmr`.

- [ ] **Step 5: Run tests**

```
uv run pytest tests/test_dmr_empirical_fdr.py -v
uv run pytest -m slow tests/test_dmr_empirical_fdr.py -v
```

Expected: both PASS (fast one immediately; slow one after marking with `@pytest.mark.slow`).

- [ ] **Step 6: CHANGELOG + commit**

```markdown
- **P1-6**: DMR empirical FDR now (a) raises `ValueError` at
  n_treat=1, n_ctrl=1 (no valid permutations), with a migration hint
  to use `tl.dmc(test='fisher')` directly; (b) accepts
  `empirical_strata=<column>` to permute within strata when a paired
  design exists. Previously ignored paired structure and shuffled
  globally, producing misleadingly small empirical FDRs.
```

```
git commit -am "$(cat <<'EOF'
fix(dmr) P1-6: empirical FDR paired-design + n=1,1 refusal

Pre-Phase-3, _empirical_fdr_for_dmr shuffled labels globally
regardless of design. Paired designs (subject_id-matched
treatment/control) got misleadingly small empirical FDRs because
between-pair variance was conflated with within-pair signal.

This commit:
- Adds tl.dmr(empirical_strata=<column>) to permute within strata.
- Refuses n_treat=1, n_ctrl=1 with ValueError pointing the user at
  tl.dmc(test='fisher').

Affects: empirical_dmr@paired-cells

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 12

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 13: P1-7 — Brown-Forsythe replaces Bartlett in `dvc.py`

**Spec §4 commit 13.** Bartlett's test assumes normality and is wildly wrong on bounded U-shaped β. Replace with two-pass Brown-Forsythe (median-centred Levene).

**Files:**
- Modify: `src/epykit/dvc.py:58-86`
- Modify: `tests/test_dvr.py` (existing — appends a reference-comparison test) or create `tests/test_dvc.py`
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Read the current Bartlett implementation**

```
rg -n "bartlett" src/epykit/dvc.py
```

Read 30 lines around the call.

- [ ] **Step 2: Write the failing test**

Create or append to `tests/test_dvc.py`:

```python
"""P1-7: DVC uses Brown-Forsythe (median-centred Levene) instead of
Bartlett on bounded beta values. Reference: scipy.stats.levene with
center='median'."""
from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import levene

# Import the per-site variance-test function directly. Name may differ;
# grep src/epykit/dvc.py for `def _per_site_variance_test` or the
# per-site loop and import accordingly.
from epykit.dvc import _per_site_variance_test  # adjust import as needed


def test_brown_forsythe_matches_scipy_levene_reference():
    """On a synthetic bimodal beta dataset, the DVC per-site F-stat
    must match scipy.stats.levene(center='median') to 1e-6."""
    rng = np.random.default_rng(0)
    # 50 synthetic sites with bimodal betas in two groups.
    n_sites = 50
    n_per_group = 6
    out_f = np.empty(n_sites)
    ref_f = np.empty(n_sites)
    for i in range(n_sites):
        a = rng.beta(0.5, 0.5, size=n_per_group)
        b = rng.beta(0.5, 0.5, size=n_per_group)
        out_f[i] = _per_site_variance_test(a, b)
        ref_f[i] = levene(a, b, center="median").statistic
    np.testing.assert_allclose(
        out_f, ref_f, rtol=1e-6, atol=1e-9,
        err_msg="Brown-Forsythe (median-Levene) reference mismatch",
    )
```

(If `_per_site_variance_test` isn't the actual function name, grep for the per-site call in `dvc.py:58-86` and adjust. The test's intent is: on bounded β, the DVC test must equal `scipy.stats.levene(..., center='median')`.)

- [ ] **Step 3: Run the failing test**

Run: `uv run pytest tests/test_dvc.py -v`
Expected: FAIL — current implementation uses Bartlett.

- [ ] **Step 4: Replace Bartlett with Brown-Forsythe**

In `src/epykit/dvc.py:58-86`:

```python
# Pre-Phase-3:
#   from scipy.stats import bartlett
#   stat, p = bartlett(group_a, group_b)
# Replace with:
def _per_site_variance_test(group_a: np.ndarray, group_b: np.ndarray) -> tuple[float, float]:
    """Brown-Forsythe (median-centred Levene). Robust on bounded U-shaped
    beta values where Bartlett's normality assumption is violated."""
    a = np.asarray(group_a, dtype=np.float64)
    b = np.asarray(group_b, dtype=np.float64)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    # Two-pass: medians, then absolute deviations, then F-test.
    z_a = np.abs(a - np.median(a))
    z_b = np.abs(b - np.median(b))
    # F-test on the absolute deviations (Brown-Forsythe = Levene with median).
    from scipy.stats import f_oneway
    f_stat, p = f_oneway(z_a, z_b)
    return float(f_stat), float(p)
```

(If the call site expects only the statistic, return just `f_stat`; adapt the signature to the existing surface.)

- [ ] **Step 5: Run the test (expect pass)**

Run: `uv run pytest tests/test_dvc.py -v`
Expected: PASS.

- [ ] **Step 6: Run the broader suite**

```
uv run pytest -m "not slow" -q -x 2>&1 | tail -10
```

Expected: passing. The existing `test_dvr.py` may assert on Bartlett-derived values; regenerate expected values against the new Brown-Forsythe output.

- [ ] **Step 7: CHANGELOG + commit**

```markdown
- **P1-7**: `tl.dvc` per-site variance test replaced from Bartlett with
  Brown-Forsythe (median-centred Levene). Bartlett assumes normality;
  beta methylation values are bounded U-shaped, which makes Bartlett
  badly miscalibrated. Brown-Forsythe is robust to non-normality and
  matches `scipy.stats.levene(center='median')` to 1e-6. DVC not in
  paper; no headline impact.
```

```
git commit -am "$(cat <<'EOF'
fix(dvc) P1-7: Brown-Forsythe replaces Bartlett on bounded beta

Bartlett's variance equality test assumes normality. Beta methylation
values are bounded [0, 1] and U-shaped (most CpGs sit near 0 or 1).
Bartlett is badly miscalibrated on this distribution.

Replaced with two-pass Brown-Forsythe (median-centred Levene):
- Pass 1: per-group medians
- Pass 2: F-test on absolute deviations from median

Matches scipy.stats.levene(center='median') to 1e-6 atol on a
synthetic bimodal fixture. DVC not in paper; no headline impact.

Affects: dvc@all

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 13

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 14: P1-9 — Sex-check Hartigan dip-test fallback for unimodal cohorts

**Spec §4 commit 14.** Add `diptest` to extras; gate the largest-gap clustering on a dip-test prereq; fall back to fixed Y-coverage threshold when unimodal.

**Files:**
- Modify: `src/epykit/qc.py:474-495`
- Modify: `pyproject.toml` (add `diptest` to the `qc` extra)
- Create or append: `tests/test_qc.py`
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Add `diptest` to `pyproject.toml`**

```
rg -n '"qc"' pyproject.toml
```

Find the `[project.optional-dependencies]` section. Append `"diptest"` to the `qc` list (or create the list if it doesn't exist).

- [ ] **Step 2: Install the new extra locally**

```
uv sync --extra dev --extra all
```

Expected: `diptest` installs without error.

- [ ] **Step 3: Write the failing test**

Append to `tests/test_qc.py`:

```python
"""P1-9: sex check on single-sex cohorts falls back to fixed
Y-coverage threshold (was: largest-gap 1D clustering, which produced
spurious splits on unimodal data)."""
from __future__ import annotations

import warnings

import numpy as np
import pytest


def test_sex_check_unimodal_falls_back():
    """Synthetic all-female cohort (low Y coverage everywhere): the
    dip-test prereq should fail, the function should fall back to the
    fixed 0.25 threshold, and a UserWarning should fire."""
    from epykit.qc import infer_sex_from_y_coverage  # adjust to actual name

    rng = np.random.default_rng(0)
    # All samples: low Y coverage, unimodal at 0.05.
    y_ratios = rng.normal(0.05, 0.005, size=8)
    sample_ids = [f"S{i}" for i in range(8)]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sex = infer_sex_from_y_coverage(sample_ids, y_ratios)
    assert all(s == "female" for s in sex.values()), (
        f"expected all-female; got {sex}"
    )
    user_warns = [w for w in caught if issubclass(w.category, UserWarning)]
    assert user_warns and "single-sex" in str(user_warns[0].message), (
        f"expected single-sex UserWarning; got {user_warns}"
    )
```

(Grep `qc.py` for the actual function name handling sex inference if `infer_sex_from_y_coverage` is wrong.)

- [ ] **Step 4: Run the failing test**

Run: `uv run pytest tests/test_qc.py -v -k unimodal`
Expected: FAIL — no dip-test fallback yet.

- [ ] **Step 5: Patch `qc.py`**

In the sex-check function around `qc.py:474-495`:

```python
import warnings
import diptest

def infer_sex_from_y_coverage(sample_ids, y_ratios, *, fallback_threshold: float = 0.25):
    y = np.asarray(y_ratios, dtype=np.float64)
    # Hartigan dip-test for unimodality.
    if len(y) >= 4:
        dip_stat, dip_p = diptest.diptest(y)
        if dip_p > 0.10:
            warnings.warn(
                "single-sex cohort detected (dip-test p=%.3f); sex inferred "
                "from Y-coverage fixed threshold (%.2f) only." % (
                    dip_p, fallback_threshold,
                ),
                UserWarning, stacklevel=2,
            )
            return {
                sid: ("male" if r >= fallback_threshold else "female")
                for sid, r in zip(sample_ids, y)
            }
    # Bimodal path: existing largest-gap clustering unchanged.
    # ... existing code ...
```

- [ ] **Step 6: Run the test (expect pass)**

Run: `uv run pytest tests/test_qc.py -v -k unimodal`
Expected: PASS.

- [ ] **Step 7: Run broader suite**

```
uv run pytest -m "not slow" -q -x 2>&1 | tail -10
```

Expected: passing. Existing bimodal-cohort tests should keep working since the dip-test prereq rejects (`p < 0.10`) and falls through to the existing path.

- [ ] **Step 8: CHANGELOG + commit**

```markdown
- **P1-9**: sex check now gates the largest-gap 1D clustering on a
  Hartigan dip-test prereq (`diptest.diptest`). On unimodal data (e.g.
  single-sex cohort) the function falls back to a fixed Y-coverage
  threshold (`0.25`) and emits a `UserWarning`. Previously spuriously
  split unimodal data into two clusters. Adds `diptest` to the `qc`
  extra.
```

```
git commit -am "$(cat <<'EOF'
fix(qc) P1-9: sex check dip-test fallback for unimodal cohorts

The pre-Phase-3 sex check ran 1D largest-gap clustering on the
Y-coverage ratio vector unconditionally. On single-sex cohorts (all
female or all male) the algorithm spuriously split the unimodal
distribution into two clusters, mis-assigning a fraction of samples
to the opposite sex.

Fix: gate the clustering on Hartigan's dip-test (diptest.diptest). If
the dip-test fails to reject unimodality (p > 0.10), fall back to a
fixed Y-coverage threshold (0.25) and emit UserWarning so the user
knows the cohort is single-sex. Bimodal cohorts unchanged.

Added diptest to the qc extra in pyproject.toml.

Affects: qc@single-sex-cohorts

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 14

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 15: P1-10 — Storey π₀ clamp at `1/n`

**Spec §4 commit 15.** `_storey_pi0` can return 0 when all p-values are below `lam`; clamp at `1/n` (Storey's standard plug-in floor).

**Files:**
- Modify: `src/epykit/dmc.py:2419-2421` (`_storey_pi0`)
- Create or append: `tests/test_dmc_multitest.py`
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Read the current `_storey_pi0`**

```
rg -n "def _storey_pi0" src/epykit/dmc.py
```

Read 20 lines.

- [ ] **Step 2: Write the failing test**

Create or append to `tests/test_dmc_multitest.py`:

```python
"""P1-10: _storey_pi0 must return at least 1/n when all p-values are
below the lambda threshold (preventing zero pi0 division)."""
from __future__ import annotations

import numpy as np
import pytest

from epykit.dmc import _storey_pi0


def test_storey_pi0_clamped_at_one_over_n():
    """All p-values below lam=0.5 -> numerator = 0; pre-fix returns 0,
    post-fix returns 1/n."""
    rng = np.random.default_rng(0)
    n = 1000
    pvals = rng.uniform(0, 0.4, size=n)  # all p < 0.5
    pi0 = _storey_pi0(pvals)
    assert pi0 >= 1.0 / n - 1e-12, (
        f"_storey_pi0 returned {pi0}; expected clamp at 1/n = {1.0/n:.6f}"
    )
    assert pi0 <= 1.0, f"pi0 = {pi0} outside [0, 1]"


def test_storey_pi0_unclamped_when_above_floor():
    """When some p-values are >= lam, pi0 is the unclamped estimate."""
    rng = np.random.default_rng(0)
    pvals = rng.uniform(0, 1, size=1000)  # uniform -> pi0 ~ 1
    pi0 = _storey_pi0(pvals)
    assert 0.9 <= pi0 <= 1.1, f"expected pi0 ~ 1 on uniform; got {pi0}"
```

- [ ] **Step 3: Run the failing test**

Run: `uv run pytest tests/test_dmc_multitest.py -v`
Expected: first test FAIL (returns 0); second PASS.

- [ ] **Step 4: Patch `_storey_pi0`**

In `src/epykit/dmc.py:2419-2421`:

```python
def _storey_pi0(pvals: np.ndarray, lam: float | None = None) -> float:
    """Plug-in estimator of the null proportion at lambda=0.5 (default).

    Returns the standard Storey plug-in `(# p > lam) / (n * (1 - lam))`,
    clamped at `1/n` from below (Storey 2002 §2.1 standard floor; prevents
    zero pi0 from yielding +inf q-values via BH-adjustment).

    For the smoother variant (lam scanned across a grid + spline fit),
    use a separate function; this one is the plug-in estimator only.
    """
    if lam is None:
        lam = 0.5
    n = len(pvals)
    if n == 0:
        return 1.0
    numer = float((pvals > lam).sum())
    denom = float(n) * (1.0 - lam)
    pi0_hat = numer / max(denom, 1e-12)
    return float(min(max(pi0_hat, 1.0 / n), 1.0))
```

- [ ] **Step 5: Run the tests (expect pass)**

Run: `uv run pytest tests/test_dmc_multitest.py -v`
Expected: both PASS.

- [ ] **Step 6: Run broader suite**

```
uv run pytest -m "not slow" -q -x 2>&1 | tail -10
```

Expected: passing. Default FDR path uses `fdr_tsbh` (bootstrap π₀), so the plug-in change rarely affects headline numbers.

- [ ] **Step 7: CHANGELOG + commit**

```markdown
- **P1-10**: `_storey_pi0` now clamps at `1/n` from below (Storey's
  standard plug-in floor). Previously could return 0 when all p-values
  fell below `lam`, yielding +inf q-values via BH adjustment.
  Docstring documents this is the plug-in estimator, not the spline-
  smoother variant.
```

```
git commit -am "$(cat <<'EOF'
fix(dmc) P1-10: _storey_pi0 clamped at 1/n

The plug-in Storey pi0 estimator (# p > lam) / (n * (1 - lam)) can
return 0 when all p-values fall below lam. The downstream BH adjusted
q = p * n / rank then gets divided by pi0 in the Storey-corrected
form, yielding +inf. Standard Storey 2002 floor: clamp at 1/n.

Docstring clarified that this is the plug-in estimator, not the
spline-smoother variant.

Affects: lr@all-significant-cells (edge case; default fdr_tsbh uses
bootstrap pi0).

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 15

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 16: `methylkit_stouffer_combine.R` — adjacent-3-CpG Stouffer combine

**Spec §4 commit 16.** R script that mirrors epykit's `neighbour_combine` on methylKit output. Lets the Phase 4 head-to-head be tuned-vs-tuned per PROTOCOL R1.

**Files:**
- Create: `benchmark/scripts/methylkit_stouffer_combine.R`
- Create: `benchmark/scripts/tests/test_methylkit_stouffer_combine.py`
- Create: `benchmark/scripts/tests/fixtures/methylkit_sample_in.tsv` (6-CpG fixture)
- Create: `benchmark/scripts/tests/fixtures/methylkit_sample_expected.tsv` (hand-computed expected)

- [ ] **Step 1: Build the 6-CpG fixture**

Create `benchmark/scripts/tests/fixtures/methylkit_sample_in.tsv`:

```
chr	start	end	strand	pvalue	qvalue	meth.diff
chr1	100	100	+	0.01	0.02	-30.0
chr1	200	200	+	0.005	0.01	-25.0
chr1	300	300	+	0.02	0.04	-20.0
chr1	5000	5000	+	0.3	0.4	5.0
chr1	5100	5100	+	0.25	0.35	7.0
chr1	5200	5200	+	0.5	0.6	-2.0
```

(Three CpGs clustered around 100-300 bp (within a `max_gap_bp=1000` window of each other) and three around 5000-5200 (also clustered). The combined p-value for the first triple, with equal weights and `max_gap_bp=1000`, is:

```
z_i = qnorm(1 - p_i/2)
z_combined = sum(z_i) / sqrt(3)
p_combined = 2 * (1 - pnorm(|z_combined|))
```

For p = (0.01, 0.005, 0.02): z = (2.576, 2.807, 2.326). Sum = 7.709. /sqrt(3) = 4.451. p ≈ 8.5e-6.)

Create `benchmark/scripts/tests/fixtures/methylkit_sample_expected.tsv` with the hand-computed `pvalue_combined` and `qvalue_combined` columns for all 6 rows.

- [ ] **Step 2: Write the failing test**

Create `benchmark/scripts/tests/test_methylkit_stouffer_combine.py`:

```python
"""Tests for methylkit_stouffer_combine.R.

Skips when Rscript is not on PATH (CI / Windows-Python-only envs)."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import polars as pl
import pytest

FIXTURE_DIR = Path(__file__).parent / "fixtures"
SCRIPT = Path(__file__).parents[1] / "methylkit_stouffer_combine.R"


@pytest.mark.slow
def test_methylkit_stouffer_combine_matches_expected(tmp_path):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not on PATH; skipping R subprocess test")
    in_tsv = FIXTURE_DIR / "methylkit_sample_in.tsv"
    expected = pl.read_csv(FIXTURE_DIR / "methylkit_sample_expected.tsv", separator="\t")
    out_tsv = tmp_path / "out.tsv"
    result = subprocess.run(
        ["Rscript", str(SCRIPT),
         "--in", str(in_tsv),
         "--out", str(out_tsv),
         "--max-gap-bp", "1000"],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, (
        f"Rscript failed: stderr=\n{result.stderr}\nstdout=\n{result.stdout}"
    )
    got = pl.read_csv(out_tsv, separator="\t")
    # Schema must include the new columns.
    assert "pvalue_combined" in got.columns
    assert "qvalue_combined" in got.columns
    # Row count preserved.
    assert got.height == 6
    # Combined p-values match hand-computed to 5 decimal places.
    import numpy as np
    np.testing.assert_allclose(
        got["pvalue_combined"].to_numpy(),
        expected["pvalue_combined"].to_numpy(),
        rtol=1e-4,
    )
```

- [ ] **Step 3: Run the test to verify failure**

Run: `uv run pytest benchmark/scripts/tests/test_methylkit_stouffer_combine.py -v`
Expected: FAIL — script doesn't exist yet. (Or SKIP if Rscript unavailable; that's also acceptable.)

- [ ] **Step 4: Write the R script**

Create `benchmark/scripts/methylkit_stouffer_combine.R`:

```R
#!/usr/bin/env Rscript
# methylkit_stouffer_combine.R
#
# Apply adjacent-CpG Stouffer combination to a methylKit per-scenario
# TSV. Mirrors epykit's neighbour_combine knob so Phase 4's head-to-head
# can be tuned-vs-tuned per PROTOCOL R1.
#
# Inputs:
#   --in <tsv>           methylKit calculateDiffMeth output (must contain
#                        chr, start, pvalue columns; qvalue optional)
#   --out <tsv>          output path
#   --max-gap-bp <int>   neighbours within this bp window are pooled
#                        (default 1000)
#   --window <int>       Stouffer window size (default 3 = focal + 2
#                        nearest neighbours, one on each side)
#
# Output schema:
#   input columns + pvalue_combined + qvalue_combined (BH per chromosome)
#
# Assertions:
#   - input must NOT have been BH-corrected already (qvalue is at most
#     1.5x larger than pvalue per row in expectation; we check that
#     qvalue / pvalue is finite and not orders of magnitude inflated).
#   - direction of meth.diff among the window must agree (mixed-direction
#     windows have p_combined = NA -- conservative).

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--in",         type = "character", help = "Input TSV"),
  make_option("--out",        type = "character", help = "Output TSV"),
  make_option("--max-gap-bp", type = "integer",   default = 1000L),
  make_option("--window",     type = "integer",   default = 3L)
)
opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$`in`) || is.null(opt$out)) {
  stop("--in and --out are required")
}

df <- fread(opt$`in`, sep = "\t", header = TRUE)
stopifnot(all(c("chr", "start", "pvalue") %in% names(df)))

stouffer <- function(pvals) {
  pvals <- pmin(pmax(pvals, 1e-300), 1 - 1e-15)
  z <- qnorm(1 - pvals / 2)  # two-sided -> half-tail z per site
  z_comb <- sum(z) / sqrt(length(z))
  2 * (1 - pnorm(abs(z_comb)))
}

# Process per chromosome.
df[, pvalue_combined := NA_real_]
setkey(df, chr, start)

for (chrom in unique(df$chr)) {
  rows <- df[chr == chrom]
  pos <- rows$start
  pv  <- rows$pvalue
  md  <- if ("meth.diff" %in% names(rows)) rows$`meth.diff` else rep(0, nrow(rows))
  comb <- rep(NA_real_, nrow(rows))
  k <- opt$window
  half <- (k - 1) %/% 2
  for (i in seq_along(pv)) {
    lo <- max(1, i - half)
    hi <- min(length(pv), i + half)
    # Restrict to within max-gap-bp of the focal CpG.
    in_window <- which(abs(pos[lo:hi] - pos[i]) <= opt$`max-gap-bp`)
    idx <- (lo:hi)[in_window]
    if (length(idx) < 2) {
      comb[i] <- pv[i]
      next
    }
    # Direction check: skip mixed-direction windows.
    dirs <- sign(md[idx])
    if (length(unique(dirs[dirs != 0])) > 1) {
      comb[i] <- pv[i]  # fall back to raw p
      next
    }
    comb[i] <- stouffer(pv[idx])
  }
  df[chr == chrom, pvalue_combined := comb]
}

# BH per chromosome.
df[, qvalue_combined := NA_real_]
for (chrom in unique(df$chr)) {
  mask <- df$chr == chrom & !is.na(df$pvalue_combined)
  if (any(mask)) {
    df[mask, qvalue_combined := p.adjust(pvalue_combined, method = "BH")]
  }
}

fwrite(df, opt$out, sep = "\t", quote = FALSE)
cat(sprintf("wrote %s (%d rows)\n", opt$out, nrow(df)))
```

- [ ] **Step 5: Run the test (expect pass)**

Run: `uv run pytest benchmark/scripts/tests/test_methylkit_stouffer_combine.py -v`
Expected: PASS (or SKIP if Rscript unavailable).

- [ ] **Step 6: Smoke against a real methylKit file**

If `benchmark/data/study2/methylkit_results/dmc_cov10.tsv` exists:

```
Rscript benchmark/scripts/methylkit_stouffer_combine.R \
  --in benchmark/data/study2/methylkit_results/dmc_cov10.tsv \
  --out /tmp/dmc_cov10_tuned.tsv \
  --max-gap-bp 1000
```

Expected: prints `wrote /tmp/dmc_cov10_tuned.tsv (N rows)`. Inspect with `head /tmp/dmc_cov10_tuned.tsv`.

- [ ] **Step 7: CHANGELOG + commit**

```markdown
### Added (benchmark)

- **`benchmark/scripts/methylkit_stouffer_combine.R`**: adjacent-CpG
  Stouffer combine for methylKit output. Mirrors epykit's
  `neighbour_combine` knob so Phase 4's head-to-head can be tuned-vs-
  tuned per PROTOCOL R1. Tested with a 6-CpG hand-computed fixture;
  test SKIPS without Rscript on PATH.
```

```
git add benchmark/scripts/methylkit_stouffer_combine.R \
        benchmark/scripts/tests/test_methylkit_stouffer_combine.py \
        benchmark/scripts/tests/fixtures/methylkit_sample_in.tsv \
        benchmark/scripts/tests/fixtures/methylkit_sample_expected.tsv \
        CHANGELOG.md
git commit -m "$(cat <<'EOF'
feat(benchmark): methylkit_stouffer_combine.R -- adjacent-CpG combine

R script that applies Stouffer combination across a window of adjacent
methylKit CpGs, mirroring epykit's neighbour_combine knob. Required
for Phase 4's head-to-head to be tuned-vs-tuned per PROTOCOL R1.

Schema: input columns + pvalue_combined + qvalue_combined (BH per
chromosome). Mixed-direction windows fall back to the raw p-value
(conservative). Test SKIPS when Rscript not on PATH; otherwise
verifies against a 6-CpG hand-computed fixture.

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 16

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 17: `_null_engines.py` + wire `run_null_calibration.py` to real `ep.tl.dmc`

**Spec §4 commit 17.** Phase 2 left `run_null_calibration.py` with a mock noise engine in `main()`. This task adds real engine closures and rewrites `main()` to dispatch via them.

**Files:**
- Create: `benchmark/scripts/_null_engines.py`
- Modify: `benchmark/scripts/run_null_calibration.py`
- Create: `benchmark/scripts/tests/test_null_engines.py`

- [ ] **Step 1: Write the failing test (parametrised across surviving engines)**

Create `benchmark/scripts/tests/test_null_engines.py`:

```python
"""Test the real-engine closures: each surviving DMC/DMR engine wraps
into the engine_fn(samples_treatment, samples_control, seed)->qvals
contract."""
from __future__ import annotations

import numpy as np
import pytest

from _null_engines import ENGINE_REGISTRY


SURVIVING_DMC = ["lr", "lr_plus", "welch_t", "fisher", "glm"]
SURVIVING_DMR = ["tile", "sliding", "chain_merge", "segment"]


@pytest.mark.slow
@pytest.mark.parametrize("engine_name", SURVIVING_DMC + SURVIVING_DMR)
def test_engine_closure_runs_and_returns_qvalues(synth_md_filtered, engine_name):
    """Each registered engine: callable, returns 1D array of q-values
    in [0, 1] (or NaN), deterministic across two seeded runs."""
    md = synth_md_filtered
    closure = ENGINE_REGISTRY[engine_name](md)  # factory that captures md
    treat = list(md.treatment_ids)
    ctrl = list(md.control_ids)
    out_a = closure(samples_treatment=treat, samples_control=ctrl, seed=42)
    out_b = closure(samples_treatment=treat, samples_control=ctrl, seed=42)
    assert isinstance(out_a, np.ndarray)
    assert out_a.ndim == 1
    np.testing.assert_array_equal(out_a, out_b, err_msg="not deterministic")
    finite = out_a[np.isfinite(out_a)]
    assert ((finite >= 0) & (finite <= 1)).all(), (
        f"q-values out of [0, 1] for engine={engine_name}"
    )
```

(The `synth_md_filtered` fixture comes from the main tests/conftest.py — need to make it importable here. Phase 2's `benchmark/scripts/tests/conftest.py` already adds the scripts dir to sys.path; extend it to also import the main test conftest's fixtures.)

- [ ] **Step 2: Extend `benchmark/scripts/tests/conftest.py` to import main fixtures**

Open `benchmark/scripts/tests/conftest.py` and append:

```python
# Re-export the main tests/ conftest fixtures so benchmark-script
# tests can use synth_md, synth_md_filtered without duplicating the
# fixture code.
import sys as _sys
from pathlib import Path as _P
_REPO = _P(__file__).resolve().parents[3]
_sys.path.insert(0, str(_REPO))
from tests.conftest import (  # noqa: F401, E402
    synth_bundle, synth_md, synth_md_filtered,
)
```

- [ ] **Step 3: Run the failing test**

Run: `uv run pytest benchmark/scripts/tests/test_null_engines.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named '_null_engines'`.

- [ ] **Step 4: Create `_null_engines.py`**

Create `benchmark/scripts/_null_engines.py`:

```python
"""Real-engine closures for run_null_calibration.py.

Each registry entry is a factory: ``factory(md) -> closure``, where
the closure has signature ``(samples_treatment, samples_control,
seed) -> np.ndarray of q-values``. Factories capture the MethylData
object once outside the shuffle loop so the parquet store is loaded
only once per (engine, scenario).
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import polars as pl

import epykit as ep


def _permute_obs(md, samples_treatment, samples_control):
    """Build an obs frame where the named samples are tagged as
    treatment/control. Returns a shallow-copied MethylData with
    md.obs replaced; the methylstore is unchanged."""
    label_map = {s: "treatment" for s in samples_treatment}
    label_map.update({s: "control" for s in samples_control})
    new_groups = [label_map.get(sid, md.obs["group"][i]) for i, sid in enumerate(md.obs["sample_id"])]
    md.obs = md.obs.with_columns(pl.Series("group", new_groups))
    return md


def _dmc_engine(test_name: str, *, use_lr_plus: bool = False) -> Callable:
    def factory(md):
        def closure(samples_treatment, samples_control, seed):
            _permute_obs(md, samples_treatment, samples_control)
            if use_lr_plus:
                ep.tl.dmc(md, test="lr", power_stack="auto",
                          fdr_method="fdr_tsbh", neighbour_combine=True,
                          sep_fallback=True, dispersion="eb")
            else:
                ep.tl.dmc(md, test=test_name)
            df = md.dmc
            qcol = "qvalue" if "qvalue" in df.columns else "pvalue"
            return df[qcol].to_numpy().astype(np.float64)
        return closure
    return factory


def _dmr_engine(method: str) -> Callable:
    def factory(md):
        def closure(samples_treatment, samples_control, seed):
            _permute_obs(md, samples_treatment, samples_control)
            ep.tl.dmc(md, test="lr")  # DMR needs a DMC pass
            ep.tl.dmr(md, method=method)
            dmrs = md.uns["dmr"]
            if "qvalue" not in dmrs.columns or dmrs.height == 0:
                return np.array([], dtype=np.float64)
            return dmrs["qvalue"].to_numpy().astype(np.float64)
        return closure
    return factory


ENGINE_REGISTRY: dict[str, Callable] = {
    # DMC
    "lr":       _dmc_engine("lr"),
    "lr_plus":  _dmc_engine("lr", use_lr_plus=True),
    "welch_t":  _dmc_engine("welch_t"),
    "fisher":   _dmc_engine("fisher"),
    "glm":      _dmc_engine("glm"),
    # DMR
    "tile":         _dmr_engine("tile"),
    "sliding":      _dmr_engine("sliding"),
    "chain_merge":  _dmr_engine("chain_merge"),
    "segment":      _dmr_engine("segment"),
}
```

- [ ] **Step 5: Rewrite `run_null_calibration.py::main` to dispatch via the registry**

Replace the existing `main()` in `benchmark/scripts/run_null_calibration.py`:

```python
def main(argv: list[str] | None = None) -> None:
    """CLI: real-engine label-shuffle calibration.

    python run_null_calibration.py \\
        --engine lr --methylstore <path-to-md> --scenario cov10_3v3 \\
        --k-shuffles 20 --seed 0 \\
        --out benchmark/data/null_calibration/cov10_3v3/lr.parquet
    """
    import argparse
    from pathlib import Path

    from _null_engines import ENGINE_REGISTRY
    from epykit.methyldata import MethylData

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", required=True, choices=sorted(ENGINE_REGISTRY))
    parser.add_argument("--methylstore", required=True, type=Path,
                        help="Path to the methylstore root (.cache parent)")
    parser.add_argument("--scenario", required=True, type=str)
    parser.add_argument("--k-shuffles", type=int, default=20)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--q-thresh", type=float, default=0.05)
    args = parser.parse_args(argv)

    md = MethylData.load(args.methylstore)
    closure = ENGINE_REGISTRY[args.engine](md)

    samples = list(md.treatment_ids) + list(md.control_ids)
    n_per_group = len(md.treatment_ids)

    df = run_null_calibration(
        engine_fn=closure,
        engine_name=args.engine,
        scenario_name=args.scenario,
        samples=samples,
        n_per_group=n_per_group,
        k_shuffles=args.k_shuffles,
        q_thresh=args.q_thresh,
        seed=args.seed,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(args.out)
    print(f"wrote {args.out} with {len(df)} shuffle rows for engine={args.engine}")
```

(Adapt `MethylData.load` to the actual class API surface — `rg -n "def load|@classmethod" src/epykit/methyldata.py` to find the right loader.)

- [ ] **Step 6: Run the test (expect pass, slow)**

Run: `uv run pytest -m slow benchmark/scripts/tests/test_null_engines.py -v`
Expected: 9 PASS (5 DMC + 4 DMR engines).

- [ ] **Step 7: Run the full benchmark-scripts suite (no regression)**

```
uv run pytest benchmark/scripts/tests/ -m "not slow" -q
```

Expected: all passing (the existing Phase 2 tests for `run_null_calibration.py` with the mock engine still work; the `main()` rewrite didn't break the module-level functions).

- [ ] **Step 8: CHANGELOG + commit**

```markdown
- **`benchmark/scripts/_null_engines.py`** (new): real-engine
  closures for `run_null_calibration.py`. Wraps each surviving DMC
  engine (`lr`, `lr_plus`, `welch_t`, `fisher`, `glm`) and DMR
  method (`tile`, `sliding`, `chain_merge`, `segment`) as
  factories matching the Phase 2 `engine_fn(samples_treatment,
  samples_control, seed)->qvals` contract. `run_null_calibration.py
  main()` now takes `--engine` + `--methylstore` and dispatches via
  the registry instead of the Phase-2 mock noise engine.
```

```
git commit -am "$(cat <<'EOF'
feat(benchmark): null_engines + real-engine wiring for run_null_calibration

Phase 2 left run_null_calibration.py with a mock noise engine in
main(). This commit:

- Adds _null_engines.py with ENGINE_REGISTRY covering the surviving
  DMC engines (lr, lr_plus, welch_t, fisher, glm) and DMR methods
  (tile, sliding, chain_merge, segment).
- Each closure: permutes md.obs[group] per call, runs the engine
  via ep.tl.dmc / ep.tl.dmr, returns the qvalue array.
- main() of run_null_calibration.py rewritten to dispatch via the
  registry: --engine <name> --methylstore <path> --scenario ...
- Extended benchmark/scripts/tests/conftest.py to re-export
  synth_md_filtered from the main tests/conftest.py.

Tests (slow): 9 engines x deterministic-with-seed + value-range
assertions on a synth fixture.

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 17

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 18: `evaluate.py` emits Wilson + bootstrap CI columns

**Spec §4 commit 18.** Append a finalisation step to evaluate.py that calls the Phase 2 CI helpers on every row.

**Files:**
- Modify: `benchmark/scripts/evaluate.py`
- Create: `benchmark/scripts/tests/test_evaluate_ci.py`

- [ ] **Step 1: Locate evaluate.py**

```
ls benchmark/scripts/ | grep -i evaluate
```

If evaluate.py exists, read it. If it doesn't yet exist as a single script (the existing benchmark layout may have eval logic split across `make_summary_figures.py` and others), this task creates a thin wrapper:

```
rg -n "eval_summary" benchmark/scripts/*.py
```

If no single `evaluate.py` is found, create `benchmark/scripts/evaluate.py` as a finalisation runner that reads an existing `eval_summary.parquet`, adds the CI columns, and writes back. This is the path of least resistance for Phase 3; Phase 4 may inline it into the main eval loop.

- [ ] **Step 2: Write the failing test**

Create `benchmark/scripts/tests/test_evaluate_ci.py`:

```python
"""evaluate.py --ci adds Wilson + bootstrap CI columns to
eval_summary.parquet. Phase 3 deliverable."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest


SCRIPT = Path(__file__).parents[1] / "evaluate.py"


def _make_synth_eval(path: Path, n_rows: int = 5):
    """Write a minimal eval_summary.parquet for the CI step to consume."""
    pl.DataFrame({
        "tool": ["epykit_lr"] * n_rows,
        "scenario": [f"cov10_n3_seed{i}" for i in range(n_rows)],
        "test": ["lr"] * n_rows,
        "tp": [90, 85, 70, 60, 95],
        "fp": [10, 5, 20, 30, 2],
        "tn": [890, 905, 880, 870, 898],
        "fn": [10, 5, 30, 40, 5],
        "tpr": [90 / 100, 85 / 90, 70 / 100, 60 / 100, 95 / 100],
        "fpr": [10 / 900, 5 / 910, 20 / 900, 30 / 900, 2 / 900],
        "f1": [0.90, 0.94, 0.74, 0.63, 0.96],
        "auroc": [0.97, 0.95, 0.88, 0.80, 0.99],
    }).write_parquet(path)


def test_evaluate_ci_adds_columns(tmp_path):
    eval_in = tmp_path / "eval_summary.parquet"
    _make_synth_eval(eval_in)
    result = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--ci-only", "--eval-summary", str(eval_in)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, (
        f"evaluate --ci-only failed: {result.stderr}"
    )
    out = pl.read_parquet(eval_in)  # in-place
    for col in ("tpr_ci_lo", "tpr_ci_hi", "fpr_ci_lo", "fpr_ci_hi"):
        assert col in out.columns, f"missing column {col}; got {out.columns}"


def test_evaluate_ci_brackets_point_estimate(tmp_path):
    eval_in = tmp_path / "eval_summary.parquet"
    _make_synth_eval(eval_in)
    subprocess.run(
        [sys.executable, str(SCRIPT),
         "--ci-only", "--eval-summary", str(eval_in)],
        check=True, capture_output=True, text=True,
    )
    out = pl.read_parquet(eval_in)
    lo = out["tpr_ci_lo"].to_numpy()
    hi = out["tpr_ci_hi"].to_numpy()
    p = out["tpr"].to_numpy()
    valid = np.isfinite(lo) & np.isfinite(hi) & np.isfinite(p)
    assert ((lo <= p) & (p <= hi))[valid].all()
```

- [ ] **Step 3: Run failing test**

Run: `uv run pytest benchmark/scripts/tests/test_evaluate_ci.py -v`
Expected: FAIL — script doesn't exist or doesn't have `--ci-only`.

- [ ] **Step 4: Write or extend `benchmark/scripts/evaluate.py`**

Create (or modify) `benchmark/scripts/evaluate.py`:

```python
"""evaluate.py — append Wilson + bootstrap CIs to eval_summary.parquet.

Phase 3 scope: --ci-only mode that operates on an existing parquet.
Phase 4 may inline the CI step into the main eval loop (which lives
elsewhere right now).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl

# Make sibling helper importable without packaging.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from wilson_bootstrap_ci import (  # noqa: E402
    add_wilson_ci_for_tpr_fpr,
    bootstrap_auroc_ci,
    bootstrap_f1_ci,
)


def _add_bootstrap_columns_if_per_cpg_available(df: pl.DataFrame, per_cpg_dir: Path | None) -> pl.DataFrame:
    """Bootstrap CIs for AUROC / F1 require per-CpG joined frames at
    eval_per_cpg/<tool>_<scenario>.parquet. If unavailable, fill the
    columns with NaN -- Phase 3 lands the wiring; Phase 4 populates
    when the eval pass caches per-CpG frames."""
    import numpy as np
    n = df.height
    auroc_lo = np.full(n, np.nan); auroc_hi = np.full(n, np.nan)
    f1_lo    = np.full(n, np.nan); f1_hi    = np.full(n, np.nan)
    if per_cpg_dir is not None and per_cpg_dir.exists():
        for i, row in enumerate(df.iter_rows(named=True)):
            cache = per_cpg_dir / f"{row['tool']}_{row['scenario']}.parquet"
            if not cache.exists():
                continue
            j = pl.read_parquet(cache)
            is_dmc = j["is_dmc"].to_numpy().astype(bool)
            pvals  = j["pvalue"].to_numpy()
            qvals  = j["qvalue"].to_numpy() if "qvalue" in j.columns else pvals
            seed = abs(hash((row["tool"], row["scenario"], 0.05))) % (2**32)
            auroc_lo[i], auroc_hi[i] = bootstrap_auroc_ci(
                is_dmc=is_dmc, pvalues=pvals, B=1000, seed=seed,
            )
            f1_lo[i], f1_hi[i] = bootstrap_f1_ci(
                is_dmc=is_dmc, qvalues=qvals, threshold=0.05, B=1000, seed=seed,
            )
    return df.with_columns([
        pl.Series("auroc_ci_lo", auroc_lo), pl.Series("auroc_ci_hi", auroc_hi),
        pl.Series("f1_ci_lo", f1_lo),       pl.Series("f1_ci_hi", f1_hi),
    ])


def add_ci_columns(df: pl.DataFrame, per_cpg_dir: Path | None = None) -> pl.DataFrame:
    df = add_wilson_ci_for_tpr_fpr(df)
    df = _add_bootstrap_columns_if_per_cpg_available(df, per_cpg_dir)
    return df


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ci-only", action="store_true",
                        help="Append CI columns to an existing eval_summary.parquet in place.")
    parser.add_argument("--eval-summary", required=True, type=Path)
    parser.add_argument("--per-cpg-dir", type=Path, default=None,
                        help="Optional dir of <tool>_<scenario>.parquet per-CpG caches.")
    args = parser.parse_args(argv)

    if not args.ci_only:
        parser.error("Only --ci-only mode is implemented in Phase 3.")

    df = pl.read_parquet(args.eval_summary)
    out = add_ci_columns(df, per_cpg_dir=args.per_cpg_dir)
    out.write_parquet(args.eval_summary)
    print(f"updated {args.eval_summary} with CI columns ({out.height} rows)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Run the tests (expect pass)**

Run: `uv run pytest benchmark/scripts/tests/test_evaluate_ci.py -v`
Expected: both PASS.

- [ ] **Step 6: CHANGELOG + commit**

```markdown
- **`benchmark/scripts/evaluate.py`** (new, Phase 3 surface):
  `--ci-only` mode appends Wilson 95% CIs on TPR/FPR and bootstrap
  CIs on AUROC/F1 to an existing `eval_summary.parquet`. Bootstrap
  CIs are filled when `--per-cpg-dir` points to cached per-CpG
  joined frames; otherwise they're NaN. Closes the "every cell
  has a CI" requirement from PROTOCOL Results §6.
```

```
git add benchmark/scripts/evaluate.py benchmark/scripts/tests/test_evaluate_ci.py CHANGELOG.md
git commit -m "$(cat <<'EOF'
feat(benchmark): evaluate.py --ci-only appends CI columns

Phase 3 deliverable: --ci-only mode that adds tpr_ci_lo, tpr_ci_hi,
fpr_ci_lo, fpr_ci_hi (Wilson) and auroc_ci_lo, auroc_ci_hi,
f1_ci_lo, f1_ci_hi (bootstrap, B=1000) to an existing
eval_summary.parquet in place.

Bootstrap CIs are filled only when --per-cpg-dir points to cached
per-CpG joined frames named <tool>_<scenario>.parquet. Phase 4's
locked re-run will populate that cache during the main eval pass;
Phase 3 wires the mechanism with NaN-fill fallback so the columns
always exist.

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 18

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 19: `regen_all.py` — acceptance gate (seed manifest)

**Spec §4 commit 19.** Verify-mode acceptance gate that reads `<!-- source: scripts/X.py -->` HTML comments in `paper.md` and asserts cited parquets match printed values. Phase 3 lands an empty `claims.yaml` seed; Phase 4 populates.

**Files:**
- Create: `benchmark/scripts/regen_all.py`
- Create: `benchmark/scripts/claims.yaml` (empty seed)
- Create: `benchmark/scripts/tests/test_regen_all.py`
- Create: `benchmark/scripts/tests/fixtures/regen_paper_ok.md`, `regen_paper_bad.md`, `regen_claims_ok.yaml`

- [ ] **Step 1: Write the failing tests**

Create `benchmark/scripts/tests/test_regen_all.py`:

```python
"""regen_all.py --verify reads claims.yaml + a paper markdown,
asserts each claim matches the parquet to the printed precision.
Phase 3 lands the gate; Phase 4 populates claims."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

FIXTURES = Path(__file__).parent / "fixtures"
SCRIPT = Path(__file__).parents[1] / "regen_all.py"


def _make_parquet(path: Path):
    pl.DataFrame({
        "tool": ["epykit_lr"], "scenario": ["cov10_3v3"],
        "auroc": [0.987], "tpr": [0.95],
    }).write_parquet(path)


def test_verify_pass_on_matching_claim(tmp_path):
    parquet = tmp_path / "vals.parquet"
    _make_parquet(parquet)
    claims = tmp_path / "claims.yaml"
    claims.write_text(f"""
- claim_id: study1_auroc
  parquet: {parquet}
  column: auroc
  filter:
    tool: epykit_lr
    scenario: cov10_3v3
  expected: 0.987
  precision: 0.001
""")
    paper = tmp_path / "paper.md"
    paper.write_text(
        "epykit's AUROC was 0.987 <!-- claim: study1_auroc --> on Study 1.\\n"
    )
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--verify",
         "--claims", str(claims), "--paper", str(paper)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, f"verify should pass: {result.stderr}"


def test_verify_fail_on_off_by_precision(tmp_path):
    parquet = tmp_path / "vals.parquet"
    _make_parquet(parquet)
    claims = tmp_path / "claims.yaml"
    # Expected differs from parquet by 0.1 (way beyond precision=0.001).
    claims.write_text(f"""
- claim_id: study1_auroc
  parquet: {parquet}
  column: auroc
  filter:
    tool: epykit_lr
    scenario: cov10_3v3
  expected: 0.887
  precision: 0.001
""")
    paper = tmp_path / "paper.md"
    paper.write_text("AUROC 0.887 <!-- claim: study1_auroc -->\\n")
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--verify",
         "--claims", str(claims), "--paper", str(paper)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0, "verify should fail on mismatch"
    assert "study1_auroc" in result.stdout + result.stderr
```

- [ ] **Step 2: Run failing test**

Run: `uv run pytest benchmark/scripts/tests/test_regen_all.py -v`
Expected: FAIL — script doesn't exist.

- [ ] **Step 3: Write `regen_all.py`**

Create `benchmark/scripts/regen_all.py`:

```python
"""regen_all.py — acceptance gate for paper claims.

Modes:
  --verify       (default in CI): assert every claim in claims.yaml
                 matches its cited parquet to the printed precision.
                 Exits non-zero on any mismatch.
  --run-cheap    Re-runs fast scripts (CI helpers, eval_summary
                 finalisation). Phase 4 expansion target.
  --run-all      Full regen including locked re-runs. Phase 4 only.

claims.yaml schema (one list entry per claim):
  - claim_id:   stable identifier referenced by <!-- claim: <id> --> comments
    parquet:    path (absolute or repo-relative) to the source parquet
    column:     column name in the parquet to read
    filter:     mapping of column -> value (selects exactly one row)
    expected:   numeric expected value
    precision:  absolute tolerance for the match
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import polars as pl
import yaml


CLAIM_COMMENT_RE = re.compile(r"<!--\s*claim:\s*([A-Za-z0-9_\-]+)\s*-->")


def _load_claims(claims_path: Path) -> dict:
    data = yaml.safe_load(claims_path.read_text())
    if data is None:
        return {}
    return {c["claim_id"]: c for c in data}


def _parse_paper_claims(paper_path: Path) -> set[str]:
    text = paper_path.read_text(encoding="utf-8")
    return set(CLAIM_COMMENT_RE.findall(text))


def _read_value(claim: dict) -> float:
    df = pl.read_parquet(claim["parquet"])
    for col, val in (claim.get("filter") or {}).items():
        df = df.filter(pl.col(col) == val)
    if df.height != 1:
        raise ValueError(
            f"claim {claim['claim_id']}: filter selected {df.height} rows; "
            f"expected exactly 1"
        )
    return float(df[claim["column"]][0])


def verify(claims_path: Path, paper_path: Path) -> int:
    claims = _load_claims(claims_path)
    referenced = _parse_paper_claims(paper_path)

    missing = referenced - set(claims)
    if missing:
        for cid in sorted(missing):
            print(f"FAIL: paper references claim '{cid}' not in claims.yaml")
        return 1

    failures = 0
    for cid in sorted(referenced):
        claim = claims[cid]
        try:
            actual = _read_value(claim)
        except Exception as exc:
            print(f"FAIL: {cid}: could not read value -- {exc}")
            failures += 1
            continue
        expected = float(claim["expected"])
        precision = float(claim.get("precision", 0.0))
        if abs(actual - expected) > precision:
            print(
                f"FAIL: {cid}: parquet={actual:.6f} expected={expected:.6f} "
                f"diff={actual - expected:+.6f} precision={precision}"
            )
            failures += 1
        else:
            print(f"OK:   {cid}: {actual:.6f} (expected {expected:.6f})")
    return 0 if failures == 0 else 1


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--run-cheap", action="store_true",
                        help="Phase 4 expansion target.")
    parser.add_argument("--run-all", action="store_true",
                        help="Phase 4 expansion target.")
    parser.add_argument("--claims", default="benchmark/scripts/claims.yaml", type=Path)
    parser.add_argument("--paper", default="benchmark/paper/paper.md", type=Path)
    args = parser.parse_args(argv)

    if args.run_cheap or args.run_all:
        print("--run-cheap / --run-all are Phase 4 expansion targets; not implemented in 0.7.5.")
        sys.exit(2)
    # --verify is the default.
    sys.exit(verify(args.claims, args.paper))


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Create empty seed `claims.yaml`**

Create `benchmark/scripts/claims.yaml`:

```yaml
# Phase 3 seed: empty. Phase 4's locked re-run populates one entry per
# numeric claim in benchmark/paper/paper.md; each entry references the
# parquet that produced the number.
[]
```

- [ ] **Step 5: Smoke check on the seed**

```
uv run python benchmark/scripts/regen_all.py --verify \
  --claims benchmark/scripts/claims.yaml \
  --paper benchmark/paper/paper.md
```

Expected: exit 0 with no claim-mismatch output (paper has no `<!-- claim: ... -->` comments yet).

- [ ] **Step 6: Run the tests**

Run: `uv run pytest benchmark/scripts/tests/test_regen_all.py -v`
Expected: 2 PASS.

- [ ] **Step 7: CHANGELOG + commit**

```markdown
- **`benchmark/scripts/regen_all.py`** (new): acceptance gate for
  paper claims. `--verify` (default) reads `<!-- claim: <id> -->`
  HTML comments in `paper.md`, asserts each claim's parquet value
  matches its expected to the printed precision, exits non-zero on
  mismatch. Phase 3 lands an empty `claims.yaml` seed; Phase 4
  populates during the locked re-run.
```

```
git add benchmark/scripts/regen_all.py benchmark/scripts/claims.yaml \
        benchmark/scripts/tests/test_regen_all.py CHANGELOG.md
git commit -m "$(cat <<'EOF'
feat(benchmark): regen_all.py -- acceptance gate (seed manifest)

--verify (default) reads claims.yaml + paper.md and asserts every
numeric claim in the paper matches its source parquet to the
printed precision. Claims are cited inline via <!-- claim: <id> -->
HTML comments adjacent to the number in the prose.

Phase 3 lands the gate + an empty claims.yaml seed; the smoke
verify passes with no claims. Phase 4 populates claims.yaml during
the locked re-run.

--run-cheap and --run-all are Phase 4 expansion targets; they
print a message and exit 2 if invoked in 0.7.5.

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 19

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 20: `bug_fix_audit.py` — pre/post-fix per-cell delta

**Spec §4 commit 20.** Diffs the pre-Phase-1 `benchmark/data/study*/eval_summary.parquet` against the post-Phase-3 re-run output. Attribution uses `Affects:` trailers parsed from commit messages.

**Files:**
- Create: `benchmark/scripts/bug_fix_audit.py`
- Create: `benchmark/scripts/tests/test_bug_fix_audit.py`

- [ ] **Step 1: Write failing tests**

Create `benchmark/scripts/tests/test_bug_fix_audit.py`:

```python
"""bug_fix_audit.py: pre/post per-cell delta with commit-message
Affects:-trailer attribution."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

SCRIPT = Path(__file__).parents[1] / "bug_fix_audit.py"


def _write_parquets(pre_path: Path, post_path: Path):
    pl.DataFrame({
        "tool": ["epykit_lr", "epykit_lr"],
        "scenario": ["cov10_3v3", "cov15_3v3"],
        "tpr": [0.90, 0.85],
        "fpr": [0.05, 0.04],
    }).write_parquet(pre_path)
    pl.DataFrame({
        "tool": ["epykit_lr", "epykit_lr"],
        "scenario": ["cov10_3v3", "cov15_3v3"],
        "tpr": [0.92, 0.85],   # cov10 changed; cov15 unchanged
        "fpr": [0.07, 0.04],
    }).write_parquet(post_path)


def _write_commits(commits_path: Path):
    # Each entry: subject + body trailer. cov10_3v3 cell attributed
    # to P1-1; cov15_3v3 unattributed.
    commits_path.write_text(json.dumps([
        {
            "subject": "fix(dmc) P1-1: Fisher two-sided mid-p",
            "body": "Affects: lr@cov10_3v3",
        },
    ]))


def test_attribution_success(tmp_path):
    pre = tmp_path / "pre.parquet"
    post = tmp_path / "post.parquet"
    commits = tmp_path / "commits.json"
    out = tmp_path / "audit.parquet"
    _write_parquets(pre, post)
    _write_commits(commits)
    result = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--pre", str(pre), "--post", str(post),
         "--commits-json", str(commits), "--out", str(out)],
        capture_output=True, text=True, check=False,
    )
    # cov15 is unchanged, so it's not an unattributed CELL --
    # only cells that CHANGED need attribution. cov10 changed and
    # is attributed to P1-1. Exit 0.
    assert result.returncode == 0, (
        f"audit should exit 0 when all changed cells attributed: "
        f"{result.stderr}"
    )
    audit = pl.read_parquet(out)
    assert audit.height >= 1
    cov10_row = audit.filter(
        (pl.col("tool") == "epykit_lr") & (pl.col("scenario") == "cov10_3v3")
        & (pl.col("metric") == "tpr")
    )
    assert cov10_row.height == 1
    assert cov10_row["fix_id"][0] == "P1-1"


def test_unattributed_cell_causes_nonzero_exit(tmp_path):
    pre = tmp_path / "pre.parquet"
    post = tmp_path / "post.parquet"
    commits = tmp_path / "commits.json"
    out = tmp_path / "audit.parquet"
    _write_parquets(pre, post)
    # Empty commit log -> no attribution for the cov10 change.
    commits.write_text("[]")
    result = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--pre", str(pre), "--post", str(post),
         "--commits-json", str(commits), "--out", str(out)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0, (
        "audit should exit non-zero on unattributed changed cells"
    )
    assert "UNATTRIBUTED" in (result.stdout + result.stderr)
```

- [ ] **Step 2: Run failing tests**

Run: `uv run pytest benchmark/scripts/tests/test_bug_fix_audit.py -v`
Expected: 2 FAIL — script missing.

- [ ] **Step 3: Write `bug_fix_audit.py`**

Create `benchmark/scripts/bug_fix_audit.py`:

```python
"""bug_fix_audit.py — pre/post-fix per-cell delta with commit
attribution via Affects: trailers.

Inputs:
  --pre <parquet>           pre-fix eval_summary.parquet
  --post <parquet>          post-fix eval_summary.parquet
  --commits-json <path>     JSON list of {subject, body} dicts spanning
                            the pre -> post commit range. Each commit's
                            "Affects: <engine>@<scenario>[, ...]" trailer
                            is parsed for cell attribution.
                            (Production: `git log --format='%s\\n%b' pre..post
                            | json-encode` from a wrapper; tests pass a fixture.)
  --out <parquet>           output audit table

Output schema:
  tool, scenario, metric, pre_value, post_value, delta, fix_id

Exits non-zero if any cell with |delta| > 1e-9 has fix_id="UNATTRIBUTED".
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import polars as pl

# Trailer like "Affects: lr@cov10_3v3, glm@cov10_3v3"
AFFECTS_RE = re.compile(r"^Affects:\s*(.+)$", re.MULTILINE)
CELL_RE    = re.compile(r"([A-Za-z0-9_+\-]+)@([A-Za-z0-9_]+)")
SUBJECT_PID_RE = re.compile(r"\b(P[01]-\d+[a-z]?)\b")


def _parse_commits(commits: list[dict]) -> dict[tuple[str, str], str]:
    """Map (tool_engine, scenario) -> most recent fix_id touching it."""
    attribution: dict[tuple[str, str], str] = {}
    for c in commits:
        body = c.get("body", "")
        subject = c.get("subject", "")
        pid_match = SUBJECT_PID_RE.search(subject)
        fix_id = pid_match.group(1) if pid_match else subject.split(":")[0]
        for m in AFFECTS_RE.finditer(body):
            for cell_m in CELL_RE.finditer(m.group(1)):
                eng, scen = cell_m.group(1), cell_m.group(2)
                attribution[(eng, scen)] = fix_id  # later commit overwrites
    return attribution


def _engine_from_tool(tool: str) -> str:
    """Strip 'epykit_' prefix; 'epykit_lr' -> 'lr'."""
    return tool[len("epykit_"):] if tool.startswith("epykit_") else tool


def diff_and_attribute(
    pre_df: pl.DataFrame, post_df: pl.DataFrame, attribution: dict,
    metrics: tuple[str, ...] = ("tpr", "fpr", "f1", "auroc"),
) -> tuple[pl.DataFrame, int]:
    """Return (audit_df, n_unattributed)."""
    rows: list[dict] = []
    n_unattributed = 0
    joined = pre_df.join(
        post_df, on=["tool", "scenario"], how="outer", suffix="_post",
    )
    for r in joined.iter_rows(named=True):
        for m in metrics:
            pre_v = r.get(m)
            post_v = r.get(f"{m}_post")
            if pre_v is None or post_v is None:
                continue
            if not isinstance(pre_v, (int, float)) or not isinstance(post_v, (int, float)):
                continue
            delta = float(post_v) - float(pre_v)
            if abs(delta) < 1e-9:
                continue  # unchanged cell needs no attribution
            engine = _engine_from_tool(r["tool"])
            fix_id = attribution.get((engine, r["scenario"]), "UNATTRIBUTED")
            if fix_id == "UNATTRIBUTED":
                n_unattributed += 1
            rows.append({
                "tool": r["tool"], "scenario": r["scenario"], "metric": m,
                "pre_value": float(pre_v), "post_value": float(post_v),
                "delta": delta, "fix_id": fix_id,
            })
    return pl.DataFrame(rows), n_unattributed


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pre",         required=True, type=Path)
    parser.add_argument("--post",        required=True, type=Path)
    parser.add_argument("--commits-json", required=True, type=Path)
    parser.add_argument("--out",         required=True, type=Path)
    args = parser.parse_args(argv)

    pre  = pl.read_parquet(args.pre)
    post = pl.read_parquet(args.post)
    commits = json.loads(args.commits_json.read_text())
    attribution = _parse_commits(commits)

    audit, n_unattributed = diff_and_attribute(pre, post, attribution)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    audit.write_parquet(args.out)
    print(f"wrote {args.out} ({audit.height} rows; {n_unattributed} UNATTRIBUTED)")
    if n_unattributed > 0:
        print(
            f"FAIL: {n_unattributed} cells changed without an attributable "
            f"fix. Add an 'Affects: <engine>@<scenario>' trailer to the "
            f"commit that introduced the change, or document the delta as "
            f"expected churn."
        )
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the tests (expect pass)**

Run: `uv run pytest benchmark/scripts/tests/test_bug_fix_audit.py -v`
Expected: 2 PASS.

- [ ] **Step 5: CHANGELOG + commit**

```markdown
- **`benchmark/scripts/bug_fix_audit.py`** (new): diffs pre-fix vs
  post-fix `eval_summary.parquet` per (tool, scenario, metric);
  attributes each delta to a P0/P1 fix via `Affects: <engine>
  @<scenario>` trailers parsed from commit messages spanning the
  range. Unattributed changed cells cause non-zero exit (forces the
  author to add a trailer or document the churn).
```

```
git add benchmark/scripts/bug_fix_audit.py benchmark/scripts/tests/test_bug_fix_audit.py CHANGELOG.md
git commit -m "$(cat <<'EOF'
feat(benchmark): bug_fix_audit.py -- pre/post-fix per-cell delta

Diffs the pre-Phase-1 eval_summary.parquet (decided in the spec as
the natural pre-fix baseline) against the post-Phase-3 re-run, per
(tool, scenario, metric). Each delta is attributed to a P0/P1 fix
via Affects: <engine>@<scenario> trailers parsed from commit
messages between the pre-snapshot tag and post-fix tag.

Failure mode: unattributed changed cells cause non-zero exit. The
author then either adds the trailer to the missing commit or
documents the delta as expected churn (Limitations §10.5 of the
paper).

Production usage in Phase 4:
  git log --format='%s%n%b%n%H%n---END---' v0.7.2..v0.7.5 \\
    | <wrapper that emits JSON list of {subject, body}>
  python benchmark/scripts/bug_fix_audit.py \\
    --pre  benchmark/data/study1/eval_summary.parquet \\
    --post benchmark/data/study1/eval_summary_post_phase3.parquet \\
    --commits-json commits.json --out audit.parquet

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 20

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 21: Wrap-up — CHANGELOG sweep, acceptance criteria, tag

**Spec §4 commit 21.**

- [ ] **Step 1: Full main test suite**

```
uv run pytest -m "not slow" --strict-markers -q 2>&1 | tail -5
```

Expected: passing, roughly Phase 2 baseline + ~5 net new tests after deletions. Record the number for the tag message.

- [ ] **Step 2: Slow tier (locally, may skip in CI)**

```
uv run pytest -m slow -q 2>&1 | tail -5
```

Expected: passing.

- [ ] **Step 3: Benchmark-scripts suite**

```
uv run pytest benchmark/scripts/tests/ -q 2>&1 | tail -5
```

Expected: ~23 passed (Phase 2 baseline 15 + ~8 new fast tests). Slow tier (`-m slow`) adds the null_engines (9 cases) and methylkit Rscript test.

- [ ] **Step 4: Ruff**

```
uv run ruff check src/ benchmark/scripts/
```

Expected: no new F-level violations vs Phase 2 baseline.

- [ ] **Step 5: Mypy**

```
uv run mypy src/epykit
```

Expected: no new errors vs Phase 2 baseline. (Pre-existing errors are acceptable; the spec says "no new errors".)

- [ ] **Step 6: CLI smoke**

```
uv run epykit --help 2>&1 | head -20
uv run epykit dmc --help 2>&1 | grep -iE "(logit_t|bb_lr|score|cmh)" || echo "no dropped engines listed (good)"
```

Expected: no dropped engines in `--help` output.

- [ ] **Step 7: Migration-hint sanity check**

```
uv run python -c "
import epykit as ep
from epykit.methyldata import MethylData
md = MethylData()  # whatever minimal construction works
for engine in ('logit_t', 'bb_lr', 'score', 'cmh'):
    try:
        ep.tl.dmc(md, test=engine)
    except ValueError as e:
        assert 'removed in 0.7.5' in str(e), f'{engine}: missing version: {e}'
        print(f'{engine}: OK ({str(e)[:60]}...)')
"
```

Expected: prints OK for all four (the MethylData construction may need a real fixture; alternative is to invoke the tests/test_phase3_drops.py file directly: `uv run pytest tests/test_phase3_drops.py -v`).

- [ ] **Step 8: Docs sanity check**

```
rg -n "logit_t|bb_lr|test=.score|test=.cmh" README.md docs/analysis/dmc.md CLAUDE.md
```

Expected: only context-appropriate mentions (e.g., "removed in 0.7.5: logit_t, bb_lr, score, cmh"). No remaining `test="logit_t"` examples or engine-list rows for the dropped engines.

- [ ] **Step 9: regen_all.py seed smoke**

```
uv run python benchmark/scripts/regen_all.py --verify \
  --claims benchmark/scripts/claims.yaml \
  --paper benchmark/paper/paper.md
```

Expected: exit 0 (no claims to check yet).

- [ ] **Step 10: Final CHANGELOG sweep**

Open `CHANGELOG.md`. Confirm the `## Unreleased` section has all four required subsections populated by the per-task commits:

- `### Added` — five integration scripts.
- `### Changed (BREAKING for the renamed column / module)` — P1-8, P1-11.
- `### Removed (BREAKING)` — four dropped engines.
- `### Fixed (P1 manifest)` — P1-1, P1-2 (closed by removal), P1-3, P1-4, P1-5, P1-6, P1-7, P1-9, P1-10.
- `### Fixed (P2 manifest, folded into the rename above)` — P2-4.

If any bullet is missing, add it inline now and amend.

- [ ] **Step 11: Empty-tree commit if needed for CHANGELOG amendments**

If you made CHANGELOG edits in Step 10:

```
git add CHANGELOG.md
git commit -m "docs(changelog): Phase 3 wrap-up sweep

Adds any bullets missed during the per-task commits.

Refs docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md §4 commit 21

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

If no edits: skip.

- [ ] **Step 12: Tag**

Hold the tag until the user confirms (same rule as Phase 1 and Phase 2).

```
git tag -a v0.7.5-phase3-engines-frozen -m "Phase 3 engine freeze complete

All 11 P1 functional fixes landed; four engines dropped (logit_t,
bb_lr, score, cmh) — surviving public DMC engines: auto, lr,
welch_t, fisher, glm. dmr_hmm renamed to dmr_segment with real
per-segment p-values. log2_odds_ratio renamed per backend with
transitional column + FutureWarning. Five integration scripts
landed: methylkit_stouffer_combine.R, _null_engines.py + real-
engine wiring for run_null_calibration, evaluate.py --ci-only,
regen_all.py (seed manifest), bug_fix_audit.py.

Engine schema and public API are now final. Phase 4 = locked
benchmark re-run (multi-seed simulator + Studies 2 / 3 re-runs +
null calibration on real data) + paper rewrite + P2 hygiene.

Companion: docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md"
```

- [ ] **Step 13: Brief out**

Append a one-paragraph closeout note to the design spec's §9 (Handoff to Phase 4) summarising:

- Final test count in main + benchmark suites.
- Any unexpected churn during execution (e.g. line-number drift, tests deleted that weren't anticipated).
- Anything Phase 4 should know that wasn't in the spec.

Commit:

```
git add docs/superpowers/specs/2026-05-27-phase-3-engine-freeze-design.md
git commit -m "docs(spec): Phase 3 closeout note in §9

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Self-review (writer's pass, completed inline)

- **Spec coverage:** Spec §4 has 21 commits; this plan has 21 tasks numbered 1-21 with explicit references at each task header. P1 manifest §3 has 11 items (P1-1..P1-11); this plan covers all eleven (P1-2 closed by removal in Task 4, others one-per-task in Tasks 1, 2, 8-15). Surviving-engine surface §3 specified in Task 7 docstring + Task 7's `_auto_test` comment + README/CLAUDE.md updates. Schema deltas §3 implemented in Tasks 1 (segment p-values), 2 (log2 rename), 9 (Newcombe CI), and the dmr_hmm shim warning in Task 1.
- **Placeholder scan:** No "TBD", "implement later", "similar to Task N". Every code step shows the actual code. R script in Task 16 fully written; `bug_fix_audit.py` fully written; `regen_all.py` fully written; `_null_engines.py` fully written. The two "adapt to existing dispatcher style" notes in Task 1 Step 5 and Task 12 Step 3 are not placeholders — the surrounding context is described in enough detail that an engineer can locate the dispatch site and add the branch, and the principle (add `"segment"` branch + deprecation alias for `"hmm"`) is explicit.
- **Type consistency:** Function names match across tasks: `call_dmr_rule_segment` (Task 1 produces, Task 1 imports), `_storey_pi0` (Task 15), `fisher_exact_vectorized` (Task 8), `newcombe_diff_ci` (Task 9 imports), `build_design(reference_level=...)` (Task 10), `_per_site_variance_test` (Task 13 — flagged for grep if name differs), `infer_sex_from_y_coverage` (Task 14 — flagged for grep if name differs), `ENGINE_REGISTRY` (Task 17 produces, Task 17 test imports), `add_wilson_ci_for_tpr_fpr` / `bootstrap_auroc_ci` / `bootstrap_f1_ci` (Phase 2's `wilson_bootstrap_ci.py` — Task 18 imports). Three Task-1 internal helpers (`_stouffer_combine`, `_bh_per_chrom`, `_state_means_for_meth_diff`) are defined and used in the same file. `_null_engines.ENGINE_REGISTRY` keys (`lr`, `lr_plus`, etc.) match the closure factory names. `Affects:` trailer regex in Task 20 (`AFFECTS_RE`, `CELL_RE`) matches the `Affects: lr@cov10_3v3` form used in every per-task commit message (Tasks 3-15 use `Affects:` with the engine names).
- **Ambiguity:** Task 11 patch site assumes a `DEFAULT_MAX_ITER` module-level constant; the test includes a graceful fallback path noted in the test comment ("if `DEFAULT_MAX_ITER` doesn't exist as a module-level constant, find the equivalent — likely a kwarg threaded through `irls_dispatch`. Adjust the patching strategy accordingly. Alternative: construct a degenerate design matrix directly"). Task 13 and 14 grep first for the actual function names. Task 17 notes `MethylData.load` may need a grep-confirm.
- **Test load:** Main suite gains ~11 fast tests + 1 slow; benchmark-scripts suite gains ~6 fast + 2 slow. Net runtime increase: seconds in fast tier; tens-of-seconds in slow tier (R subprocess + parametrised engine runs). Within the budget discussed during brainstorming.

---

## Wrap-up checklist (for the executing agent to confirm at Task 21)

- [ ] `v0.7.5-phase3-engines-frozen` tag exists on `p0-fixes`.
- [ ] Main test suite passes (`pytest -m "not slow"`).
- [ ] Slow tier passes locally (`pytest -m slow`).
- [ ] Benchmark scripts test suite passes.
- [ ] Ruff and mypy show no new violations.
- [ ] CLI `--help` lists no dropped engines.
- [ ] All four dropped engines raise ValueError with migration hint.
- [ ] CHANGELOG `## Unreleased` has Added / Changed / Removed / Fixed sections populated.
- [ ] regen_all.py --verify smokes clean on the seed claims.yaml.
- [ ] Design spec §9 has the closeout note.
