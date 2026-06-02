# CSV / TSV Output Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give every tabular epykit result a human-readable TSV/CSV escape hatch so users can `head`, `awk`, or open in Excel without writing Python.

**Architecture:** Add per-result-type writer functions to `src/epykit/export.py` (delegate to a shared private helper that does the Polars `write_csv`). Wire a `csv=` kwarg into `tl.dmc` / `tl.dmr` / `tl.dvc` / `tl.qc` that calls the corresponding writer when set. CLI subcommands (`dmc`, `dmr`, `annotate`, `qc-report`) auto-emit a sibling TSV next to their parquet output unless `--no-csv` or `EPYKIT_NO_AUTO_CSV=1` suppresses it. Delimiter is derived from the path suffix (`.csv` → comma, otherwise tab); no `--csv-format` flag.

**Tech Stack:** Polars `write_csv`, argparse, pytest. The `synth_md_filtered` fixture in `tests/conftest.py` provides a ready-to-DMC MethylData.

**Spec:** [docs/superpowers/specs/2026-06-02-csv-output-design.md](../specs/2026-06-02-csv-output-design.md)

### Spec refinements applied here

While drafting this plan, three small spec inaccuracies surfaced (the CLI handlers don't quite match what the spec assumed). The plan implements the **corrected** behavior; the spec doc can be updated post-hoc:

1. **`annotate`** writes a single parquet at `--output X.parquet`, not a saved MethylData. Auto-emit is therefore a single sibling: `X.parquet` → `X.tsv`. (Spec said "walks the annotated frames inside the saved MethylData" — wrong.)
2. **`qc-report`** writes parquets into `--output-dir`, not an HTML file. Auto-emit writes sibling TSVs into the same `--output-dir`: `global_methylation.tsv` and `coverage_uniformity.tsv`. (Spec said "qc.html → qc.tsv" — wrong; that's the `report` subcommand, which produces HTML only and gets no CSV.)
3. **`tl.qc`** merges per-sample metrics into `md.obs` (not `md.uns['qc']`). `qc_to_tsv(md, path)` therefore writes `md.obs`. Side tables like `qc_sex_check` stay in `md.uns` but are out of scope for the canonical "QC table" export.
4. **DVC has no CLI subcommand today.** `tl.dvc(md, csv=...)` is API-only; no CLI changes needed.
5. **`aggregate-regions`** mutates a saved MethylData (no single tabular `--output`). Deferred — not in this plan.

---

## File Structure

**Created:**
- `tests/test_export_csv.py` — every test for the CSV/TSV writers and the CLI auto-emit.

**Modified:**
- `src/epykit/export.py` — add `_write_table`, `dmc_to_tsv`, `dmr_to_tsv`, `dvc_to_tsv`, `qc_to_tsv` plus update `__all__`.
- `src/epykit/tl.py` — add `csv=`/`csv_full=`/`csv_alpha=` kwargs to `dmc`, `dmr`, `dvc`, `qc`.
- `src/epykit/cli.py` — add `--no-csv` / `--csv PATH` / `--csv-alpha` / `--csv-full` flags to `dmc`, `dmr`, `annotate`, `qc-report` subcommands, and shared auto-emit helper at the bottom of each handler.

---

## Task 1: Shared `_write_table` helper + `dmr_to_tsv`

The DMR writer is the simplest case (full table, no filter, no engine variants), so it anchors the helper.

**Files:**
- Modify: `src/epykit/export.py` (add helper + `dmr_to_tsv` near the bottom, before `__all__`)
- Create: `tests/test_export_csv.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_export_csv.py` with the first two tests:

```python
"""Tests for CSV/TSV export of epykit result tables."""
from __future__ import annotations

import csv
from pathlib import Path

import polars as pl
import pytest

import epykit as ep
from epykit import export


def _make_dmr_frame() -> pl.DataFrame:
    return pl.DataFrame({
        "chrom": ["chr1", "chr1", "chr2"],
        "start": [200, 100, 50],
        "end":   [300, 200, 150],
        "meth_diff": [0.3, -0.4, 0.2],
        "qvalue":  [0.01, 0.02, 0.5],
        "dmr_type": ["hyper", "hypo", "hyper"],
    })


def _stub_md_with_dmr(tmp_path: Path):
    """Build a minimal MethylData carrying a DMR table in md.uns['dmr']."""
    from epykit.methyldata import MethylData
    obs = pl.DataFrame({
        "sample_id": ["s1", "s2"],
        "group": ["case", "ctrl"],
    })
    store = tmp_path / "store"
    store.mkdir()
    md = MethylData(obs=obs, store=str(store))
    md.uns["dmr"] = _make_dmr_frame()
    return md


def test_dmr_to_tsv_writes_full_table_chrom_start_sorted(tmp_path):
    md = _stub_md_with_dmr(tmp_path)
    out = tmp_path / "dmr.tsv"
    export.dmr_to_tsv(md, str(out))

    text = out.read_text(encoding="utf-8")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    assert len(rows) == 3
    assert [(r["chrom"], int(r["start"])) for r in rows] == [
        ("chr1", 100), ("chr1", 200), ("chr2", 50),
    ]


def test_dmr_to_csv_uses_comma_for_csv_suffix(tmp_path):
    md = _stub_md_with_dmr(tmp_path)
    out = tmp_path / "dmr.csv"
    export.dmr_to_tsv(md, str(out))

    text = out.read_text(encoding="utf-8")
    # Header line must contain commas, no tabs.
    header = text.splitlines()[0]
    assert "," in header
    assert "\t" not in header
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: FAIL with `AttributeError: module 'epykit.export' has no attribute 'dmr_to_tsv'`

- [ ] **Step 3: Implement the helper + `dmr_to_tsv`**

In `src/epykit/export.py`, add **before** the existing `__all__` block:

```python
def _separator_for(path: str) -> str:
    """Tab unless the path ends in .csv (case-insensitive)."""
    return "," if str(path).lower().endswith(".csv") else "\t"


def _write_table(df: pl.DataFrame, path: str) -> str:
    """Write `df` to `path` using the suffix-derived delimiter.

    Returns the resolved absolute path. Creates parent directories.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(str(out), separator=_separator_for(path))
    logger.info("Wrote table: %s (%d rows)", out, len(df))
    return str(out.resolve())


def dmr_to_tsv(md: MethylData, path: str) -> str:
    """Write the DMR table (md.uns['dmr']) as TSV/CSV.

    Full table, sorted by (chrom, start). Delimiter is derived from the
    path suffix (.csv -> comma, otherwise tab).
    """
    df = _resolve_dmr_table(md)
    sort_cols = ["chrom", "start"] if "start" in df.columns else ["chrom"]
    return _write_table(df.sort(sort_cols), path)
```

Then update the `__all__` tuple at the bottom of the file:

```python
__all__ = [
    "to_bedgraph",
    "to_bigwig",
    "dmcs_to_bed",
    "dmrs_to_bed",
    "dmr_to_tsv",
]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add src/epykit/export.py tests/test_export_csv.py
git commit -m "feat(export): add dmr_to_tsv writer with suffix-driven delimiter"
```

---

## Task 2: `dmc_to_tsv` with significant/full/alpha + `pvalue_combined` handling

**Files:**
- Modify: `src/epykit/export.py`
- Modify: `tests/test_export_csv.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_export_csv.py`:

```python
def _make_dmc_frame(with_combined: bool = False) -> pl.DataFrame:
    rows = {
        "chrom": ["chr1", "chr1", "chr1", "chr2", "chr2"],
        "pos":   [100, 200, 300, 50, 150],
        "meth_diff": [0.3, -0.4, 0.05, 0.2, -0.1],
        "pvalue": [1e-6, 1e-5, 0.5, 1e-3, 0.4],
        "qvalue": [1e-5, 1e-4, 0.6, 1e-2, 0.5],
    }
    if with_combined:
        # Flip combined values so significance differs from raw qvalue.
        rows["pvalue_combined"] = [0.5, 1e-7, 0.5, 0.5, 1e-9]
        rows["qvalue_combined"] = [0.6, 1e-6, 0.6, 0.6, 1e-8]
    return pl.DataFrame(rows)


def _stub_md_with_dmc(tmp_path: Path, *, key: str = "dmc_lr",
                       with_combined: bool = False):
    from epykit.methyldata import MethylData
    obs = pl.DataFrame({
        "sample_id": ["s1", "s2"], "group": ["case", "ctrl"],
    })
    store = tmp_path / "store"
    store.mkdir()
    md = MethylData(obs=obs, store=str(store))
    md.varm[key] = _make_dmc_frame(with_combined=with_combined)
    md.uns["dmc"] = {"last_key": key}
    return md


def test_dmc_to_tsv_significant_only_qvalue_asc(tmp_path):
    md = _stub_md_with_dmc(tmp_path)
    out = tmp_path / "dmc.significant.tsv"
    export.dmc_to_tsv(md, str(out))  # default alpha=0.05, full=False

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    # 4 rows are q<0.05 (qvalue 1e-5, 1e-4, 1e-2, plus row with q=0.5 / 0.6 dropped)
    assert len(rows) == 3
    # qvalue ascending
    qvalues = [float(r["qvalue"]) for r in rows]
    assert qvalues == sorted(qvalues)


def test_dmc_to_tsv_full_writes_all_rows_genomic_order(tmp_path):
    md = _stub_md_with_dmc(tmp_path)
    out = tmp_path / "dmc.tsv"
    export.dmc_to_tsv(md, str(out), full=True)

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    assert len(rows) == 5
    assert [(r["chrom"], int(r["pos"])) for r in rows] == [
        ("chr1", 100), ("chr1", 200), ("chr1", 300),
        ("chr2", 50),  ("chr2", 150),
    ]


def test_dmc_to_tsv_alpha_override(tmp_path):
    md = _stub_md_with_dmc(tmp_path)
    out = tmp_path / "dmc.strict.tsv"
    export.dmc_to_tsv(md, str(out), alpha=1e-3)

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    # Only rows with qvalue < 1e-3 -> qvalue 1e-5 and 1e-4
    assert len(rows) == 2
    assert all(float(r["qvalue"]) < 1e-3 for r in rows)


def test_dmc_to_tsv_uses_qvalue_combined_when_present(tmp_path):
    md = _stub_md_with_dmc(tmp_path, with_combined=True)
    out = tmp_path / "dmc.combined.tsv"
    export.dmc_to_tsv(md, str(out))

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    # qvalue_combined < 0.05 for rows with combined values 1e-6 and 1e-8.
    assert len(rows) == 2
    qc = sorted(float(r["qvalue_combined"]) for r in rows)
    assert qc == [1e-8, 1e-6]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: 4 new tests FAIL with `AttributeError: ... 'dmc_to_tsv'`

- [ ] **Step 3: Implement `dmc_to_tsv`**

In `src/epykit/export.py`, add **after** `dmr_to_tsv`:

```python
def dmc_to_tsv(
    md: MethylData,
    path: str,
    *,
    alpha: float = 0.05,
    full: bool = False,
    test: str | None = None,
) -> str:
    """Write the DMC table as TSV/CSV.

    Default: significant-only (qvalue < alpha) sorted by qvalue ascending.
    full=True: every row, sorted by (chrom, pos). When the frame carries
    `qvalue_combined` (from the lr+ neighbour-combine knob), that column
    drives the significance filter and the sort; otherwise `qvalue`.

    Delimiter is derived from the path suffix (.csv -> comma, else tab).
    """
    df = _resolve_dmc_table(md, test)
    q_col = "qvalue_combined" if "qvalue_combined" in df.columns else "qvalue"
    p_col = "pvalue_combined" if "pvalue_combined" in df.columns else "pvalue"

    if full:
        out_df = df.sort(["chrom", "pos"])
    else:
        # Significance gate: prefer qvalue (combined or raw); fall back to
        # pvalue when no qvalue column is present at all.
        gate_col = q_col if q_col in df.columns else p_col
        out_df = (
            df.filter(
                pl.col(gate_col).is_not_null()
                & pl.col(gate_col).is_not_nan()
                & (pl.col(gate_col) < alpha)
            )
            .sort(gate_col)
        )
    return _write_table(out_df, path)
```

Update `__all__`:

```python
__all__ = [
    "to_bedgraph",
    "to_bigwig",
    "dmcs_to_bed",
    "dmrs_to_bed",
    "dmr_to_tsv",
    "dmc_to_tsv",
]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: PASS (6 passed)

- [ ] **Step 5: Commit**

```bash
git add src/epykit/export.py tests/test_export_csv.py
git commit -m "feat(export): add dmc_to_tsv with significant/full/alpha + pvalue_combined"
```

---

## Task 3: `dvc_to_tsv` and `qc_to_tsv`

DVC mirrors DMC's significance logic on its own qvalue columns; QC just writes `md.obs`.

**Files:**
- Modify: `src/epykit/export.py`
- Modify: `tests/test_export_csv.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_export_csv.py`:

```python
def _stub_md_with_dvc(tmp_path: Path):
    from epykit.methyldata import MethylData
    obs = pl.DataFrame({
        "sample_id": ["s1", "s2"], "group": ["case", "ctrl"],
    })
    store = tmp_path / "store"
    store.mkdir()
    md = MethylData(obs=obs, store=str(store))
    md.varm["dvc"] = pl.DataFrame({
        "chrom": ["chr1", "chr1", "chr2"],
        "pos":   [100, 200, 50],
        "var_log_ratio": [1.2, 0.1, 0.9],
        "p_variance": [1e-5, 0.6, 1e-3],
        "q_variance": [1e-4, 0.7, 1e-2],
        "is_dvc": [True, False, True],
    })
    return md


def test_dvc_to_tsv_significant_only(tmp_path):
    md = _stub_md_with_dvc(tmp_path)
    out = tmp_path / "dvc.tsv"
    export.dvc_to_tsv(md, str(out))

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    assert len(rows) == 2
    qs = [float(r["q_variance"]) for r in rows]
    assert qs == sorted(qs)


def test_dvc_to_tsv_full(tmp_path):
    md = _stub_md_with_dvc(tmp_path)
    out = tmp_path / "dvc.full.tsv"
    export.dvc_to_tsv(md, str(out), full=True)

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    assert len(rows) == 3


def test_qc_to_tsv_writes_md_obs(tmp_path):
    from epykit.methyldata import MethylData
    obs = pl.DataFrame({
        "sample_id": ["s1", "s2"],
        "group": ["case", "ctrl"],
        "mean_coverage": [12.3, 8.1],
    })
    store = tmp_path / "store"
    store.mkdir()
    md = MethylData(obs=obs, store=str(store))

    out = tmp_path / "qc.tsv"
    export.qc_to_tsv(md, str(out))

    text = out.read_text(encoding="utf-8")
    header = text.splitlines()[0].split("\t")
    assert set(header) == {"sample_id", "group", "mean_coverage"}
    assert "s1\t" in text and "s2\t" in text
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: 3 new tests FAIL.

- [ ] **Step 3: Implement both writers**

In `src/epykit/export.py`, add **after** `dmc_to_tsv`:

```python
def dvc_to_tsv(
    md: MethylData,
    path: str,
    *,
    alpha: float = 0.05,
    full: bool = False,
) -> str:
    """Write the DVC table (md.varm['dvc']) as TSV/CSV.

    Default: significant-only (q_variance < alpha) sorted by q_variance
    ascending. full=True keeps every row in (chrom, pos) order.
    Delimiter is derived from the path suffix.
    """
    df = md.varm.get("dvc")
    if df is None or len(df) == 0:
        raise ValueError(
            "No DVC results on this MethylData. Run ep.tl.dvc(md) first."
        )
    if full:
        sort_cols = ["chrom", "pos"] if "pos" in df.columns else ["chrom"]
        out_df = df.sort(sort_cols)
    else:
        gate = "q_variance" if "q_variance" in df.columns else "p_variance"
        out_df = (
            df.filter(
                pl.col(gate).is_not_null()
                & pl.col(gate).is_not_nan()
                & (pl.col(gate) < alpha)
            )
            .sort(gate)
        )
    return _write_table(out_df, path)


def qc_to_tsv(md: MethylData, path: str) -> str:
    """Write the per-sample QC summary (md.obs) as TSV/CSV.

    After ep.tl.qc(md), md.obs carries the per-sample metrics joined onto
    the existing samplesheet columns. This writer dumps it verbatim.
    Delimiter is derived from the path suffix.
    """
    return _write_table(md.obs, path)
```

Update `__all__`:

```python
__all__ = [
    "to_bedgraph",
    "to_bigwig",
    "dmcs_to_bed",
    "dmrs_to_bed",
    "dmr_to_tsv",
    "dmc_to_tsv",
    "dvc_to_tsv",
    "qc_to_tsv",
]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: PASS (9 passed)

- [ ] **Step 5: Commit**

```bash
git add src/epykit/export.py tests/test_export_csv.py
git commit -m "feat(export): add dvc_to_tsv and qc_to_tsv writers"
```

---

## Task 4: Wire `csv=` kwarg into `tl.dmc` and `tl.dvc`

Both use the same `csv` / `csv_full` / `csv_alpha` triple.

**Files:**
- Modify: `src/epykit/tl.py` (function signatures + bodies of `dmc` and `dvc`)
- Modify: `tests/test_export_csv.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_export_csv.py`:

```python
def test_tl_dmc_csv_kwarg_writes_file(tmp_path, synth_md_filtered):
    out = tmp_path / "dmc.significant.tsv"
    ep.tl.dmc(synth_md_filtered, test="lr", csv=str(out))

    assert out.exists()
    text = out.read_text(encoding="utf-8")
    # Header from md.varm['dmc_lr']: chrom, pos, ... qvalue ...
    header = text.splitlines()[0].split("\t")
    assert "chrom" in header and "pos" in header and "qvalue" in header


def test_tl_dmc_csv_full_writes_every_row(tmp_path, synth_md_filtered):
    full_out = tmp_path / "dmc.tsv"
    ep.tl.dmc(synth_md_filtered, test="lr", csv=str(full_out), csv_full=True)

    n_full = len(full_out.read_text(encoding="utf-8").splitlines()) - 1
    n_varm = len(synth_md_filtered.varm["dmc_lr"])
    assert n_full == n_varm


def test_tl_dvc_csv_kwarg_writes_file(tmp_path, synth_md_filtered):
    ep.tl.dvc(synth_md_filtered, test="bartlett")
    out = tmp_path / "dvc.significant.tsv"
    ep.tl.dvc(synth_md_filtered, test="bartlett", csv=str(out))
    # File should exist; may be empty (no significant DVCs in the fixture)
    assert out.exists()
    header = out.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "chrom" in header and "pos" in header
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_export_csv.py -v -k "tl_dmc or tl_dvc"`
Expected: FAIL with `TypeError: dmc() got an unexpected keyword argument 'csv'`

- [ ] **Step 3: Modify `tl.dmc` to accept the new kwargs**

In `src/epykit/tl.py`, modify the `def dmc(` signature at line 252. Add the three new kwargs as the last keyword-only parameters (after `reference_level`):

```python
    reference_level: str | None = None,
    csv: str | None = None,
    csv_full: bool = False,
    csv_alpha: float = 0.05,
) -> None:
```

At the end of `tl.dmc`'s body (just before the function returns — `tl.dmc` returns `None` implicitly after writing varm), add:

```python
    if csv is not None:
        from .export import dmc_to_tsv
        dmc_to_tsv(md, csv, alpha=csv_alpha, full=csv_full)
```

Find the right insertion point: scroll to the end of the `dmc` function body (it ends when the next `def ` starts — likely around line 790 before `_run_dmc_contrast`). The insertion goes inside `tl.dmc`, just before its closing implicit `return`.

- [ ] **Step 4: Modify `tl.dvc` the same way**

In `src/epykit/tl.py`, modify the `def dvc(` signature at line 1528 — add three kwargs at the end:

```python
    backend: str = "sequential",
    n_workers: int | None = None,
    csv: str | None = None,
    csv_full: bool = False,
    csv_alpha: float = 0.05,
) -> None:
```

At the end of `tl.dvc`'s body (after the `md.uns["dvc"] = { ... }` block around line 1586), add:

```python
    if csv is not None:
        from .export import dvc_to_tsv
        dvc_to_tsv(md, csv, alpha=csv_alpha, full=csv_full)
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/test_export_csv.py -v -k "tl_dmc or tl_dvc"`
Expected: PASS (3 passed)

Also run the full file to confirm no regressions:

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: PASS (12 passed)

- [ ] **Step 6: Commit**

```bash
git add src/epykit/tl.py tests/test_export_csv.py
git commit -m "feat(tl): csv=/csv_full=/csv_alpha= kwargs on tl.dmc and tl.dvc"
```

---

## Task 5: Wire `csv=` kwarg into `tl.dmr` and `tl.qc`

Simpler shape — only the `csv` path (no full/alpha knobs).

**Files:**
- Modify: `src/epykit/tl.py`
- Modify: `tests/test_export_csv.py`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_export_csv.py`:

```python
def test_tl_dmr_csv_kwarg_writes_file(tmp_path, synth_md_filtered):
    ep.tl.dmc(synth_md_filtered, test="lr")
    out = tmp_path / "dmr.tsv"
    ep.tl.dmr(synth_md_filtered, method="chain_merge", csv=str(out))

    assert out.exists()
    header = out.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "chrom" in header and "start" in header and "end" in header


def test_tl_qc_csv_kwarg_writes_md_obs(tmp_path, synth_md_filtered):
    out = tmp_path / "qc.tsv"
    ep.tl.qc(synth_md_filtered, csv=str(out))

    assert out.exists()
    text = out.read_text(encoding="utf-8")
    header = text.splitlines()[0].split("\t")
    assert "sample_id" in header
    # Body has one row per sample
    n_rows = len(text.splitlines()) - 1
    assert n_rows == len(synth_md_filtered.obs)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_export_csv.py -v -k "tl_dmr or tl_qc"`
Expected: FAIL with `TypeError: ... unexpected keyword argument 'csv'`

- [ ] **Step 3: Modify `tl.dmr`**

In `src/epykit/tl.py`, modify the `def dmr(` signature at line 927 — add `csv` after the last existing kwarg (`merge_adjacent: bool = True`):

```python
    merge_adjacent: bool = True,
    csv: str | None = None,
) -> None:
```

At the end of `tl.dmr`'s body (after `md.uns["dmr"] = ...` is set — search for that assignment within the function; it's the place where the result lands), add:

```python
    if csv is not None:
        from .export import dmr_to_tsv
        dmr_to_tsv(md, csv)
```

- [ ] **Step 4: Modify `tl.qc`**

In `src/epykit/tl.py`, modify the `def qc(` signature at line 146 — add `csv` after `expected_sex_col`:

```python
    expected_sex_col: str | None = None,
    csv: str | None = None,
) -> None:
```

At the very end of `tl.qc`'s body (after the final `md.obs = obs` reassignment / return — find the last statement of the function before the next `def `), add:

```python
    if csv is not None:
        from .export import qc_to_tsv
        qc_to_tsv(md, csv)
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/test_export_csv.py -v -k "tl_dmr or tl_qc"`
Expected: PASS (2 passed)

Also run the full file:

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: PASS (14 passed)

- [ ] **Step 6: Commit**

```bash
git add src/epykit/tl.py tests/test_export_csv.py
git commit -m "feat(tl): csv= kwarg on tl.dmr and tl.qc"
```

---

## Task 6: CLI auto-emit for `epykit dmc`

Adds `--no-csv`, `--csv PATH`, `--csv-alpha`, `--csv-full` flags + sibling-TSV auto-emit at end of `_cmd_dmc`.

**Files:**
- Modify: `src/epykit/cli.py`
- Modify: `tests/test_export_csv.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_export_csv.py`:

```python
def test_cli_dmc_auto_emits_sibling_significant_tsv(tmp_path, synth_bundle, monkeypatch):
    """`epykit dmc --output X.parquet` writes X.significant.tsv next to it."""
    import sys
    from epykit.cli import main

    out_parquet = tmp_path / "dmc.parquet"
    sibling_tsv = tmp_path / "dmc.significant.tsv"

    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", synth_bundle.store_root,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(out_parquet),
        "--test", "lr",
    ])
    main()

    assert out_parquet.exists()
    assert sibling_tsv.exists()
    # Should contain a header line and some data rows
    lines = sibling_tsv.read_text(encoding="utf-8").splitlines()
    assert len(lines) >= 1
    assert "chrom" in lines[0]


def test_cli_dmc_no_csv_suppresses_sibling(tmp_path, synth_bundle, monkeypatch):
    import sys
    from epykit.cli import main

    out_parquet = tmp_path / "dmc.parquet"
    sibling_tsv = tmp_path / "dmc.significant.tsv"

    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", synth_bundle.store_root,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(out_parquet),
        "--test", "lr",
        "--no-csv",
    ])
    main()

    assert out_parquet.exists()
    assert not sibling_tsv.exists()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_export_csv.py -v -k "cli_dmc"`
Expected: FAIL — first test fails because the sibling TSV is not created, second fails on `--no-csv` being an unknown flag.

- [ ] **Step 3: Add the CLI flags + auto-emit helper**

In `src/epykit/cli.py`, add a small helper near the top of the file (after the existing `_add_min_samples_args` helper around line 28) so every subcommand can share it:

```python
def _auto_csv_path(parquet_path: str, *, suffix: str = "") -> str:
    """Derive a sibling .tsv path from a --output parquet path.

    `dmc.parquet` -> `dmc.significant.tsv` (suffix=".significant")
    `dmr.parquet` -> `dmr.tsv`             (suffix="")
    Strips a `.parquet` extension if present; otherwise appends.
    """
    p = Path(parquet_path)
    stem = p.stem if p.suffix.lower() == ".parquet" else p.name
    return str(p.with_name(f"{stem}{suffix}.tsv"))


def _csv_suppressed(args) -> bool:
    """True if the user opted out of the CLI auto-emit."""
    if getattr(args, "no_csv", False):
        return True
    if os.environ.get("EPYKIT_NO_AUTO_CSV") in ("1", "true", "True"):
        return True
    return False
```

Also add `import os` at the top of `cli.py` if it's not already imported (it likely is — check first).

Find the `# dmc` argparse block around line 540 in `_cmd_dmc`'s parser setup (`p_dmc.add_argument(...)` cluster) and add these flags after the existing `--allow-n1`:

```python
    p_dmc.add_argument(
        "--no-csv", action="store_true", dest="no_csv", default=False,
        help="Suppress the sibling .significant.tsv auto-emit.",
    )
    p_dmc.add_argument(
        "--csv", dest="csv_path", default=None,
        help=(
            "Override sibling TSV/CSV path. Suffix .csv selects comma "
            "delimiter; otherwise tab. Implies the file is written."
        ),
    )
    p_dmc.add_argument(
        "--csv-alpha", dest="csv_alpha", type=float, default=0.05,
        help="qvalue threshold for significant-only CSV. Default 0.05.",
    )
    p_dmc.add_argument(
        "--csv-full", dest="csv_full", action="store_true", default=False,
        help="Also write the full (unfiltered) TSV next to the parquet.",
    )
```

In `_cmd_dmc` (around line 124), at the very end (after the existing `print(f"  Significant (q<0.05): ...")` line), add the auto-emit logic:

```python
    if not _csv_suppressed(args):
        from .methyldata import MethylData
        from .export import dmc_to_tsv
        # Build a transient MethylData carrying just the dmc result so the
        # writer can re-use the same delegation path the API uses.
        obs = pl.DataFrame({"sample_id": treatment_samples + control_samples})
        md_tmp = MethylData(obs=obs, store=str(args.methylstore))
        md_tmp.varm["dmc_lr"] = results
        md_tmp.uns["dmc"] = {"last_key": "dmc_lr"}

        sig_path = args.csv_path or _auto_csv_path(
            args.output, suffix=".significant"
        )
        dmc_to_tsv(md_tmp, sig_path, alpha=args.csv_alpha)
        print(f"  Significant CSV:    {sig_path}")
        if args.csv_full:
            full_path = _auto_csv_path(args.output)
            dmc_to_tsv(md_tmp, full_path, full=True)
            print(f"  Full CSV:           {full_path}")
```

(The `import polars as pl` at top of cli.py should already be present from other handlers; if not, add it.)

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_export_csv.py -v -k "cli_dmc"`
Expected: PASS (2 passed)

Run the full test file to confirm no regressions:

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: PASS (16 passed)

- [ ] **Step 5: Commit**

```bash
git add src/epykit/cli.py tests/test_export_csv.py
git commit -m "feat(cli): auto-emit sibling significant.tsv from epykit dmc"
```

---

## Task 7: CLI auto-emit for `epykit dmr`

Simpler than dmc — no significance filter (full table only) and no `--csv-full`.

**Files:**
- Modify: `src/epykit/cli.py`
- Modify: `tests/test_export_csv.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_export_csv.py`:

```python
def test_cli_dmr_auto_emits_sibling_tsv(tmp_path, synth_bundle, monkeypatch):
    import sys
    from epykit.cli import main

    # First make a DMC parquet, since chain_merge consumes one.
    dmc_parquet = tmp_path / "dmc.parquet"
    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", synth_bundle.store_root,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(dmc_parquet),
        "--test", "lr",
        "--no-csv",  # don't pollute tmp_path with the dmc sibling
    ])
    main()

    dmr_parquet = tmp_path / "dmr.parquet"
    sibling = tmp_path / "dmr.tsv"
    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmr",
        "--method", "chain_merge",
        "--dmc-results", str(dmc_parquet),
        "--output", str(dmr_parquet),
        "--preset", "permissive",
    ])
    main()

    assert dmr_parquet.exists()
    assert sibling.exists()
    header = sibling.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "chrom" in header and "start" in header
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_export_csv.py -v -k "cli_dmr"`
Expected: FAIL — sibling TSV is not created.

- [ ] **Step 3: Add the CLI flags + auto-emit**

Find the `# dmr` argparse block around line 603 (the `--alpha`/`--min-abs-meth-diff` cluster ends at `p_dmr.set_defaults(func=_cmd_dmr)`). Add **before** the `set_defaults` line:

```python
    p_dmr.add_argument(
        "--no-csv", action="store_true", dest="no_csv", default=False,
        help="Suppress the sibling .tsv auto-emit.",
    )
    p_dmr.add_argument(
        "--csv", dest="csv_path", default=None,
        help="Override sibling TSV/CSV path. .csv suffix -> comma delim.",
    )
```

In `_cmd_dmr` (around line 204), at the very end (after the existing `print(dmr_results.head(10))` line — or after the final `print(f"DMR results written to {args.output}")` block), add:

```python
    if not _csv_suppressed(args):
        from .methyldata import MethylData
        from .export import dmr_to_tsv
        obs = pl.DataFrame({"sample_id": []})
        md_tmp = MethylData(obs=obs, store="")
        md_tmp.uns["dmr"] = dmr_results

        tsv_path = args.csv_path or _auto_csv_path(args.output)
        dmr_to_tsv(md_tmp, tsv_path)
        print(f"DMR CSV: {tsv_path}")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_export_csv.py -v -k "cli_dmr"`
Expected: PASS (1 passed)

Run the full file:

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: PASS (17 passed)

- [ ] **Step 5: Commit**

```bash
git add src/epykit/cli.py tests/test_export_csv.py
git commit -m "feat(cli): auto-emit sibling .tsv from epykit dmr"
```

---

## Task 8: CLI auto-emit for `epykit annotate`

`annotate --output X.parquet` reads a DMC/DMR parquet, annotates it, writes a parquet — single sibling: `X.tsv` (full table, no filter — keep all annotated rows so the user can see context).

**Files:**
- Modify: `src/epykit/cli.py`
- Modify: `tests/test_export_csv.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_export_csv.py`:

```python
def test_cli_annotate_auto_emits_sibling_tsv(tmp_path, synth_bundle, monkeypatch):
    import sys
    from epykit.cli import main

    # Build a tiny DMC parquet to feed into annotate.
    dmc_parquet = tmp_path / "dmc.parquet"
    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", synth_bundle.store_root,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(dmc_parquet),
        "--test", "lr",
        "--no-csv",
    ])
    main()

    annotated = tmp_path / "annotated.parquet"
    sibling   = tmp_path / "annotated.tsv"
    # No --gtf / --cpg-islands -> annotate is a no-op pass-through, but the
    # sibling TSV must still be written.
    monkeypatch.setattr(sys, "argv", [
        "epykit", "annotate",
        "--input", str(dmc_parquet),
        "--output", str(annotated),
    ])
    main()

    assert annotated.exists()
    assert sibling.exists()
    header = sibling.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "chrom" in header and "pos" in header
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_export_csv.py -v -k "cli_annotate"`
Expected: FAIL — sibling TSV is not created.

- [ ] **Step 3: Add the CLI flags + auto-emit**

Find the `# annotate` argparse block around line 707. Add **before** `p_ann.set_defaults(func=_cmd_annotate)`:

```python
    p_ann.add_argument(
        "--no-csv", action="store_true", dest="no_csv", default=False,
        help="Suppress the sibling .tsv auto-emit.",
    )
    p_ann.add_argument(
        "--csv", dest="csv_path", default=None,
        help="Override sibling TSV/CSV path. .csv suffix -> comma delim.",
    )
```

In `_cmd_annotate` (around line 318), at the very end (after `print(f"Annotated results written to {args.output}")`), add:

```python
    if not _csv_suppressed(args):
        tsv_path = args.csv_path or _auto_csv_path(args.output)
        _write_table_local(sites, tsv_path)
        print(f"Annotated CSV: {tsv_path}")
```

And add a tiny local helper near the top of `cli.py` (next to `_auto_csv_path`), since the annotate handler holds a raw Polars frame rather than a MethylData:

```python
def _write_table_local(df, path: str) -> str:
    """Same as export._write_table but operates on a raw Polars frame.

    Kept local to cli.py so handlers that hold a frame directly don't have
    to wrap it in a stub MethylData just to write a TSV.
    """
    from pathlib import Path as _P
    out = _P(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    sep = "," if str(path).lower().endswith(".csv") else "\t"
    df.write_csv(str(out), separator=sep)
    return str(out.resolve())
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_export_csv.py -v -k "cli_annotate"`
Expected: PASS (1 passed)

Run the full file:

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: PASS (18 passed)

- [ ] **Step 5: Commit**

```bash
git add src/epykit/cli.py tests/test_export_csv.py
git commit -m "feat(cli): auto-emit sibling .tsv from epykit annotate"
```

---

## Task 9: CLI auto-emit for `epykit qc-report`

`qc-report --output-dir DIR` writes parquets into DIR. Add sibling TSVs in the same DIR.

**Files:**
- Modify: `src/epykit/cli.py`
- Modify: `tests/test_export_csv.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_export_csv.py`:

```python
def test_cli_qc_report_auto_emits_sibling_tsvs(tmp_path, synth_bundle, monkeypatch):
    import sys
    from epykit.cli import main

    out_dir = tmp_path / "qc_out"
    samples = ",".join(synth_bundle.treatment_ids + synth_bundle.control_ids)

    monkeypatch.setattr(sys, "argv", [
        "epykit", "qc-report",
        "--methylstore", synth_bundle.store_root,
        "--samples", samples,
        "--output-dir", str(out_dir),
    ])
    main()

    assert (out_dir / "global_methylation.parquet").exists()
    assert (out_dir / "global_methylation.tsv").exists()
    assert (out_dir / "coverage_uniformity.parquet").exists()
    assert (out_dir / "coverage_uniformity.tsv").exists()

    header = (out_dir / "global_methylation.tsv").read_text(
        encoding="utf-8"
    ).splitlines()[0]
    assert "sample" in header or "context" in header
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_export_csv.py -v -k "cli_qc_report"`
Expected: FAIL — sibling TSVs are not created.

- [ ] **Step 3: Add the CLI flags + auto-emit**

Find the `# qc-report` argparse block around line 726. Add **before** `p_qc.set_defaults(func=_cmd_qc_report)`:

```python
    p_qc.add_argument(
        "--no-csv", action="store_true", dest="no_csv", default=False,
        help="Suppress the sibling .tsv auto-emit alongside the parquets.",
    )
```

In `_cmd_qc_report` (around line 344), modify the two write blocks. The first block currently:

```python
    if args.output_dir:
        out = Path(args.output_dir)
        out.mkdir(parents=True, exist_ok=True)
        meth_report.write_parquet(str(out / "global_methylation.parquet"))
```

becomes:

```python
    if args.output_dir:
        out = Path(args.output_dir)
        out.mkdir(parents=True, exist_ok=True)
        meth_report.write_parquet(str(out / "global_methylation.parquet"))
        if not _csv_suppressed(args):
            _write_table_local(meth_report, str(out / "global_methylation.tsv"))
```

The second block:

```python
    if cov_frames and args.output_dir:
        combined = pl.concat(cov_frames)
        combined.write_parquet(str(Path(args.output_dir) / "coverage_uniformity.parquet"))
        print(f"\nQC reports written to {args.output_dir}")
```

becomes:

```python
    if cov_frames and args.output_dir:
        combined = pl.concat(cov_frames)
        combined.write_parquet(str(Path(args.output_dir) / "coverage_uniformity.parquet"))
        if not _csv_suppressed(args):
            _write_table_local(
                combined,
                str(Path(args.output_dir) / "coverage_uniformity.tsv"),
            )
        print(f"\nQC reports written to {args.output_dir}")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_export_csv.py -v -k "cli_qc_report"`
Expected: PASS (1 passed)

Run the full file:

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: PASS (19 passed)

- [ ] **Step 5: Commit**

```bash
git add src/epykit/cli.py tests/test_export_csv.py
git commit -m "feat(cli): auto-emit sibling .tsv from epykit qc-report"
```

---

## Task 10: `EPYKIT_NO_AUTO_CSV` env-var opt-out + full-suite regression check

The env var was already wired into `_csv_suppressed` in Task 6, but we never tested it end-to-end. This task locks the contract and runs the full epykit test suite to confirm no regressions.

**Files:**
- Modify: `tests/test_export_csv.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_export_csv.py`:

```python
def test_env_var_suppresses_cli_auto_emit(tmp_path, synth_bundle, monkeypatch):
    import sys
    from epykit.cli import main

    out_parquet = tmp_path / "dmc.parquet"
    sibling = tmp_path / "dmc.significant.tsv"

    monkeypatch.setenv("EPYKIT_NO_AUTO_CSV", "1")
    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", synth_bundle.store_root,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(out_parquet),
        "--test", "lr",
    ])
    main()

    assert out_parquet.exists()
    assert not sibling.exists(), (
        "EPYKIT_NO_AUTO_CSV=1 must suppress the sibling write"
    )


def test_explicit_csv_path_wins_over_auto_emit_name(tmp_path, synth_bundle, monkeypatch):
    """`--csv` flag overrides the derived `<stem>.significant.tsv` name."""
    import sys
    from epykit.cli import main

    out_parquet = tmp_path / "dmc.parquet"
    explicit = tmp_path / "my_hits.csv"   # .csv suffix -> comma delim
    default_sibling = tmp_path / "dmc.significant.tsv"

    monkeypatch.setattr(sys, "argv", [
        "epykit", "dmc",
        "--methylstore", synth_bundle.store_root,
        "--samplesheet", synth_bundle.samplesheet,
        "--treatment-group", "treatment",
        "--control-group", "control",
        "--output", str(out_parquet),
        "--test", "lr",
        "--csv", str(explicit),
    ])
    main()

    assert explicit.exists()
    assert not default_sibling.exists()
    # Comma delimiter because suffix is .csv
    header = explicit.read_text(encoding="utf-8").splitlines()[0]
    assert "," in header and "\t" not in header
```

- [ ] **Step 2: Run tests to verify they pass**

Run: `uv run pytest tests/test_export_csv.py -v -k "env_var or explicit_csv_path"`
Expected: PASS (2 passed) — these should already pass because Task 6 wired both code paths. If either fails, the bug is in `_csv_suppressed` (env path) or in `_cmd_dmc`'s `args.csv_path or _auto_csv_path(...)` precedence.

- [ ] **Step 3: Run the full export test file**

Run: `uv run pytest tests/test_export_csv.py -v`
Expected: PASS (21 passed)

- [ ] **Step 4: Run the full epykit test suite (non-slow) for regressions**

Run: `uv run pytest -m "not slow" --strict-markers -ra`
Expected: All previously passing tests continue to pass. If something fails, the most likely culprits are: (a) the new kwargs on tl.* broke a caller that uses `**kwargs` introspection, or (b) the CLI handlers' end-of-function additions broke the existing stdout-line assertions in other CLI tests. Inspect and fix before committing.

- [ ] **Step 5: Commit**

```bash
git add tests/test_export_csv.py
git commit -m "test(export): EPYKIT_NO_AUTO_CSV env var + explicit --csv override"
```

---

## Done

Final sanity:

- [ ] `uv run pytest -m "not slow" --strict-markers -ra` is green.
- [ ] `uv run ruff check src/` produces no new findings.
- [ ] `uv run mypy src/epykit` produces no new errors.
- [ ] `git log --oneline -10` shows the 10 commits from this plan in order.
