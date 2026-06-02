# CSV / TSV output for human-readable inspection

**Date:** 2026-06-02
**Status:** Approved (brainstorming)
**Owner:** Deniz

## Motivation

Today every tabular result `epykit` produces lands in Parquet (or stays in
`md.varm` / `md.uns` as a Polars frame). That is fine for downstream
Python code, but it forces a user who just wants to eyeball "which CpGs
came back significant?" to open a notebook, write `polars.read_parquet`,
`.filter`, `.head`. The request is to give every tabular result a
human-readable escape hatch — by default for CLI users, opt-in for the
Python API — so users can `head`, `awk`, or open in Excel without
writing code.

## Scope

CSV/TSV output is added for the following result surfaces:

- **DMC** (`md.varm['dmc_<test>']`) — per-CpG calls. Significant-only by
  default (huge table at whole-genome scale); full table opt-in.
- **DMR** (`md.uns['dmr']`) — region calls. Full table always
  (small enough).
- **Annotated DMC / DMR** — same tables after `tl.annotate`, with the
  added gene / CpG-island context columns.
- **DVC** (`md.varm['dvc']`) — iEVORA differentially-variable CpG calls.
  Same significant-only-by-default rule as DMC.
- **QC** (`md.uns['qc']`) — per-sample QC metrics. Full table.
- **Aggregated regions** — per-region beta matrix produced by
  `aggregate-regions`. Full table.

Out of scope: BED / BigWig / BedGraph stay where they are (already
human-readable enough, and they have BED-specific semantics that don't
map onto a generic CSV writer).

## Format defaults

- **Delimiter:** TSV (`\t`) by default — genomics convention, and safe
  against gene-symbol / contrast strings that may contain commas. CSV
  is opt-in via the file suffix (`.csv`) on the path.
- **Sort order (DMC / DVC, significant-only mode):** `qvalue`
  ascending — most significant hits at the top, matching the
  "eyeball-the-first-50-rows" use case. Full mode and DMR keep genomic
  order (chrom, pos / chrom, start).
- **Columns:** every column from the underlying Polars frame, verbatim
  — no dropping, no renaming. Users can hide columns in Excel.
- **`lr+` neighbour-combined p-values:** when `pvalue_combined` is
  present on the DMC frame, the significant-only filter uses
  `qvalue_combined` (falling back to `qvalue` when absent), to match the
  semantics documented in `CLAUDE.md`.

## Implementation

### `src/epykit/export.py`

Add one writer per result type, all delegating to a private helper that
performs the actual `polars.DataFrame.write_csv(separator=...)`. The
separator is derived from the path suffix (`.csv` → `,`, otherwise
`\t`); writers do not take a separate `format` argument.

```python
def dmc_to_tsv(
    md: MethylData,
    path: str,
    *,
    alpha: float = 0.05,
    full: bool = False,
    test: str | None = None,
) -> str: ...

def dvc_to_tsv(
    md: MethylData,
    path: str,
    *,
    alpha: float = 0.05,
    full: bool = False,
) -> str: ...

def dmr_to_tsv(md: MethylData, path: str) -> str: ...
def qc_to_tsv(md: MethylData, path: str) -> str: ...
def region_matrix_to_tsv(md: MethylData, path: str) -> str: ...
```

All return the resolved output path (`str`), matching the existing
`to_bedgraph` / `dmcs_to_bed` style.

### Python API on `tl.*`

The user-facing entry points grow a single `csv` path parameter plus,
where significance filtering is meaningful, optional `csv_full` and
`csv_alpha` knobs. When `csv` is `None` (the default), nothing is
written — preserving existing behavior.

```python
ep.tl.dmc(md, csv="dmc.significant.tsv")                          # default: q<0.05, qvalue asc
ep.tl.dmc(md, csv="dmc.tsv", csv_full=True)                       # full table, genomic order
ep.tl.dmc(md, csv="dmc.significant.csv", csv_alpha=0.01)          # stricter, CSV format from suffix
ep.tl.dvc(md, csv="dvc.significant.tsv")                          # same shape as dmc
ep.tl.dmr(md, csv="dmr.tsv")                                      # full only
ep.tl.qc(md, csv="qc.tsv")
```

Signatures:

- `tl.dmc(..., csv: str | None = None, csv_full: bool = False, csv_alpha: float = 0.05)`
- `tl.dvc(..., csv, csv_full, csv_alpha)` — same shape
- `tl.dmr(..., csv: str | None = None)` — no filter knobs
- `tl.qc(..., csv: str | None = None)`

Each `tl.*` function delegates to its `export.*_to_tsv` counterpart
after the computation finishes; the underlying writer stays the single
source of truth.

### CLI

Each CLI subcommand that produces a tabular result gains an automatic
sibling TSV write derived from `--output`, plus flags to tune or
suppress it.

| Subcommand          | Auto-emit path                                | Knobs                                                       |
| ------------------- | --------------------------------------------- | ----------------------------------------------------------- |
| `dmc`               | `<stem>.significant.tsv` (+ `<stem>.tsv` if `--csv-full`) | `--no-csv`, `--csv PATH`, `--csv-alpha`, `--csv-full`       |
| `dmr`               | `<stem>.tsv`                                  | `--no-csv`, `--csv PATH`                                    |
| `qc-report`         | `<stem>.tsv` (from `--output qc.html`)        | `--no-csv`, `--csv PATH`                                    |
| `annotate`          | `<save_stem>.<varm_key>.tsv` per annotated frame | `--no-csv`, `--csv-dir DIR`                              |
| `aggregate-regions` | `<stem>.tsv`                                  | `--no-csv`, `--csv PATH`                                    |

The auto-emit synthesizes the TSV path from the `--output` stem; if
`--csv PATH` is supplied, that path wins. Format (TSV vs CSV) is
derived from the suffix on the resolved path.

`annotate` is the odd one out — it mutates a saved `MethylData` rather
than emitting a single result file, so its auto-emit walks the
annotated `varm` / `uns` keys (DMC tables, DMR, DVC) and writes one
TSV per table under `<save_stem>.<key>.tsv`. `--csv-dir DIR` overrides
the destination directory.

### Global opt-out

`EPYKIT_NO_AUTO_CSV=1` env var suppresses the automatic sibling write
across every CLI subcommand. The per-command `--no-csv` flag continues
to work. The env var has no effect on explicit `--csv PATH` or on the
Python API (those are already opt-in).

## Tests (`tests/test_export_csv.py`)

- DMC TSV with significant-only filter (synthetic md with ~10 rows,
  known qvalues — assert correct subset + qvalue-asc sort).
- DMC TSV `full=True` — no filter, genomic order.
- DMR TSV — full table, chrom/start sort.
- Path suffix `.csv` → comma-delimited; `.tsv` → tab.
- Embedded commas in annotated gene-symbol columns survive a CSV
  round-trip (proper quoting via Polars).
- `EPYKIT_NO_AUTO_CSV=1` suppresses the CLI sibling write
  (subprocess-style test invoking the CLI handler).
- `--no-csv` flag suppresses the sibling write.
- `lr+` `pvalue_combined` / `qvalue_combined` columns: the filter
  picks `qvalue_combined` when present and falls back to `qvalue`
  otherwise.
- `tl.dmc(md, csv="...")` produces the same file content as a direct
  `export.dmc_to_tsv(md, "...")` call (delegation contract).

## Non-goals / explicitly deferred

- **No `csv_format` / `--csv-format` parameter.** The path suffix is
  the single source of truth for delimiter choice. Keeps the surface
  small; users who want CSV write `.csv`, users who want TSV write
  `.tsv`.
- **No streaming write.** The DMC significant-only filter cuts the
  table to a Polars-friendly size (typically << 1M rows even at
  whole-genome scale at q<0.05), and the full-table opt-in is a power
  user choice — single `write_csv` is fine. Streaming would only
  matter at hundreds of millions of rows.
- **No new `epykit export csv` subcommand.** The auto-emit at
  `dmc`/`dmr`/`qc-report` time plus the API entry points cover the
  use case; a separate post-hoc CSV exporter would only matter for the
  "I forgot the flag and don't want to re-run" case, which is solved
  by `export.dmc_to_tsv(md, ...)` from a notebook.
- **No CLI flags for the `lr+` knobs.** Already deferred to 1.1 per
  `CLAUDE.md`; CSV output is orthogonal.
