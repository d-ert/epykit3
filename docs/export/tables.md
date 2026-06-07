# Tabular Exports (TSV / CSV)

Every result table epykit produces — DMCs, DMRs, DVCs, and the per-sample QC
summary — can be written to a plain-text **TSV** (or CSV) for opening in Excel /
R / pandas, sharing with collaborators, or feeding a downstream script. There
are three ways to get them, from most convenient to most explicit.

!!! tip "Why TSV by default?"
    epykit writes **tab-delimited** files by default. Tab is the genomics
    convention (BED, GTF, VCF, bedGraph, methylKit are all tab-delimited)
    because gene names and annotation fields routinely contain commas, which
    would corrupt a CSV. Pass a path ending in `.csv` to get comma-delimited
    output instead — the delimiter is always chosen from the file suffix.

## 1. Auto-emit from the analysis call

The **main analyses — `tl.dmc`, `tl.dmr`, `tl.annotate` — write a human-readable
TSV by default.** With no extra arguments, the table lands in
`<analysis_root>/results/` (the same folder `md.save()` uses) as soon as the
analysis finishes:

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38",
                     store_dir="methylstore")   # sets the analysis_root

ep.tl.dmc(md, test="lr")          # -> methylstore/results/dmc.significant.tsv
ep.tl.dmr(md)                     # -> methylstore/results/dmr.tsv
ep.tl.annotate(md, gtf="...")     # -> methylstore/results/dmc_annotated.tsv
```

You can redirect or turn it off per call:

| `tsv=` value | Effect |
|--------------|--------|
| *(omitted)* / `True` | Auto-emit to `<analysis_root>/results/<name>.tsv` (the default) |
| `False` | Don't write anything |
| `"path.tsv"` | Write to this exact path (`.csv` suffix → comma-delimited) |

```python
ep.tl.dmc(md, test="lr", tsv="custom/dmc.tsv")  # explicit path
ep.tl.dmc(md, test="lr", tsv=False)             # skip the auto-emit
```

Auto-emit is also suppressed globally by `EPYKIT_NO_AUTO_TSV=1`, and is skipped
silently when there's no `analysis_root` to anchor on (an in-memory
`MethylData` built without `store_dir`). The auto write is **best-effort** — if
it can't be written for some reason it logs a warning rather than failing the
analysis (the results are already on `md`).

For DMC you also have:

| Argument | Effect |
|----------|--------|
| `tsv_full=True` | Write **every** tested CpG instead of just the significant ones |
| `tsv_alpha=0.01` | Significance threshold for the significant table (default `0.05`) |

!!! note "DVC and QC are opt-in"
    `tl.dvc` and `tl.qc` do **not** auto-emit — they only write when you pass a
    path: `ep.tl.dvc(md, tsv="dvc.tsv")`, `ep.tl.qc(md, tsv="qc.tsv")`. Or grab
    everything at once with `md.export_tables()` (below).

!!! warning "`csv=` is deprecated"
    Earlier versions used `csv=` / `csv_full=` / `csv_alpha=`. Those still work
    but emit a `DeprecationWarning` and will be removed in a future release —
    they were misleading, since the output was always TSV. Switch to `tsv=` /
    `tsv_full=` / `tsv_alpha=`.

## 2. `md.export_tables()` — dump everything at once

After a full analysis, one call writes whichever result tables exist and skips
the rest (nothing is raised for a missing table):

```python
md.export_tables("results/tables")             # significant tables, as .tsv
md.export_tables("results/tables", full=True)  # + full DMC / DVC tables
md.export_tables("results/tables", fmt="csv")  # comma-delimited .csv files
md.export_tables("results/tables", dvc=False)  # per-table switches
```

Files written (only for tables present on `md`):

| File | Source |
|------|--------|
| `<dmc_key>.significant.tsv` | significant DMCs (`+ <dmc_key>.tsv` with `full=True`) |
| `dmr.tsv` | the DMR table (`md.uns["dmr"]`) |
| `dvc.significant.tsv` | significant DVCs (`+ dvc.tsv` with `full=True`) |
| `qc_summary.tsv` | per-sample QC summary (`md.obs`) |

It returns a `dict` of logical table name → path written, so you know exactly
what landed:

```python
written = md.export_tables("results/tables")
# {'dmc_significant': '.../dmc_lr.significant.tsv', 'dmr': '.../dmr.tsv', ...}
```

`ep.export_tables(md, "results/tables")` is the same call as a top-level
function if you prefer.

!!! note "`full=True` can be large"
    The full DMC / DVC tables hold *every* tested CpG — on whole-genome data
    that is tens of millions of rows (a multi-GB file for DVC). The default
    (significant-only) tables are small. Reach for `full=True` deliberately.

## 3. Individual writers — one table, fine control

For a single table with explicit options, the writers live in `ep.export`:

```python
from epykit import export

export.dmc_to_tsv(md, "dmc.tsv", alpha=0.05, full=False)  # significant DMCs
export.dmc_to_tsv(md, "dmc_all.tsv", full=True)           # every tested CpG
export.dmr_to_tsv(md, "dmr.tsv")                          # full DMR table
export.dvc_to_tsv(md, "dvc.tsv", alpha=0.05)              # significant DVCs
export.qc_to_tsv(md, "qc.tsv")                            # per-sample QC (md.obs)
```

These are exactly what `tsv=` and `export_tables()` call under the hood, so the
output is identical. A few details worth knowing:

- **`dmc_to_tsv`** resolves the *annotated* DMC table when annotation has run,
  so gene / feature columns are included. List columns
  (`all_overlapping_genes`, `all_overlapping_features`) are flattened to
  `; `-joined strings so they survive both tab and comma delimiters.
- **`dmc_to_tsv`** understands the `lr+` neighbour-combine output: when
  `qvalue_combined` is present it drives the significance gate and sort.
- The delimiter is derived from the path suffix in every writer
  (`.csv` → comma, otherwise tab).

## From the CLI

The file-producing CLI commands auto-emit a sibling TSV next to their parquet
output by default:

```bash
epykit dmc --output dmc.parquet ...        # also writes dmc.significant.tsv
epykit dmr --output dmr.parquet ...        # also writes dmr.tsv
epykit annotate --output ann.parquet ...   # also writes ann.tsv
epykit qc-report --output-dir qc/ ...      # writes *.tsv next to the parquets
```

Control it with:

| Flag / env | Effect |
|------------|--------|
| `--tsv PATH` | Override the sibling path (`.csv` suffix → comma) |
| `--tsv-full` | (dmc) also write the full, unfiltered table |
| `--tsv-alpha` | (dmc) significance threshold for the significant table |
| `--no-tsv` | Suppress the auto-emit |
| `EPYKIT_NO_AUTO_TSV=1` | Suppress the auto-emit globally |

The legacy `--csv` / `--csv-full` / `--csv-alpha` / `--no-csv` flags and
`EPYKIT_NO_AUTO_CSV` still work but are deprecated (they print a one-line
deprecation notice). Use the `--tsv*` / `EPYKIT_NO_AUTO_TSV` forms.
