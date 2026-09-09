# Changelog

The full version history is maintained in the project's
[`CHANGELOG.md`](https://github.com/d-ert/epykit3/blob/main/CHANGELOG.md)
at the repository root. Deprecated names that still work, and what replaces
them, are listed on the [Deprecations](deprecations.md) page.

## 1.1.0 highlights (2026-09-05)

1.1 is a maintenance release: the pre-submission review fixes, a redesigned
HTML report, and a tighter build. No public API was removed.

- **Python >= 3.10.** Python 3.9 reached end of life; the CI matrix is
  `{ubuntu, windows} × {3.10, 3.12, 3.13}`.
- **The `dev` extra is gone.** Contributor tooling lives in
  `[dependency-groups]`: `uv sync --group dev` or
  `pip install -e . --group dev` (pip >= 25.1). The user-facing extras
  (`all`, `report`, `export`, ...) are unchanged, and `uv.lock` is enforced
  in CI.
- **Bismark `.cov` coordinates.** Standard 1-based Bismark coverage is now
  shifted onto the 0-based store correctly. The raw-store manifest moved to
  version 2, so a store built from real `.cov` files is rebuilt on the next
  `read_bismark`.
- **Redesigned HTML report.** `md.report()` / `epykit report` render a
  MultiQC-style dashboard; `self_contained=True` (the default) embeds Plotly
  so the file works offline.
- **Region-level annotation.** `tl.annotate()` applies gene features and
  CpG-island context to the DMR table; `pl.genomic_context_bar()` and
  `pl.cpg_island_pie()` accept `level="dmr"`.
- **Deprecations deferred.** The `log2_odds_ratio` column and the
  `epykit.dmr_hmm` shim, announced for removal in 1.1, are retained and now
  scheduled for 1.2.

## 1.0.0 highlights (2026-06-02)

1.0 is the SemVer-stable release. Three targeted breaking changes shipped at
the major-version cutover; each carries a deprecation shim so 0.7.6 code
continues to run with warnings.

- **`tl.dmc(power_stack="auto")` now engages the full lr+ stack at any
  sample size.** `"auto"`, `"lr+"`, and `True` are aliases that flip all
  four knobs (`neighbour_combine`, `fdr_method="fdr_tsbh"`, `sep_fallback`,
  `dispersion="eb"`). The pre-1.0 conservative behaviour is reachable via
  `power_stack="conservative"`. `"off"` / `False` is the unchanged default.
  See [lr+ Power Stack](../analysis/lr-plus.md).
- **`pp.unite()` renamed to `pp.set_unite_type()`.** The old name suggested
  a verb performing a union; the function only writes a state marker. The
  old name continues to work as a deprecation wrapper through 1.x.
- **`method="hmm"` removed from `tl.dmr`** (deprecated in 0.7.5). Use
  `method="segment"` -- same rule-based 3-state segmenter, honest name.
- **`process_chromosomes_dmc`, `apply_multiple_testing_correction`,
  `empirical_fdr_for_dmc`, `fisher_exact_vectorized`, `shrink_meth_diff`
  removed from the top-level `epykit.*` namespace.** Use the recommended
  `tl.dmc` wrapper, or import explicitly via `from epykit.dmc import ...`.
  A `__getattr__` shim continues to accept the old top-level access pattern
  with a `DeprecationWarning`; removed in 1.2.
- **CLI `dmr --method` default changed from `tile` to `chain_merge`**,
  matching the `tl.dmr` library default. Scripts that relied on the
  implicit tile path must now pass `--method tile` explicitly.
- **CLI `chain_merge` surface.** `epykit dmr --method chain_merge
  --dmc-results <dmc.parquet> [--preset strict|default|permissive]`
  closes the API/CLI parity gap.
- **CLI auto-emits a sibling `.tsv` from `dmc`, `dmr`, `annotate`, and
  `qc-report`.** Override with `--csv PATH`, suppress with `--no-csv`, or
  set `EPYKIT_NO_AUTO_CSV=1` globally. See [CLI Reference](../cli/index.md).
- **`MethylData.analysis_root`** is the new public name for the
  analysis-root attribute; `_analysis_root` is a deprecated alias.

For the full Added / Changed / Fixed / Internal breakdown -- and the
release history before 1.0 -- read the root
[CHANGELOG.md](https://github.com/d-ert/epykit3/blob/main/CHANGELOG.md).
