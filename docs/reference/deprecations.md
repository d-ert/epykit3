# Deprecations

epykit 1.2 removes no public API. Every name below still works and emits a
warning that names its replacement. A deprecated name stays until the
maintainer decides on a removal version and a migration path. The warning
text is not that decision.

## Retained names and their replacements

| Deprecated | Use instead | Warning | Schedule |
|---|---|---|---|
| `ep.pp.unite(md, type=...)` | `ep.pp.set_unite_type(md, type=...)`. Same semantics; the name no longer suggests that a join is materialised. | `DeprecationWarning` on call | Announced for 2.0 |
| `import epykit.dmr_hmm` and `call_dmr_hmm` | `from epykit.dmr_segment import call_dmr_rule_segment`. The shim re-exports the same function. | `DeprecationWarning` on import | A future major release |
| `log2_odds_ratio` column in `md.varm["dmc_<test>"]` (NaN-filled since 0.7.5) | `log2_odds_ratio_pooled` for `lr`, `fisher` and `welch_t`; `coef_treatment_log2` for `glm`. | `FutureWarning` on every `tl.dmc` call | A future major release |
| `csv`, `csv_full`, `csv_alpha` keywords on `tl.qc`, `tl.dmc`, `tl.dmr`, `tl.dvc` and `tl.annotate` | `tsv`, `tsv_full`, `tsv_alpha`. When both are given the `tsv*` value wins. | `DeprecationWarning` when any `csv*` value is passed | A future release |
| `--csv`, `--no-csv`, `--csv-alpha`, `--csv-full` CLI flags and `EPYKIT_NO_AUTO_CSV` | `--tsv`, `--no-tsv`, `--tsv-alpha`, `--tsv-full` and `EPYKIT_NO_AUTO_TSV`; see the [CLI reference](../cli/index.md#sibling-tsv-auto-emit). | Logged warning when any `csv` flag is used | A future release |
| `epykit.process_chromosomes_dmc`, `epykit.apply_multiple_testing_correction`, `epykit.empirical_fdr_for_dmc`, `epykit.fisher_exact_vectorized`, `epykit.shrink_meth_diff` (top-level access) | `from epykit.dmc import ...`, or the `ep.tl.dmc` wrapper. | `DeprecationWarning` on attribute access | A future major release |
| `MethylData._analysis_root` | `MethylData.analysis_root` | `DeprecationWarning` on access | Announced for 2.0 |
| `tl.dmc(use_smoothed=True)` pseudo-count path | `tl.dmc(smoothing=True)`, the DSS-style moving average over raw counts. | `DeprecationWarning` on call | A future minor release |
| `tl.dvc(test="bartlett")` | `tl.dvc(test="brown_forsythe")`, which is what runs either way. | `UserWarning` on call | None announced |

The `log2_odds_ratio` column is the only one of these that changes what a
result table looks like. It is always NaN; read the engine's effect column
listed in the [DMC output columns](../analysis/dmc.md#output-columns) instead.

## Keyword names and file formats are separate

`tsv` names the keyword, not the delimiter. The path suffix selects the
format: `tsv="out.tsv"` writes tab-delimited text and `tsv="out.csv"` writes
comma-delimited text. The deprecated `csv` keyword follows the same rule, so
`csv="out.csv"` and `tsv="out.csv"` both request commas, and `csv="out.tsv"`
still writes tabs. See [Tabular exports](../export/tables.md).

## History of the schedules

- `epykit.dmr_hmm` and the `log2_odds_ratio` column were announced for
  removal in 0.8, then 1.1, then 1.2. 1.2 retains both. Their warnings now
  say "a future major release" because no removal version has been decided.
- The top-level DMC names demoted in 1.0 were announced for removal in 1.2.
  1.2 retains the shim; its warning now says "a future major release".
- `pp.unite()` has named 2.0 since 1.0. That schedule is unchanged.
- The `csv*` keyword aliases have said "a future release" since 1.0. That
  wording is unchanged.

The release notes that announced each earlier schedule are preserved as
written in the root
[`CHANGELOG.md`](https://github.com/d-ert/epykit3/blob/main/CHANGELOG.md).
