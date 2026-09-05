# Contributing to epykit

Thanks for helping. This page is the short version of what a change needs
before it can merge. `CLAUDE.md` holds the architecture rules in more depth;
`docs/advanced/architecture.md` is the canonical engine map.

## Set up

```bash
uv sync --frozen --extra dev --extra all     # matches CI; --frozen keeps uv.lock untouched
uv run --frozen pytest -m "not slow" --strict-markers -ra
```

Plain `pip install -e ".[dev,all]"` works too. Python 3.9 to 3.12 on Linux,
macOS and Windows are all CI targets; `pysam`-based extras (`bam`,
`methylkit`) and `pyBigWig` (`export`) have no Windows wheels and are gated
in `pyproject.toml`.

## Gates

Every pull request must pass all of these locally before review:

| Gate | Command | Why it exists |
|---|---|---|
| Fast tier | `uv run --frozen pytest -m "not slow" --strict-markers -ra` | The CI matrix. About two minutes. |
| Slow tier | `uv run --frozen pytest -m slow --strict-markers -ra` | Null calibration and accuracy tests against a synthetic truth table. About one minute. |
| Lint and types | `uv run --frozen ruff check src/` and `uv run --frozen mypy src/epykit` | Baseline is pyflakes plus mypy; see `pyproject.toml` for the ratchet plan. |
| Engine hashes | `uv run --frozen python benchmark/scripts/regen_small.py` | Hashes `lr` and `lr+` output on a fixed simulator slice and diffs against the committed reference. Any change to a p-value, q-value or `meth_diff` at eight decimals fails it. |

If a change is *meant* to move engine output, re-snapshot with
`regen_small.py --update`, commit the hash file in the same change, and say
so in the commit message and in `CHANGELOG.md`.

Tests slower than about five seconds carry `@pytest.mark.slow`.
`--strict-markers` is on, so unregistered marks fail collection.

## Rules that are load-bearing

- **Never load the whole genome into one frame.** Iterate per chromosome or
  use `pl.scan_parquet` on a partition glob. Anything under `tl.dmc` must
  keep peak memory at O(largest chromosome).
- **Library code never calls `print()`.** Use
  `logger = logging.getLogger(__name__)`. Only `cli.py` prints, and only the
  final result line. A test enforces this.
- **Preprocessing state is derived.** A new `pp.*` step appends to
  `md.uns["_store_history"]` and repoints `md.store`. Do not add boolean
  flags.
- **Every DMC engine emits the canonical schema.** Extra columns are fine;
  renaming or dropping a canonical column is a breaking change.
- **Coordinates are 0-based in the store.** GTF is 1-based closed and is
  converted on parse; BED is 0-based half-open; Bismark `.cov` is 1-based
  and is shifted on ingest.
- **Deprecate, then remove.** Emit `DeprecationWarning` or `FutureWarning`
  with the removal version in the message. Tests assert on these warnings,
  so do not silence them globally.
- **Temporary files go through `tempfile`.** `ep.set_tmp_dir()` redirects
  them; never hardcode a path.
- **Do not change `lr+` defaults** (`power_stack`, `neighbour_combine`,
  `fdr_method`, `sep_fallback`, `dispersion`) without re-running the
  ablations under `benchmark/`.

## Commits and pull requests

Commit messages follow the existing `type(scope): summary` style, for
example `fix(dmr): honor --min-cpgs on chain_merge`. Keep one logical change
per commit. Update `CHANGELOG.md` under `[Unreleased]` for anything a user
can notice. Do not add generated-by or co-author trailers for tools.

## Where things live

| Directory | Contents |
|---|---|
| `src/epykit/` | The package. `tl.py` orchestrates; `dmc.py`, `dmr.py`, `annotate.py`, `qc.py` are the engines. |
| `tests/` | pytest suite. `tests/fixtures/synth.py` generates the synthetic dataset with known truth. |
| `docs/` | MkDocs site. `docs/history/` is frozen planning material, excluded from the build. |
| `benchmark/` | Head-to-head benchmark scripts, frozen data, and the manuscript. |
