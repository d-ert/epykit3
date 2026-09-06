# R4: the CLI package

Branch `refactor-4-cli-package` from `refactor-3-engine-runners`. PR against it. Decision: map ticket 16. Read `src/epykit/cli.py` (1298 lines) and `tests/test_cli_api_parity.py` first; check the `[project.scripts]` entry in `pyproject.toml`.

You own `cli.py` and the new `cli/` package, and the CLI tests.

## Commits, in this order

1. `refactor(cli): move cli.py to the cli package`
   - `git mv src/epykit/cli.py src/epykit/cli/__init__.py`. Nothing else in this commit, so the move is tracked.
2. `refactor(cli): shared helpers into _common.py`
   - `_configure_logging`, `_auto_tsv_path`, `_cli_tsv_opts`, `_write_table_local`, `_add_min_samples_args`, `_read_samplesheet_groups`, `_cli_n1_and_footgun_checks`.
3. `refactor(cli): one module per command family`
   - `_ingest.py` (convert, filter, summary, smooth), `_dmc.py`, `_dmr.py`, `_downstream.py` (annotate, qc-report, report, aggregate-regions, export with its sub-parsers). Each exposes `register(sub: argparse._SubParsersAction) -> None` that adds its parsers and sets `func=`. `__init__.py` keeps `build_parser` (calling the four registrars in today's order) and `main`, so `epykit.cli:main` and `python -m epykit` are unchanged.

## Contract

`epykit --help` and every subcommand's `--help` are byte-identical before and after (capture both and diff them in the PR body). Option names, defaults and outputs unchanged; the parity test is the gate. No `lr+` flags. Regen hashes unchanged.

## Deliver

PR title: `Split the CLI into a package, one module per command family`. Then `worker_done`.
