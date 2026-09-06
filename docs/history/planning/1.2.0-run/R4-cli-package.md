# R4: split the CLI without changing its interface

Start from `refactor-3-engine-runners` after its implementation and gates pass.
Use `refactor-4-cli-package` and target R3.
Read the [run rules](README.md), issue [16](https://github.com/d-ert/epykit3/issues/16), and all imports and parser registrations in `src/epykit/cli.py`.

Own the CLI module and new package, CLI tests, and `tests/test_no_print_outside_cli.py`.

## Implement

1. Capture top-level and all nested command help on the base with a fixed terminal width and program name.
2. Move `cli.py` to `cli/__init__.py`. In that same commit, change imports of epykit siblings from one dot to two dots, including lazy imports. A bare file move is not a working commit.
3. In the move commit, update the stdout rule to allow only files under `src/epykit/cli/` by relative path. Update the CLI-presence assertion. Do not allow every `__init__.py` or a broad filename pattern.
4. Move shared helpers into `cli/_common.py`: logging, TSV helpers, minimum-sample arguments, sample-group reading, and the n=1 checks. Keep imports directed from command modules to common helpers.
5. Use four registrars called in this order: `_ingest.register`, `_dmc.register`, `_dmr.register`, `_downstream.register`.
6. Put convert, filter, and summary in `_ingest.py`. Put dmc and dmr in their respective modules. Put annotate, qc-report, smooth, report, aggregate-regions, and export in `_downstream.py`, in that order. Moving smooth into the ingest registrar would change top-level help order.
7. Keep `build_parser` and `main` in `cli/__init__.py`. Preserve the `epykit.cli:main` entry point and `python -m epykit`. Avoid an import from a command module back into the package initializer.
8. Update tests that patch private CLI helpers to patch their new defining module. Preserve behavior assertions.

## Accept when

- The top-level command order remains convert, filter, summary, dmc, dmr, annotate, qc-report, smooth, report, aggregate-regions, export.
- Before and after help output matches byte for byte under the same environment, including export subcommands.
- Existing CLI integration, API parity, TSV export, n=1, and error-path tests pass.
- The stdout test still rejects print calls outside the CLI package.
- Both entry points work from an installed package.
- No option, default, output format, or numerical result changes. New options wait for R6.
- All code-layer gates pass.

PR title: `Split the CLI into command modules without changing its interface`.
