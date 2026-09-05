## What and why

<!-- One paragraph. Link the issue if there is one. -->

## Does this change engine output?

- [ ] No. `benchmark/scripts/regen_small.py` passes unchanged.
- [ ] Yes. The hash file is re-snapshotted in this PR, and `CHANGELOG.md` explains the change to results.

## Checklist

- [ ] `uv run --frozen pytest -m "not slow" --strict-markers -ra` passes
- [ ] `uv run --frozen pytest -m slow --strict-markers -ra` passes
- [ ] `uv run --frozen ruff check src/` and `uv run --frozen mypy src/epykit` pass
- [ ] New behaviour has a test; a fixed bug has a test that failed before the fix
- [ ] `CHANGELOG.md` updated under `[Unreleased]` if a user can notice the change
- [ ] No `print()` in library code; new preprocessing steps append to `_store_history`
