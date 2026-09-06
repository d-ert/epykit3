# R1: the tl.dmc stage split

Branch `refactor-1-stage-split` from `main` after L merged. Bottom of the refactor stack; PR against `main`. Decisions: map ticket 13 (the split) and ticket 26 (the CI install change and the smoothed-store test). Read `README.md` here, `CLAUDE.md`, `src/epykit/_dmc_config.py`, and the stub at `src/epykit/_dmc_stages.py` on branch `prototype/dmc-stage-split` (commit `3836a6d`): it carries every stage signature and the `tl.py` line ranges each absorbs.

You own `tl.dmc` and `_run_dmc_contrast` in `tl.py`, `_dmc_config.py`, the new `_dmc_stages.py`, `.github/workflows/test.yml`, and the tests you add. You do not touch `tl.dmr`, `dmr.py`, `dmc.py`, `cli.py`.

## Commits, in this order

1. `ci: install the bam extra on the ubuntu legs`
   - In `test.yml`, give each matrix include entry an `extras` field: `"--extra bam"` for ubuntu, `""` for windows (pysam has no Windows wheel). The install line becomes `uv sync --locked --python ... --group dev --extra all ${{ matrix.extras }}`. The `slow` job (ubuntu) adds `--extra bam` too. Validate the YAML with PyYAML; this PR's own checks must show `test_asm.py`, `test_bam_io.py` and `test_entropy.py` running, not skipping, on ubuntu.
2. `refactor(tl): split tl.dmc into stages`
   - Create `src/epykit/_dmc_stages.py` from the stub: keep its records (`TsvPlan`, `DMCPlan`, `ContrastDesign`, `ResumeTicket`, `DMCOutcome`) and its nine stages, and move the bodies from `tl.py` into them line range by line range. `tl.dmc` keeps its signature, docstring, the `DMCConfig` build, and becomes the orchestrator body of `dmc_sketch`. `_run_dmc_contrast` becomes `run_contrast`. `publish` is the only writer of `md.uns["dmc"]`, always through `cfg.to_uns`.
   - Keep today's error order inside `plan_run`: TSV resolution, then `cfg.validate()`, then the contrast dispatch. No resume support on the contrast path. `persist_resume` after `publish`. Warnings emitted from inside a stage move from `stacklevel=2` to `3`.
   - The one observable difference, decided on the ticket: the `log2_odds_ratio` FutureWarning is emitted by `finish` on the resume cache hit too. Say so in the PR body. If a test asserts the hit is silent, update that test only.
3. `test(smoothed_store): the pseudo-count store through open_input_store`
   - One focused test: build a small smoothed store with `pp.smooth`, enter `open_input_store` with `use_smoothed=True`, and assert on the yielded store that `coverage` equals the raw coverage and `N_meth` equals `round(beta_smooth * coverage)` per site; assert the temporary directory is gone after exit. Fast tier.

## Contract

No results change. `tests/test_dmc_metadata.py` (the uns contract) and `tests/test_cli_api_parity.py` must pass untouched. Regen hashes unchanged. Report `tl.dmc`'s complexity before and after (`ruff check --select C901 --config 'lint.mccabe.max-complexity=1'`).

## Deliver

PR title: `Split tl.dmc into stages with one writer for the uns record`. Body per the run README, naming the resume-hit warning as the one difference. Then `worker_done`.
