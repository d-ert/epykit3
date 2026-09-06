# Plan review, 2026-09-06

All 18 GitHub issues were closed when this review started. No implementation PR was open.
The review covered every issue body and comment, all 12 implementation briefs, the run document, and their referenced prototype and donor commits.
The plans now specify start conditions, file ownership, behavior, checks, and delivery.
Closed decision records remain closed. Implementation starts in the order in [README](README.md).

## Evidence and corrections

Source baseline: main commit `900ea6cfc22f521b1a708ae43c664c09c563f318`.
Original planning commit: `fe3307990f324bb6af88b9699c1a070a73a17561`.

| Plan | Finding in the original plan | Correction |
|---|---|---|
| Run | Old and new tip trees were required to match after rebasing onto changed main | Compare each layer's patch series, inspect the integrated diff, and rerun affected gates |
| Run | Shared docs and CHANGELOG ownership conflicted; R5 could edit tl.dmr while salvage was active | Keep concurrent edits by function, serialize shared documentation work, and start R5 after S3 merges |
| Run | Plans required an Orca coordinator and authorized tags and branch deletion by implication | Use the active session's handoff mechanism; leave publishing and destructive operations to separate authorization |
| L | The field guide depended on an untracked file supplied later | Archive an exact input copy with its hash on the planning branch |
| L | The plan added a redundant MkDocs exclusion and copied default-on helper text | Keep the existing history exclusion and correct helper wording for opt-in use |
| S1 | Whole donor commits edited files outside the assigned scope and referenced folded tests | Port owned changes and adapt tests to their current files |
| S1 | The donor chain_merge permutation helper materialized DMC results and did not apply the observed DMC correction | Require a streamed replay of the observed engine, chromosome universe, correction, and region filters |
| S1 | Existing max-T behavior and new count-ratio diagnostics were not distinguished in sufficient detail | Define both modes, additive columns, failed and empty permutations, and explicit nonfinite behavior |
| S2 | Filter scope for non-tile DMR methods was unspecified | Make the option tile-only and direct DMC-derived callers to upstream filtering |
| S2 | Cache invalidation did not cover stale partitions or legacy manifests | Add false-to-true-to-false cache checks and explicit legacy behavior |
| S3 | Proposed DMR smoothing metadata was actually a DMC donor change already present on main | Keep only the current permissive-preset doc correction |
| R1 | The prototype used fixed warning stacklevels and left shared helper imports unresolved | Verify caller locations and preserve DMR helpers without circular imports |
| R2 | Validation only covered the high-level config, although the CLI calls the engine directly | Validate low-level entry points too and remove unused registry fields |
| R3 | EngineInput omitted strand, the selected engine, and minimum-sample thresholds | Carry the full canonical frame and all finalization inputs; preserve GLM adjusted effects and standard errors |
| R4 | A pure move broke relative imports and the stdout guard; grouping smooth with ingestion changed help order | Fix imports and the test with the move, and preserve parser registration order |
| R5 | The issue claimed all four removals were promised for 1.2 | Keep APIs; pp.unite actually names 2.0 and CSV aliases name no version; reserve removal decisions for the maintainer |
| R6 | The plan assumed CLI handlers always used high-level APIs and omitted resume invalidation | Name the real dispatch paths, update both cache levels and metadata checks, and expose only supported flags |
| C4 | A strict xfail asserted the old bug and would unexpectedly pass | Assert desired behavior, record failure on old code, then fix it in a passing commit |
| C4 | The output contract claimed only wrong results changed | Document conservative loss of anchors and test safe signal retention |
| REL | A project-version-only bump was claimed to preserve lock validity | Permit the matching root epykit lock entry update with all dependencies pinned |

## Primary code references

All links below use the reviewed source commit.

- [Deprecation promises in pp.py](https://github.com/d-ert/epykit3/blob/900ea6cfc22f521b1a708ae43c664c09c563f318/src/epykit/pp.py) and [CSV helpers and public signatures in tl.py](https://github.com/d-ert/epykit3/blob/900ea6cfc22f521b1a708ae43c664c09c563f318/src/epykit/tl.py).
- [DMC config, recorded smoothing, and resume signature](https://github.com/d-ert/epykit3/blob/900ea6cfc22f521b1a708ae43c664c09c563f318/src/epykit/_dmc_config.py).
- [Engine inputs, effect finalization, and cache logic](https://github.com/d-ert/epykit3/blob/900ea6cfc22f521b1a708ae43c664c09c563f318/src/epykit/dmc.py).
- [CLI imports and command order](https://github.com/d-ert/epykit3/blob/900ea6cfc22f521b1a708ae43c664c09c563f318/src/epykit/cli.py) and [stdout guard](https://github.com/d-ert/epykit3/blob/900ea6cfc22f521b1a708ae43c664c09c563f318/tests/test_no_print_outside_cli.py).
- [Root package version in uv.lock](https://github.com/d-ert/epykit3/blob/900ea6cfc22f521b1a708ae43c664c09c563f318/uv.lock) and [CI matrix and uv version](https://github.com/d-ert/epykit3/blob/900ea6cfc22f521b1a708ae43c664c09c563f318/.github/workflows/test.yml).
- [Donor region FDR and chain_merge helper](https://github.com/d-ert/epykit3/blob/d5cb25448e5067ac3256abb9ce1c542a20ee68ba/src/epykit/dmr.py) and [bench-tune changes](https://github.com/d-ert/epykit3/commit/47d787823f7c2d4a509493a621d44823f38591de).
- [ASM implementation](https://github.com/d-ert/epykit3/blob/900ea6cfc22f521b1a708ae43c664c09c563f318/src/epykit/asm.py), [synthetic ASM test](https://github.com/d-ert/epykit3/blob/900ea6cfc22f521b1a708ae43c664c09c563f318/tests/test_asm.py), and [Bismark's XG and XM tag documentation](https://felixkrueger.github.io/Bismark/usage/alignment/).

## Verification performed

- Queried GitHub through both the issue list and REST API. Both reported zero open issues. Reviewed all 18 closed planning issues and their comments.
- Parsed source signatures and parser registrations. Confirmed no DMR smoothing parameter, all five public CSV-alias callers, and the exact CLI command order.
- Copied pyproject and the lockfile into an isolated temporary project. With uv 0.12.10, a version-only bump to 1.2.0 made `uv lock --check --offline` fail as stale. Running `uv lock --offline` changed only the root epykit version line.
- Ran the strict MkDocs build successfully. Planning files are excluded from the public site, so their relative links need a separate check.
- Validated 14 Markdown files, 28 relative links, nine commit references, the input hash, whitespace, and unresolved placeholder tokens before commit.
- Checked the proposed C4 confound table with SciPy: Fisher's exact p-value for 10/0 versus 10/20 is 0.0004359197907585005. This confirms that the planned null fixture can expose a false positive.

This is a planning review. Product fixes and their numerical tests remain work for the implementation layers.
The source tree, dependency lockfile, main branch, and original untracked field guide are unchanged by this review.
