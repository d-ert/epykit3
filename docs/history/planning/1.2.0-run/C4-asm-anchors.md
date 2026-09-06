# C4: bisulfite-safe phasing anchors for ASM

Branch `fix-asm-anchor-classes` from `main` after L merged. Standalone PR against `main`. Decision: map ticket 24; evidence: the research note `docs/review/2026-09-06-c4-asm-research.md` on branch `research/c4-asm` (read it first) and the peer review, `docs/review/2026-06-06-epykit-peer-review.md` lines 53 to 57.

Environment: Linux or macOS with `uv sync --locked --group dev --extra all --extra bam`. Confirm `tests/test_asm.py` runs (not skips) before you start.

You own `src/epykit/asm.py` and `tests/test_asm.py`.

## Commits, in this order

1. `test(asm): a planted C/T anchor fabricates ASM today`
   - Extend `_write_synth_bam_and_vcf` so it can plant a C/T het SNV whose C-allele reads split into methylated (base C at the SNV, `Z` at the CpG) and unmethylated (base T, `z`), and T-allele reads with the same split, with an `XG` tag per read. Add the test that asserts today's code reports a significant ASM call on that data (the confound), marked `xfail(strict=True)` with the reason. Tag the existing A/G test's reads `XG:Z:CT` (A/G is safe on CT reads) so it keeps passing after the fix.
2. `fix(asm): filter phasing anchors by SNV class and read strand`
   - In `_call_asm_one_sample`, between the het check and the fetch: drop C/G anchors always; for the other classes decide per read where the base is read: with `XG:Z:CT` accept A/T, G/A, G/T; with `XG:Z:GA` accept A/T, C/T, C/A; skip the read otherwise. Reads without an `XG` tag: accept A/T anchors only. Do not use `is_reverse` as a strand proxy. No new public parameter. Log one INFO line per sample: anchors kept, anchors dropped by class, reads skipped by the strand rule. Update the `call_asm` docstring with the rule and keep the research-grade warning.
   - Turn commit 1's `xfail` into a passing test: with CT-only reads the anchor is skipped and no ASM is called; with GA reads added the call is correct.
3. `docs(changelog): C4`
   - Under `[Unreleased]` / `### Fixed`, saying what was wrong and which inputs were affected (every Bismark or MethylDackel BAM with C/T or G/A anchors).

## Contract

ASM output changes, by design, only where the old code was wrong. No other module is touched. Regen hashes unchanged (the benchmark slice does not run ASM).

## Deliver

PR title: `ASM: exclude bisulfite-confounded phasing anchors`. Body per the run README, with the strand table from the research note. Then `worker_done`.
