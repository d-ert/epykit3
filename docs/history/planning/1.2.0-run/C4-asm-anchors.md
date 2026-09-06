# C4: reject bisulfite-confounded ASM anchors

Start from main after R1 merges, so this PR's CI runs the BAM tests.
Use `fix-asm-anchor-classes` and target main.
Read the [run rules](README.md), issues [23](https://github.com/d-ert/epykit3/issues/23) and [24](https://github.com/d-ert/epykit3/issues/24), and the research note on `research/c4-asm`.
The note is [docs/review/2026-09-06-c4-asm-research.md](https://github.com/d-ert/epykit3/blob/6825b1b2583793c0ab59e0e0536b3ef69eb572b7/docs/review/2026-09-06-c4-asm-research.md). Treat its findings as context, not approval to expand this fix.

Own `src/epykit/asm.py`, `tests/test_asm.py`, `docs/analysis/asm.md`, and the C4 changelog entry.
Use Linux or macOS with `uv sync --locked --group dev --extra all --extra bam`.
Run the current ASM test first and confirm it executes.

## Implement the conservative rule

Bismark records the genome conversion state in XG, distinct from the ordinary alignment flag.
See the [Bismark alignment documentation](https://felixkrueger.github.io/Bismark/usage/alignment/).
Use that tag on each read, not `is_reverse`.

For each heterozygous biallelic SNV, normalize REF and ALT to uppercase and compare their unordered class.

| Unordered SNV class | XG=CT | XG=GA | Missing or unrecognized XG |
|---|---|---|---|
| A/T | Accept | Accept | Accept |
| A/G | Accept | Reject | Reject |
| G/T | Accept | Reject | Reject |
| C/T | Reject | Accept | Reject |
| A/C | Reject | Accept | Reject |
| C/G | Reject | Reject | Reject |

The accepted classes exclude a potentially converted allele from the selected conversion strand.
This deliberately conservative policy can drop usable anchors. Do not claim that all changes are limited to false positives.

1. Apply the class check after the heterozygous SNV check. Skip invalid bases and C/G anchors before fetching reads.
2. Apply the XG rule before assigning a read to an allele. Preserve the existing mapping, sequence, and phasing checks.
3. Keep the public API unchanged. Add no fallback switch that restores unsafe phasing.
4. Log one summary per sample, including no-phaseable-read results. Define counts as anchors that assigned at least one accepted read, anchors rejected before fetch by class, and read-anchor observations rejected by the XG rule.
5. Document the conservative untagged-read fallback and research limitations. The current `call_asm` docstring has no existing research warning to preserve; add plain scope text without inventing a new runtime warning.

## Prove the confound and the retained signal

Extend the synthetic BAM and VCF builder rather than adding a large binary fixture.

Create a null C/T anchor with 20 reads per true allele.
For the C allele, 10 reads have C at the anchor and Z at the measured CpG. The other 10 have converted T at the anchor and z at the CpG.
For the T allele, all 20 reads have T at the anchor. Ten have Z at the CpG and ten have z.
True methylation is 50% for both alleles. CT-only reads make the old phaser split the data into 10/0 versus 10/20 counts.

Write the desired assertion: the unsafe CT-only anchor contributes no ASM result.
Run it against the old code and record the failure, then implement the fix in the same passing commit.
Do not mark an assertion of the old false positive as strict xfail; it would unexpectedly pass.

Also verify:

- GA-tagged C/T reads with balanced methylation produce the expected balanced counts and no significant difference.
- A genuine GA C/T signal remains detectable.
- The existing A/G signal remains detectable after its reads receive XG=CT.
- A parameterized test covers the six classes, both REF/ALT orders, both XG values, and the missing or invalid tag fallback.
- Untagged A/T signal remains usable. C/G and unsafe reads never contribute to phasing.

## Accept when

The focused tests and all code-layer gates pass with BAM support installed.
Describe changed ASM results and possible loss of anchors in the changelog.
Other numerical results remain unchanged. LR hash parity is necessary but does not test ASM.

PR title: `Exclude bisulfite-confounded ASM phasing anchors`.
