# What the C4 ASM fix needs

Research note for d-ert/epykit3 issue #23. Repository state examined: `main` at
`900ea6cfc22f521b1a708ae43c664c09c563f318` (2026-09-06). Read-only; nothing was
installed or edited. Each section separates what was verified (read in a file
or on a primary web page) from what is inferred (my reasoning from those facts).

## 1. Which SNV classes are bisulfite-safe on each strand, and why

### The confound in plain terms

Bisulfite converts every unmethylated cytosine to uracil, which PCR reads as
thymine; a methylated cytosine stays cytosine. Bis-SNP, Background: "Bisulfite
treatment of DNA converts unmethylated cytosines to uracils, which are replaced
by thymines during amplification." Bismark's docs index: "unmethylated C reads
as T, methylated C stays C."

The two strands convert separately. A read from the original top strand
(Bismark "OT") shows an unmethylated top-strand C as T. A read from the original
bottom strand ("OB") has its own C converted, which in reference coordinates is
a G read as A. Bismark `--help`: reads are "transformed into a bisulfite
converted forward strand version (C->T conversion) or into a bisulfite treated
reverse strand (G->A conversion of the forward strand)". So at a heterozygous
SNV, a T on a top-strand read is either a true T allele or an unmethylated C
allele, and an A on a bottom-strand read is either a true A allele or a G allele
whose bottom-strand C was unmethylated. Code that assigns reads to alleles by
comparing the raw base to REF and ALT makes allele membership a function of
methylation state. That is the fabricated ASM.

BS-SNPer, Introduction: "true C/T SNPs in samples cannot be distinguished from
C/T substitutions caused by bisulfite conversion"; "a C/T SNP in Watson strand
and G/A SNP in Crick strand should be interpreted as unmethylated cytosine
variation. Thus, Watson and Crick strands should be independently treated".
BISCUIT (Zhou et al. 2024) names "the ambiguity of whether C to T or G to A
conversions are due to the presence of a SNP or an unmethylated cytosine" and
keeps ambiguous calls as IUPAC Y (C or T) and R (G or A).

### Which strand rescues which class

Bis-SNP, Background: in the directional protocol "approximately half the reads
at a given cytosine position (those mapping to the 'C-strand') can be used for
methylation quantification but cannot distinguish C>T SNPs. The other half
(those mapping to the 'G-strand') ... yield no methylation information but can
be used to identify C>T SNPs." And: "When these C>T SNPs are heterozygous, they
can be used in the analysis of allele specific methylation." Its genotyper
treats "Reads on the G-strand opposite the cytosine" with the normal GATK model
and C-strand reads with "an alternate model that considers C>T substitutions as
either potential errors or bisulfite conversions."

In Bismark BAMs the per-read strand is the `XG:Z` tag (genome conversion
state). Verified in the Bismark source: `$genome_conversion` is `'CT'` or
`'GA'`, written as `XG:Z:$genome_conversion`; `XR:Z` is the read conversion
state; (read, genome) pairs map to strands CT/CT OT, CT/GA OB, GA/CT CTOT,
GA/GA CTOB. Directional libraries yield only OT and OB (help text: directional
mode "will only run 2 alignment threads to the original top (OT) or bottom (OB)
strands").

Inferred from those rules (the table is my derivation): `XG:Z:CT` reads (OT,
CTOT) may show a top-strand C as T, while C, G, A are literal; `XG:Z:GA` reads
(OB, CTOB) may show a top-strand G as A, while C, G, T are literal. Per
unordered class {REF, ALT}, with plain base matching as asm.py does it:

| Class | XG:CT reads | XG:GA reads | Failure mode if used anyway |
| --- | --- | --- | --- |
| A/T | safe | safe | none |
| C/T | not usable | safe | CT reads: unmethylated C allele counted as T allele (fabricated ASM) |
| G/A | safe | not usable | GA reads: unmethylated G allele counted as A allele (fabricated ASM) |
| C/A | not usable | safe | CT reads: unmethylated C-allele reads match neither base, dropped |
| G/T | safe | not usable | GA reads: unmethylated G-allele reads dropped |
| C/G | not usable | not usable | each strand drops one allele when unmethylated |

With literal matching and no strand information only A/T is safe, as the review
says. Strand-aware exclusion (skip a read when its strand cannot resolve the
class) makes every class but C/G usable from one strand. Conversion-aware
folding goes further: on a CT read an observed T at a C/G or C/A site can only
be the C allele (T is not an allele), so fold T to C; on a GA read fold A to G
at G/C and G/T sites. C/T on CT reads and G/A on GA reads can never be folded,
because the ambiguous base is itself an allele. Folding is the idea behind
BISCUIT's Y/R alphabet and Bis-SNP's alternate C-strand model. dnmtools
(MethPipe's successor) sidesteps SNVs: `amrfinder` fits a "single-allele model"
against a "two-allele model" on epiread patterns; a different design, not a fix.

`read.is_reverse` is not a substitute for `XG` (inferred from Bismark source
comments): in paired-end directional data read 2 of an OT pair aligns to the
minus strand yet carries `XG:Z:CT` ("CT read1/GA read2/CT genome, original top
strand"). `is_reverse` matches the strand only for single-end directional data.

## 2. What the peer review specified as the fix

Verified in `docs/review/2026-06-06-epykit-peer-review.md`:

- Line 53, heading C4: "ASM phasing reads the bisulfite-converted base at the
  SNV; C/T & G/A SNVs confound allele with methylation".
- Line 54, location: "`asm.py:167-176` (`base = read.query_sequence[qry]; if
  base==ref → hap1 elif base==alt → hap2 else drop`); `query_sequence` holds
  converted bases in a Bismark BAM (`bam_io.py:204`)". The asm.py lines have
  drifted to 174 to 183 (format commits b6757c2, c3e6d9b); `bam_io.py:204` is
  still `allele = seq[query_idx] if seq else ""`.
- Line 55, mechanism: "an unmethylated ref-C is miscalled as the alt allele, so
  haplotype membership becomes a direct function of methylation state"; for
  other C/G-containing SNVs "the unmethylated C reads as T, matches neither
  allele, and the read is dropped"; "Only A/T SNVs are safe, and there is no
  SNV-type filtering, strand-aware exclusion, or methylation-aware genotyping
  anywhere."
- Line 56, fix: "Restrict anchors to bisulfite-safe SNV classes (exclude C/T,
  A/G and strand-context-confounded types), or use a methylation-aware genotyper
  (Bis-SNP / biscuit). At minimum document the restriction and warn."
- Line 57: "certain (agent) on the mechanism". Line 158 lists C4 among the
  research-grade surface "shipped as stable 1.0 with real scientific bugs";
  line 114 (M-PKG5) proposes moving `tl.asm` behind `epykit.experimental.*`.

Verified in `docs/review/2026-06-06-remediation-summary.md`:

- Line 68, under "Cannot be validated in this environment": C4 is "the
  highest-severity deferred item. The fix (restrict phasing anchors to
  bisulfite-safe SNV classes, or gate behind `experimental`) is well-specified
  in the review, but `asm.py` is pysam-gated and its tests
  `importorskip("pysam")` on Windows, so it cannot be exercised here. Do this on
  Linux with the `bam` extra installed."
- Line 83: "On a Linux box with the `bam` extra, do C4 (ASM)"; "it's the
  highest-severity deferred item and silently fabricates ASM on WGBS."
  Line 65 restates M-PKG5 (experimental gating) as unfinished.

Inferred: the review allows three tiers. Minimal: A/T anchors only, plus a
warning and a docstring note. Middle: strand-aware exclusion using `XG`. Full:
methylation-aware genotyping or folding. The remediation summary names only
the minimal tier and the experimental gate.

## 3. What asm.py does today at the anchor-selection step

Verified in `src/epykit/asm.py` (336 lines). Entry points: `call_asm` (line
58) and the `MethylData` wrapper `asm` (line 296), re-exported as `tl.asm` in
`src/epykit/tl.py` lines 2210 to 2246. The docstring (lines 3 to 8, 23 to 24)
names heterozygous biallelic SNVs as anchors and says nothing about bisulfite
classes. Flow in `_call_asm_one_sample` (line 125):

- Line 140 opens the VCF and BAM. Lines 141 to 147 iterate VCF records, skip
  non-biallelic (142 to 143) and indels (146 to 147: `if len(ref) != 1 or
  len(alt) != 1`). Lines 150 to 151 keep heterozygous records via `_is_het`
  (283 to 293).
- Lines 153 to 158 convert to 0-based `pos` and `bam.fetch(chrom, pos, pos +
  1)`. Lines 159 to 163 drop unmapped, low-MAPQ, secondary, supplementary reads.
  Lines 165 to 176 find the query index via
  `get_aligned_pairs(matches_only=True)` and take `read.query_sequence`.
- Lines 177 to 183 assign the allele: `base = seq[qry].upper()`; `base == ref`
  is haplotype 1, `base == alt` is haplotype 2, else drop. No strand or class
  check exists.
- Lines 184 to 193 keep one haplotype per `query_name` and drop reads whose
  SNVs disagree; lines 199 to 200 apply `min_phased_snvs`; lines 207 to 215
  join onto `meth_df` by `read_id`; lines 220 to 253 build the per-CpG 2x2
  table; line 266 runs `fisher_exact_vectorized`.

`meth_df` comes from `bam_io.read_methylation_calls` (line 83). That frame has
a `strand` column from `read.is_reverse` (bam_io.py line 155) and `allele_base`
(line 204), but asm.py uses neither; it re-reads bases from the BAM. Nothing in
`src/epykit` reads the Bismark `XG` or `XR` tags (grep verified).

Where a fix goes (inferred):

- Class filter: between lines 147 and 150, on `{ref.upper(), alt.upper()}`.
  Minimal fix: `continue` unless the set is `{"A", "T"}`; log skipped anchors.
- Strand-aware exclusion: between lines 176 and 177, once `read` is in hand.
  Read `read.get_tag("XG")` and skip the read when the class is not resolvable
  on that strand per the section 1 table. Folding, if wanted, replaces the
  comparisons at lines 178 to 181.
- The `methyldackel` path (MM/ML BAMs, usually bwa-meth) has no `XG` tag;
  strand needs another source there, or the minimal A/T rule.

Gating: `call_asm` calls `_require_pysam()` (asm.py line 75; bam_io.py lines
54 to 65), which raises `ImportError` saying "Install with: pip install
'epykit[bam]' (Linux/macOS only ...)". `pyproject.toml` lines 95 to 98 define
`bam = ["pysam>=0.22"]`; lines 91 to 94 the twin `methylkit` extra. The `all`
extra (lines 123 to 134) does not include pysam.

## 4. What the tests need to run

Verified:

- `tests/test_asm.py` is the only ASM test file. Line 17:
  `pysam = pytest.importorskip("pysam", ...)`. One test,
  `test_asm_recovers_planted_signal` (line 78).
- The fixture is built at run time, not checked in. `_write_synth_bam_and_vcf`
  (lines 22 to 75) writes a 20-read BAM into `tmp_path` with pysam: contig
  `chr_asm`, reads at 0-based start 45, `flag = 0`, CIGAR `30M`, a het A/G SNV at
  0-based 50 (VCF line 69: `chr_asm 51 . A G ... GT 0/1`), a CpG at 60 encoded
  in `XM` as `Z` (reads 0 to 9, A allele) or `z` (reads 10 to 19, G allele).
  Reads carry no `XG` or `XR` tag. BAM sorted and indexed; VCF bgzipped and
  tabix-indexed.
- No BAM, BAI, CRAM, SAM or FASTA exists under `tests/` (find verified).
  `tests/fixtures/synth.py` makes Bismark `.cov.gz` files only; `conftest.py`
  has no pysam fixture. `tests/test_bam_io.py` line 20 and
  `tests/test_entropy.py` line 19 skip the same way; `test_bam_io.py` builds
  its BAM in-test too.
- CI (`.github/workflows/test.yml` line 52) installs `--group dev --extra all`,
  never `--extra bam`, so all three files skip on every CI leg. The remediation
  summary's "5 platform-skips" (line 76) agrees. The local `.venv` on this Mac
  (Python 3.12) has no pysam, so `test_asm.py` skips here as well.

Inferred:

- No committed BAM fixture is needed. A C4 regression test extends the in-test
  builder: plant a C/T het SNV with C-allele reads split between methylated
  (base C at the SNV, `Z` at the CpG) and unmethylated (base T, `z`), and
  T-allele reads with the same split. Today's code reports a strong fake ASM;
  fixed code must skip the anchor (no `XG`, or all reads `XG:Z:CT`) or, with
  `XG:Z:GA` reads added, report no ASM.
- The existing A/G anchor is safe on CT reads under a strand-aware rule but
  fails an A/T-only rule, so the test changes with whichever rule lands. If the
  rule reads `XG`, the fixture needs `r.set_tag("XG", "CT")` or a defined
  fallback for untagged reads.
- Run: `uv sync --locked --group dev --extra all --extra bam`, then
  `uv run pytest tests/test_asm.py tests/test_bam_io.py tests/test_entropy.py`.
  Adding `--extra bam` to the CI install line would end the silent skips.

## 5. Whether pysam installs cleanly on Arch Linux for Python 3.10 to 3.13

Verified in `uv.lock`: the block at lines 3001 to 3037 pins `pysam 0.24.0`
(uploaded 2026-04-27) and lists, for each of cp310 through cp314,
`manylinux_2_27_x86_64.manylinux_2_28_x86_64`,
`manylinux_2_27_aarch64.manylinux_2_28_aarch64`, `musllinux_1_2_{x86_64,aarch64}`,
`macosx_*_x86_64` and `macosx_11_0_arm64` wheels. cp38/cp39 wheels are absent
because of the lock's `requires-python = ">=3.10"` (line 3).

Verified on PyPI (`https://pypi.org/pypi/pysam/0.24.0/json`): 43 files,
`requires_python >=3.8`, the same manylinux x86_64 and aarch64 wheels for cp310
through cp314; 0.24.0 is the latest release. No Windows wheel exists, which
matches the `bam_io.py` error text and the README platform note. Arch package
page: `core/glibc` is `2.44+r24+g16be1518495f-1` (updated 2026-08-11);
manylinux_2_28 needs glibc 2.28 or newer. The pysam install docs say `pip
install pysam` takes a "pre-built wheel" and compiles htslib only under
`--no-binary`.

Inferred, not executed: on an x86_64 Arch box `uv sync --locked --extra bam`
(CI pins uv 0.12.10, workflow line 14) resolves from the existing lock without
rewriting it and installs the manylinux wheel for 3.10, 3.11, 3.12 and 3.13
with no compiler or htslib. uv-managed CPython builds are manylinux-compatible,
so the Arch system Python version does not matter. macOS arm64 wheels exist for
every version too, so the "Linux only" phrasing in the remediation summary
reflects the Windows box, not a Mac limit; this Mac could run the tests after
`--extra bam`.

## Sources

Repository files (relative to the repository root):

- `docs/review/2026-06-06-epykit-peer-review.md` lines 53 to 57, 114, 158.
- `docs/review/2026-06-06-remediation-summary.md` lines 65, 68, 76, 83.
- `src/epykit/asm.py` lines 3 to 24, 58, 75, 125 to 200, 207 to 266, 283 to 333.
- `src/epykit/bam_io.py` lines 54 to 65, 155, 181, 204; `src/epykit/tl.py`
  lines 2210 to 2246.
- `tests/test_asm.py` lines 17, 22 to 75, 78 to 100; `tests/test_bam_io.py`
  lines 20 to 71; `tests/test_entropy.py` line 19.
- `pyproject.toml` lines 11, 31 to 34, 91 to 98, 123 to 134; `uv.lock` lines
  3, 857, 869, 932 to 933, 3001 to 3037.
- `.github/workflows/test.yml` lines 14, 24 to 35, 52; `CHANGELOG.md` lines
  998 to 1007; `.python-version`; `.venv/lib/python3.12/site-packages`.

Web sources:

- Bis-SNP, Liu et al. 2012, Genome Biology 13:R61:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3491382/ (also
  https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3491382/fullTextXML).
- BS-SNPer, Gao et al. 2015, Bioinformatics 31:4006:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4673977/ ; repo
  https://github.com/hellbelly/BS-Snper (README has usage only).
- BISCUIT, Zhou et al. 2024, NAR 52:e32:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11014253/ ; docs
  https://huishenlab.github.io/biscuit/docs/pileup.html and
  https://huishenlab.github.io/biscuit/docs/epiread.html
- dnmtools amrfinder: https://dnmtools.readthedocs.io/en/latest/amrfinder/ and
  https://dnmtools.readthedocs.io/en/latest/allelic/
- Bismark docs https://felixkrueger.github.io/Bismark/ and
  https://felixkrueger.github.io/Bismark/usage/alignment/ ; `bismark --help`
  mirror https://hpc.nih.gov/apps/bismark/bismark_help.txt ; source
  https://raw.githubusercontent.com/FelixKrueger/Bismark/master/bismark (lines
  4405 to 4444, 7113 to 7120, 8611, 8615, 8746 to 8760, 9605 to 9609).
- PyPI JSON https://pypi.org/pypi/pysam/0.24.0/json and
  https://pypi.org/pypi/pysam/json ; pysam install docs
  https://pysam.readthedocs.io/en/latest/installation.html
- Arch glibc package https://archlinux.org/packages/core/x86_64/glibc/
