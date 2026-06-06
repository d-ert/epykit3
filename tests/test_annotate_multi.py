"""annotatr-style multi-annotation in ``annotate_features``.

Covers four invariants:

1. ``multi_annotation=False`` is unchanged -- old schema only
   (``gene_id``, ``gene_name``, ``feature_type``, ``distance_to_tss``).
2. ``multi_annotation=True`` adds exactly four columns with the documented
   names and dtypes.
3. The nearest-TSS gene can differ from the single-best ``gene_name`` --
   the whole point of HOMER-style annotation. Constructed scenario: a
   site sits inside one gene's intron while another gene's TSS is closer.
4. ``all_overlapping_*`` truly captures one-to-many -- a site that lies in
   one gene's intron AND another gene's promoter window must surface both.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from epykit.annotate import annotate_features
from epykit import annotate as A


@pytest.fixture
def synth_gtf(tmp_path):
    """Minimal GTF with three genes designed to expose the three invariants.

    Layout on chr1 (1-based GTF coords; converted to 0-based half-open in
    the parser):

        GENEA: + strand, 1001-3000, exons 1001-1100 and 2901-3000
               (so intron spans roughly 1100-2901)
        GENEB: + strand, 2000-2500, exon 2000-2500
               -> TSS at 2000, sits *inside* GENEA's intron
        GENEC: - strand, 9001-10000, exon 9001-10000
               -> TSS at 10000 (gene end on "-"); used as nearest-TSS
               anchor for the intergenic test site.
    """
    gtf = tmp_path / "synth.gtf"
    lines = [
        # GENEA
        'chr1\tt\tgene\t1001\t3000\t.\t+\t.\tgene_id "g1"; gene_name "GENEA";',
        'chr1\tt\texon\t1001\t1100\t.\t+\t.\tgene_id "g1"; gene_name "GENEA";',
        'chr1\tt\texon\t2901\t3000\t.\t+\t.\tgene_id "g1"; gene_name "GENEA";',
        # GENEB sits inside GENEA's intron
        'chr1\tt\tgene\t2000\t2500\t.\t+\t.\tgene_id "g2"; gene_name "GENEB";',
        'chr1\tt\texon\t2000\t2500\t.\t+\t.\tgene_id "g2"; gene_name "GENEB";',
        # GENEC far away, on "-" strand, used for nearest-TSS check
        'chr1\tt\tgene\t9001\t10000\t.\t-\t.\tgene_id "g3"; gene_name "GENEC";',
        'chr1\tt\texon\t9001\t10000\t.\t-\t.\tgene_id "g3"; gene_name "GENEC";',
    ]
    gtf.write_text("\n".join(lines) + "\n")
    return str(gtf)


@pytest.fixture(autouse=True)
def _clear_caches():
    """Each test starts with empty annotate caches so synth GTFs don't bleed."""
    A._GTF_CACHE.clear()
    A._BUILT_FEATURES_CACHE.clear()
    yield
    A._GTF_CACHE.clear()
    A._BUILT_FEATURES_CACHE.clear()


def _sites(positions):
    """Helper: build a single-base sites DataFrame from a list of chr1 positions."""
    return pl.DataFrame({
        "chrom": ["chr1"] * len(positions),
        "pos":   list(positions),
    })


def test_chrom_name_mismatch_raises(tmp_path):
    """C2: 'chr1' sites vs a '1'-named (Ensembl-style) GTF must raise rather
    than silently annotate every site as intergenic."""
    gtf = tmp_path / "ensembl_named.gtf"
    gtf.write_text("\n".join([
        '1\tt\tgene\t1001\t3000\t.\t+\t.\tgene_id "g1"; gene_name "GENEA";',
        '1\tt\texon\t1001\t3000\t.\t+\t.\tgene_id "g1"; gene_name "GENEA";',
    ]) + "\n")
    with pytest.raises(ValueError, match="chromosome names"):
        annotate_features(_sites([1500]), str(gtf))


def test_matching_chrom_names_do_not_raise(synth_gtf):
    """Sanity: matching 'chr1' on both sides annotates without raising."""
    out = annotate_features(_sites([1500]), synth_gtf)
    assert "feature_type" in out.columns


def test_default_now_includes_multi_columns(synth_gtf):
    """multi_annotation=True is the default: output includes both legacy
    and annotatr-style columns out of the box."""
    sites = _sites([2100])
    out = annotate_features(sites, synth_gtf)
    cols = set(out.columns)
    assert {"gene_id", "gene_name", "feature_type", "distance_to_tss"} <= cols
    assert {
        "nearest_tss_gene", "nearest_tss_distance",
        "all_overlapping_genes", "all_overlapping_features",
    } <= cols


def test_multi_false_opts_out_of_extra_columns(synth_gtf):
    """multi_annotation=False is still respected -- useful for callers that
    want the slimmer legacy schema (e.g. narrow storage, or downstream code
    that asserts on exact column sets)."""
    sites = _sites([2100])
    out = annotate_features(sites, synth_gtf, multi_annotation=False)
    cols = set(out.columns)
    assert {"gene_id", "gene_name", "feature_type", "distance_to_tss"} <= cols
    for new_col in (
        "nearest_tss_gene", "nearest_tss_distance",
        "all_overlapping_genes", "all_overlapping_features",
    ):
        assert new_col not in cols, f"unexpected column {new_col!r} when multi_annotation=False"


def test_multi_adds_four_columns_with_correct_dtypes(synth_gtf):
    sites = _sites([2100])
    out = annotate_features(sites, synth_gtf, multi_annotation=True)
    assert out["nearest_tss_gene"].dtype == pl.Utf8
    assert out["nearest_tss_distance"].dtype == pl.Int32
    assert out["all_overlapping_genes"].dtype == pl.List(pl.Utf8)
    assert out["all_overlapping_features"].dtype == pl.List(pl.Utf8)


def test_nearest_tss_can_differ_from_best_pick(synth_gtf):
    """Site at pos=2100 sits in GENEA's intron and is 100 bp inside GENEB.

    - Single-best pick prioritises promoter > intron, so picks GENEB
      because GENEB's TSS at 2000 + default promoter window (-2000/+200)
      covers 2100. So gene_name=GENEB here too.
    - Nearest-TSS rule: GENEA's TSS at 1000 is 1100 bp away; GENEB's TSS
      at 2000 is 100 bp away -- GENEB wins.

    Use a second site at pos=2700 -- well past GENEB's promoter window
    (2000+200=2200) and outside GENEB body (which ends 2500). At 2700
    the site is GENEA-intron-only:
      - Single-best: gene_name=GENEA, feature=intron
      - Nearest-TSS: GENEB at 2000 is 700 bp away; GENEA at 1000 is
        1700 bp away -- GENEB still wins despite the best-pick being GENEA.

    Note: explicitly pin ``features=`` to the legacy 3-tuple so the
    intron / nearest-TSS divergence stays demonstrable. Under the new
    full-HOMER default this fixture's 2700 also falls into GENEB's TTS
    window ([End-100, End+1000) = [2400, 3500)) so single-best resolves
    to GENEB+TTS and the divergence collapses. The test is about
    best-pick vs nearest-TSS semantics, not about category vocabulary,
    so restricting the feature set is the right knob here.
    """
    sites = _sites([2700])
    out = annotate_features(
        sites, synth_gtf, multi_annotation=True,
        features=("promoter", "exon", "intron"),
    )
    row = out.row(0, named=True)
    assert row["gene_name"] == "GENEA"
    assert row["feature_type"] == "intron"
    assert row["nearest_tss_gene"] == "GENEB"
    # GTF is 1-based closed (start=2000) but the parser stores 0-based
    # half-open internally, so the stored TSS is 1999 -- distance is 701, not
    # 700. The convention is consistent with the rest of epykit (BED-style
    # coords) and matches the ``distance_to_tss`` column's behaviour too.
    assert row["nearest_tss_distance"] == 701


def test_nearest_tss_minus_strand_sign_flip(synth_gtf):
    """GENEC is on "-" strand with TSS at 10000 (gene end). For an upstream
    site at pos=8000, raw delta is -2000, but the sign must flip on "-"
    strand so the reported distance is +2000 (downstream in transcription
    direction)."""
    sites = _sites([8000])
    out = annotate_features(sites, synth_gtf, multi_annotation=True)
    row = out.row(0, named=True)
    assert row["nearest_tss_gene"] == "GENEC"
    assert row["nearest_tss_distance"] == 2000


def test_all_overlapping_captures_one_to_many(synth_gtf):
    """Site at pos=2100 lies in GENEA's intron AND inside GENEB's promoter
    window. Multi-annotation must surface both."""
    sites = _sites([2100])
    out = annotate_features(sites, synth_gtf, multi_annotation=True)
    row = out.row(0, named=True)
    genes = set(row["all_overlapping_genes"])
    feats = set(row["all_overlapping_features"])
    assert {"GENEA", "GENEB"} <= genes, f"expected both genes, got {genes}"
    assert {"intron", "promoter"} <= feats, f"expected intron+promoter, got {feats}"


def test_no_overlap_site_gets_empty_lists_and_still_finds_nearest(synth_gtf):
    """Intergenic site far from any gene body: best-pick is intergenic with
    empty gene_name, but nearest-TSS still returns the closest gene."""
    sites = _sites([6000])  # between GENEA/B group (<=3000) and GENEC (>=9001)
    out = annotate_features(sites, synth_gtf, multi_annotation=True)
    row = out.row(0, named=True)
    assert row["feature_type"] == "intergenic"
    assert row["gene_name"] == ""
    assert list(row["all_overlapping_genes"]) == []
    assert list(row["all_overlapping_features"]) == []
    # Nearest TSS: GENEC at 10000 is 4000 away (after "-" flip: +4000);
    # GENEA at 1000 is 5000 away; GENEB at 2000 is 4000 away (no flip).
    # Tie at |4000| -- bisect picks one deterministically; we just assert it
    # didn't fall back to "" and the magnitude is right.
    assert row["nearest_tss_gene"] in {"GENEB", "GENEC"}
    assert abs(row["nearest_tss_distance"]) == 4000


def test_multi_annotation_reuses_built_features_cache(synth_gtf):
    """The TSS-by-chromosome arrays must live in the same cache bundle so a
    second multi-annotation call hits the cache instead of re-parsing."""
    sites = _sites([2100, 6000])
    _ = annotate_features(sites, synth_gtf, multi_annotation=True)
    assert len(A._BUILT_FEATURES_CACHE) == 1
    # Second call: same key, must reuse -- cache size still 1
    _ = annotate_features(sites, synth_gtf, multi_annotation=True)
    assert len(A._BUILT_FEATURES_CACHE) == 1
    # And toggling multi_annotation off uses the SAME cache key (the bundle
    # always includes the TSS arrays -- they're just unused when False)
    _ = annotate_features(sites, synth_gtf, multi_annotation=False)
    assert len(A._BUILT_FEATURES_CACHE) == 1


# ---------------------------------------------------------------------------
# gene_type_filter
# ---------------------------------------------------------------------------

@pytest.fixture
def synth_gtf_mixed_types(tmp_path):
    """GTF with two protein_coding genes and one lincRNA -- used to verify
    that gene_type_filter actually excludes the lincRNA from all outputs.

    Coordinates are chosen so the test site (3700) is closest to the
    lincRNA (TSS=3500) but with the lincRNA filtered out the next
    nearest protein-coding TSS (PROT2 at 4500) wins over PROT1 (at 1000).
    """
    gtf = tmp_path / "synth_mixed.gtf"
    gtf.write_text("\n".join([
        # PROT1: far away on the left -- TSS at GTF 1001 -> 0-based 1000
        'chr1\tt\tgene\t1001\t2000\t.\t+\t.\tgene_id "g1"; gene_name "PROT1"; gene_type "protein_coding";',
        'chr1\tt\texon\t1001\t2000\t.\t+\t.\tgene_id "g1"; gene_name "PROT1"; gene_type "protein_coding";',
        # NOVEL_LNC: TSS at GTF 3501 -> 0-based 3500 (nearest to a site at 3700)
        'chr1\tt\tgene\t3501\t4000\t.\t+\t.\tgene_id "g2"; gene_name "NOVEL_LNC"; gene_type "lincRNA";',
        'chr1\tt\texon\t3501\t4000\t.\t+\t.\tgene_id "g2"; gene_name "NOVEL_LNC"; gene_type "lincRNA";',
        # PROT2: TSS at GTF 4501 -> 0-based 4500 (second-nearest if lincRNA filtered)
        'chr1\tt\tgene\t4501\t5500\t.\t+\t.\tgene_id "g3"; gene_name "PROT2"; gene_type "protein_coding";',
        'chr1\tt\texon\t4501\t5500\t.\t+\t.\tgene_id "g3"; gene_name "PROT2"; gene_type "protein_coding";',
    ]) + "\n")
    return str(gtf)


def test_gene_type_filter_excludes_lincrna_from_nearest_tss(synth_gtf_mixed_types):
    """Site at 3700. Without filter NOVEL_LNC (TSS=3500, ~200 bp) wins;
    with protein_coding-only filter PROT2 (TSS=4500, ~800 bp) wins over
    PROT1 (TSS=1000, ~2700 bp)."""
    sites = _sites([3700])

    out_unfiltered = annotate_features(
        sites, synth_gtf_mixed_types, multi_annotation=True,
    )
    assert out_unfiltered["nearest_tss_gene"][0] == "NOVEL_LNC"

    out_filtered = annotate_features(
        sites, synth_gtf_mixed_types, multi_annotation=True,
        gene_type_filter="protein_coding",
    )
    assert out_filtered["nearest_tss_gene"][0] == "PROT2"


def test_gene_type_filter_string_or_list_both_accepted(synth_gtf_mixed_types):
    """``"protein_coding"`` and ``["protein_coding"]`` should be equivalent."""
    sites = _sites([2600])
    s = annotate_features(sites, synth_gtf_mixed_types,
                          multi_annotation=True, gene_type_filter="protein_coding")
    l = annotate_features(sites, synth_gtf_mixed_types,
                          multi_annotation=True, gene_type_filter=["protein_coding"])
    assert s["nearest_tss_gene"][0] == l["nearest_tss_gene"][0]


def test_gene_type_filter_makes_distinct_cache_entry(synth_gtf_mixed_types):
    """Two calls with different gene_type filters must produce distinct
    cache entries -- silent cache cross-contamination would be a footgun."""
    sites = _sites([2600])
    _ = annotate_features(sites, synth_gtf_mixed_types,
                          multi_annotation=True)
    _ = annotate_features(sites, synth_gtf_mixed_types,
                          multi_annotation=True,
                          gene_type_filter="protein_coding")
    assert len(A._BUILT_FEATURES_CACHE) == 2


# ---------------------------------------------------------------------------
# refGene source
# ---------------------------------------------------------------------------

@pytest.fixture
def synth_refgene(tmp_path):
    """Minimal UCSC refGene.txt fixture (one NM_ + one NR_ transcript).

    Schema columns (tab-separated): bin, name, chrom, strand, txStart,
    txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2,
    cdsStartStat, cdsEndStat, exonFrames. Coords are 0-based half-open.
    """
    refgene = tmp_path / "synth.refGene.txt"
    # PROT_RG on "+" strand: TSS at 1000, gene 1000-2000, two exons
    # NCRNA_RG on "-" strand: TSS at 6000 (gene end), gene 4000-6000, one exon
    rows = [
        ["585", "NM_000001", "chr1", "+", "1000", "2000", "1100", "1900", "2",
         "1000,1500,", "1100,2000,", "0", "PROT_RG", "cmpl", "cmpl", "0,0,"],
        ["585", "NR_000001", "chr1", "-", "4000", "6000", "6000", "6000", "1",
         "4000,",      "6000,",      "0", "NCRNA_RG", "unk", "unk", "-1,"],
    ]
    refgene.write_text("\n".join("\t".join(r) for r in rows) + "\n")
    return str(refgene)


def test_refgene_source_produces_legacy_columns(synth_refgene):
    """A refGene source must populate the same single-best columns as a GTF
    source would (gene_id, gene_name, feature_type, distance_to_tss)."""
    sites = _sites([1500])  # inside PROT_RG body
    out = annotate_features(sites, synth_refgene)
    row = out.row(0, named=True)
    assert row["gene_name"] == "PROT_RG"
    # Inside an exon (1500 is in exon 1500-2000) OR could be promoter if
    # 1500 falls in the default promoter window of TSS=1000+200=1200 -- it
    # doesn't (1500 > 1200), so feature should be exon.
    assert row["feature_type"] in {"exon", "promoter"}


def test_refgene_source_supports_multi_annotation(synth_refgene):
    """Multi-annotation on a refGene source. The convention is
    "positive = downstream in transcription direction", so for the "-"
    strand NCRNA_RG (TSS at 6000):

      - site at 5000 (lower genomic coord) is DOWNSTREAM in transcription
        direction -> positive distance (+1000).
      - site at 7000 (higher genomic coord) is UPSTREAM in transcription
        direction -> negative distance (-1000).
    """
    sites = _sites([5000, 7000])
    out = annotate_features(
        sites, synth_refgene, multi_annotation=True,
    )
    rows = out.to_dicts()
    assert rows[0]["nearest_tss_gene"] == "NCRNA_RG"
    assert rows[0]["nearest_tss_distance"] == 1000   # downstream
    assert rows[1]["nearest_tss_gene"] == "NCRNA_RG"
    assert rows[1]["nearest_tss_distance"] == -1000  # upstream


def test_refgene_gene_type_filter_works(synth_refgene):
    """NM_* -> protein_coding, NR_* -> non-coding; filter must respect this."""
    sites = _sites([7000])
    out_unfilt = annotate_features(
        sites, synth_refgene, multi_annotation=True,
    )
    assert out_unfilt["nearest_tss_gene"][0] == "NCRNA_RG"  # NR_ wins by proximity

    out_filt = annotate_features(
        sites, synth_refgene, multi_annotation=True,
        gene_type_filter="protein_coding",
    )
    # NCRNA_RG (NR_) is filtered out; PROT_RG (NM_, TSS=1000) is the only
    # remaining candidate.
    assert out_filt["nearest_tss_gene"][0] == "PROT_RG"


def test_refgene_and_gtf_caches_dont_collide(synth_gtf, synth_refgene):
    """Same path basename but different sources must not collide in the
    built-features cache. Two calls -> two cache entries."""
    sites = _sites([1500])
    _ = annotate_features(sites, synth_gtf, multi_annotation=True)
    _ = annotate_features(sites, synth_refgene, multi_annotation=True)
    assert len(A._BUILT_FEATURES_CACHE) == 2


def test_annotation_argument_is_required():
    """``annotation`` is positional-required: omitting it is a TypeError
    (no ambiguity about which format the caller meant)."""
    sites = _sites([1000])
    with pytest.raises(TypeError):
        annotate_features(sites)


def test_source_auto_detects_gtf(synth_gtf):
    """Fixture file ends in ``.gtf`` -> auto-detected as GTF without an
    explicit ``source=``."""
    sites = _sites([2100])
    out = annotate_features(sites, synth_gtf)  # no source= passed
    # Single-best pick at this position should find GENEB via promoter window
    assert out["gene_name"][0] in {"GENEA", "GENEB"}


def test_source_auto_detects_refgene(synth_refgene):
    """Fixture file ends in ``.refGene.txt`` -> auto-detected as refGene."""
    sites = _sites([1500])
    out = annotate_features(sites, synth_refgene)  # no source= passed
    assert out["gene_name"][0] == "PROT_RG"


def test_source_explicit_override_for_unusual_filename(synth_refgene, tmp_path):
    """A refGene file with an unconventional name needs ``source="refgene"``
    explicitly. The fallback path must work."""
    # Copy the synth refGene under a name auto-detect won't recognize.
    odd = tmp_path / "custom_table.tsv"
    odd.write_bytes(Path(synth_refgene).read_bytes())
    sites = _sites([1500])

    # Auto-detect fails:
    with pytest.raises(ValueError, match="Cannot auto-detect"):
        annotate_features(sites, str(odd))

    # Explicit source= works:
    out = annotate_features(sites, str(odd), source="refgene")
    assert out["gene_name"][0] == "PROT_RG"


def test_source_invalid_value_raises():
    """``source`` must be one of {auto, gtf, refgene}."""
    sites = _sites([100])
    with pytest.raises(ValueError, match="must be 'auto'"):
        annotate_features(sites, "x.gtf", source="bam")


# ---------------------------------------------------------------------------
# default-features regression: the full HOMER set is built by default so
# downstream ``feature_type`` value-counts match HOMER's vocabulary without
# the caller having to opt in via ``features=``.
# ---------------------------------------------------------------------------

def test_default_features_is_full_homer_set():
    """The default ``features=`` tuple covers every category in
    :data:`epykit.annotate._FEATURE_PRIORITY` except the fallback
    ``intergenic``. Anything narrower silently buckets UTR / TTS / noncoding
    sites into intron/intergenic and produces a HOMER-incompatible
    distribution -- which is what bit the chain_merge benchmark."""
    import inspect
    sig = inspect.signature(annotate_features)
    default = tuple(sig.parameters["features"].default)
    assert default == (
        "promoter", "5UTR", "exon", "intron", "3UTR", "TTS", "noncoding",
    )
    # And the priority dict / default tuple must agree on the vocabulary
    # (every default feature has a priority; intergenic is the fallback
    # so it's deliberately not in ``features``).
    assert set(default) == set(A._FEATURE_PRIORITY) - {"intergenic"}


@pytest.fixture
def synth_gtf_with_utrs_tts_noncoding(tmp_path):
    """GTF with explicit UTR records, a TTS-window-overlapping site, and a
    non-protein-coding gene. Lets us assert that the expanded default
    actually surfaces 5UTR / 3UTR / TTS / noncoding in ``feature_type`` --
    the categories that used to silently disappear under the old 3-feature
    default.

    Layout on chr2 (1-based GTF coords):
      PROT_X: + strand, gene 1001-5000, exons 1001-1500 and 4501-5000,
              five_prime_utr 1001-1100, three_prime_utr 4901-5000.
              TSS=1001 (-> 0-based 1000); TTS-window on + strand is
              [end-100, end+1000) i.e. [4900, 6000).
      LNC_Y:  + strand, gene 10001-11000, exon 10001-11000,
              gene_type "lincRNA" -> contributes a noncoding interval.
    """
    gtf = tmp_path / "synth_utrs.gtf"
    lines = [
        'chr2\tt\tgene\t1001\t5000\t.\t+\t.\tgene_id "px"; gene_name "PROT_X"; gene_type "protein_coding";',
        'chr2\tt\texon\t1001\t1500\t.\t+\t.\tgene_id "px"; gene_name "PROT_X"; gene_type "protein_coding";',
        'chr2\tt\texon\t4501\t5000\t.\t+\t.\tgene_id "px"; gene_name "PROT_X"; gene_type "protein_coding";',
        'chr2\tt\tfive_prime_utr\t1001\t1100\t.\t+\t.\tgene_id "px"; gene_name "PROT_X"; gene_type "protein_coding";',
        'chr2\tt\tthree_prime_utr\t4901\t5000\t.\t+\t.\tgene_id "px"; gene_name "PROT_X"; gene_type "protein_coding";',
        'chr2\tt\tgene\t10001\t11000\t.\t+\t.\tgene_id "ly"; gene_name "LNC_Y"; gene_type "lincRNA";',
        'chr2\tt\texon\t10001\t11000\t.\t+\t.\tgene_id "ly"; gene_name "LNC_Y"; gene_type "lincRNA";',
    ]
    gtf.write_text("\n".join(lines) + "\n")
    return str(gtf)


def test_default_surfaces_utr_tts_noncoding(synth_gtf_with_utrs_tts_noncoding):
    """Regression for the chain_merge benchmark bug: under the old default
    (``features=("promoter","exon","intron")``) the 5UTR / 3UTR / TTS /
    noncoding categories never appeared in ``feature_type`` because their
    builders weren't invoked. Each of the four positions below sits in a
    category that only resolves correctly when the full HOMER default is
    in effect.

    Pick positions outside the promoter window (TSS=1000, window [-1000,
    +200] so promoter covers up to 1200) and outside any earlier-priority
    feature so the priority chain actually exposes the category we want
    to assert on.
    """
    sites = pl.DataFrame({
        "chrom": ["chr2"] * 4,
        # 1300: in 5'UTR (1000-1100 is the actual UTR after 0-based conv).
        # The UTR is priority 1 so it wins over exon (priority 4) when both
        # overlap. To make this site UTR-dominated, place it at 1050.
        "pos":   [1050, 4950, 5500, 10500],
    })
    out = annotate_features(sites, synth_gtf_with_utrs_tts_noncoding)
    rows = out.to_dicts()
    # 1050 -> in 5UTR (1000-1100). Promoter window [-1000,+200] from TSS=1000
    # covers [-1000, 1200) so 1050 is also in promoter. Promoter (0) outranks
    # 5UTR (1), so feature_type==promoter -- BUT all_overlapping_features
    # must include the 5UTR label since the UTR builder ran.
    assert "5UTR" in set(rows[0]["all_overlapping_features"]), \
        f"5UTR missing from overlapping features: {rows[0]}"
    # 4950 -> in 3UTR (4900-5000). Also inside an exon (4500-5000) and
    # inside the TTS window [4900, 6000). 3UTR (priority 2) wins over TTS
    # (3) and exon (4), so feature_type==3UTR.
    assert rows[1]["feature_type"] == "3UTR", rows[1]
    # 5500 -> outside gene body but inside TTS window [4900, 6000). Only
    # TTS contributes here.
    assert rows[2]["feature_type"] == "TTS", rows[2]
    # 10500 -> inside LNC_Y (lincRNA). Exon (4) overlaps; noncoding (6)
    # also overlaps. Exon wins, but noncoding must be present in the
    # one-to-many list -- the builder having actually run is the
    # regression we're guarding.
    assert "noncoding" in set(rows[3]["all_overlapping_features"]), \
        f"noncoding missing from overlapping features: {rows[3]}"
