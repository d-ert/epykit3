"""Unit tests for the shared canonical-chromosome helper (``epykit._chroms``).

The helper is the single definition of "main-assembly chromosome" behind the
opt-in ``canonical_only`` options. These tests pin its naming rules, the
fixed human-style set, order-preserving filtering, and the INFO audit line.
"""

import logging

import pytest

from epykit._chroms import (
    CANONICAL_CHROM_CORES,
    CANONICAL_CHROMS_UCSC,
    filter_canonical,
    filter_canonical_logged,
    is_canonical_chrom,
)

_LOGGER = "epykit._chroms"


@pytest.mark.parametrize(
    "name",
    [
        # UCSC convention
        "chr1",
        "chr9",
        "chr22",
        "chrX",
        "chrY",
        # Ensembl convention (no prefix)
        "1",
        "9",
        "22",
        "X",
        "Y",
    ],
)
def test_canonical_kept_under_both_naming_styles(name):
    assert is_canonical_chrom(name)


@pytest.mark.parametrize("name", ["chrM", "chrMT", "M", "MT", "chrm", "mt"])
def test_mitochondrion_kept_under_every_spelling(name):
    assert is_canonical_chrom(name)


@pytest.mark.parametrize("name", ["CHR1", "Chr1", "cHrX", "chrx", "x"])
def test_prefix_and_core_are_case_insensitive(name):
    assert is_canonical_chrom(name)


@pytest.mark.parametrize(
    "name",
    [
        # unplaced / unlocalised / alt contigs (the real GSE263850 offenders)
        "chr14_KI270722v1_random",
        "chr1_KI270706v1_random",
        "chr14_GL000194v1_random",
        "chrUn_KI270742v1",
        "chrUn_GL000216v2",
        "GL000216v2",
        "KI270722.1",
        # out-of-range / malformed
        "chr23",
        "chr0",
        "chr01",
        "chr",
        "chrUn",
        "",
        # not a human-style name: the helper is not species-aware
        "chrI",
        "2L",
    ],
)
def test_noncanonical_dropped(name):
    assert not is_canonical_chrom(name)


def test_cores_are_the_fixed_human_style_set():
    assert CANONICAL_CHROM_CORES == frozenset(
        [str(i) for i in range(1, 23)] + ["X", "Y", "M", "MT"]
    )


def test_ucsc_order_is_genome_order_and_all_canonical():
    assert CANONICAL_CHROMS_UCSC[:3] == ("chr1", "chr2", "chr3")
    assert CANONICAL_CHROMS_UCSC[-3:] == ("chrX", "chrY", "chrM")
    assert len(CANONICAL_CHROMS_UCSC) == 25
    assert all(is_canonical_chrom(c) for c in CANONICAL_CHROMS_UCSC)


def test_filter_canonical_preserves_input_order():
    chroms = [
        "chr3",
        "chrUn_KI270742v1",
        "chr1",
        "chr14_GL000194v1_random",
        "chrX",
        "MT",
    ]
    assert filter_canonical(chroms) == ["chr3", "chr1", "chrX", "MT"]


def test_filter_canonical_accepts_any_iterable_and_empty_input():
    assert filter_canonical([]) == []
    assert filter_canonical(iter(("chr2", "GL000216v2"))) == ["chr2"]


def test_filter_canonical_logged_names_dropped_contigs_and_the_opt_out(caplog):
    chroms = ["chr1", "chrUn_KI270742v1", "chr2", "chr14_GL000194v1_random"]
    with caplog.at_level(logging.INFO, logger=_LOGGER):
        kept = filter_canonical_logged(chroms, context="dmc")

    assert kept == ["chr1", "chr2"]
    records = [r for r in caplog.records if r.name == _LOGGER]
    assert len(records) == 1, caplog.text
    msg = records[0].getMessage()
    assert records[0].levelno == logging.INFO
    assert "[dmc]" in msg
    assert "keeping 2 canonical chromosome(s)" in msg
    assert "dropping 2 contig(s)" in msg
    assert "chrUn_KI270742v1, chr14_GL000194v1_random" in msg
    assert "Omit canonical_only=True to retain them" in msg
    assert "--all-contigs" not in msg


def test_filter_canonical_logged_summarises_a_long_drop_list(caplog):
    dropped = [f"chrUn_GL{i:06d}v1" for i in range(8)]
    with caplog.at_level(logging.INFO, logger=_LOGGER):
        kept = filter_canonical_logged(["chr1", *dropped])

    assert kept == ["chr1"]
    msg = next(r.getMessage() for r in caplog.records if r.name == _LOGGER)
    assert msg.startswith("canonical_only:")  # no context tag without ``context``
    assert ", ".join(dropped[:5]) in msg
    assert dropped[5] not in msg
    assert "and 3 more" in msg


def test_filter_canonical_logged_is_silent_when_nothing_is_dropped(caplog):
    with caplog.at_level(logging.INFO, logger=_LOGGER):
        kept = filter_canonical_logged(["chr1", "X", "MT"], context="dmr/tile")

    assert kept == ["chr1", "X", "MT"]
    assert [r for r in caplog.records if r.name == _LOGGER] == []
