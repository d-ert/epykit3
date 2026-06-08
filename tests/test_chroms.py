"""Unit tests for the shared canonical-chromosome helper (epykit._chroms)."""

import pytest

from epykit._chroms import (
    CANONICAL_CHROM_CORES,
    filter_canonical,
    is_canonical_chrom,
)


@pytest.mark.parametrize(
    "name",
    [
        # UCSC convention
        "chr1", "chr9", "chr22", "chrX", "chrY", "chrM", "chrMT",
        # Ensembl convention (no prefix)
        "1", "9", "22", "X", "Y", "M", "MT",
        # prefix case-insensitivity
        "CHR1", "Chr1",
    ],
)
def test_canonical_kept(name):
    assert is_canonical_chrom(name)


@pytest.mark.parametrize(
    "name",
    [
        # unplaced / unlocalized / alt contigs (the real GSE263850 offenders)
        "chr14_KI270722v1_random",
        "chr1_KI270706v1_random",
        "chr14_GL000194v1_random",
        "chrUn_KI270742v1",
        "chrUn_GL000216v2",
        "GL000216v2",
        "KI270722.1",
        # out-of-range / malformed
        "chr23", "chr0", "chr01", "chr", "chrUn", "",
    ],
)
def test_noncanonical_dropped(name):
    assert not is_canonical_chrom(name)


def test_filter_canonical_preserves_order():
    chroms = [
        "chr3",
        "chrUn_KI270742v1",
        "chr1",
        "chr14_GL000194v1_random",
        "chrX",
        "MT",
    ]
    assert filter_canonical(chroms) == ["chr3", "chr1", "chrX", "MT"]


def test_filter_canonical_empty():
    assert filter_canonical([]) == []


def test_cores_cover_human_autosomes_plus_sex_and_mito():
    # 22 autosomes + X, Y, M, MT
    assert CANONICAL_CHROM_CORES == frozenset(
        [str(i) for i in range(1, 23)] + ["X", "Y", "M", "MT"]
    )
