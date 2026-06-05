"""Tests for caller-agnostic ingest helpers."""

import pytest

from meiva.io.base import MAX_PLAUSIBLE_INSERTION_BP, normalize_contig, plausible_length


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("chr1", "chr1"),
        ("1", "chr1"),
        ("chrX", "chrX"),
        ("X", "chrX"),
        ("MT", "chrM"),
        ("M", "chrM"),
        ("chrM", "chrM"),
    ],
)
def test_normalize_contig(raw, expected):
    assert normalize_contig(raw) == expected


def test_normalize_contig_empty():
    with pytest.raises(ValueError):
        normalize_contig("   ")


@pytest.mark.parametrize(
    "svlen,length,unreliable",
    [
        (None, None, False),
        (244, 244, False),
        (2584, 2584, False),
        (MAX_PLAUSIBLE_INSERTION_BP, MAX_PLAUSIBLE_INSERTION_BP, False),
        (MAX_PLAUSIBLE_INSERTION_BP + 1, None, True),
        (35_139_274, None, True),  # the real xTEA transduction value
        (0, None, True),
        (-5, None, True),
    ],
)
def test_plausible_length(svlen, length, unreliable):
    assert plausible_length(svlen) == (length, unreliable)
