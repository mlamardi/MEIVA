"""Tests for the canonical MEI data model.

These pin the invariants the rest of the codebase will rely on. The
coordinate-arithmetic tests matter most: an off-by-one here propagates into
every annotation layer.
"""

from dataclasses import FrozenInstanceError

import pytest

from meiva.model import MEIFamily, MEISite, SampleGenotype, Strand


# --------------------------------------------------------------------------- #
# MEIFamily.from_raw                                                          #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "raw,expected",
    [
        ("ALU", MEIFamily.ALU),
        ("Alu", MEIFamily.ALU),
        ("AluYa5", MEIFamily.ALU),
        ("L1", MEIFamily.L1),
        ("LINE1", MEIFamily.L1),
        ("LINE/L1", MEIFamily.L1),
        ("L1HS", MEIFamily.L1),
        ("SVA", MEIFamily.SVA),
        ("SVA_E", MEIFamily.SVA),
        ("HERVK", MEIFamily.HERVK),
        ("HERV-K", MEIFamily.HERVK),
        ("nonsense", MEIFamily.OTHER),
        ("", MEIFamily.OTHER),
    ],
)
def test_family_from_raw(raw, expected):
    assert MEIFamily.from_raw(raw) is expected


# --------------------------------------------------------------------------- #
# Construction & validation                                                   #
# --------------------------------------------------------------------------- #
def test_minimal_valid_site():
    s = MEISite(chrom="chr6", pos=31_000_000, family=MEIFamily.ALU)
    assert s.is_precise
    assert s.strand is Strand.UNKNOWN


def test_chrom_must_be_prefixed():
    with pytest.raises(ValueError, match="chr"):
        MEISite(chrom="6", pos=100, family=MEIFamily.ALU)


def test_empty_chrom_rejected():
    with pytest.raises(ValueError):
        MEISite(chrom="", pos=100, family=MEIFamily.ALU)


def test_pos_must_be_one_based():
    with pytest.raises(ValueError, match="pos"):
        MEISite(chrom="chr1", pos=0, family=MEIFamily.L1)


@pytest.mark.parametrize("lo,hi", [(1, 0), (0, -1), (3, 5)])
def test_ci_sign_convention_enforced(lo, hi):
    with pytest.raises(ValueError):
        MEISite(chrom="chr1", pos=100, family=MEIFamily.L1, ci_lower=lo, ci_upper=hi)


def test_nonpositive_length_rejected():
    with pytest.raises(ValueError):
        MEISite(chrom="chr1", pos=100, family=MEIFamily.L1, length=0)


@pytest.mark.parametrize("af", [-0.1, 1.5, 2.0])
def test_allele_freq_out_of_range_rejected(af):
    with pytest.raises(ValueError, match="allele_freq"):
        MEISite(chrom="chr1", pos=100, family=MEIFamily.ALU, allele_freq=af)


@pytest.mark.parametrize("af", [0.0, 0.5, 1.0])
def test_allele_freq_in_range_ok(af):
    s = MEISite(chrom="chr1", pos=100, family=MEIFamily.ALU, allele_freq=af)
    assert s.allele_freq == af


def test_wrong_family_type_rejected():
    with pytest.raises(TypeError):
        MEISite(chrom="chr1", pos=100, family="ALU")  # type: ignore[arg-type]


# --------------------------------------------------------------------------- #
# search_interval: the arithmetic that everything downstream trusts           #
# --------------------------------------------------------------------------- #
def test_search_interval_precise():
    # pos=100, exact -> single base, 0-based [99, 100)
    s = MEISite(chrom="chr1", pos=100, family=MEIFamily.ALU)
    assert s.search_interval() == (99, 100)


def test_search_interval_with_ci():
    # pos=100, [-5,+7] -> 1-based 95..107 -> 0-based [94, 107)
    s = MEISite(chrom="chr1", pos=100, family=MEIFamily.ALU, ci_lower=-5, ci_upper=7)
    assert s.search_interval() == (94, 107)


def test_search_interval_padding():
    s = MEISite(chrom="chr1", pos=100, family=MEIFamily.ALU)
    assert s.search_interval(padding=50) == (49, 150)


def test_search_interval_clamps_lower_bound():
    s = MEISite(chrom="chr1", pos=3, family=MEIFamily.ALU)
    assert s.search_interval(padding=50)[0] == 0


def test_negative_padding_rejected():
    s = MEISite(chrom="chr1", pos=100, family=MEIFamily.ALU)
    with pytest.raises(ValueError):
        s.search_interval(padding=-1)


# --------------------------------------------------------------------------- #
# Immutability                                                                #
# --------------------------------------------------------------------------- #
def test_frozen_scalar():
    s = MEISite(chrom="chr1", pos=100, family=MEIFamily.ALU)
    with pytest.raises(FrozenInstanceError):
        s.pos = 200  # type: ignore[misc]


def test_raw_info_is_immutable_and_copied():
    src = {"FOO": "bar"}
    s = MEISite(chrom="chr1", pos=100, family=MEIFamily.ALU, raw_info=src)
    # external mutation of the source dict must not leak into the record
    src["FOO"] = "mutated"
    assert s.raw_info["FOO"] == "bar"
    with pytest.raises(TypeError):
        s.raw_info["FOO"] = "x"  # type: ignore[index]


# --------------------------------------------------------------------------- #
# Identity / equality keys on locus, not payload                              #
# --------------------------------------------------------------------------- #
def test_equality_ignores_payload():
    a = MEISite(chrom="chr1", pos=100, family=MEIFamily.ALU, qual=99.0)
    b = MEISite(
        chrom="chr1",
        pos=100,
        family=MEIFamily.ALU,
        qual=10.0,
        genotypes={"s1": SampleGenotype("s1", (0, 1))},
    )
    assert a == b
    assert hash(a) == hash(b)
    assert len({a, b}) == 1


def test_different_locus_not_equal():
    a = MEISite(chrom="chr1", pos=100, family=MEIFamily.ALU)
    b = MEISite(chrom="chr1", pos=100, family=MEIFamily.L1)
    assert a != b


# --------------------------------------------------------------------------- #
# SampleGenotype                                                              #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "alleles,missing,carrier,dosage",
    [
        ((0, 0), False, False, 0),
        ((0, 1), False, True, 1),
        ((1, 1), False, True, 2),
        ((None, None), True, False, None),
        ((None, 1), False, True, 1),  # partial-missing counts observed alt
        ((), True, False, None),
    ],
)
def test_sample_genotype(alleles, missing, carrier, dosage):
    gt = SampleGenotype("s1", alleles)
    assert gt.is_missing is missing
    assert gt.is_carrier is carrier
    assert gt.dosage == dosage
