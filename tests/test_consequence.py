"""Tests for Layer 2 consequence classification.

These drive ``classify_consequence`` directly with hand-built ``GenicContext``
objects so each rule -- and each graceful-degradation path -- is exercised
deterministically, independent of the GENCODE machinery.
"""

from meiva.annotate.consequence import (
    FULL_LENGTH_L1_BP,
    IMPACT_SEVERITY,
    Consequence,
    ConsequenceResult,
    Impact,
    classify_consequence,
)
from meiva.annotate.genic import GenicContext, GenicRegion, InsertionOrientation
from meiva.model import MEIFamily, MEISite, Strand


def _site(family: MEIFamily, length: int | None = None) -> MEISite:
    # strand on the site is irrelevant to Layer 2 (it reads ctx.orientation);
    # family and length are what the rules consume.
    return MEISite(chrom="chr1", pos=1000, family=family, strand=Strand.UNKNOWN, length=length)


def _ctx(
    region: GenicRegion,
    orientation: InsertionOrientation = InsertionOrientation.UNKNOWN,
    biotype: str | None = None,
) -> GenicContext:
    return GenicContext(region=region, orientation=orientation, gene_biotype=biotype)


def _classify(
    family: MEIFamily,
    region: GenicRegion,
    orientation: InsertionOrientation = InsertionOrientation.UNKNOWN,
    length: int | None = None,
    biotype: str | None = None,
) -> ConsequenceResult:
    return classify_consequence(_site(family, length), _ctx(region, orientation, biotype))


# --------------------------------------------------------------------------- #
# Region ceilings that don't depend on family/orientation                     #
# --------------------------------------------------------------------------- #
def test_cds_is_coding_disruption_high():
    r = _classify(MEIFamily.ALU, GenicRegion.CDS)
    assert r.consequence is Consequence.CODING_DISRUPTION
    assert r.impact is Impact.HIGH
    assert r.flags == ()


def test_splice_region_is_splice_disruption_high():
    r = _classify(MEIFamily.L1, GenicRegion.SPLICE_REGION)
    assert r.consequence is Consequence.SPLICE_DISRUPTION
    assert r.impact is Impact.HIGH


def test_promoter_is_moderate():
    assert _classify(MEIFamily.ALU, GenicRegion.PROMOTER).impact is Impact.MODERATE


def test_utr5_outranks_utr3():
    utr5 = _classify(MEIFamily.ALU, GenicRegion.UTR5)
    utr3 = _classify(MEIFamily.ALU, GenicRegion.UTR3)
    assert utr5.impact is Impact.MODERATE
    assert utr3.impact is Impact.LOW
    assert IMPACT_SEVERITY[utr5.impact] > IMPACT_SEVERITY[utr3.impact]


def test_flanking_and_intergenic_are_modifiers():
    assert _classify(MEIFamily.ALU, GenicRegion.UPSTREAM).impact is Impact.MODIFIER
    assert _classify(MEIFamily.ALU, GenicRegion.DOWNSTREAM).impact is Impact.MODIFIER
    inter = _classify(MEIFamily.ALU, GenicRegion.INTERGENIC)
    assert inter.consequence is Consequence.INTERGENIC
    assert inter.impact is Impact.MODIFIER


# --------------------------------------------------------------------------- #
# Alu intronic: antisense -> exonization; otherwise plain intronic            #
# --------------------------------------------------------------------------- #
def test_alu_intron_antisense_is_exonization():
    r = _classify(MEIFamily.ALU, GenicRegion.INTRON, InsertionOrientation.ANTISENSE)
    assert r.consequence is Consequence.EXONIZATION_CANDIDATE
    assert r.impact is Impact.MODERATE
    assert r.flags == ()


def test_alu_intron_sense_is_plain_intronic():
    r = _classify(MEIFamily.ALU, GenicRegion.INTRON, InsertionOrientation.SENSE)
    assert r.consequence is Consequence.INTRONIC
    assert r.flags == ()  # sense is a definite call, no degradation


def test_alu_intron_unknown_orientation_degrades_with_flag():
    r = _classify(MEIFamily.ALU, GenicRegion.INTRON, InsertionOrientation.UNKNOWN)
    assert r.consequence is Consequence.INTRONIC
    assert "orientation_unknown" in r.flags


# --------------------------------------------------------------------------- #
# L1 intronic: full-length sense -> polyA interference; else degrade          #
# --------------------------------------------------------------------------- #
def test_l1_intron_sense_full_length_is_polya_interference():
    r = _classify(MEIFamily.L1, GenicRegion.INTRON, InsertionOrientation.SENSE, length=6000)
    assert r.consequence is Consequence.POLYA_INTERFERENCE
    assert r.impact is Impact.MODERATE
    assert "full_length_l1" in r.flags


def test_l1_intron_sense_truncated_is_intronic():
    r = _classify(
        MEIFamily.L1, GenicRegion.INTRON, InsertionOrientation.SENSE, length=FULL_LENGTH_L1_BP - 1
    )
    assert r.consequence is Consequence.INTRONIC
    assert "full_length_l1" not in r.flags


def test_l1_intron_sense_unknown_length_degrades_with_flag():
    r = _classify(MEIFamily.L1, GenicRegion.INTRON, InsertionOrientation.SENSE, length=None)
    assert r.consequence is Consequence.INTRONIC
    assert "length_unknown" in r.flags


def test_l1_intron_antisense_is_intronic_no_flag():
    r = _classify(MEIFamily.L1, GenicRegion.INTRON, InsertionOrientation.ANTISENSE, length=6000)
    assert r.consequence is Consequence.INTRONIC
    assert r.flags == ()  # antisense L1 is better tolerated; a definite call


def test_l1_intron_unknown_orientation_degrades_with_flag():
    r = _classify(MEIFamily.L1, GenicRegion.INTRON, InsertionOrientation.UNKNOWN, length=6000)
    assert r.consequence is Consequence.INTRONIC
    assert "orientation_unknown" in r.flags


# --------------------------------------------------------------------------- #
# Non-coding exon: elevated for a lncRNA host (ties in the biotype work)       #
# --------------------------------------------------------------------------- #
def test_noncoding_exon_in_lncrna_is_elevated():
    r = _classify(MEIFamily.ALU, GenicRegion.EXON_NONCODING, biotype="lncRNA")
    assert r.consequence is Consequence.NONCODING_EXON_INSERTION
    assert r.impact is Impact.MODERATE


def test_noncoding_exon_non_lncrna_stays_low():
    r = _classify(MEIFamily.ALU, GenicRegion.EXON_NONCODING, biotype="processed_pseudogene")
    assert r.consequence is Consequence.NONCODING_EXON_INSERTION
    assert r.impact is Impact.LOW


def test_noncoding_exon_unknown_biotype_stays_low():
    assert _classify(MEIFamily.ALU, GenicRegion.EXON_NONCODING).impact is Impact.LOW


# --------------------------------------------------------------------------- #
# SVA regulatory flag                                                         #
# --------------------------------------------------------------------------- #
def test_sva_in_promoter_gets_regulatory_flag():
    r = _classify(MEIFamily.SVA, GenicRegion.PROMOTER)
    assert r.consequence is Consequence.PROMOTER_INSERTION
    assert "sva_regulatory" in r.flags


def test_non_sva_in_promoter_has_no_regulatory_flag():
    r = _classify(MEIFamily.ALU, GenicRegion.PROMOTER)
    assert "sva_regulatory" not in r.flags


def test_sva_in_intron_has_no_regulatory_flag():
    # intron is not a regulatory context for the SVA flag
    r = _classify(MEIFamily.SVA, GenicRegion.INTRON)
    assert r.consequence is Consequence.INTRONIC
    assert "sva_regulatory" not in r.flags


# --------------------------------------------------------------------------- #
# Severity ordering sanity                                                    #
# --------------------------------------------------------------------------- #
def test_impact_severity_is_strictly_ordered():
    assert (
        IMPACT_SEVERITY[Impact.HIGH]
        > IMPACT_SEVERITY[Impact.MODERATE]
        > IMPACT_SEVERITY[Impact.LOW]
        > IMPACT_SEVERITY[Impact.MODIFIER]
    )
