"""Tests for the xTEA parser, run against real VDA xTEA output.

The fixtures in tests/data are genuine xTEA single-sample VCFs, so these tests
double as a regression guard on the exact format quirks we designed around.
"""

from pathlib import Path

import pytest

from meiva.io import XteaParser, detect_parser
from meiva.model import MEIFamily, Strand

DATA = Path(__file__).parent / "data"
SAMPLE_VCF = DATA / "SAMPLE_REDACTED.vcf"

# The real VDA VCFs are gitignored (real cohort data), so they're absent in CI
# and fresh clones. Skip rather than error; re-activates once a committable
# fixture (e.g. an HGDP/TSI xTEA VCF) is placed in tests/data.
pytestmark = pytest.mark.skipif(
    not SAMPLE_VCF.exists(), reason="xTEA VCF fixture not present (gitignored real data)"
)


@pytest.fixture(scope="module")
def sites_by_pos():
    """All sites from the sample VCF, keyed by genomic position."""
    parser = XteaParser()
    sites = list(parser.parse(SAMPLE_VCF))
    return {s.pos: s for s in sites}


# --------------------------------------------------------------------------- #
# Detection                                                                   #
# --------------------------------------------------------------------------- #
def test_sniff_recognises_xtea():
    assert XteaParser.sniff(SAMPLE_VCF) is True


def test_detect_parser_picks_xtea():
    assert isinstance(detect_parser(SAMPLE_VCF), XteaParser)


# --------------------------------------------------------------------------- #
# Whole-file invariants                                                       #
# --------------------------------------------------------------------------- #
def test_parses_every_record(sites_by_pos):
    n_records = sum(
        1 for line in SAMPLE_VCF.read_text().splitlines() if line and not line.startswith("#")
    )
    assert len(sites_by_pos) == n_records


def test_all_sites_chr_prefixed_and_attributed(sites_by_pos):
    for site in sites_by_pos.values():
        assert site.chrom.startswith("chr")
        assert site.source_caller == "xtea"
        assert site.filters == ("PASS",)
        assert site.ci_lower == 0 and site.ci_upper == 0  # xTEA reports no CIPOS


# --------------------------------------------------------------------------- #
# Specific records (values read directly from the VCF)                        #
# --------------------------------------------------------------------------- #
def test_alu_het_record(sites_by_pos):
    s = sites_by_pos[6_188_790]
    assert s.family is MEIFamily.ALU
    assert s.chrom == "chr1"
    assert s.strand is Strand.MINUS
    assert s.length == 244
    assert s.tsd == "CATTTTTTTTTTTTTT"  # leading '+' stripped
    assert s.allele_freq == pytest.approx(0.3333333, abs=1e-5)
    gt = s.genotypes["SAMPLE_REDACTED"]
    assert gt.is_carrier and gt.dosage == 1


def test_sva_unknown_strand_hom(sites_by_pos):
    s = sites_by_pos[13_396_328]
    assert s.family is MEIFamily.SVA
    assert s.strand is Strand.UNKNOWN  # STRAND=.
    assert s.length == 283
    assert s.allele_freq == pytest.approx(1.0)
    assert s.genotypes["SAMPLE_REDACTED"].dosage == 2


def test_transduction_length_is_nulled_and_flagged(sites_by_pos):
    # orphan transduction: SVLEN=35139274 is a coordinate, not a length
    s = sites_by_pos[22_950_447]
    assert s.family is MEIFamily.L1
    assert s.length is None
    assert s.raw_info.get("MEIVA_LENGTH_UNRELIABLE") == "1"
    # the raw (bogus) value is preserved, never silently discarded
    assert s.raw_info["SVLEN"] == "35139274"


def test_null_tsd_becomes_none(sites_by_pos):
    s = sites_by_pos[27_878_440]
    assert s.tsd is None  # TSD=NULL
    assert s.strand is Strand.UNKNOWN


def test_raw_info_is_lossless(sites_by_pos):
    # a caller-specific field we don't promote must still survive in raw_info
    s = sites_by_pos[6_188_790]
    assert "GENE_INFO" in s.raw_info
    assert "RPL22" in s.raw_info["GENE_INFO"]
