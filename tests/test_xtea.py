"""Tests for the xTEA parser, run against committed public HGDP/TSI xTEA output.

``tests/data/HGDP01162.tsi.vcf`` is a genuine xTEA single-sample VCF from an
openly consented HGDP Tuscan (TSI) sample, trimmed to chr1. Because it is
committed -- unlike the gitignored VDA cohort data -- these tests run in CI and
double as a regression guard on the exact xTEA format quirks we designed around.
"""

from pathlib import Path

import pytest

from meiva.io import XteaParser, detect_parser
from meiva.model import MEIFamily, Strand

DATA = Path(__file__).parent / "data"
SAMPLE_VCF = DATA / "HGDP01162.tsi.vcf"
SAMPLE = "HGDP01162.alt_bwamem_GRCh38DH.20181023.Tuscan"


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
    assert len(sites_by_pos) == n_records == 27


def test_all_sites_chr_prefixed_and_attributed(sites_by_pos):
    for site in sites_by_pos.values():
        assert site.chrom.startswith("chr")
        assert site.source_caller == "xtea"
        assert site.filters == ("PASS",)
        assert site.ci_lower == 0 and site.ci_upper == 0  # xTEA reports no CIPOS


# --------------------------------------------------------------------------- #
# Specific records (values read directly from the VCF)                        #
# --------------------------------------------------------------------------- #
def test_alu_het_minus_strand_record(sites_by_pos):
    s = sites_by_pos[72_173_322]
    assert s.family is MEIFamily.ALU
    assert s.chrom == "chr1"
    assert s.strand is Strand.MINUS
    assert s.length == 237
    assert s.tsd == "AGCAATCTTATTTTC"  # leading '+' stripped
    assert s.allele_freq == pytest.approx(0.6666667, abs=1e-5)
    gt = s.genotypes[SAMPLE]
    assert gt.is_carrier and gt.dosage == 1


def test_l1_unknown_strand_hom_null_tsd(sites_by_pos):
    s = sites_by_pos[72_339_812]
    assert s.family is MEIFamily.L1
    assert s.strand is Strand.UNKNOWN  # STRAND=.
    assert s.length == 371
    assert s.tsd is None  # TSD=NULL
    assert s.allele_freq == pytest.approx(1.0)
    assert s.genotypes[SAMPLE].dosage == 2


def test_alu_null_tsd_unknown_strand(sites_by_pos):
    s = sites_by_pos[27_878_440]
    assert s.family is MEIFamily.ALU
    assert s.tsd is None  # TSD=NULL
    assert s.strand is Strand.UNKNOWN


def test_raw_info_is_lossless(sites_by_pos):
    # a caller-specific field we don't promote must still survive in raw_info
    s = sites_by_pos[72_173_322]
    assert "GENE_INFO" in s.raw_info
    assert "NEGR1" in s.raw_info["GENE_INFO"]


# --------------------------------------------------------------------------- #
# Transduction length guard -- synthetic & deterministic (no fixture needed).  #
# xTEA can emit a genomic coordinate in SVLEN for orphan/transduction calls;   #
# we null the implausible value and flag it rather than report a 35 Mb insert. #
# --------------------------------------------------------------------------- #
_TRANSDUCTION_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##ALT=<ID=INS:ME:LINE1,Description="L1">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="">
##INFO=<ID=AF_FMAP,Number=1,Type=Integer,Description="xTEA signature">
##INFO=<ID=LC_CLUSTER,Number=1,Type=String,Description="xTEA signature">
##INFO=<ID=RD_AKR_NRC,Number=1,Type=Integer,Description="xTEA signature">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
chr1\t22950447\t.\tA\t<INS:ME:LINE1>\t.\tPASS\tSVTYPE=INS:ME:LINE1;SVLEN=35139274;AF_FMAP=0;LC_CLUSTER=.;RD_AKR_NRC=0\tGT\t0/1
"""


def test_transduction_length_is_nulled_and_flagged(tmp_path):
    vcf = tmp_path / "transduction.vcf"
    vcf.write_text(_TRANSDUCTION_VCF)
    assert XteaParser.sniff(vcf) is True
    (s,) = list(XteaParser().parse(vcf))
    assert s.family is MEIFamily.L1
    assert s.length is None  # 35 Mb is implausible -> nulled
    assert s.raw_info.get("MEIVA_LENGTH_UNRELIABLE") == "1"
    assert s.raw_info["SVLEN"] == "35139274"  # bogus value preserved, never discarded
