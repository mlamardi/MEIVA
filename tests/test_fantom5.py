"""Tests for FANTOM5 enhancer and CAGE-promoter annotation.

Fixtures reproduce the real file formats, including the peak identifiers that
encode hg19 coordinates and the ``p1@GENE`` promoter-rank naming.
"""

import gzip

import pytest

from meiva.annotate.fantom5 import (
    Fantom5Model,
    RegulatoryContext,
    load_fantom5,
)

# BED12 enhancers, as in F5.hg38.enhancers.bed.gz (unstranded, ~300 bp)
_ENH = (
    "chr1\t1000\t1300\tchr1:1000-1300\t35\t.\t1150\t1151\t0,0,0\t2\t100,50\t0,250\n"
    "chr1\t5000\t5400\tchr1:5000-5400\t101\t.\t5200\t5201\t0,0,0\t2\t40,60\t0,340\n"
    "chr2\t2000\t2300\tchr2:2000-2300\t12\t.\t2150\t2151\t0,0,0\t2\t30,20\t0,280\n"
)

# BED9 CAGE peaks; col 7 (thickStart) is the representative TSS. Deliberately
# narrow, as the real median width is 14 bp.
_PEAKS = (
    "chr1\t8000\t8020\thg19::chr1:7000..7020,+;hg_1.1\t2972\t+\t8010\t8011\t255,0,0\n"
    "chr1\t9000\t9014\thg19::chr1:8000..8014,-;hg_2.1\t397\t-\t9007\t9008\t0,0,255\n"
    "chr2\t4000\t4030\thg19::chr2:3000..3030,+;hg_3.1\t587\t+\t4015\t4016\t255,0,0\n"
)

_NAMES = (
    "#TSS_ID\tNew_TSS_name\tOld_TSS_name\n"
    "hg19::chr1:7000..7020,+;hg_1.1\tp1@GENEA\n"
    "hg19::chr1:8000..8014,-;hg_2.1\tp3@GENEB\n"
    "hg19::chr2:3000..3030,+;hg_3.1\tp@GENEC\n"  # no rank
)


@pytest.fixture
def files(tmp_path):
    e = tmp_path / "enh.bed"
    p = tmp_path / "peaks.bed"
    n = tmp_path / "names.txt"
    e.write_text(_ENH)
    p.write_text(_PEAKS)
    n.write_text(_NAMES)
    return e, p, n


@pytest.fixture
def model(files):
    e, p, n = files
    return load_fantom5(enhancers=e, cage_peaks=p, peak_names=n)


# --------------------------------------------------------------------------- #
# Loading                                                                      #
# --------------------------------------------------------------------------- #
def test_counts(model):
    assert model.n_enhancers == 3
    assert model.n_peaks == 3


def test_requires_at_least_one_track():
    with pytest.raises(ValueError, match="at least one"):
        load_fantom5()


def test_peak_names_require_peaks(files):
    e, _, n = files
    with pytest.raises(ValueError, match="requires cage_peaks"):
        load_fantom5(enhancers=e, peak_names=n)


def test_reads_gzip(tmp_path):
    p = tmp_path / "enh.bed.gz"
    p.write_bytes(gzip.compress(_ENH.encode()))
    assert load_fantom5(enhancers=p).n_enhancers == 3


def test_rejects_empty_bed(tmp_path):
    p = tmp_path / "empty.bed"
    p.write_text("")
    with pytest.raises(ValueError, match="no BED records"):
        load_fantom5(enhancers=p)


def test_rejects_truncated_bed(tmp_path):
    p = tmp_path / "bad.bed"
    p.write_text("chr1\t100\n")
    with pytest.raises(ValueError, match="at least 4 BED columns"):
        load_fantom5(enhancers=p)


# --------------------------------------------------------------------------- #
# Enhancers                                                                    #
# --------------------------------------------------------------------------- #
def test_inside_enhancer(model):
    # BED is 0-based half-open; MEIVA pos is 1-based, so pos 1001 is offset 1000
    ctx = model.annotate("chr1", 1001)
    assert ctx.enhancer_id == "chr1:1000-1300"
    assert ctx.enhancer_distance == 0
    assert ctx.in_enhancer


def test_enhancer_boundaries_are_half_open(model):
    assert model.annotate("chr1", 1300).in_enhancer is True  # offset 1299, last base
    assert model.annotate("chr1", 1301).in_enhancer is False  # offset 1300, past the end


def test_outside_enhancer_reports_distance(model):
    ctx = model.annotate("chr1", 2001)  # offset 2000, enhancer midpoint is 1150
    assert ctx.enhancer_id is None
    assert ctx.enhancer_distance == 850


def test_distance_beyond_cap_is_none(files):
    e, p, _ = files
    m = load_fantom5(enhancers=e, cage_peaks=p, max_distance=100)
    ctx = m.annotate("chr1", 3001)
    assert ctx.enhancer_id is None
    assert ctx.enhancer_distance is None


def test_unknown_contig_is_empty(model):
    ctx = model.annotate("chrUn_gl000220", 5000)
    assert ctx == RegulatoryContext()


# --------------------------------------------------------------------------- #
# CAGE promoters                                                               #
# --------------------------------------------------------------------------- #
def test_inside_peak_carries_gene_and_rank(model):
    ctx = model.annotate("chr1", 8011)  # offset 8010 == the representative TSS
    assert ctx.promoter_id == "hg19::chr1:7000..7020,+;hg_1.1"
    assert ctx.promoter_gene == "GENEA"
    assert ctx.promoter_rank == 1
    assert ctx.promoter_distance == 0
    assert ctx.in_promoter


def test_peak_without_rank(model):
    ctx = model.annotate("chr2", 4016)
    assert ctx.promoter_gene == "GENEC"
    assert ctx.promoter_rank is None


def test_promoter_distance_is_signed_on_plus_strand(model):
    # offset 8015 is 5 bp downstream of the TSS at 8010
    assert model.annotate("chr1", 8016).promoter_distance == 5
    # offset 8005 is 5 bp upstream
    assert model.annotate("chr1", 8006).promoter_distance == -5


def test_promoter_distance_flips_on_minus_strand(model):
    # peak 2 is on the minus strand with its TSS at 9007; a higher coordinate is
    # UPSTREAM in the peak's own orientation, so the sign must invert
    assert model.annotate("chr1", 9013).promoter_distance == -5
    assert model.annotate("chr1", 9003).promoter_distance == 5


def test_near_but_not_in_a_peak(model):
    """The point of reporting distance: peaks are ~14 bp, so near misses matter."""
    ctx = model.annotate("chr1", 8100)  # 90 bp past a 20 bp peak
    assert ctx.promoter_id is None  # not inside
    assert ctx.promoter_gene is None
    assert ctx.promoter_distance == 89  # but the evidence is still reported
    assert ctx.in_promoter is False


def test_peak_names_optional(files):
    e, p, _ = files
    m = load_fantom5(enhancers=e, cage_peaks=p)
    ctx = m.annotate("chr1", 8011)
    assert ctx.promoter_id is not None
    assert ctx.promoter_gene is None  # no name file supplied


def test_enhancers_only(files):
    _, p, _ = files
    m = load_fantom5(cage_peaks=p)
    assert m.n_enhancers == 0
    assert m.annotate("chr1", 1001).enhancer_id is None
    assert m.annotate("chr1", 8011).promoter_id is not None


def test_empty_model_annotates_to_nothing():
    assert Fantom5Model().annotate("chr1", 1000) == RegulatoryContext()
