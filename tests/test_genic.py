"""Tests for Layer 1 genic classification, using a synthetic gene model."""

from meiva.annotate import (
    Gene,
    GenicRegion,
    InMemoryGeneModel,
    InsertionOrientation,
    Transcript,
    annotate_genic,
)
from meiva.model import MEIFamily, MEISite, Strand

# A coding gene on the + strand. All coords 0-based half-open.
#   exons:    [1000,1200) [1500,1700) [2000,2200)
#   CDS span: [1100,2100)
#   => 5'UTR 1000-1100 (in exon1), CDS 1100-2100, 3'UTR 2100-2200 (in exon3)
PLUS_GENE = Gene(
    gene_id="ENSG_PLUS",
    gene_name="PLUSGENE",
    chrom="chr1",
    start=1000,
    end=2200,
    strand=Strand.PLUS,
    transcripts=(
        Transcript(
            transcript_id="ENST_PLUS",
            exons=((1000, 1200), (1500, 1700), (2000, 2200)),
            cds_start=1100,
            cds_end=2100,
        ),
    ),
)

# Same structure, minus strand: UTR sides flip.
MINUS_GENE = Gene(
    gene_id="ENSG_MINUS",
    gene_name="MINUSGENE",
    chrom="chr2",
    start=1000,
    end=2200,
    strand=Strand.MINUS,
    transcripts=(
        Transcript(
            transcript_id="ENST_MINUS",
            exons=((1000, 1200), (1500, 1700), (2000, 2200)),
            cds_start=1100,
            cds_end=2100,
        ),
    ),
)

NONCODING_GENE = Gene(
    gene_id="ENSG_NC",
    gene_name="LINCSOMETHING",
    chrom="chr3",
    start=1000,
    end=1200,
    strand=Strand.PLUS,
    transcripts=(Transcript(transcript_id="ENST_NC", exons=((1000, 1200),)),),
)

MODEL = InMemoryGeneModel([PLUS_GENE, MINUS_GENE, NONCODING_GENE])


def _site(chrom: str, pos_1based: int, strand: Strand = Strand.UNKNOWN) -> MEISite:
    return MEISite(chrom=chrom, pos=pos_1based, family=MEIFamily.ALU, strand=strand)


def _region(chrom: str, pos: int, strand: Strand = Strand.UNKNOWN) -> GenicRegion:
    return annotate_genic(_site(chrom, pos, strand), MODEL).region


# --------------------------------------------------------------------------- #
# Plus-strand coding gene                                                     #
# --------------------------------------------------------------------------- #
def test_cds():
    assert _region("chr1", 1151) is GenicRegion.CDS  # exon1, within CDS


def test_five_prime_utr_plus():
    assert _region("chr1", 1051) is GenicRegion.UTR5  # exon1, before CDS start


def test_three_prime_utr_plus():
    assert _region("chr1", 2151) is GenicRegion.UTR3  # exon3, after CDS end


def test_intron():
    assert _region("chr1", 1351) is GenicRegion.INTRON  # between exon1 and exon2


def test_splice_region():
    # 1205 (1-based) -> interval [1204,1205); 4 bp from the exon1/intron boundary at 1200
    assert _region("chr1", 1205) is GenicRegion.SPLICE_REGION


def test_intergenic():
    assert _region("chr1", 50_000) is GenicRegion.INTERGENIC


# --------------------------------------------------------------------------- #
# Minus-strand gene: UTR sides flip                                           #
# --------------------------------------------------------------------------- #
def test_utr_sides_flip_on_minus_strand():
    assert _region("chr2", 1051) is GenicRegion.UTR3  # genomic-left of CDS, minus => 3'
    assert _region("chr2", 2151) is GenicRegion.UTR5  # genomic-right of CDS, minus => 5'


# --------------------------------------------------------------------------- #
# Non-coding exon                                                             #
# --------------------------------------------------------------------------- #
def test_noncoding_exon():
    assert _region("chr3", 1100) is GenicRegion.EXON_NONCODING


# --------------------------------------------------------------------------- #
# Orientation relative to the gene                                            #
# --------------------------------------------------------------------------- #
def test_orientation_sense():
    ctx = annotate_genic(_site("chr1", 1151, Strand.PLUS), MODEL)
    assert ctx.orientation is InsertionOrientation.SENSE


def test_orientation_antisense():
    ctx = annotate_genic(_site("chr1", 1151, Strand.MINUS), MODEL)
    assert ctx.orientation is InsertionOrientation.ANTISENSE


def test_orientation_unknown_when_site_strand_unknown():
    ctx = annotate_genic(_site("chr1", 1151, Strand.UNKNOWN), MODEL)
    assert ctx.orientation is InsertionOrientation.UNKNOWN


# --------------------------------------------------------------------------- #
# Multi-gene precedence: most severe region wins                              #
# --------------------------------------------------------------------------- #
def test_most_severe_region_across_overlapping_genes():
    # An intronic-only gene overlapping the same locus as PLUS_GENE's CDS.
    intronic_gene = Gene(
        gene_id="ENSG_INTRONIC",
        gene_name="INTRONIC",
        chrom="chr1",
        start=900,
        end=2300,
        strand=Strand.PLUS,
        transcripts=(Transcript(transcript_id="ENST_I", exons=((900, 1000), (2200, 2300))),),
    )
    model = InMemoryGeneModel([PLUS_GENE, intronic_gene])
    ctx = annotate_genic(_site("chr1", 1151), model)
    assert ctx.region is GenicRegion.CDS
    assert ctx.gene_id == "ENSG_PLUS"


# --------------------------------------------------------------------------- #
# Flanking: promoter / upstream / downstream + distance                       #
# --------------------------------------------------------------------------- #
def test_distance_zero_when_genic():
    ctx = annotate_genic(_site("chr1", 1151), MODEL)
    assert ctx.region is GenicRegion.CDS
    assert ctx.distance == 0


def test_promoter_then_upstream_plus():
    # PLUS_GENE TSS is gene.start (1000). Use explicit small windows to separate.
    near = annotate_genic(_site("chr1", 960), MODEL, promoter_window=100, upstream_window=1000)
    far = annotate_genic(_site("chr1", 700), MODEL, promoter_window=100, upstream_window=1000)
    assert near.region is GenicRegion.PROMOTER
    assert near.distance == 40
    assert far.region is GenicRegion.UPSTREAM
    assert far.distance == 300


def test_downstream_plus():
    ctx = annotate_genic(_site("chr1", 2500), MODEL)  # past gene.end (2200)
    assert ctx.region is GenicRegion.DOWNSTREAM
    assert ctx.distance == 299


def test_flank_sides_flip_on_minus():
    # MINUS_GENE (chr2): TSS at gene.end (2200), so the high-coordinate side is upstream.
    up = annotate_genic(_site("chr2", 2260), MODEL, promoter_window=100)
    down = annotate_genic(_site("chr2", 700), MODEL)
    assert up.region is GenicRegion.PROMOTER
    assert down.region is GenicRegion.DOWNSTREAM


def test_intergenic_far_has_no_gene():
    ctx = annotate_genic(_site("chr1", 500_000), MODEL)
    assert ctx.region is GenicRegion.INTERGENIC
    assert ctx.gene_id is None
    assert ctx.distance is None
