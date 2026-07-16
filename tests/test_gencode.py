"""Tests for the GENCODE GTF loader and the indexed gene model.

Uses a tiny synthetic GTF in real GENCODE format (1-based, tab-separated,
key "value" attributes). The same gene geometry mirrors the synthetic gene used
in test_genic so the classifier behaviour is directly comparable.
"""

import gzip

import pytest

from meiva.annotate import (
    IndexedGeneModel,
    InMemoryGeneModel,
    annotate_genic,
    load_gencode,
    parse_gencode_gtf,
)
from meiva.annotate.genic import GenicRegion
from meiva.model import MEIFamily, MEISite, Strand

# 1-based inclusive coords. ENSGTEST: exons 1001-1200, 1501-1700, 2001-2200;
# CDS 1101-1200, 1501-1700, 2001-2100 => CDS span 1101..2100. MANE transcript.
# Plus a second, non-coding transcript ENSTN.
GTF = "\t".join
_LINES = [
    "#description: synthetic",
    GTF(
        [
            "chr1",
            "HAVANA",
            "gene",
            "1001",
            "2200",
            ".",
            "+",
            ".",
            'gene_id "ENSGTEST.1"; gene_name "TESTGENE"; gene_type "protein_coding";',
        ]
    ),
    GTF(
        [
            "chr1",
            "HAVANA",
            "transcript",
            "1001",
            "2200",
            ".",
            "+",
            ".",
            'gene_id "ENSGTEST.1"; transcript_id "ENSTM.1"; tag "MANE_Select";',
        ]
    ),
    GTF(
        [
            "chr1",
            "HAVANA",
            "exon",
            "1001",
            "1200",
            ".",
            "+",
            ".",
            'gene_id "ENSGTEST.1"; transcript_id "ENSTM.1";',
        ]
    ),
    GTF(
        [
            "chr1",
            "HAVANA",
            "exon",
            "1501",
            "1700",
            ".",
            "+",
            ".",
            'gene_id "ENSGTEST.1"; transcript_id "ENSTM.1";',
        ]
    ),
    GTF(
        [
            "chr1",
            "HAVANA",
            "exon",
            "2001",
            "2200",
            ".",
            "+",
            ".",
            'gene_id "ENSGTEST.1"; transcript_id "ENSTM.1";',
        ]
    ),
    GTF(
        [
            "chr1",
            "HAVANA",
            "CDS",
            "1101",
            "1200",
            ".",
            "+",
            "0",
            'gene_id "ENSGTEST.1"; transcript_id "ENSTM.1";',
        ]
    ),
    GTF(
        [
            "chr1",
            "HAVANA",
            "CDS",
            "1501",
            "1700",
            ".",
            "+",
            "0",
            'gene_id "ENSGTEST.1"; transcript_id "ENSTM.1";',
        ]
    ),
    GTF(
        [
            "chr1",
            "HAVANA",
            "CDS",
            "2001",
            "2100",
            ".",
            "+",
            "0",
            'gene_id "ENSGTEST.1"; transcript_id "ENSTM.1";',
        ]
    ),
    GTF(
        [
            "chr1",
            "HAVANA",
            "transcript",
            "1001",
            "1700",
            ".",
            "+",
            ".",
            'gene_id "ENSGTEST.1"; transcript_id "ENSTN.1";',
        ]
    ),
    GTF(
        [
            "chr1",
            "HAVANA",
            "exon",
            "1001",
            "1700",
            ".",
            "+",
            ".",
            'gene_id "ENSGTEST.1"; transcript_id "ENSTN.1";',
        ]
    ),
    # a second gene on the minus strand, far away
    GTF(
        [
            "chr1",
            "HAVANA",
            "gene",
            "50001",
            "60000",
            ".",
            "-",
            ".",
            'gene_id "ENSG2.1"; gene_name "GENE2"; gene_type "lncRNA";',
        ]
    ),
    GTF(
        [
            "chr1",
            "HAVANA",
            "transcript",
            "50001",
            "60000",
            ".",
            "-",
            ".",
            'gene_id "ENSG2.1"; transcript_id "ENST2.1";',
        ]
    ),
    GTF(
        [
            "chr1",
            "HAVANA",
            "exon",
            "50001",
            "60000",
            ".",
            "-",
            ".",
            'gene_id "ENSG2.1"; transcript_id "ENST2.1";',
        ]
    ),
]


@pytest.fixture
def gtf_file(tmp_path):
    p = tmp_path / "synthetic.gtf"
    p.write_text("\n".join(_LINES) + "\n")
    return p


# --------------------------------------------------------------------------- #
# Parsing                                                                     #
# --------------------------------------------------------------------------- #
def test_parses_genes_and_transcripts(gtf_file):
    genes = {g.gene_id: g for g in parse_gencode_gtf(gtf_file)}
    assert set(genes) == {"ENSGTEST.1", "ENSG2.1"}
    g = genes["ENSGTEST.1"]
    assert g.gene_name == "TESTGENE"
    assert g.chrom == "chr1"
    assert g.strand is Strand.PLUS
    assert g.biotype == "protein_coding"
    assert genes["ENSG2.1"].biotype == "lncRNA"
    assert (g.start, g.end) == (1000, 2200)  # 1-based 1001..2200 -> 0-based half-open
    assert len(g.transcripts) == 2


def test_coding_transcript_cds_and_mane(gtf_file):
    g = next(x for x in parse_gencode_gtf(gtf_file) if x.gene_id == "ENSGTEST.1")
    mane = next(t for t in g.transcripts if t.transcript_id == "ENSTM.1")
    assert mane.exons == ((1000, 1200), (1500, 1700), (2000, 2200))
    assert mane.cds_start == 1100 and mane.cds_end == 2100  # min/max over CDS features
    assert mane.is_mane_select is True


def test_noncoding_transcript(gtf_file):
    g = next(x for x in parse_gencode_gtf(gtf_file) if x.gene_id == "ENSGTEST.1")
    nc = next(t for t in g.transcripts if t.transcript_id == "ENSTN.1")
    assert nc.cds_start is None and nc.cds_end is None
    assert nc.is_mane_select is False


def test_minus_strand_gene(gtf_file):
    g = next(x for x in parse_gencode_gtf(gtf_file) if x.gene_id == "ENSG2.1")
    assert g.strand is Strand.MINUS


def test_handles_gzipped_input(tmp_path):
    p = tmp_path / "synthetic.gtf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("\n".join(_LINES) + "\n")
    genes = parse_gencode_gtf(p)
    assert {g.gene_id for g in genes} == {"ENSGTEST.1", "ENSG2.1"}


# --------------------------------------------------------------------------- #
# Classification through the loaded model + MANE reporting                     #
# --------------------------------------------------------------------------- #
def test_load_gencode_and_classify(gtf_file):
    model = load_gencode(gtf_file)
    ctx = annotate_genic(MEISite(chrom="chr1", pos=1151, family=MEIFamily.ALU), model)
    assert ctx.region is GenicRegion.CDS
    assert ctx.gene_id == "ENSGTEST.1"
    assert ctx.transcript_id == "ENSTM.1"
    assert ctx.is_mane_select is True


# --------------------------------------------------------------------------- #
# IndexedGeneModel matches the reference InMemory implementation               #
# --------------------------------------------------------------------------- #
def test_indexed_model_matches_in_memory(gtf_file):
    genes = parse_gencode_gtf(gtf_file)
    indexed = IndexedGeneModel(genes, bin_size=1000)  # small bins to exercise multi-bin spans
    in_mem = InMemoryGeneModel(genes)
    for start, end in [(1149, 1150), (0, 1), (55_000, 55_001), (1500, 2300), (59_999, 70_000)]:
        a = {g.gene_id for g in indexed.genes_overlapping("chr1", start, end)}
        b = {g.gene_id for g in in_mem.genes_overlapping("chr1", start, end)}
        assert a == b
        a2 = {g.gene_id for g in indexed.nearby_genes("chr1", start, end, 5000)}
        b2 = {g.gene_id for g in in_mem.nearby_genes("chr1", start, end, 5000)}
        assert a2 == b2


def test_indexed_bin_size_validation():
    with pytest.raises(ValueError):
        IndexedGeneModel([], bin_size=0)
