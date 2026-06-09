"""Layer 1 — genic context.

Given a (merged) :class:`~meiva.model.MEISite` and a gene model, work out where
the insertion sits relative to gene structure: coding exon, UTR, non-coding
exon, splice region, intron, or intergenic — plus the insertion's orientation
*relative to the gene*, which Layer 2's consequence model depends on.

This module is deliberately decoupled from any particular annotation source:
the classifier talks to a :class:`GeneModel` protocol. A GENCODE-backed
implementation (and the reference-data cache that feeds it) is a later step;
:class:`InMemoryGeneModel` here is enough to develop and test the logic.

Conventions:

* All feature coordinates are **0-based half-open**, matching
  :meth:`MEISite.search_interval`. This is the same convention used throughout
  the codebase; keeping the gene model in it avoids re-deriving off-by-ones.
* "Orientation" is the insertion strand *relative to the gene's* strand:
  ``SENSE`` (same), ``ANTISENSE`` (opposite), or ``UNKNOWN`` if either strand
  is unknown. A sense vs antisense intronic L1 behaves very differently, so
  this distinction is first-class.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from typing import Protocol

from meiva.model import MEISite, Strand

__all__ = [
    "Gene",
    "GeneModel",
    "GenicContext",
    "GenicRegion",
    "InMemoryGeneModel",
    "InsertionOrientation",
    "Transcript",
    "annotate_genic",
]

#: default splice-region half-width (bp) either side of an exon/intron boundary
DEFAULT_SPLICE_WINDOW = 8


class GenicRegion(Enum):
    """Where an insertion sits relative to gene structure.

    Intergenic sub-classes (promoter / upstream / downstream) and the
    nearest-gene distance are added in the next step, alongside the GENCODE
    loader; only the members the classifier can currently assign are defined.
    """

    CDS = "CDS"
    UTR5 = "5_PRIME_UTR"
    UTR3 = "3_PRIME_UTR"
    EXON_NONCODING = "NONCODING_EXON"
    SPLICE_REGION = "SPLICE_REGION"
    INTRON = "INTRON"
    INTERGENIC = "INTERGENIC"


# Severity ordering, used only to pick the primary region when a site overlaps
# multiple transcripts/genes. Functional *impact* is Layer 2's job, not this.
_SEVERITY: dict[GenicRegion, int] = {
    GenicRegion.CDS: 8,
    GenicRegion.SPLICE_REGION: 7,
    GenicRegion.UTR5: 6,
    GenicRegion.UTR3: 5,
    GenicRegion.EXON_NONCODING: 4,
    GenicRegion.INTRON: 3,
    GenicRegion.INTERGENIC: 0,
}


class InsertionOrientation(Enum):
    SENSE = "SENSE"
    ANTISENSE = "ANTISENSE"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True, slots=True)
class Transcript:
    """A transcript as exon intervals plus an optional CDS span (genomic, 0-based)."""

    transcript_id: str
    exons: tuple[tuple[int, int], ...]  # sorted, 0-based half-open
    cds_start: int | None = None  # None => non-coding transcript
    cds_end: int | None = None

    @property
    def span(self) -> tuple[int, int]:
        return (self.exons[0][0], self.exons[-1][1])


@dataclass(frozen=True, slots=True)
class Gene:
    """A gene: identity, span, strand, and its transcripts."""

    gene_id: str
    gene_name: str
    chrom: str
    start: int  # 0-based half-open span of the gene
    end: int
    strand: Strand
    transcripts: tuple[Transcript, ...] = ()


@dataclass(frozen=True, slots=True)
class GenicContext:
    """Result of genic classification for one site."""

    region: GenicRegion
    gene_id: str | None = None
    gene_name: str | None = None
    gene_strand: Strand | None = None
    transcript_id: str | None = None
    orientation: InsertionOrientation = InsertionOrientation.UNKNOWN


class GeneModel(Protocol):
    """Minimal query surface the classifier needs from a gene annotation source."""

    def genes_overlapping(self, chrom: str, start: int, end: int) -> list[Gene]:
        """Return all genes whose span overlaps the 0-based half-open ``[start, end)``."""
        ...


# --------------------------------------------------------------------------- #
# Geometry helpers                                                            #
# --------------------------------------------------------------------------- #
def _overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start < b_end and b_start < a_end


def _point_distance(start: int, end: int, point: int) -> int:
    """Distance from a 0-based half-open interval to a boundary coordinate."""
    if point < start:
        return start - point
    if point >= end:
        return point - end + 1
    return 0


def _orientation(site_strand: Strand, gene_strand: Strand) -> InsertionOrientation:
    if site_strand is Strand.UNKNOWN or gene_strand is Strand.UNKNOWN:
        return InsertionOrientation.UNKNOWN
    return (
        InsertionOrientation.SENSE if site_strand == gene_strand else InsertionOrientation.ANTISENSE
    )


# --------------------------------------------------------------------------- #
# Classification                                                              #
# --------------------------------------------------------------------------- #
def _exon_region(tx: Transcript, gene_strand: Strand, s: int, e: int) -> GenicRegion:
    """Classify an exonic hit as CDS, UTR (strand-aware), or non-coding exon."""
    if tx.cds_start is None or tx.cds_end is None:
        return GenicRegion.EXON_NONCODING
    if _overlaps(s, e, tx.cds_start, tx.cds_end):
        return GenicRegion.CDS
    # Outside the CDS but in an exon => UTR. Genomic side maps to 5'/3' by strand.
    before_cds = e <= tx.cds_start
    if gene_strand is Strand.MINUS:
        return GenicRegion.UTR3 if before_cds else GenicRegion.UTR5
    return GenicRegion.UTR5 if before_cds else GenicRegion.UTR3


def _internal_boundaries(tx: Transcript) -> list[int]:
    """Exon/intron boundary coordinates (excludes the transcript's outer ends)."""
    tx_start, tx_end = tx.span
    boundaries: list[int] = []
    for es, ee in tx.exons:
        if es != tx_start:
            boundaries.append(es)
        if ee != tx_end:
            boundaries.append(ee)
    return boundaries


def _region_in_transcript(
    tx: Transcript, gene_strand: Strand, s: int, e: int, splice_window: int
) -> GenicRegion | None:
    for es, ee in tx.exons:
        if _overlaps(s, e, es, ee):
            return _exon_region(tx, gene_strand, s, e)
    for boundary in _internal_boundaries(tx):
        if _point_distance(s, e, boundary) <= splice_window:
            return GenicRegion.SPLICE_REGION
    tx_start, tx_end = tx.span
    if _overlaps(s, e, tx_start, tx_end):
        return GenicRegion.INTRON
    return None


def _classify_within_gene(
    gene: Gene, s: int, e: int, splice_window: int
) -> tuple[GenicRegion, str | None] | None:
    """Most severe region across the gene's transcripts, with the transcript id."""
    best: tuple[GenicRegion, str | None] | None = None
    for tx in gene.transcripts:
        region = _region_in_transcript(tx, gene.strand, s, e, splice_window)
        if region is None:
            continue
        if best is None or _SEVERITY[region] > _SEVERITY[best[0]]:
            best = (region, tx.transcript_id)
    # A gene with no transcript detail still counts as intronic when overlapped.
    if best is None and gene.transcripts == ():
        return (GenicRegion.INTRON, None)
    return best


def annotate_genic(
    site: MEISite, model: GeneModel, *, splice_window: int = DEFAULT_SPLICE_WINDOW
) -> GenicContext:
    """Classify a site against ``model``.

    Returns the most severe genic context across all overlapping genes, or
    ``INTERGENIC`` when nothing overlaps. Ties are broken deterministically by
    gene id.
    """
    s, e = site.search_interval()
    genes = model.genes_overlapping(site.chrom, s, e)

    best: GenicContext | None = None
    for gene in sorted(genes, key=lambda g: g.gene_id):
        result = _classify_within_gene(gene, s, e, splice_window)
        if result is None:
            continue
        region, tx_id = result
        if best is None or _SEVERITY[region] > _SEVERITY[best.region]:
            best = GenicContext(
                region=region,
                gene_id=gene.gene_id,
                gene_name=gene.gene_name,
                gene_strand=gene.strand,
                transcript_id=tx_id,
                orientation=_orientation(site.strand, gene.strand),
            )
    return best if best is not None else GenicContext(region=GenicRegion.INTERGENIC)


class InMemoryGeneModel:
    """A simple :class:`GeneModel` backed by an in-memory list of genes.

    Adequate for tests and small inputs. The GENCODE-backed model (next step)
    will use an interval index; this one does a per-contig linear scan.
    """

    def __init__(self, genes: Iterable[Gene]) -> None:
        self._by_chrom: dict[str, list[Gene]] = {}
        for gene in genes:
            self._by_chrom.setdefault(gene.chrom, []).append(gene)

    def genes_overlapping(self, chrom: str, start: int, end: int) -> list[Gene]:
        return [g for g in self._by_chrom.get(chrom, []) if _overlaps(start, end, g.start, g.end)]
