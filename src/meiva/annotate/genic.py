"""Layer 1 — genic context.

Given a (merged) :class:`~meiva.model.MEISite` and a gene model, work out where
the insertion sits relative to gene structure: coding exon, UTR, non-coding
exon, splice region, intron, promoter, up/downstream, or intergenic — plus the
insertion's orientation *relative to the gene*, which Layer 2's consequence
model depends on, and the distance to the relevant gene.

The classifier is decoupled from any annotation source: it talks to a
:class:`GeneModel` protocol. A GENCODE-backed implementation (and the
reference-data cache that feeds it) is a later step; :class:`InMemoryGeneModel`
here is enough to develop and test the logic.

Conventions:

* All feature coordinates are **0-based half-open**, matching
  :meth:`MEISite.search_interval` — the convention used throughout the codebase.
* "Orientation" is the insertion strand *relative to the gene's* strand:
  ``SENSE`` (same), ``ANTISENSE`` (opposite), or ``UNKNOWN`` if either is
  unknown. A sense vs antisense intronic L1 behaves very differently.
* Up/downstream are defined in *transcription* terms: the promoter/upstream side
  is the TSS side, which is the gene's low coordinate on the + strand and its
  high coordinate on the - strand.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from typing import Protocol

from meiva.model import MEISite, Strand

__all__ = [
    "DEFAULT_DOWNSTREAM_WINDOW",
    "DEFAULT_PROMOTER_WINDOW",
    "DEFAULT_SPLICE_WINDOW",
    "DEFAULT_UPSTREAM_WINDOW",
    "Gene",
    "GeneModel",
    "GenicContext",
    "GenicRegion",
    "InMemoryGeneModel",
    "IndexedGeneModel",
    "InsertionOrientation",
    "Transcript",
    "annotate_genic",
]

#: default splice-region half-width (bp) either side of an exon/intron boundary
DEFAULT_SPLICE_WINDOW = 8
#: default promoter window (bp) upstream of the TSS
DEFAULT_PROMOTER_WINDOW = 1000
#: default upstream window (bp) beyond the promoter, on the TSS side
DEFAULT_UPSTREAM_WINDOW = 5000
#: default downstream window (bp) past the gene 3' end
DEFAULT_DOWNSTREAM_WINDOW = 5000


class GenicRegion(Enum):
    """Where an insertion sits relative to gene structure."""

    CDS = "CDS"
    UTR5 = "5_PRIME_UTR"
    UTR3 = "3_PRIME_UTR"
    EXON_NONCODING = "NONCODING_EXON"
    SPLICE_REGION = "SPLICE_REGION"
    INTRON = "INTRON"
    PROMOTER = "PROMOTER"
    UPSTREAM = "UPSTREAM"
    DOWNSTREAM = "DOWNSTREAM"
    INTERGENIC = "INTERGENIC"


# Severity ordering, used only to pick the primary region when a site relates to
# multiple transcripts/genes. Functional *impact* is Layer 2's job, not this.
_SEVERITY: dict[GenicRegion, int] = {
    GenicRegion.CDS: 10,
    GenicRegion.SPLICE_REGION: 9,
    GenicRegion.UTR5: 8,
    GenicRegion.UTR3: 7,
    GenicRegion.EXON_NONCODING: 6,
    GenicRegion.INTRON: 5,
    GenicRegion.PROMOTER: 4,
    GenicRegion.UPSTREAM: 3,
    GenicRegion.DOWNSTREAM: 2,
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
    is_mane_select: bool = False

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
    biotype: str | None = None  # GENCODE gene_type, e.g. "protein_coding", "lncRNA"


@dataclass(frozen=True, slots=True)
class GenicContext:
    """Result of genic classification for one site."""

    region: GenicRegion
    gene_id: str | None = None
    gene_name: str | None = None
    gene_biotype: str | None = None  # GENCODE gene_type of the reported gene
    gene_strand: Strand | None = None
    transcript_id: str | None = None
    orientation: InsertionOrientation = InsertionOrientation.UNKNOWN
    distance: int | None = None  # bp to gene; 0 if within gene; None if far intergenic
    is_mane_select: bool = False  # whether the reported transcript is MANE Select


class GeneModel(Protocol):
    """Minimal query surface the classifier needs from a gene annotation source."""

    def genes_overlapping(self, chrom: str, start: int, end: int) -> list[Gene]:
        """Genes whose span overlaps the 0-based half-open ``[start, end)``."""
        ...

    def nearby_genes(self, chrom: str, start: int, end: int, max_distance: int) -> list[Gene]:
        """Genes whose span is within ``max_distance`` bp of ``[start, end)``."""
        ...


# --------------------------------------------------------------------------- #
# Geometry helpers                                                            #
# --------------------------------------------------------------------------- #
def _overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start < b_end and b_start < a_end


def _point_distance(start: int, end: int, point: int) -> int:
    if point < start:
        return start - point
    if point >= end:
        return point - end + 1
    return 0


def _gap(s: int, e: int, gene: Gene) -> int:
    """Distance between a site interval and a gene span (0 if overlapping)."""
    if e <= gene.start:
        return gene.start - e
    if s >= gene.end:
        return s - gene.end
    return 0


def _orientation(site_strand: Strand, gene_strand: Strand) -> InsertionOrientation:
    if site_strand is Strand.UNKNOWN or gene_strand is Strand.UNKNOWN:
        return InsertionOrientation.UNKNOWN
    return (
        InsertionOrientation.SENSE if site_strand == gene_strand else InsertionOrientation.ANTISENSE
    )


# --------------------------------------------------------------------------- #
# Within-gene classification                                                  #
# --------------------------------------------------------------------------- #
def _exon_region(tx: Transcript, gene_strand: Strand, s: int, e: int) -> GenicRegion:
    if tx.cds_start is None or tx.cds_end is None:
        return GenicRegion.EXON_NONCODING
    if _overlaps(s, e, tx.cds_start, tx.cds_end):
        return GenicRegion.CDS
    before_cds = e <= tx.cds_start
    if gene_strand is Strand.MINUS:
        return GenicRegion.UTR3 if before_cds else GenicRegion.UTR5
    return GenicRegion.UTR5 if before_cds else GenicRegion.UTR3


def _internal_boundaries(tx: Transcript) -> list[int]:
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
) -> tuple[GenicRegion, str | None, bool] | None:
    """Most severe region across transcripts, with the transcript id and MANE flag.

    Region selection is by severity; ties are broken by preferring a MANE Select
    transcript, then by transcript id for determinism. MANE preference only
    decides *which* transcript is reported, never the region itself.
    """
    candidates: list[tuple[GenicRegion, Transcript]] = []
    for tx in gene.transcripts:
        region = _region_in_transcript(tx, gene.strand, s, e, splice_window)
        if region is not None:
            candidates.append((region, tx))
    if not candidates:
        if gene.transcripts == ():
            return (GenicRegion.INTRON, None, False)
        return None
    candidates.sort(
        key=lambda rt: (-_SEVERITY[rt[0]], not rt[1].is_mane_select, rt[1].transcript_id)
    )
    region, tx = candidates[0]
    return (region, tx.transcript_id, tx.is_mane_select)


# --------------------------------------------------------------------------- #
# Flanking classification (no overlap)                                        #
# --------------------------------------------------------------------------- #
def _flank_region(
    s: int, e: int, gene: Gene, promoter_w: int, upstream_w: int, downstream_w: int
) -> tuple[GenicRegion, int] | None:
    if gene.strand is Strand.UNKNOWN:
        return None  # can't orient promoter/up/down without a strand
    if e <= gene.start:
        gap, left = gene.start - e, True
    elif s >= gene.end:
        gap, left = s - gene.end, False
    else:
        return None  # overlapping; handled elsewhere
    upstream_side = (left and gene.strand is Strand.PLUS) or (
        not left and gene.strand is Strand.MINUS
    )
    if upstream_side:
        if gap <= promoter_w:
            return (GenicRegion.PROMOTER, gap)
        if gap <= upstream_w:
            return (GenicRegion.UPSTREAM, gap)
        return None
    if gap <= downstream_w:
        return (GenicRegion.DOWNSTREAM, gap)
    return None


# --------------------------------------------------------------------------- #
# Entry point                                                                 #
# --------------------------------------------------------------------------- #
def annotate_genic(
    site: MEISite,
    model: GeneModel,
    *,
    splice_window: int = DEFAULT_SPLICE_WINDOW,
    promoter_window: int = DEFAULT_PROMOTER_WINDOW,
    upstream_window: int = DEFAULT_UPSTREAM_WINDOW,
    downstream_window: int = DEFAULT_DOWNSTREAM_WINDOW,
) -> GenicContext:
    """Classify a site against ``model``; ties broken deterministically by gene id."""
    s, e = site.search_interval()

    # 1) genic overlap takes precedence
    best: GenicContext | None = None
    for gene in sorted(model.genes_overlapping(site.chrom, s, e), key=lambda g: g.gene_id):
        result = _classify_within_gene(gene, s, e, splice_window)
        if result is None:
            continue
        region, tx_id, is_mane = result
        if best is None or _SEVERITY[region] > _SEVERITY[best.region]:
            best = GenicContext(
                region=region,
                gene_id=gene.gene_id,
                gene_name=gene.gene_name,
                gene_biotype=gene.biotype,
                gene_strand=gene.strand,
                transcript_id=tx_id,
                orientation=_orientation(site.strand, gene.strand),
                distance=0,
                is_mane_select=is_mane,
            )
    if best is not None:
        return best

    # 2) flanking: promoter / upstream / downstream, else nearest-gene intergenic
    max_w = max(promoter_window, upstream_window, downstream_window)
    best_flank: GenicContext | None = None
    best_key: tuple[int, int] | None = None
    nearest_gene: Gene | None = None
    nearest_gap = 0
    for gene in sorted(model.nearby_genes(site.chrom, s, e, max_w), key=lambda g: g.gene_id):
        gap = _gap(s, e, gene)
        if nearest_gene is None or gap < nearest_gap:
            nearest_gene, nearest_gap = gene, gap
        rel = _flank_region(s, e, gene, promoter_window, upstream_window, downstream_window)
        if rel is None:
            continue
        region, dist = rel
        key = (_SEVERITY[region], -dist)
        if best_key is None or key > best_key:
            best_key, best_flank = (
                key,
                GenicContext(
                    region=region,
                    gene_id=gene.gene_id,
                    gene_name=gene.gene_name,
                    gene_biotype=gene.biotype,
                    gene_strand=gene.strand,
                    orientation=_orientation(site.strand, gene.strand),
                    distance=dist,
                ),
            )
    if best_flank is not None:
        return best_flank
    if nearest_gene is not None:
        return GenicContext(
            region=GenicRegion.INTERGENIC,
            gene_id=nearest_gene.gene_id,
            gene_name=nearest_gene.gene_name,
            gene_biotype=nearest_gene.biotype,
            gene_strand=nearest_gene.strand,
            orientation=_orientation(site.strand, nearest_gene.strand),
            distance=nearest_gap,
        )
    return GenicContext(region=GenicRegion.INTERGENIC)


class InMemoryGeneModel:
    """A simple :class:`GeneModel` backed by an in-memory list of genes.

    Adequate for tests and small inputs (per-contig linear scan). The
    GENCODE-backed model will use an interval index.
    """

    def __init__(self, genes: Iterable[Gene]) -> None:
        self._by_chrom: dict[str, list[Gene]] = {}
        for gene in genes:
            self._by_chrom.setdefault(gene.chrom, []).append(gene)

    def genes_overlapping(self, chrom: str, start: int, end: int) -> list[Gene]:
        return [g for g in self._by_chrom.get(chrom, []) if _overlaps(start, end, g.start, g.end)]

    def nearby_genes(self, chrom: str, start: int, end: int, max_distance: int) -> list[Gene]:
        return [g for g in self._by_chrom.get(chrom, []) if _gap(start, end, g) <= max_distance]


class IndexedGeneModel:
    """A :class:`GeneModel` with a binned index, for genome-scale gene sets.

    Each gene is placed in every fixed-width bin its span covers; a query scans
    only the bins overlapping the (optionally padded) query interval, then
    filters precisely. Dependency-free and fast enough for whole-genome GENCODE
    against a cohort's worth of sites; genes rarely span more than a few bins.
    """

    def __init__(self, genes: Iterable[Gene], *, bin_size: int = 1_000_000) -> None:
        if bin_size <= 0:
            raise ValueError("bin_size must be > 0")
        self._bin_size = bin_size
        self._bins: dict[tuple[str, int], list[Gene]] = {}
        for gene in genes:
            last = max(gene.end - 1, gene.start)
            for b in range(gene.start // bin_size, last // bin_size + 1):
                self._bins.setdefault((gene.chrom, b), []).append(gene)

    def _candidates(self, chrom: str, start: int, end: int) -> list[Gene]:
        bs = self._bin_size
        seen: set[int] = set()
        out: list[Gene] = []
        last = max(end - 1, start)
        for b in range(start // bs, last // bs + 1):
            for gene in self._bins.get((chrom, b), ()):
                if id(gene) not in seen:
                    seen.add(id(gene))
                    out.append(gene)
        return out

    def genes_overlapping(self, chrom: str, start: int, end: int) -> list[Gene]:
        return [
            g for g in self._candidates(chrom, start, end) if _overlaps(start, end, g.start, g.end)
        ]

    def nearby_genes(self, chrom: str, start: int, end: int, max_distance: int) -> list[Gene]:
        qs, qe = max(start - max_distance, 0), end + max_distance
        return [g for g in self._candidates(chrom, qs, qe) if _gap(start, end, g) <= max_distance]
