"""GENCODE GTF loader.

Parses a GENCODE GTF (GRCh38, e.g. release 46) into the :class:`~meiva.annotate.genic.Gene`
model the classifier consumes. We keep only what genic classification needs:
gene span/strand, transcript exons, the CDS span per transcript (min..max over
CDS features), and the MANE Select tag. Explicit UTR features are ignored —
the classifier derives 5'/3' UTR from exon-minus-CDS using gene strand.

Coordinate handling: GTF is 1-based inclusive; we convert to the codebase-wide
0-based half-open convention (``start - 1`` .. ``end``).

The file is large (~1.5 GB uncompressed for comprehensive GENCODE), so parsing
streams line by line and accepts a gzipped or plain path.
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

from meiva.annotate.genic import Gene, IndexedGeneModel, Transcript
from meiva.model import Strand

__all__ = ["load_gencode", "parse_gencode_gtf"]

_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_GENE_NAME = re.compile(r'gene_name "([^"]+)"')
_GENE_TYPE = re.compile(r'gene_type "([^"]+)"')
_TX_ID = re.compile(r'transcript_id "([^"]+)"')
_MANE_SELECT = re.compile(r'tag "MANE_Select"')

_STRAND = {"+": Strand.PLUS, "-": Strand.MINUS}


@dataclass
class _TxBuilder:
    gene_id: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds_start: int | None = None
    cds_end: int | None = None
    is_mane: bool = False

    def add_cds(self, start: int, end: int) -> None:
        self.cds_start = start if self.cds_start is None else min(self.cds_start, start)
        self.cds_end = end if self.cds_end is None else max(self.cds_end, end)


@dataclass
class _GeneBuilder:
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: Strand
    biotype: str | None = None


def _open(path: str | Path) -> Iterator[str]:
    p = Path(path)
    if p.suffix == ".gz":
        with gzip.open(p, "rt") as fh:
            yield from fh
    else:
        with open(p) as fh:
            yield from fh


def _search(pattern: re.Pattern[str], text: str) -> str | None:
    m = pattern.search(text)
    return m.group(1) if m else None


def parse_gencode_gtf(path: str | Path) -> list[Gene]:
    """Parse a GENCODE GTF into a list of :class:`Gene` objects."""
    genes: dict[str, _GeneBuilder] = {}
    txs: dict[str, _TxBuilder] = {}

    for line in _open(path):
        if not line or line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) != 9:
            continue
        chrom, _src, feature, start_s, end_s, _score, strand_s, _frame, attrs = cols
        # GTF 1-based inclusive -> 0-based half-open
        start, end = int(start_s) - 1, int(end_s)
        strand = _STRAND.get(strand_s, Strand.UNKNOWN)

        if feature == "gene":
            gene_id = _search(_GENE_ID, attrs)
            if gene_id is None:
                continue
            genes[gene_id] = _GeneBuilder(
                gene_name=_search(_GENE_NAME, attrs) or gene_id,
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
                biotype=_search(_GENE_TYPE, attrs),
            )
        elif feature == "transcript":
            tx_id = _search(_TX_ID, attrs)
            gene_id = _search(_GENE_ID, attrs)
            if tx_id is None or gene_id is None:
                continue
            txs[tx_id] = _TxBuilder(gene_id=gene_id, is_mane=bool(_MANE_SELECT.search(attrs)))
        elif feature == "exon":
            tx_id = _search(_TX_ID, attrs)
            if tx_id is not None and tx_id in txs:
                txs[tx_id].exons.append((start, end))
        elif feature == "CDS":
            tx_id = _search(_TX_ID, attrs)
            if tx_id is not None and tx_id in txs:
                txs[tx_id].add_cds(start, end)

    # assemble transcripts under their genes
    by_gene: dict[str, list[Transcript]] = {}
    for tx_id, tb in txs.items():
        if not tb.exons:
            continue  # a transcript with no exons can't be placed
        by_gene.setdefault(tb.gene_id, []).append(
            Transcript(
                transcript_id=tx_id,
                exons=tuple(sorted(tb.exons)),
                cds_start=tb.cds_start,
                cds_end=tb.cds_end,
                is_mane_select=tb.is_mane,
            )
        )

    out: list[Gene] = []
    for gene_id, gb in genes.items():
        out.append(
            Gene(
                gene_id=gene_id,
                gene_name=gb.gene_name,
                chrom=gb.chrom,
                start=gb.start,
                end=gb.end,
                strand=gb.strand,
                transcripts=tuple(by_gene.get(gene_id, ())),
                biotype=gb.biotype,
            )
        )
    return out


def load_gencode(path: str | Path, *, bin_size: int = 1_000_000) -> IndexedGeneModel:
    """Parse a GENCODE GTF and return an indexed gene model ready for queries."""
    return IndexedGeneModel(parse_gencode_gtf(path), bin_size=bin_size)
