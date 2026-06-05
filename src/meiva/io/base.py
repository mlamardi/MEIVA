"""Ingest base layer: the parser contract and shared, caller-agnostic plumbing.

Design split (mirrors the model's): **parsers own normalisation**, the model
only validates. Everything in this module that touches a caller's quirks lives
*here or in the concrete subclass*, never in :mod:`meiva.model`.

Three tiers:

* :class:`MEIParser` -- the abstract contract every caller parser satisfies
  (``parse`` + ``sniff``). Pure; imports nothing heavy.
* :class:`Cyvcf2Parser` -- a concrete base that handles opening a VCF with
  cyvcf2, iterating records, and extracting the *VCF-standard* bits (FILTER,
  GT, raw INFO). Subclasses implement only the caller-specific field mapping in
  ``_build_site``.
* Free helpers (:func:`normalize_contig`, :func:`plausible_length`, …) shared
  across callers and unit-testable in isolation.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import cyvcf2

from meiva.model import MEISite, SampleGenotype

__all__ = [
    "MAX_PLAUSIBLE_INSERTION_BP",
    "Cyvcf2Parser",
    "MEIParser",
    "genotypes_from_cyvcf2",
    "normalize_contig",
    "plausible_length",
    "raw_info_from_cyvcf2",
]

# Comfortably above a full-length HERV-K (~9.5 kb) and any Alu/L1/SVA, but far
# below the 20-35 Mb coordinate values xTEA leaks into SVLEN for orphan
# transductions. Anything above this is not a real inserted-element length.
MAX_PLAUSIBLE_INSERTION_BP = 50_000


def normalize_contig(contig: str) -> str:
    """Normalise a contig name to GRCh38 'chr'-prefixed space.

    MEIVA works in 'chr' space (GENCODE / gnomAD-SV v4). This is the single
    place that decision is enforced, so a GRCh37-style ``"6"`` or ``"MT"``
    becomes ``"chr6"`` / ``"chrM"`` here rather than silently failing to
    overlap anything later. It is intentionally conservative: it does not try
    to remap alt/decoy contigs, only to add the prefix and fix the
    mitochondrion alias.
    """
    c = contig.strip()
    if not c:
        raise ValueError("empty contig name")
    if c.startswith("chr"):
        return c
    if c in {"MT", "M"}:
        return "chrM"
    return "chr" + c


def plausible_length(svlen: int | None) -> tuple[int | None, bool]:
    """Validate a caller-reported insertion length.

    Returns ``(length, unreliable)``. ``length`` is ``None`` when the value
    cannot be a real inserted-element length -- chiefly the transduction case
    where xTEA writes a genomic coordinate (tens of Mb) into ``SVLEN``. The
    caller should preserve the raw value and flag the site rather than drop it.
    """
    if svlen is None:
        return None, False
    if svlen <= 0 or svlen > MAX_PLAUSIBLE_INSERTION_BP:
        return None, True
    return svlen, False


def _stringify(value: Any) -> str:
    """Render a cyvcf2 INFO value as a string for lossless ``raw_info`` storage."""
    if isinstance(value, bool):  # INFO flags; must precede the int check
        return "1" if value else "0"
    if isinstance(value, (list, tuple)):
        return ",".join(str(x) for x in value)
    return str(value)


def raw_info_from_cyvcf2(variant: Any) -> dict[str, str]:
    """Capture the complete INFO field as strings (nothing dropped)."""
    return {key: _stringify(value) for key, value in variant.INFO}


def genotypes_from_cyvcf2(variant: Any, samples: list[str]) -> dict[str, SampleGenotype]:
    """Build per-sample genotypes from cyvcf2's representation.

    cyvcf2 yields one ``[allele_1, ..., phased]`` row per sample, using ``-1``
    for a missing allele. We translate ``-1`` to ``None`` (matching the model's
    convention) and read the trailing element as the phase flag.
    """
    out: dict[str, SampleGenotype] = {}
    gts = variant.genotypes
    if gts is None:
        return out
    for sample, row in zip(samples, gts, strict=True):
        *alleles, phased = row
        norm = tuple(None if a is None or a < 0 else int(a) for a in alleles)
        out[sample] = SampleGenotype(sample=sample, alleles=norm, phased=bool(phased))
    return out


class MEIParser(ABC):
    """Abstract contract for a caller-specific MEI parser."""

    #: short, lowercase caller identifier recorded on each site (e.g. "xtea")
    caller: str = "unknown"

    @abstractmethod
    def parse(self, path: str | Path) -> Iterator[MEISite]:
        """Yield normalised :class:`MEISite` objects from ``path``."""
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def sniff(cls, path: str | Path) -> bool:
        """Return True if this parser recognises ``path`` from its header."""
        raise NotImplementedError


class Cyvcf2Parser(MEIParser):
    """Concrete base handling cyvcf2 iteration and VCF-standard extraction.

    Subclasses implement :meth:`_build_site` (caller-specific INFO mapping) and
    :meth:`sniff` (header detection). Everything generic lives here.
    """

    def parse(self, path: str | Path) -> Iterator[MEISite]:
        vcf = cyvcf2.VCF(str(path))
        try:
            samples = list(vcf.samples)
            for variant in vcf:
                site = self._build_site(variant, samples)
                if site is not None:
                    yield site
        finally:
            vcf.close()

    @abstractmethod
    def _build_site(self, variant: Any, samples: list[str]) -> MEISite | None:
        """Map one cyvcf2 variant to a :class:`MEISite`, or None to skip it."""
        raise NotImplementedError

    @staticmethod
    def _filters(variant: Any) -> tuple[str, ...]:
        """Extract the FILTER column.

        Note the cyvcf2 quirk: ``variant.FILTER`` is ``None`` for *both* PASS
        and an unset ('.') filter -- the two cannot be distinguished. We treat
        ``None`` as PASS, which is correct for callers (like xTEA) that always
        emit an explicit filter value.
        """
        f = variant.FILTER
        if f is None:
            return ("PASS",)
        return tuple(f.split(";"))
