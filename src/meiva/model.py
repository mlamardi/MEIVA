"""Canonical data model for mobile element insertions (MEIs).

This module defines the single source of truth that every other layer in
MEIVA operates on: :class:`MEISite`. The division of responsibility is
deliberate and strict:

* **Parsers** (``meiva.io``) translate caller-specific VCFs into *already
  normalised* :class:`MEISite` objects. All coercion lives there.
* **This module** only *validates* invariants. It never silently fixes,
  coerces, or guesses. A malformed input raises at construction time rather
  than producing a subtly wrong annotation 200 lines downstream.

Keeping normalisation out of the record means the normalisation contract
exists in exactly one place, and a record that exists is, by construction, a
record we can trust.

Critical conventions -- read these before touching anything downstream:

* ``pos`` is **1-based**, matching the VCF POS column. It is the single most
  likely breakpoint as reported by the caller.
* ``ci_lower`` / ``ci_upper`` encode breakpoint *imprecision* as signed
  offsets from ``pos`` (VCF ``CIPOS`` semantics: ``ci_lower <= 0 <=
  ci_upper``). MEI callers almost never pin a base-exact breakpoint, so
  essentially every downstream overlap or cross-callset match MUST use the
  interval (see :meth:`MEISite.search_interval`), not ``pos`` alone. Treating
  an MEI as a point is the most common correctness bug in this domain.
* ``chrom`` is stored **with** the ``chr`` prefix (GRCh38 / GENCODE /
  gnomAD-SV v4 convention). Mismatched contig naming is a silent killer: it
  yields zero overlaps and *no error*. We therefore make it loud -- the
  record refuses to exist without the prefix, and the parser is responsible
  for normalising to it.
* An MEI is an *insertion*; it occupies ~zero reference bases. We model its
  reference footprint as the imprecision interval, NOT the length of the
  inserted element (``length`` is a separate, optional attribute of the
  inserted sequence).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType

__all__ = ["MEIFamily", "MEISite", "SampleGenotype", "Strand"]


class MEIFamily(str, Enum):
    """Recently-active human MEI families (plus escape hatches).

    Subclassing ``str`` keeps members JSON- and VCF-friendly while still
    giving typo-proof identity comparisons. ``OTHER`` is intentional: an
    unexpected caller label should be *flagged*, never crash ingest or get
    silently dropped. The original string is preserved on the record via
    ``raw_info`` so nothing is lost.
    """

    ALU = "ALU"
    L1 = "L1"
    SVA = "SVA"
    HERVK = "HERVK"
    OTHER = "OTHER"

    @classmethod
    def from_raw(cls, raw: str) -> MEIFamily:
        """Map a caller's family label onto the canonical vocabulary.

        This is a *parser-facing* helper: callers disagree wildly on
        spelling (``ALU`` / ``AluY`` / ``Alu``, ``L1`` / ``LINE1`` /
        ``LINE/L1``, ``SVA``, ``HERVK`` / ``HERV-K``). We normalise to
        alphanumerics, uppercase, then match by prefix. Unknown labels
        resolve to :attr:`OTHER` rather than raising, so one odd record can
        never abort a whole-cohort run.
        """
        key = "".join(ch for ch in raw.upper() if ch.isalnum())
        if key.startswith("ALU"):
            return cls.ALU
        if key.startswith("L1") or key.startswith("LINE"):
            return cls.L1
        if key.startswith("SVA"):
            return cls.SVA
        if key.startswith("HERV"):
            return cls.HERVK
        return cls.OTHER


class Strand(str, Enum):
    """Insertion orientation relative to the reference plus strand.

    ``UNKNOWN`` is the default on purpose. Orientation drives consequence
    prediction (a sense-oriented intronic L1 behaves very differently from an
    antisense one), so the consequence layer MUST handle "we don't know"
    explicitly rather than us silently defaulting to ``PLUS`` and inventing a
    mechanism that isn't supported by the data.
    """

    PLUS = "+"
    MINUS = "-"
    UNKNOWN = "."


@dataclass(frozen=True, slots=True)
class SampleGenotype:
    """A single sample's genotype at an MEI site.

    ``alleles`` mirrors VCF ``GT``: a tuple of allele indices, with ``None``
    for a missing allele (``.``). For a biallelic insertion the alt index is
    ``1``; any index ``>= 1`` is treated as carrying the insertion so that
    multiallelic edge cases degrade gracefully rather than being miscounted
    as reference.
    """

    sample: str
    alleles: tuple[int | None, ...] = ()
    phased: bool = False

    @property
    def is_missing(self) -> bool:
        """True if no allele was observed (empty or all ``.``)."""
        return len(self.alleles) == 0 or all(a is None for a in self.alleles)

    @property
    def is_carrier(self) -> bool:
        """True if at least one observed allele is the insertion."""
        return any(a is not None and a >= 1 for a in self.alleles)

    @property
    def dosage(self) -> int | None:
        """Count of insertion alleles, or ``None`` if fully missing.

        Partial-missing genotypes (e.g. ``./1``) return the count of the
        *observed* alt alleles; this is documented rather than hidden because
        partial calls are common in low-coverage MEI data and silently
        treating them as ``0`` or ``2`` would bias frequency estimates.
        """
        if self.is_missing:
            return None
        return sum(1 for a in self.alleles if a is not None and a >= 1)


@dataclass(frozen=True, slots=True)
class MEISite:
    """A single mobile-element-insertion site on the reference genome.

    The record is a *site*, not a per-sample call: genomic and functional
    annotation depend only on the locus, so we annotate once and carry
    per-sample genotypes alongside. This is both a correctness decision
    (annotation can't differ between carriers of the same site) and a
    performance one (annotate N sites, not N x S sample-calls).

    Identity / hashing keys on the locus-defining fields only (chrom, pos,
    imprecision, family, subfamily, strand, source). Payload fields
    (genotypes, raw_info, qual, filters, etc.) are excluded from equality so
    two records describing the same locus compare equal regardless of which
    samples carry it -- which is what cross-callset matching needs.
    """

    # --- locus identity (compared / hashed) ---
    chrom: str
    pos: int  # 1-based breakpoint (VCF POS)
    family: MEIFamily
    ci_lower: int = 0  # offset from pos, <= 0  (CIPOS semantics)
    ci_upper: int = 0  # offset from pos, >= 0
    strand: Strand = Strand.UNKNOWN
    subfamily: str | None = None  # e.g. "AluYa5", "L1HS"
    source_caller: str | None = None  # e.g. "xtea", "melt"

    # --- payload / evidence (excluded from equality & hash) ---
    length: int | None = field(default=None, compare=False)  # inserted element bp
    tsd: str | None = field(default=None, compare=False)  # target site duplication
    site_id: str | None = field(default=None, compare=False)  # caller-assigned ID
    qual: float | None = field(default=None, compare=False)
    allele_freq: float | None = field(default=None, compare=False)  # caller-estimated VAF in [0,1]
    filters: tuple[str, ...] = field(default=(), compare=False)  # VCF FILTER column
    genotypes: Mapping[str, SampleGenotype] = field(default_factory=dict, compare=False)
    raw_info: Mapping[str, str] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        # VALIDATION ONLY. No coercion -- the parser must pre-normalise.
        if not self.chrom:
            raise ValueError("chrom must be a non-empty string")
        if not self.chrom.startswith("chr"):
            raise ValueError(
                f"chrom {self.chrom!r} is not 'chr'-prefixed. MEIVA works in "
                "GRCh38 'chr' space; normalise contig names in the parser. "
                "(Mismatched contig naming silently produces zero overlaps.)"
            )
        if self.pos < 1:
            raise ValueError(f"pos must be >= 1 (1-based), got {self.pos}")
        if self.ci_lower > 0:
            raise ValueError(f"ci_lower must be <= 0 (offset from pos), got {self.ci_lower}")
        if self.ci_upper < 0:
            raise ValueError(f"ci_upper must be >= 0 (offset from pos), got {self.ci_upper}")
        if self.length is not None and self.length <= 0:
            raise ValueError(f"length, if given, must be > 0, got {self.length}")
        if self.allele_freq is not None and not 0.0 <= self.allele_freq <= 1.0:
            raise ValueError(f"allele_freq, if given, must be in [0, 1], got {self.allele_freq}")
        if not isinstance(self.family, MEIFamily):
            raise TypeError(f"family must be a MEIFamily, got {type(self.family).__name__}")
        if not isinstance(self.strand, Strand):
            raise TypeError(f"strand must be a Strand, got {type(self.strand).__name__}")

        # Install immutable, defensively-copied views of the mappings. A frozen
        # value object whose contents are mutable is a foot-gun; the
        # object.__setattr__ dance is the accepted way to do this under frozen.
        object.__setattr__(self, "genotypes", MappingProxyType(dict(self.genotypes)))
        object.__setattr__(self, "raw_info", MappingProxyType(dict(self.raw_info)))
        object.__setattr__(self, "filters", tuple(self.filters))

    # ------------------------------------------------------------------ #
    # Derived geometry                                                    #
    # ------------------------------------------------------------------ #
    @property
    def is_precise(self) -> bool:
        """True when the caller reported a base-exact breakpoint."""
        return self.ci_lower == 0 and self.ci_upper == 0

    def search_interval(self, padding: int = 0) -> tuple[int, int]:
        """Return the breakpoint footprint as a **0-based half-open** interval.

        This is the canonical conversion from the record's 1-based ``pos`` +
        signed CI into the ``[start, end)`` convention used by every interval
        library and by BED. Centralising the arithmetic here means the
        off-by-one lives in one audited place instead of being re-derived
        (and re-broken) in every annotation layer.

        ``padding`` widens the interval symmetrically -- use it for tolerant
        cross-callset matching (e.g. the +/- 50 bp window MEI callers
        conventionally use to call two breakpoints "the same site").

        The lower bound is clamped at 0; the upper bound is **not** clamped to
        the contig length here, because the record has no knowledge of contig
        sizes -- the annotator clamps against the reference it loads.
        """
        if padding < 0:
            raise ValueError(f"padding must be >= 0, got {padding}")
        start_1based = self.pos + self.ci_lower  # ci_lower <= 0
        end_1based = self.pos + self.ci_upper  # ci_upper >= 0
        # 1-based inclusive -> 0-based half-open: subtract 1 from start; the
        # inclusive 1-based end equals the exclusive 0-based end.
        start = start_1based - 1 - padding
        end = end_1based + padding
        return (max(start, 0), end)

    def __str__(self) -> str:
        sub = f"({self.subfamily})" if self.subfamily else ""
        ci = "" if self.is_precise else f" [{self.ci_lower:+d},{self.ci_upper:+d}]"
        return f"{self.chrom}:{self.pos:,} {self.family.value}{sub}{ci} {self.strand.value}"
