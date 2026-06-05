"""xTEA VCF parser.

Maps xTEA's symbolic-insertion records onto :class:`~meiva.model.MEISite`.
The non-obvious decisions, all justified by inspecting real xTEA output:

* **Family** comes from the third ``SVTYPE`` field (``INS:ME:ALU`` -> ALU,
  ``INS:ME:LINE1`` -> L1, ``INS:ME:HERV-K`` -> HERVK). Only ``INS:ME:*``
  records are MEIs; ``INS:MT`` / ``INS:PSDGN`` are out of scope and skipped.
* **Length**: xTEA's ``SVLEN`` is a real element length for normal calls but a
  20-35 Mb genomic coordinate for ``orphan_or_sibling_transduction`` records.
  We never trust it blindly -- :func:`plausible_length` nulls implausible
  values and we flag the site with ``MEIVA_LENGTH_UNRELIABLE`` rather than
  emitting a fictitious multi-megabase insertion.
* **Imprecision**: xTEA reports no ``CIPOS``. ``END`` is always exactly
  ``POS + TSDLEN`` (the TSD footprint, a biological feature, not breakpoint
  uncertainty), so we keep ``ci=0`` and rely on padded matching downstream.
* **TSD**: the ``+``/``-`` prefix is an orientation marker; we store the bare
  sequence, and ``NULL`` becomes ``None``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from meiva.io.base import (
    Cyvcf2Parser,
    genotypes_from_cyvcf2,
    normalize_contig,
    plausible_length,
    raw_info_from_cyvcf2,
)
from meiva.model import MEIFamily, MEISite, Strand

_STRAND = {"+": Strand.PLUS, "-": Strand.MINUS}

# INFO IDs emitted by xTEA but not by MELT/other callers -- used for sniffing.
_XTEA_SIGNATURE = ("AF_FMAP", "LC_CLUSTER", "RD_AKR_NRC")


class XteaParser(Cyvcf2Parser):
    """Parser for VCFs produced by xTEA."""

    caller = "xtea"

    @classmethod
    def sniff(cls, path: str | Path) -> bool:
        with open(path) as fh:
            for line in fh:
                if not line.startswith("##"):
                    return False  # past the header without a match
                if any(sig in line for sig in _XTEA_SIGNATURE):
                    return True
        return False

    def _build_site(self, variant: Any, samples: list[str]) -> MEISite | None:
        svtype = variant.INFO.get("SVTYPE")
        if svtype is None or not str(svtype).startswith("INS:ME"):
            return None  # not a mobile-element insertion (e.g. INS:MT) -> skip

        family = MEIFamily.from_raw(str(svtype).split(":")[-1])
        length, length_unreliable = plausible_length(variant.INFO.get("SVLEN"))

        tsd_raw = variant.INFO.get("TSD")
        tsd: str | None = None
        if tsd_raw is not None and tsd_raw != "NULL":
            tsd = tsd_raw[1:] if tsd_raw[:1] in "+-" else tsd_raw

        allele_freq = self._allele_freq(variant.INFO.get("AF"))
        strand = _STRAND.get(str(variant.INFO.get("STRAND") or "."), Strand.UNKNOWN)

        raw_info = raw_info_from_cyvcf2(variant)
        if length_unreliable:
            raw_info["MEIVA_LENGTH_UNRELIABLE"] = "1"

        return MEISite(
            chrom=normalize_contig(variant.CHROM),
            pos=variant.POS,
            family=family,
            strand=strand,
            length=length,
            tsd=tsd,
            allele_freq=allele_freq,
            qual=variant.QUAL,
            filters=self._filters(variant),
            site_id=variant.ID,
            source_caller=self.caller,
            genotypes=genotypes_from_cyvcf2(variant, samples),
            raw_info=raw_info,
        )

    @staticmethod
    def _allele_freq(af: Any) -> float | None:
        """Coerce xTEA's ``AF`` (Number=A; a list under cyvcf2) to a float in [0,1]."""
        if af is None:
            return None
        if isinstance(af, (list, tuple)):
            if not af:
                return None
            af = af[0]
        # Clamp only against float drift; xTEA AF is already bounded in (0,1].
        return min(max(float(af), 0.0), 1.0)
