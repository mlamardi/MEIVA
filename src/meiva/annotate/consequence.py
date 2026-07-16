"""Layer 2: the MEI-aware consequence model.

Layer 1 (:mod:`meiva.annotate.genic`) answers *where* an insertion sits relative
to a gene. This module answers *what it is likely to do*, in the same spirit as
VEP/SnpEff: it maps each site to a mechanistic **consequence** term and a derived
ordinal **impact** tier, plus optional evidence flags.

Design choices (deliberate, see ROADMAP):

* **Transparent rule table, not a score or a learned model.** Every call traces
  to a known mobile-element mechanism, so the output is reproducible, testable,
  and defensible in a paper. The result schema leaves room for a future
  confidence score without breaking the columns.
* **Region sets the ceiling; family / orientation / length / biotype modify it.**
  A coding-sequence hit is disruptive regardless of family; an intronic hit's
  meaning depends on the element and its orientation relative to the gene.
* **Graceful degradation on missing inputs.** xTEA leaves ~1 in 8 *genic* calls
  with an unknown strand (``STRAND=.``), and the transduction guard can null a
  length. When orientation or length is unknown we fall back to the generic term
  and record a flag rather than guessing -- the orientation-dependent calls
  (Alu antisense exonization, sense full-length L1 polyA interference) are only
  made when the evidence is actually present.

The mechanisms encoded:

* **coding_disruption** -- insertion in a CDS; disrupts the protein. HIGH.
* **splice_disruption** -- insertion in the splice donor/acceptor core. HIGH.
* **exonization_candidate** -- Alu in an intron, *antisense*: the classic
  cryptic-splice-site exonization substrate. MODERATE.
* **polyA_interference** -- full-length L1 in an intron, *sense*: premature
  polyadenylation / transcriptional interference. MODERATE.
* **promoter_insertion** / **utr5_insertion** / **utr3_insertion** -- regulatory
  contexts; 5'UTR ranked above 3'UTR.
* **noncoding_exon_insertion** -- exon of a non-coding gene; elevated to MODERATE
  when the host is a lncRNA (a focus of the project).
* **intronic** -- intronic, not meeting the special cases above. LOW.
* **upstream_gene** / **downstream_gene** / **intergenic** -- MODIFIER.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from meiva.annotate.genic import GenicContext, GenicRegion, InsertionOrientation
from meiva.model import MEIFamily, MEISite

__all__ = [
    "FULL_LENGTH_L1_BP",
    "IMPACT_SEVERITY",
    "Consequence",
    "ConsequenceResult",
    "Impact",
    "classify_consequence",
]

# An L1 at least this long is treated as (near) full-length: it retains its 5'UTR
# internal promoter and ORFs and is the form associated with sense-orientation
# transcriptional interference. Full-length L1 is ~6 kb; a 5.5 kb floor captures
# near-complete copies while excluding the common 5'-truncated fragments.
FULL_LENGTH_L1_BP = 5500

# GENCODE gene_type values we treat as long non-coding for the purpose of
# elevating a non-coding exon hit. v47 uses "lncRNA"; the rest are legacy/related
# biotypes kept for robustness against other releases.
_LNCRNA_BIOTYPES = frozenset(
    {
        "lncRNA",
        "lincRNA",
        "antisense",
        "sense_intronic",
        "sense_overlapping",
        "macro_lncRNA",
        "bidirectional_promoter_lncRNA",
    }
)


class Impact(Enum):
    """Ordinal severity tier, mirroring VEP/SnpEff IMPACT."""

    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    MODIFIER = "MODIFIER"


# Higher = more severe; for sorting / picking the most severe of several.
IMPACT_SEVERITY: dict[Impact, int] = {
    Impact.HIGH: 3,
    Impact.MODERATE: 2,
    Impact.LOW: 1,
    Impact.MODIFIER: 0,
}


class Consequence(Enum):
    """Mechanistic consequence term for an MEI relative to the reported gene."""

    CODING_DISRUPTION = "coding_disruption"
    SPLICE_DISRUPTION = "splice_disruption"
    EXONIZATION_CANDIDATE = "exonization_candidate"
    POLYA_INTERFERENCE = "polyA_interference"
    PROMOTER_INSERTION = "promoter_insertion"
    UTR5_INSERTION = "utr5_insertion"
    UTR3_INSERTION = "utr3_insertion"
    NONCODING_EXON_INSERTION = "noncoding_exon_insertion"
    INTRONIC = "intronic"
    UPSTREAM_GENE = "upstream_gene"
    DOWNSTREAM_GENE = "downstream_gene"
    INTERGENIC = "intergenic"


_IMPACT: dict[Consequence, Impact] = {
    Consequence.CODING_DISRUPTION: Impact.HIGH,
    Consequence.SPLICE_DISRUPTION: Impact.HIGH,
    Consequence.EXONIZATION_CANDIDATE: Impact.MODERATE,
    Consequence.POLYA_INTERFERENCE: Impact.MODERATE,
    Consequence.PROMOTER_INSERTION: Impact.MODERATE,
    Consequence.UTR5_INSERTION: Impact.MODERATE,
    Consequence.UTR3_INSERTION: Impact.LOW,
    Consequence.NONCODING_EXON_INSERTION: Impact.LOW,  # elevated for lncRNA hosts
    Consequence.INTRONIC: Impact.LOW,
    Consequence.UPSTREAM_GENE: Impact.MODIFIER,
    Consequence.DOWNSTREAM_GENE: Impact.MODIFIER,
    Consequence.INTERGENIC: Impact.MODIFIER,
}

# SVA carries strong regulatory/expression-modulating potential; we flag it when
# it lands in a regulatory context (it doesn't change the term, only annotates it).
_SVA_REGULATORY_REGIONS = frozenset(
    {
        GenicRegion.PROMOTER,
        GenicRegion.UTR5,
        GenicRegion.UTR3,
        GenicRegion.EXON_NONCODING,
    }
)


@dataclass(frozen=True, slots=True)
class ConsequenceResult:
    """A consequence call: a mechanistic term, its impact tier, and evidence flags."""

    consequence: Consequence
    impact: Impact
    flags: tuple[str, ...] = field(default_factory=tuple)


def _is_lncrna(biotype: str | None) -> bool:
    return biotype is not None and biotype in _LNCRNA_BIOTYPES


def _classify_intronic(
    family: MEIFamily,
    orientation: InsertionOrientation,
    length: int | None,
    flags: list[str],
) -> Consequence:
    """Intronic insertions: orientation- and family-specific special cases.

    Mutates ``flags`` to record why a special case could not be applied.
    """
    if family is MEIFamily.ALU:
        if orientation is InsertionOrientation.ANTISENSE:
            return Consequence.EXONIZATION_CANDIDATE
        if orientation is InsertionOrientation.UNKNOWN:
            # could be the antisense exonization case -- we just can't tell
            flags.append("orientation_unknown")
        return Consequence.INTRONIC

    if family is MEIFamily.L1:
        if orientation is InsertionOrientation.SENSE:
            if length is None:
                flags.append("length_unknown")  # can't confirm full-length
                return Consequence.INTRONIC
            if length >= FULL_LENGTH_L1_BP:
                flags.append("full_length_l1")
                return Consequence.POLYA_INTERFERENCE
            return Consequence.INTRONIC  # truncated sense L1: better tolerated
        if orientation is InsertionOrientation.UNKNOWN:
            # could be the sense polyA-interference case -- can't tell
            flags.append("orientation_unknown")
        return Consequence.INTRONIC

    # SVA and anything else in an intron: no special intronic mechanism in v1
    return Consequence.INTRONIC


def classify_consequence(site: MEISite, ctx: GenicContext) -> ConsequenceResult:
    """Map a site + its Layer-1 genic context to a consequence term and impact tier.

    Pure and deterministic: depends only on the region, the element family, the
    orientation relative to the gene, the element length, and the host biotype.
    """
    region = ctx.region
    family = site.family
    flags: list[str] = []

    if region is GenicRegion.CDS:
        cons = Consequence.CODING_DISRUPTION
    elif region is GenicRegion.SPLICE_REGION:
        cons = Consequence.SPLICE_DISRUPTION
    elif region is GenicRegion.INTRON:
        cons = _classify_intronic(family, ctx.orientation, site.length, flags)
    elif region is GenicRegion.PROMOTER:
        cons = Consequence.PROMOTER_INSERTION
    elif region is GenicRegion.UTR5:
        cons = Consequence.UTR5_INSERTION
    elif region is GenicRegion.UTR3:
        cons = Consequence.UTR3_INSERTION
    elif region is GenicRegion.EXON_NONCODING:
        cons = Consequence.NONCODING_EXON_INSERTION
    elif region is GenicRegion.UPSTREAM:
        cons = Consequence.UPSTREAM_GENE
    elif region is GenicRegion.DOWNSTREAM:
        cons = Consequence.DOWNSTREAM_GENE
    else:  # GenicRegion.INTERGENIC
        cons = Consequence.INTERGENIC

    impact = _IMPACT[cons]

    # A non-coding exon hit in a lncRNA host is more interesting than a generic one.
    if cons is Consequence.NONCODING_EXON_INSERTION and _is_lncrna(ctx.gene_biotype):
        impact = Impact.MODERATE

    if family is MEIFamily.SVA and region in _SVA_REGULATORY_REGIONS:
        flags.append("sva_regulatory")

    return ConsequenceResult(consequence=cons, impact=impact, flags=tuple(flags))
