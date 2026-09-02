"""Annotation layers (Layer 1+): genic context, consequence, frequency, ...

Each layer takes MEISite / cohort objects and enriches them. Layer 1 (genic
context) lives in `genic.py`.
"""

from __future__ import annotations

from meiva.annotate.consequence import (
    Consequence,
    ConsequenceResult,
    Impact,
    classify_consequence,
)
from meiva.annotate.fantom5 import (
    Fantom5Model,
    RegulatoryContext,
    load_fantom5,
)
from meiva.annotate.fantom6 import (
    EvidenceTier,
    Fantom6Evidence,
    evidence_by_ensembl,
    load_cat_gene_ids,
    load_fantom6_evidence,
    target_to_ensembl,
)
from meiva.annotate.gencode import load_gencode, parse_gencode_gtf
from meiva.annotate.genic import (
    Gene,
    GeneModel,
    GenicContext,
    GenicRegion,
    IndexedGeneModel,
    InMemoryGeneModel,
    InsertionOrientation,
    Transcript,
    annotate_genic,
)

__all__ = [
    "Consequence",
    "ConsequenceResult",
    "EvidenceTier",
    "Fantom5Model",
    "Fantom6Evidence",
    "Gene",
    "GeneModel",
    "GenicContext",
    "GenicRegion",
    "Impact",
    "InMemoryGeneModel",
    "IndexedGeneModel",
    "InsertionOrientation",
    "RegulatoryContext",
    "Transcript",
    "annotate_genic",
    "classify_consequence",
    "evidence_by_ensembl",
    "load_cat_gene_ids",
    "load_fantom5",
    "load_fantom6_evidence",
    "load_gencode",
    "parse_gencode_gtf",
    "target_to_ensembl",
]
