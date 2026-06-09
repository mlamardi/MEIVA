"""Annotation layers (Layer 1+): genic context, consequence, frequency, ...

Each layer takes MEISite / cohort objects and enriches them. Layer 1 (genic
context) lives in `genic.py`.
"""

from __future__ import annotations

from meiva.annotate.genic import (
    Gene,
    GeneModel,
    GenicContext,
    GenicRegion,
    InMemoryGeneModel,
    InsertionOrientation,
    Transcript,
    annotate_genic,
)

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
