"""MEIVA — caller-agnostic functional annotation of mobile element insertions.

"VEP for mobile elements." See the README for the design philosophy and roadmap.
"""

from meiva.model import MEIFamily, MEISite, SampleGenotype, Strand

__version__ = "0.0.1"
__all__ = ["MEIFamily", "MEISite", "SampleGenotype", "Strand", "__version__"]
