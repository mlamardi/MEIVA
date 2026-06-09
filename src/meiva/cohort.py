"""Cohort assembly: collapse per-sample MEI calls into a unified site set.

Per-sample MEI callers (xTEA, MELT) perform *discovery*: each sample's VCF
lists only the sites that sample carries, and the same physical insertion shows
up across samples at slightly different breakpoints (we observed +/-1 bp jitter
in real xTEA output). This module groups those per-sample records into one
canonical cohort site each, with a genotype matrix.

Two things here are easy to get wrong and are handled deliberately:

* **Clustering is seed-anchored, not single-linkage.** After sorting by
  (contig, family, position), a record joins the open cluster only if it is
  within ``window`` of the *seed* (the record that opened the cluster). This
  bounds every cluster to a ``window``-wide span and avoids the classic
  single-linkage *chaining* failure, where A~B and B~C silently merge A and C
  even though they are far apart.

* **Absence is not reference.** Because the inputs are discovery-only, a sample
  with no record at a site is of *unknown* genotype -- it was never tested
  there, not confirmed homozygous reference. We therefore store a *sparse*
  genotype matrix (carriers only) and report a ``carrier_frequency`` that is
  explicitly a discovery-based estimate, NOT a true allele frequency. A real AF
  requires force-genotyping every discovered site back against every sample's
  reads -- a later step, not this one.
"""

from __future__ import annotations

import statistics
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

from meiva.model import MEIFamily, MEISite, SampleGenotype, Strand

__all__ = ["DEFAULT_MERGE_WINDOW", "Cohort", "CohortSite", "merge_sites", "merge_vcfs"]

#: default breakpoint tolerance (bp) for treating two calls as the same site.
#: Generous relative to observed jitter (~1 bp) and the TSD footprint (~<=30 bp),
#: but far tighter than general-SV mergers, which befits precise MEI breakpoints.
DEFAULT_MERGE_WINDOW = 50


@dataclass(frozen=True, slots=True)
class CohortSite:
    """One merged site plus cohort-level summary stats.

    ``site`` is a canonical :class:`~meiva.model.MEISite` whose imprecision
    interval (``ci_lower``/``ci_upper``) spans the observed breakpoints of all
    members, and whose genotype matrix holds the carrying samples only.
    """

    site: MEISite
    n_carriers: int
    cohort_size: int
    member_callers: tuple[str, ...] = ()
    flags: tuple[str, ...] = ()

    @property
    def carrier_frequency(self) -> float:
        """Carriers / cohort size.

        Discovery-based: non-carriers are samples *not observed* to carry the
        insertion, which conflates true reference with false negatives. Treat
        this as a lower-bound estimate, not a population allele frequency.
        """
        return self.n_carriers / self.cohort_size if self.cohort_size else 0.0


@dataclass(frozen=True, slots=True)
class Cohort:
    """A merged cohort: the sample roster plus the unified site list."""

    roster: tuple[str, ...]
    sites: tuple[CohortSite, ...] = field(default_factory=tuple)

    @property
    def n_samples(self) -> int:
        return len(self.roster)

    @property
    def n_sites(self) -> int:
        return len(self.sites)


def _contig_key(chrom: str) -> tuple[int, int, str]:
    """Natural-ish sort key so chr2 precedes chr10 and sex/MT sort last."""
    name = chrom[3:] if chrom.startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if name in special:
        return (0, special[name], "")
    return (1, 0, name)


def _better_genotype(a: SampleGenotype, b: SampleGenotype) -> SampleGenotype:
    """Pick the more informative of two genotypes for the same sample.

    Higher dosage wins (a carrier beats a non-call; hom beats het); ties keep
    the incumbent. Used only on collision, which is itself flagged.
    """
    da = a.dosage if a.dosage is not None else -1
    db = b.dosage if b.dosage is not None else -1
    return b if db > da else a


def _consensus_strand(members: list[MEISite]) -> Strand:
    counts = Counter(m.strand for m in members if m.strand is not Strand.UNKNOWN)
    if not counts:
        return Strand.UNKNOWN
    # most common; deterministic tie-break on the strand symbol
    return max(counts.items(), key=lambda kv: (kv[1], kv[0].value))[0]


def _build_cohort_site(members: list[MEISite], cohort_size: int) -> CohortSite:
    positions = [m.pos for m in members]
    anchor = statistics.median_low(positions)  # always a real observed coordinate
    lo, hi = min(positions), max(positions)

    lengths = [m.length for m in members if m.length is not None]
    length = statistics.median_low(lengths) if lengths else None

    tsds = [m.tsd for m in members if m.tsd is not None]
    tsd = Counter(tsds).most_common(1)[0][0] if tsds else None

    callers = tuple(sorted({m.source_caller for m in members if m.source_caller}))
    if len(callers) == 1:
        source: str | None = callers[0]
    elif callers:
        source = "merged"
    else:
        source = None

    genotypes: dict[str, SampleGenotype] = {}
    collision = False
    for m in members:
        for sample, gt in m.genotypes.items():
            if sample in genotypes:
                collision = True
                genotypes[sample] = _better_genotype(genotypes[sample], gt)
            else:
                genotypes[sample] = gt

    flags: list[str] = []
    if collision:
        flags.append("SAMPLE_COLLISION")
    if len(callers) > 1:
        flags.append("MULTI_CALLER")

    raw_info = {
        "MEIVA_N_MEMBERS": str(len(members)),
        "MEIVA_MEMBER_POS": ",".join(str(p) for p in positions),
    }
    if callers:
        raw_info["MEIVA_MEMBER_CALLERS"] = ",".join(callers)

    site = MEISite(
        chrom=members[0].chrom,
        pos=anchor,
        family=members[0].family,
        ci_lower=lo - anchor,  # <= 0
        ci_upper=hi - anchor,  # >= 0
        strand=_consensus_strand(members),
        length=length,
        tsd=tsd,
        source_caller=source,
        genotypes=genotypes,
        raw_info=raw_info,
    )
    n_carriers = sum(1 for gt in genotypes.values() if gt.is_carrier)
    return CohortSite(
        site=site,
        n_carriers=n_carriers,
        cohort_size=cohort_size,
        member_callers=callers,
        flags=tuple(flags),
    )


def merge_sites(sites: Iterable[MEISite], *, window: int = DEFAULT_MERGE_WINDOW) -> Cohort:
    """Merge per-sample :class:`MEISite` records into a :class:`Cohort`.

    Sites are grouped when they share a contig and family and their breakpoints
    fall within ``window`` of the cluster seed. The cohort roster is the union
    of all sample IDs seen across the inputs.
    """
    if window < 0:
        raise ValueError(f"window must be >= 0, got {window}")

    materialised = list(sites)
    roster = sorted({s for site in materialised for s in site.genotypes})
    cohort_size = len(roster)

    ordered = sorted(materialised, key=lambda s: (_contig_key(s.chrom), s.family.value, s.pos))

    cohort_sites: list[CohortSite] = []
    cluster: list[MEISite] = []
    seed_chrom: str | None = None
    seed_family: MEIFamily | None = None
    seed_pos = 0

    for s in ordered:
        same = (
            s.chrom == seed_chrom
            and s.family is seed_family
            and (s.pos - seed_pos) <= window  # ordered ascending, so >= 0
        )
        if cluster and same:
            cluster.append(s)
        else:
            if cluster:
                cohort_sites.append(_build_cohort_site(cluster, cohort_size))
            cluster = [s]
            seed_chrom, seed_family, seed_pos = s.chrom, s.family, s.pos
    if cluster:
        cohort_sites.append(_build_cohort_site(cluster, cohort_size))

    return Cohort(roster=tuple(roster), sites=tuple(cohort_sites))


def merge_vcfs(paths: Iterable[str | Path], *, window: int = DEFAULT_MERGE_WINDOW) -> Cohort:
    """Convenience: auto-detect parsers, parse all ``paths``, and merge.

    Imported lazily so that ``merge_sites`` stays free of the cyvcf2 dependency.
    """
    from meiva.io import detect_parser

    def all_sites() -> Iterable[MEISite]:
        for path in paths:
            yield from detect_parser(path).parse(path)

    return merge_sites(all_sites(), window=window)
