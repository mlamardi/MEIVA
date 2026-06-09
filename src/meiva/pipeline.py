"""End-to-end glue: per-sample VCFs -> merged cohort -> genic annotation -> TSV.

This wires the components built so far into one runnable path so a cohort can be
taken from raw caller output to an annotated table. It is intentionally a
*Layer-1 preview*: the TSV carries cohort genotype summaries and genic context,
but not yet the consequence tiers (Layer 2) or population frequencies (Layer 3).

The pieces are separated for testability:

* :func:`annotate_cohort` and :func:`write_tsv` are pure (cohort + gene model in,
  rows out) and unit-testable in memory.
* :func:`run` is the convenience that also parses VCFs and loads GENCODE.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from meiva.annotate.genic import GeneModel, GenicContext, annotate_genic
from meiva.cohort import Cohort, CohortSite

__all__ = ["TSV_HEADER", "AnnotatedSite", "annotate_cohort", "run", "write_tsv"]


@dataclass(frozen=True, slots=True)
class AnnotatedSite:
    """A merged cohort site paired with its genic annotation."""

    cohort_site: CohortSite
    genic: GenicContext


def annotate_cohort(cohort: Cohort, model: GeneModel) -> list[AnnotatedSite]:
    """Run genic annotation over every site in a merged cohort."""
    return [AnnotatedSite(cs, annotate_genic(cs.site, model)) for cs in cohort.sites]


TSV_HEADER = [
    "chrom",
    "pos",
    "ci_lower",
    "ci_upper",
    "family",
    "strand",
    "length",
    "tsd",
    "n_carriers",
    "cohort_size",
    "carrier_frequency",
    "flags",
    "member_callers",
    "n_members",
    "region",
    "gene_id",
    "gene_name",
    "gene_strand",
    "transcript_id",
    "is_mane_select",
    "orientation",
    "distance",
    "carriers",
]


def _dosage_str(dosage: int | None) -> str:
    return "." if dosage is None else str(dosage)


def _row(a: AnnotatedSite) -> list[str]:
    site = a.cohort_site.site
    cs = a.cohort_site
    g = a.genic
    carriers = ";".join(
        f"{sample}:{_dosage_str(gt.dosage)}" for sample, gt in sorted(site.genotypes.items())
    )
    return [
        site.chrom,
        str(site.pos),
        str(site.ci_lower),
        str(site.ci_upper),
        site.family.value,
        site.strand.value,
        "" if site.length is None else str(site.length),
        site.tsd or "",
        str(cs.n_carriers),
        str(cs.cohort_size),
        f"{cs.carrier_frequency:.4f}",
        ";".join(cs.flags),
        ";".join(cs.member_callers),
        site.raw_info.get("MEIVA_N_MEMBERS", ""),
        g.region.value,
        g.gene_id or "",
        g.gene_name or "",
        "" if g.gene_strand is None else g.gene_strand.value,
        g.transcript_id or "",
        "true" if g.is_mane_select else "false",
        g.orientation.value,
        "" if g.distance is None else str(g.distance),
        carriers,
    ]


def write_tsv(sites: Iterable[AnnotatedSite], out: TextIO) -> None:
    """Write annotated sites as a tab-separated table with a header row."""
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(TSV_HEADER)
    for site in sites:
        writer.writerow(_row(site))


def run(
    vcf_paths: Iterable[str | Path],
    gencode_gtf: str | Path,
    out_tsv: str | Path,
    *,
    window: int | None = None,
) -> int:
    """Full path: merge ``vcf_paths``, annotate against ``gencode_gtf``, write ``out_tsv``.

    Returns the number of annotated cohort sites written. Imports the VCF and
    GENCODE machinery lazily so importing this module stays light.
    """
    from meiva.annotate.gencode import load_gencode
    from meiva.cohort import DEFAULT_MERGE_WINDOW, merge_vcfs

    cohort = merge_vcfs(vcf_paths, window=window if window is not None else DEFAULT_MERGE_WINDOW)
    model = load_gencode(gencode_gtf)
    annotated = annotate_cohort(cohort, model)
    with open(out_tsv, "w", newline="") as fh:
        write_tsv(annotated, fh)
    return len(annotated)
