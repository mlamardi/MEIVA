"""End-to-end glue: per-sample VCFs -> merged cohort -> genic annotation -> TSV.

This wires the components built so far into one runnable path so a cohort can be
taken from raw caller output to an annotated table. The TSV carries cohort
genotype summaries, Layer-1 genic context, and Layer-2 consequence terms with
impact tiers. Population frequencies (Layer 3) are not yet included --
``carrier_frequency`` is a discovery-based lower bound, not a true allele frequency.

The pieces are separated for testability:

* :func:`annotate_cohort` and :func:`write_tsv` are pure (cohort + gene model in,
  rows out) and unit-testable in memory.
* :func:`run` is the convenience that also parses VCFs and loads GENCODE.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from meiva.annotate.consequence import ConsequenceResult, classify_consequence
from meiva.annotate.fantom5 import Fantom5Model, RegulatoryContext
from meiva.annotate.fantom6 import Fantom6Evidence
from meiva.annotate.genic import GeneModel, GenicContext, annotate_genic
from meiva.cohort import Cohort, CohortSite

__all__ = [
    "TSV_HEADER",
    "AnnotatedSite",
    "annotate_cohort",
    "annotate_vcfs",
    "base_gene_id",
    "run",
    "write_tsv",
]


@dataclass(frozen=True, slots=True)
class AnnotatedSite:
    """A merged cohort site paired with its genic annotation and consequence call."""

    cohort_site: CohortSite
    genic: GenicContext
    consequence: ConsequenceResult
    fantom6: Fantom6Evidence | None = None
    """FANTOM6 knockdown evidence for the reported gene, when that gene was tested."""
    regulatory: RegulatoryContext | None = None
    """FANTOM5 enhancer and CAGE-promoter context for the insertion site."""


def base_gene_id(gene_id: str | None) -> str | None:
    """Strip the version suffix from an Ensembl gene ID (``ENSG...5`` -> ``ENSG...``)."""
    if gene_id is None:
        return None
    return gene_id.split(".", 1)[0]


def annotate_cohort(
    cohort: Cohort,
    model: GeneModel,
    *,
    fantom6: Mapping[str, Fantom6Evidence] | None = None,
    fantom5: Fantom5Model | None = None,
) -> list[AnnotatedSite]:
    """Run Layer-1 genic annotation and Layer-2 consequence over every cohort site.

    ``fantom6`` maps unversioned Ensembl gene IDs to knockdown evidence (see
    :func:`meiva.annotate.fantom6.evidence_by_ensembl`). It is joined on the reported
    gene's ID -- never on symbol -- and left absent when the gene was not tested.
    """
    annotated: list[AnnotatedSite] = []
    for cs in cohort.sites:
        genic = annotate_genic(cs.site, model)
        consequence = classify_consequence(cs.site, genic)
        evidence: Fantom6Evidence | None = None
        if fantom6:
            gid = base_gene_id(genic.gene_id)
            if gid is not None:
                evidence = fantom6.get(gid)
        regulatory = fantom5.annotate(cs.site.chrom, cs.site.pos) if fantom5 else None
        annotated.append(AnnotatedSite(cs, genic, consequence, evidence, regulatory))
    return annotated


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
    "gene_biotype",
    "gene_strand",
    "transcript_id",
    "is_mane_select",
    "orientation",
    "distance",
    "consequence",
    "impact",
    "consequence_flags",
    "fantom6_evidence",
    "fantom6_cell_types",
    "fantom5_enhancer",
    "fantom5_enhancer_distance",
    "fantom5_promoter",
    "fantom5_promoter_gene",
    "fantom5_promoter_rank",
    "fantom5_promoter_distance",
    "carriers",
]


def _dosage_str(dosage: int | None) -> str:
    return "." if dosage is None else str(dosage)


def _row(a: AnnotatedSite) -> list[str]:
    site = a.cohort_site.site
    cs = a.cohort_site
    g = a.genic
    reg = a.regulatory
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
        g.gene_biotype or "",
        "" if g.gene_strand is None else g.gene_strand.value,
        g.transcript_id or "",
        "true" if g.is_mane_select else "false",
        g.orientation.value,
        "" if g.distance is None else str(g.distance),
        a.consequence.consequence.value,
        a.consequence.impact.value,
        ";".join(a.consequence.flags),
        "" if a.fantom6 is None else a.fantom6.tier.value,
        "" if a.fantom6 is None else ";".join(a.fantom6.cell_types),
        "" if reg is None or reg.enhancer_id is None else reg.enhancer_id,
        "" if reg is None or reg.enhancer_distance is None else str(reg.enhancer_distance),
        "" if reg is None or reg.promoter_id is None else reg.promoter_id,
        "" if reg is None or reg.promoter_gene is None else reg.promoter_gene,
        "" if reg is None or reg.promoter_rank is None else str(reg.promoter_rank),
        "" if reg is None or reg.promoter_distance is None else str(reg.promoter_distance),
        carriers,
    ]


def write_tsv(sites: Iterable[AnnotatedSite], out: TextIO) -> None:
    """Write annotated sites as a tab-separated table with a header row."""
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(TSV_HEADER)
    for site in sites:
        writer.writerow(_row(site))


def annotate_vcfs(
    vcf_paths: Iterable[str | Path],
    gencode_gtf: str | Path,
    *,
    window: int | None = None,
    fantom6: Mapping[str, Fantom6Evidence] | None = None,
    fantom5: Fantom5Model | None = None,
) -> list[AnnotatedSite]:
    """Merge ``vcf_paths`` and annotate against ``gencode_gtf`` -- the pure core of :func:`run`.

    Returns the annotated cohort sites without writing anything, so callers (the
    CLI, tests) can stream them wherever they like. Imports the VCF and GENCODE
    machinery lazily so importing this module stays light.
    """
    from meiva.annotate.gencode import load_gencode
    from meiva.cohort import DEFAULT_MERGE_WINDOW, merge_vcfs

    cohort = merge_vcfs(vcf_paths, window=window if window is not None else DEFAULT_MERGE_WINDOW)
    model = load_gencode(gencode_gtf)
    return annotate_cohort(cohort, model, fantom6=fantom6, fantom5=fantom5)


def run(
    vcf_paths: Iterable[str | Path],
    gencode_gtf: str | Path,
    out_tsv: str | Path,
    *,
    window: int | None = None,
    fantom6: Mapping[str, Fantom6Evidence] | None = None,
    fantom5: Fantom5Model | None = None,
) -> int:
    """Full path: merge ``vcf_paths``, annotate against ``gencode_gtf``, write ``out_tsv``.

    Returns the number of annotated cohort sites written.
    """
    annotated = annotate_vcfs(
        vcf_paths, gencode_gtf, window=window, fantom6=fantom6, fantom5=fantom5
    )
    with open(out_tsv, "w", newline="") as fh:
        write_tsv(annotated, fh)
    return len(annotated)
