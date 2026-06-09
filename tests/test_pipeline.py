"""Tests for the end-to-end pipeline (annotate_cohort, write_tsv, run)."""

import io
from pathlib import Path

import pytest

from meiva.annotate.genic import Gene, InMemoryGeneModel, Transcript
from meiva.cohort import merge_sites
from meiva.model import MEIFamily, MEISite, SampleGenotype, Strand
from meiva.pipeline import TSV_HEADER, annotate_cohort, run, write_tsv

# A coding gene so an overlapping site annotates as CDS.
GENE = Gene(
    gene_id="ENSG_X",
    gene_name="XGENE",
    chrom="chr1",
    start=1000,
    end=2200,
    strand=Strand.PLUS,
    transcripts=(
        Transcript(
            transcript_id="ENST_X",
            exons=((1000, 1200), (1500, 1700), (2000, 2200)),
            cds_start=1100,
            cds_end=2100,
            is_mane_select=True,
        ),
    ),
)
MODEL = InMemoryGeneModel([GENE])


def _site(pos: int, sample: str, alleles: tuple[int | None, ...], strand: Strand = Strand.PLUS):
    return MEISite(
        chrom="chr1",
        pos=pos,
        family=MEIFamily.ALU,
        strand=strand,
        source_caller="xtea",
        genotypes={sample: SampleGenotype(sample, alleles)},
    )


def test_annotate_cohort_and_tsv_roundtrip():
    cohort = merge_sites(
        [
            _site(1151, "A", (0, 1)),  # in CDS
            _site(1152, "B", (1, 1)),  # same site, jittered
            _site(50_000, "A", (0, 1)),  # intergenic
        ]
    )
    annotated = annotate_cohort(cohort, MODEL)
    buf = io.StringIO()
    write_tsv(annotated, buf)

    lines = buf.getvalue().splitlines()
    assert lines[0].split("\t") == TSV_HEADER
    rows = [dict(zip(TSV_HEADER, ln.split("\t"), strict=True)) for ln in lines[1:]]

    cds = next(r for r in rows if r["pos"] == "1151")  # median_low(1151,1152)=1151
    assert cds["region"] == "CDS"
    assert cds["gene_name"] == "XGENE"
    assert cds["is_mane_select"] == "true"
    assert cds["orientation"] == "SENSE"
    assert cds["n_carriers"] == "2"
    assert cds["carrier_frequency"] == "1.0000"
    assert cds["carriers"] == "A:1;B:2"

    inter = next(r for r in rows if r["pos"] == "50000")
    assert inter["region"] == "INTERGENIC"
    assert inter["gene_id"] == ""


def test_tsv_blank_for_missing_length_and_distance():
    cohort = merge_sites([_site(500_000, "A", (0, 1))])  # far intergenic, no length
    buf = io.StringIO()
    write_tsv(annotate_cohort(cohort, MODEL), buf)
    row = dict(zip(TSV_HEADER, buf.getvalue().splitlines()[1].split("\t"), strict=True))
    assert row["length"] == ""
    assert row["distance"] == ""
    assert row["region"] == "INTERGENIC"


# --------------------------------------------------------------------------- #
# Full run() over the real VDA VCFs + a synthetic GTF (local only)            #
# --------------------------------------------------------------------------- #
_REAL = sorted((Path(__file__).parent / "data").glob("ASL_VDA_*.vcf"))
_GTF_LINE = "\t".join
_SYNTH_GTF = (
    "\n".join(
        [
            _GTF_LINE(
                [
                    "chr1",
                    "HAVANA",
                    "gene",
                    "1",
                    "100",
                    ".",
                    "+",
                    ".",
                    'gene_id "G.1"; gene_name "G";',
                ]
            ),
            _GTF_LINE(
                [
                    "chr1",
                    "HAVANA",
                    "transcript",
                    "1",
                    "100",
                    ".",
                    "+",
                    ".",
                    'gene_id "G.1"; transcript_id "T.1";',
                ]
            ),
            _GTF_LINE(
                [
                    "chr1",
                    "HAVANA",
                    "exon",
                    "1",
                    "100",
                    ".",
                    "+",
                    ".",
                    'gene_id "G.1"; transcript_id "T.1";',
                ]
            ),
        ]
    )
    + "\n"
)


@pytest.mark.skipif(len(_REAL) < 2, reason="real VCF fixtures not present (gitignored)")
def test_run_end_to_end(tmp_path):
    gtf = tmp_path / "mini.gtf"
    gtf.write_text(_SYNTH_GTF)
    out = tmp_path / "out.tsv"
    n = run(_REAL, gtf, out)
    assert n > 0
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == TSV_HEADER
    assert len(lines) - 1 == n  # header + one row per site
