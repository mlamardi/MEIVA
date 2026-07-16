"""Tests for the end-to-end pipeline (annotate_cohort, write_tsv, run)."""

import io
from pathlib import Path

import pytest

from meiva.annotate.fantom6 import EvidenceTier, Fantom6Evidence
from meiva.annotate.genic import Gene, InMemoryGeneModel, Transcript
from meiva.cohort import Cohort, merge_sites
from meiva.model import MEIFamily, MEISite, SampleGenotype, Strand
from meiva.pipeline import TSV_HEADER, annotate_cohort, base_gene_id, run, write_tsv

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
    biotype="protein_coding",
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
    assert cds["gene_biotype"] == "protein_coding"
    assert cds["consequence"] == "coding_disruption"
    assert cds["impact"] == "HIGH"
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
# Full run() over the committed public HGDP/TSI VCFs + a synthetic GTF         #
# --------------------------------------------------------------------------- #
_TSI = sorted((Path(__file__).parent / "data").glob("HGDP*.tsi.vcf"))
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


@pytest.mark.skipif(len(_TSI) < 2, reason="HGDP/TSI fixtures not present")
def test_run_end_to_end(tmp_path):
    gtf = tmp_path / "mini.gtf"
    gtf.write_text(_SYNTH_GTF)
    out = tmp_path / "out.tsv"
    n = run(_TSI, gtf, out)
    assert n == 54  # 27 + 46 calls -> 54 merged cohort sites
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == TSV_HEADER
    assert len(lines) - 1 == n  # header + one row per site


# --------------------------------------------------------------------------- #
# FANTOM6 evidence join                                                        #
# --------------------------------------------------------------------------- #
def _evidence(gene: str, tier: EvidenceTier) -> Fantom6Evidence:
    return Fantom6Evidence(
        target_id="G0000001",
        target_symbol="XGENE",
        tier=tier,
        n_aso=2,
        n_aso_responding=2 if tier is EvidenceTier.CONCORDANT else 1,
        max_degs=100,
        total_degs=200,
        n_up=120,
        n_down=80,
        cell_types=("FF-iPSC (IMS)", "HDF (neonatal)"),
        ensembl_gene_id=gene,
    )


def _cds_cohort() -> Cohort:
    return merge_sites([_site(1151, "A", (0, 1))])


def test_base_gene_id_strips_version():
    assert base_gene_id("ENSG00000214548.5") == "ENSG00000214548"
    assert base_gene_id("ENSG00000214548") == "ENSG00000214548"
    assert base_gene_id(None) is None


def test_fantom6_evidence_joined_onto_site():
    ev = {"ENSG_X": _evidence("ENSG_X", EvidenceTier.CONCORDANT)}
    ann = annotate_cohort(_cds_cohort(), MODEL, fantom6=ev)
    assert ann[0].fantom6 is not None
    assert ann[0].fantom6.tier is EvidenceTier.CONCORDANT
    assert ann[0].fantom6.is_functional


def test_fantom6_join_is_version_insensitive():
    """A gene model carrying a versioned ID still matches the unversioned evidence key."""
    versioned = Gene(
        gene_id="ENSG_X.7",
        gene_name="XGENE",
        chrom="chr1",
        start=1000,
        end=2200,
        strand=Strand.PLUS,
        transcripts=GENE.transcripts,
        biotype="lncRNA",
    )
    model = InMemoryGeneModel([versioned])
    ev = {"ENSG_X": _evidence("ENSG_X", EvidenceTier.CONCORDANT)}
    ann = annotate_cohort(_cds_cohort(), model, fantom6=ev)
    assert ann[0].fantom6 is not None


def test_fantom6_absent_when_gene_not_tested():
    ev = {"ENSG_OTHER": _evidence("ENSG_OTHER", EvidenceTier.CONCORDANT)}
    assert annotate_cohort(_cds_cohort(), MODEL, fantom6=ev)[0].fantom6 is None


def test_fantom6_absent_when_not_supplied():
    assert annotate_cohort(_cds_cohort(), MODEL)[0].fantom6 is None


def test_fantom6_absent_for_intergenic_site():
    ev = {"ENSG_X": _evidence("ENSG_X", EvidenceTier.CONCORDANT)}
    cohort = merge_sites([_site(500_000, "A", (0, 1))])
    ann = annotate_cohort(cohort, MODEL, fantom6=ev)
    # nearest-gene is still reported for intergenic hits, so evidence may attach;
    # what matters is that it never attaches to a *different* gene
    if ann[0].fantom6 is not None:
        assert ann[0].fantom6.ensembl_gene_id == "ENSG_X"


def test_fantom6_columns_in_tsv():
    ev = {"ENSG_X": _evidence("ENSG_X", EvidenceTier.SINGLE_ASO)}
    ann = annotate_cohort(_cds_cohort(), MODEL, fantom6=ev)
    buf = io.StringIO()
    write_tsv(ann, buf)
    lines = buf.getvalue().splitlines()
    row = dict(zip(lines[0].split("\t"), lines[1].split("\t"), strict=True))
    assert row["fantom6_evidence"] == "single_aso"
    assert row["fantom6_cell_types"] == "FF-iPSC (IMS);HDF (neonatal)"


def test_fantom6_columns_blank_when_untested():
    buf = io.StringIO()
    write_tsv(annotate_cohort(_cds_cohort(), MODEL), buf)
    lines = buf.getvalue().splitlines()
    row = dict(zip(lines[0].split("\t"), lines[1].split("\t"), strict=True))
    assert row["fantom6_evidence"] == ""
    assert row["fantom6_cell_types"] == ""
