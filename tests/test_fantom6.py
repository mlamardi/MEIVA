"""Tests for the FANTOM6 lncRNA functional-evidence loader.

Fixtures mirror the real file formats, including the quirk that an ASO's own
identifier need not match the target it knocks down (``ASO_C013368_02`` targets
``G0253161`` in the real release) -- the regression that motivates joining on
``perturb_id`` rather than parsing column names.
"""

import bz2
import gzip

import pytest

from meiva.annotate.fantom6 import (
    EvidenceTier,
    evidence_by_ensembl,
    load_cat_gene_ids,
    load_fantom6_evidence,
    target_to_ensembl,
)

# Two ID columns, then five ASO columns:
#   G0000001: two ASOs, both responding      -> CONCORDANT
#   G0000002: two ASOs, one responding       -> SINGLE_ASO
#   G0000003: one ASO, no response           -> NO_RESPONSE
#   ASO_C013368_02 -> targets G0000001       (name/target mismatch, as in the real data)
_DEG = "\n".join(
    [
        "geneID\tgeneSymbol\tASO_G0000001_01\tASO_C013368_02\tASO_G0000002_01\t"
        "ASO_G0000002_02\tASO_G0000003_01",
        # 12 responding-gene rows; counts per column accumulate down the file
        *[f"ENSG{i:011d}\tSYM{i}\t1\t-1\t1\t0\t0" for i in range(12)],
    ]
)

_SUMMARY = "\n".join(
    [
        "sample_id\tgeneral_sample_info.cell_type_alias\tgeneral_sample_info.perturb_id\t"
        "general_sample_info.experiment_type\tgeneral_sample_info.target_gene_symbol\t"
        "general_sample_info.target_id",
        "#string\tstring\tstring\tstring\tstring\tstring",
        "S1\tHDF (neonatal)\tASO_G0000001_01\tTargeted\tMEG3\tG0000001",
        "S2\tHDF (neonatal)\tASO_C013368_02\tTargeted\tMEG3\tG0000001",
        "S3\tFF-iPSC (IMS)\tASO_G0000002_01\tTargeted\tRP11-1.1\tG0000002",
        "S4\tHDF (neonatal)\tASO_G0000002_02\tTargeted\tRP11-1.1\tG0000002",
        "S5\tHDF (neonatal)\tASO_G0000003_01\tTargeted\tDANCR\tG0000003",
        # a control row: no target, must be ignored rather than crash
        "S6\tHDF (neonatal)\tNC_A\tNegative control\t\t",
    ]
)


@pytest.fixture
def files(tmp_path):
    deg = tmp_path / "DESeq2_genes_ASO_signif.tsv"
    summary = tmp_path / "Published_sample_summary.tsv"
    deg.write_text(_DEG)
    summary.write_text(_SUMMARY)
    return deg, summary


# --------------------------------------------------------------------------- #
# Core derivation                                                              #
# --------------------------------------------------------------------------- #
def test_targets_are_keyed_by_cat_gene_id(files):
    ev = load_fantom6_evidence(*files)
    assert set(ev) == {"G0000001", "G0000002", "G0000003"}


def test_aso_name_mismatch_is_resolved_via_perturb_id(files):
    """ASO_C013368_02 must be attributed to G0000001, not to a 'C013368' target."""
    ev = load_fantom6_evidence(*files)
    assert "C013368" not in ev
    assert ev["G0000001"].n_aso == 2  # both ASOs collapsed onto the true target


def test_concordant_tier_requires_two_responding_asos(files):
    ev = load_fantom6_evidence(*files)
    g1 = ev["G0000001"]
    assert g1.tier is EvidenceTier.CONCORDANT
    assert g1.n_aso_responding == 2
    assert g1.is_functional is True


def test_single_aso_tier_is_not_functional(files):
    ev = load_fantom6_evidence(*files)
    g2 = ev["G0000002"]
    assert g2.tier is EvidenceTier.SINGLE_ASO
    assert g2.n_aso_responding == 1
    assert g2.is_functional is False


def test_no_response_tier(files):
    ev = load_fantom6_evidence(*files)
    g3 = ev["G0000003"]
    assert g3.tier is EvidenceTier.NO_RESPONSE
    assert g3.max_degs == 0
    assert g3.is_functional is False


def test_up_and_down_counted_separately(files):
    ev = load_fantom6_evidence(*files)
    g1 = ev["G0000001"]
    assert g1.n_up == 12  # ASO_G0000001_01 moved 12 genes up
    assert g1.n_down == 12  # ASO_C013368_02 moved 12 genes down
    assert g1.total_degs == 24


def test_symbol_and_cell_types_are_carried(files):
    ev = load_fantom6_evidence(*files)
    assert ev["G0000001"].target_symbol == "MEG3"
    assert ev["G0000001"].cell_types == ("HDF (neonatal)",)
    # a target tested in two cell types keeps both, sorted
    assert ev["G0000002"].cell_types == ("FF-iPSC (IMS)", "HDF (neonatal)")


# --------------------------------------------------------------------------- #
# Threshold behaviour                                                          #
# --------------------------------------------------------------------------- #
def test_min_degs_threshold_downgrades_tier(files):
    # with a threshold above every column's count, nothing responds
    ev = load_fantom6_evidence(*files, min_degs=100)
    assert all(e.tier is EvidenceTier.NO_RESPONSE for e in ev.values())


def test_lower_threshold_can_promote_tier(files):
    # G0000002's second ASO moved 0 genes, so it stays single even at threshold 1
    ev = load_fantom6_evidence(*files, min_degs=1)
    assert ev["G0000002"].tier is EvidenceTier.SINGLE_ASO


# --------------------------------------------------------------------------- #
# Robustness                                                                   #
# --------------------------------------------------------------------------- #
def test_controls_without_target_are_ignored(files):
    ev = load_fantom6_evidence(*files)
    assert all(e.target_id.startswith("G0") for e in ev.values())


def test_unresolvable_column_is_a_hard_error(tmp_path):
    deg = tmp_path / "deg.tsv"
    deg.write_text("geneID\tgeneSymbol\tASO_GHOST_09\nENSG1\tS\t1\n")
    summary = tmp_path / "s.tsv"
    summary.write_text(_SUMMARY)
    with pytest.raises(ValueError, match="absent from the sample summary"):
        load_fantom6_evidence(deg, summary)


def test_reads_bz2_transparently(tmp_path):
    deg = tmp_path / "deg.tsv.bz2"
    summary = tmp_path / "s.tsv.bz2"
    deg.write_bytes(bz2.compress(_DEG.encode()))
    summary.write_bytes(bz2.compress(_SUMMARY.encode()))
    ev = load_fantom6_evidence(deg, summary)
    assert ev["G0000001"].tier is EvidenceTier.CONCORDANT


def test_malformed_matrix_rejected(tmp_path):
    deg = tmp_path / "bad.tsv"
    deg.write_text("notGeneID\tx\ty\n")
    summary = tmp_path / "s.tsv"
    summary.write_text(_SUMMARY)
    with pytest.raises(ValueError, match="geneID"):
        load_fantom6_evidence(deg, summary)


def test_empty_matrix_rejected(tmp_path):
    deg = tmp_path / "empty.tsv"
    deg.write_text("")
    summary = tmp_path / "s.tsv"
    summary.write_text(_SUMMARY)
    with pytest.raises(ValueError, match="empty DEG matrix"):
        load_fantom6_evidence(deg, summary)


# --------------------------------------------------------------------------- #
# CAT -> Ensembl crosswalk                                                     #
# --------------------------------------------------------------------------- #
_ID_MAPPING = "\n".join(
    [
        "geneID\ttranscriptID\tCAGEClusterID",
        "ENSG00000214548.1\tENST0001.1\tchr14:100..200,+",  # MEG3, shared with Ensembl
        "ENSG00000238266.2\tENST0002.1\tchr10:100..200,-",  # LINC00707
        "CATG00000079799.1\tFTMT0001.1\tchr1:100..200,+",  # CAT-novel, not Ensembl
    ]
)


@pytest.fixture
def cat_ids(tmp_path):
    p = tmp_path / "ID_mapping.tsv"
    p.write_text(_ID_MAPPING)
    return load_cat_gene_ids(p)


def test_load_cat_gene_ids_keeps_only_unversioned_ensembl(cat_ids):
    assert cat_ids == frozenset({"ENSG00000214548", "ENSG00000238266"})


def test_load_cat_gene_ids_reads_gz(tmp_path):
    p = tmp_path / "ID_mapping.tsv.gz"
    p.write_bytes(gzip.compress(_ID_MAPPING.encode()))
    assert "ENSG00000214548" in load_cat_gene_ids(p)


def test_load_cat_gene_ids_rejects_bad_header(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("notGeneID\tx\n")
    with pytest.raises(ValueError, match="geneID"):
        load_cat_gene_ids(p)


def test_target_digits_map_to_ensembl(cat_ids):
    assert target_to_ensembl("G0214548", cat_ids) == "ENSG00000214548"  # MEG3
    assert target_to_ensembl("G0238266", cat_ids) == "ENSG00000238266"  # LINC00707


def test_target_not_in_cat_universe_is_rejected(cat_ids):
    """The digit convention is inferred, so an unvalidated ID must not be emitted.

    G0277925 is TERC, whose true accession is ENSG00000270141 -- padding its digits
    would fabricate ENSG00000277925. It is absent from CAT, so we return None.
    """
    assert target_to_ensembl("G0277925", cat_ids) is None  # TERC
    assert target_to_ensembl("G0278144", cat_ids) is None  # NEAT1_1, FANTOM-specific
    assert target_to_ensembl("G0223811", cat_ids) is None  # CAT-novel despite G0 prefix


def test_c0_targets_never_map(cat_ids):
    assert target_to_ensembl("C0008586", cat_ids) is None


def test_non_numeric_target_is_safe(cat_ids):
    assert target_to_ensembl("G0ABCDEF", cat_ids) is None
    assert target_to_ensembl("R0123456", cat_ids) is None


def test_evidence_carries_ensembl_id_when_cat_supplied(files, cat_ids):
    # G0000001 is not in the tiny CAT fixture, so it must stay unmapped
    ev = load_fantom6_evidence(*files, cat_gene_ids=cat_ids)
    assert ev["G0000001"].ensembl_gene_id is None


def test_evidence_by_ensembl_drops_unmappable(files, tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("geneID\ttranscriptID\tCAGEClusterID\nENSG00000000001.1\tT.1\tchr1:1..2,+\n")
    ev = load_fantom6_evidence(*files, cat_gene_ids=load_cat_gene_ids(p))
    by_gene = evidence_by_ensembl(ev)
    assert set(by_gene) == {"ENSG00000000001"}  # only G0000001 resolves
    assert by_gene["ENSG00000000001"].target_id == "G0000001"


def test_ensembl_id_absent_when_cat_not_supplied(files):
    ev = load_fantom6_evidence(*files)
    assert all(e.ensembl_gene_id is None for e in ev.values())
    assert evidence_by_ensembl(ev) == {}
