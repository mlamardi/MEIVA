"""Tests for the cohort merger.

Built from in-memory MEISite objects so they need no VCF fixtures. The
merge_vcfs smoke test runs only if the (gitignored) real VCFs are present
locally, and is skipped in CI.
"""

from pathlib import Path

import pytest

from meiva.cohort import Cohort, merge_sites, merge_vcfs
from meiva.model import MEIFamily, MEISite, SampleGenotype, Strand


def site(
    chrom: str,
    pos: int,
    fam: MEIFamily,
    sample: str,
    alleles: tuple[int | None, ...],
    *,
    caller: str = "xtea",
    strand: Strand = Strand.UNKNOWN,
    length: int | None = None,
    tsd: str | None = None,
) -> MEISite:
    """Build a single-sample MEISite (mimicking one discovery call)."""
    return MEISite(
        chrom=chrom,
        pos=pos,
        family=fam,
        strand=strand,
        length=length,
        tsd=tsd,
        source_caller=caller,
        genotypes={sample: SampleGenotype(sample, alleles)},
    )


# --------------------------------------------------------------------------- #
# Grouping                                                                    #
# --------------------------------------------------------------------------- #
def test_jittered_breakpoints_merge_into_one_site():
    sites = [
        site("chr1", 100, MEIFamily.ALU, "A", (0, 1)),
        site("chr1", 101, MEIFamily.ALU, "B", (1, 1)),
        site("chr1", 102, MEIFamily.ALU, "C", (0, 1)),
    ]
    cohort = merge_sites(sites)
    assert cohort.n_sites == 1
    cs = cohort.sites[0]
    assert cs.site.pos == 101  # median_low of 100/101/102
    assert cs.site.ci_lower == -1 and cs.site.ci_upper == 1  # span captured as imprecision
    assert cs.n_carriers == 3
    assert set(cs.site.genotypes) == {"A", "B", "C"}
    assert cs.site.raw_info["MEIVA_N_MEMBERS"] == "3"


def test_different_families_do_not_merge():
    sites = [
        site("chr1", 200, MEIFamily.ALU, "A", (0, 1)),
        site("chr1", 200, MEIFamily.L1, "B", (0, 1)),
    ]
    assert merge_sites(sites).n_sites == 2


def test_different_contigs_do_not_merge():
    sites = [
        site("chr1", 300, MEIFamily.ALU, "A", (0, 1)),
        site("chr2", 300, MEIFamily.ALU, "B", (0, 1)),
    ]
    assert merge_sites(sites).n_sites == 2


def test_window_boundary_inclusive():
    sites = [
        site("chr1", 100, MEIFamily.ALU, "A", (0, 1)),
        site("chr1", 150, MEIFamily.ALU, "B", (0, 1)),  # exactly window away -> merges
        site("chr1", 151, MEIFamily.ALU, "C", (0, 1)),  # one past -> new cluster
    ]
    cohort = merge_sites(sites, window=50)
    assert cohort.n_sites == 2


def test_seed_anchored_no_chaining():
    # 100 & 140 cluster; 181 is 81 from the seed (>50) so it does NOT chain in,
    # even though it is only 41 from the previous member.
    sites = [
        site("chr1", 100, MEIFamily.ALU, "A", (0, 1)),
        site("chr1", 140, MEIFamily.ALU, "B", (0, 1)),
        site("chr1", 181, MEIFamily.ALU, "C", (0, 1)),
    ]
    cohort = merge_sites(sites, window=50)
    assert cohort.n_sites == 2
    assert [cs.n_carriers for cs in cohort.sites] == [2, 1]


# --------------------------------------------------------------------------- #
# Roster, frequency, absence-is-not-reference                                 #
# --------------------------------------------------------------------------- #
def test_roster_is_union_and_frequency_uses_cohort_size():
    sites = [
        site("chr1", 100, MEIFamily.ALU, "A", (0, 1)),
        site("chr1", 101, MEIFamily.ALU, "B", (0, 1)),
        site("chr1", 5000, MEIFamily.ALU, "A", (0, 1)),  # only A carries this one
    ]
    cohort = merge_sites(sites)
    assert cohort.roster == ("A", "B")
    shared = next(cs for cs in cohort.sites if cs.site.pos in (100, 101))
    singleton = next(cs for cs in cohort.sites if cs.site.pos == 5000)
    assert shared.carrier_frequency == 1.0  # 2 carriers / 2 samples
    assert singleton.carrier_frequency == 0.5  # 1 carrier / 2 samples
    # absence is sparse: B is simply not in the singleton's genotype matrix
    assert "B" not in singleton.site.genotypes


# --------------------------------------------------------------------------- #
# Consensus fields & collisions                                               #
# --------------------------------------------------------------------------- #
def test_strand_and_length_consensus():
    sites = [
        site("chr1", 100, MEIFamily.L1, "A", (0, 1), strand=Strand.PLUS, length=6000),
        site("chr1", 100, MEIFamily.L1, "B", (1, 1), strand=Strand.PLUS, length=5900),
        site("chr1", 100, MEIFamily.L1, "C", (0, 1), strand=Strand.UNKNOWN, length=None),
    ]
    cs = merge_sites(sites).sites[0]
    assert cs.site.strand is Strand.PLUS  # UNKNOWN ignored, PLUS is consensus
    assert cs.site.length == 5950 or cs.site.length in (5900, 6000)  # median_low of [5900,6000]


def test_same_sample_collision_is_flagged_and_resolved():
    # one sample contributes two nearby calls -> collision; higher dosage wins
    sites = [
        site("chr1", 100, MEIFamily.ALU, "A", (0, 1)),
        site("chr1", 102, MEIFamily.ALU, "A", (1, 1)),
    ]
    cohort = merge_sites(sites)
    cs = cohort.sites[0]
    assert "SAMPLE_COLLISION" in cs.flags
    assert cs.site.genotypes["A"].dosage == 2  # the 1/1 call kept
    assert cs.n_carriers == 1  # still one sample


def test_empty_input():
    cohort = merge_sites([])
    assert cohort == Cohort(roster=(), sites=())


def test_negative_window_rejected():
    with pytest.raises(ValueError):
        merge_sites([], window=-1)


# --------------------------------------------------------------------------- #
# merge_vcfs against real data (local only; gitignored fixtures)              #
# --------------------------------------------------------------------------- #
_REAL = sorted((Path(__file__).parent / "data").glob("ASL_VDA_*.vcf"))


@pytest.mark.skipif(len(_REAL) < 2, reason="real VCF fixtures not present (gitignored)")
def test_merge_vcfs_smoke():
    cohort = merge_vcfs(_REAL)
    assert cohort.n_samples == len(_REAL)
    assert cohort.n_sites > 0
    # every merged site references at least one carrier
    assert all(cs.n_carriers >= 1 for cs in cohort.sites)
