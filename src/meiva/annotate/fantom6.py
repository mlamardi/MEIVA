"""FANTOM6 lncRNA functional evidence (Layer 4, part one).

FANTOM6 is *not* a coordinate track like FANTOM5. It is a perturbation screen:
antisense oligos (ASOs) knocked down individual lncRNAs, and the molecular
response was measured by CAGE. What we extract is therefore a **gene-level
functional-evidence flag** -- "is this lncRNA experimentally characterized as
having a molecular phenotype?" -- rather than an interval overlap.

Two inputs, both from the FANTOM6 ``Core_FANTOM6/RELEASE_latest`` tree:

``analysis/DEGs/02_significant/DESeq2_genes_ASO_signif.tsv``
    A sign matrix. Rows are *responding* genes (Ensembl IDs); columns are ASO
    experiments; cells are ``-1`` (down), ``0`` (unchanged), ``1`` (up), already
    filtered to significant calls. Note the rows are the downstream genes, *not*
    the knocked-down lncRNA -- the target is identified by the column.

``metadata/Published_sample_summary.tsv``
    Maps each ``perturb_id`` (an ASO, matching a matrix column name) to the
    ``target_id`` it knocked down, plus that target's symbol and cell type.

**Join on ``perturb_id``, never by parsing the column name.** The identifier
embedded in an ASO's name is not always its target: ``ASO_C013368_02`` targets
``G0253161``. Parsing the string silently mis-assigns such experiments; the
metadata join resolves every column.

Evidence is tiered by how many independent ASOs against the same target produced
a response, mirroring the consortium's own requirement of at least two successful
knockdowns before trusting a phenotype:

* :attr:`EvidenceTier.CONCORDANT` -- two or more ASOs responded. Trustworthy.
* :attr:`EvidenceTier.SINGLE_ASO` -- exactly one responded. Weak; an off-target
  effect of that single oligo cannot be excluded.
* :attr:`EvidenceTier.NO_RESPONSE` -- the lncRNA was tested and nothing moved.
  This is informative: it is evidence of *absence*, not absence of evidence.

The evidence is keyed by FANTOM6 target ID (``G0...`` / ``C0...``). Mapping those
onto the GENCODE gene IDs MEIVA annotates against is handled by
:func:`target_to_ensembl`, which requires the FANTOM CAT ID-mapping table. Do not
join by gene symbol: roughly half of the target symbols are clone-based names
(``RP11-834C11.4``) that GENCODE renames between releases.
"""

from __future__ import annotations

import bz2
import csv
import gzip
from collections import defaultdict
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import TextIO

__all__ = [
    "DEFAULT_MIN_DEGS",
    "EvidenceTier",
    "Fantom6Evidence",
    "evidence_by_ensembl",
    "load_cat_gene_ids",
    "load_fantom6_evidence",
    "target_to_ensembl",
]

# An ASO is counted as "responding" when its knockdown moved at least this many
# genes. The FANTOM6 significant matrix is already thresholded for significance,
# so this guards only against a handful of incidental calls being read as a
# phenotype. Tunable: it is the one free parameter in this module.
DEFAULT_MIN_DEGS = 10

# Column names in Published_sample_summary.tsv.
_PERTURB_ID = "general_sample_info.perturb_id"
_TARGET_ID = "general_sample_info.target_id"
_TARGET_SYMBOL = "general_sample_info.target_gene_symbol"
_CELL_TYPE = "general_sample_info.cell_type_alias"


class EvidenceTier(Enum):
    """How strongly a lncRNA is supported as functional by the FANTOM6 screen."""

    CONCORDANT = "concordant"
    SINGLE_ASO = "single_aso"
    NO_RESPONSE = "no_response"


@dataclass(frozen=True, slots=True)
class Fantom6Evidence:
    """Aggregated knockdown evidence for one target lncRNA."""

    target_id: str
    target_symbol: str
    tier: EvidenceTier
    n_aso: int
    n_aso_responding: int
    max_degs: int
    total_degs: int
    n_up: int
    n_down: int
    cell_types: tuple[str, ...] = ()
    ensembl_gene_id: str | None = None
    """Unversioned Ensembl gene ID, when the target resolves to one (see
    :func:`target_to_ensembl`). ``None`` for FANTOM-CAT-novel or FANTOM-specific
    gene models that have no Ensembl counterpart."""

    @property
    def is_functional(self) -> bool:
        """True only for the concordant tier -- the one safe to treat as evidence."""
        return self.tier is EvidenceTier.CONCORDANT


@contextmanager
def _open_text(path: str | Path) -> Iterator[TextIO]:
    """Open plain, gzip, or bzip2 text transparently, chosen by suffix."""
    p = Path(path)
    opener: Callable[..., TextIO]
    if p.suffix == ".bz2":
        opener = bz2.open
    elif p.suffix == ".gz":
        opener = gzip.open
    else:
        opener = open
    with opener(p, "rt", encoding="utf-8") as fh:
        yield fh


def _load_crosswalk(
    sample_summary: str | Path,
) -> tuple[dict[str, str], dict[str, str], dict[str, set[str]]]:
    """Read the sample summary into perturb->target, target->symbol, target->cell types."""
    pid_to_target: dict[str, str] = {}
    target_to_symbol: dict[str, str] = {}
    target_to_cells: dict[str, set[str]] = defaultdict(set)

    with _open_text(sample_summary) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id or sample_id.startswith("#"):
                continue  # the column-type line that follows the header
            perturb = (row.get(_PERTURB_ID) or "").strip()
            target = (row.get(_TARGET_ID) or "").strip()
            if not perturb or not target:
                continue  # controls and references carry no target
            pid_to_target[perturb] = target
            target_to_symbol[target] = (row.get(_TARGET_SYMBOL) or "").strip()
            cell = (row.get(_CELL_TYPE) or "").strip()
            if cell:
                target_to_cells[target].add(cell)

    return pid_to_target, target_to_symbol, target_to_cells


def _tally_degs(deg_matrix: str | Path) -> tuple[list[str], dict[str, tuple[int, int]]]:
    """Count (up, down) significant genes for each ASO column of the sign matrix."""
    with _open_text(deg_matrix) as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"{deg_matrix}: empty DEG matrix") from exc
        if len(header) < 3 or header[0] != "geneID":
            raise ValueError(f"{deg_matrix}: expected a 'geneID' column plus ASO columns")

        aso_cols = header[2:]  # geneID, geneSymbol, then one column per ASO
        up = [0] * len(aso_cols)
        down = [0] * len(aso_cols)
        for row in reader:
            for i, value in enumerate(row[2:]):
                if value == "1":
                    up[i] += 1
                elif value == "-1":
                    down[i] += 1

    return aso_cols, {c: (up[i], down[i]) for i, c in enumerate(aso_cols)}


def load_cat_gene_ids(id_mapping: str | Path) -> frozenset[str]:
    """Read the Ensembl gene IDs present in the FANTOM CAT ID-mapping table.

    ``FANTOM_CAT.lv3_robust.info_table.ID_mapping.tsv.gz`` has a ``geneID`` column
    holding either a versioned Ensembl gene ID (for genes shared with GENCODE) or a
    ``CATG...`` identifier (for CAT-novel genes). Returns the Ensembl IDs with their
    version suffix stripped; these are what :func:`target_to_ensembl` validates against.
    """
    ids: set[str] = set()
    with _open_text(id_mapping) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None or header[0] != "geneID":
            raise ValueError(f"{id_mapping}: expected a 'geneID' first column")
        for row in reader:
            if row and row[0].startswith("ENSG"):
                ids.add(row[0].split(".", 1)[0])
    return frozenset(ids)


def target_to_ensembl(target_id: str, cat_gene_ids: frozenset[str]) -> str | None:
    """Resolve a FANTOM6 target ID to an unversioned Ensembl gene ID, or ``None``.

    FANTOM6 numbers its targets ``G0<digits>``, and for genes shared with Ensembl the
    digits are the Ensembl accession: ``G0214548`` is ``ENSG00000214548`` (MEG3). This
    convention is inferred rather than documented, so a derived ID is **only accepted
    when it exists in the FANTOM CAT gene universe** (``cat_gene_ids``).

    That check is not cosmetic. Several targets break the convention: ``G0277925``
    is TERC, whose real accession is ``ENSG00000270141``; ``G0278144`` is a
    FANTOM-specific ``NEAT1_1`` model; ``G0223811`` is a CAT-novel gene. Padding their
    digits yields plausible-looking but wrong Ensembl IDs. All of them fall outside the
    CAT gene set, so validation rejects them and they resolve to ``None``.

    ``C0...`` targets are CAT-novel lncRNAs with no Ensembl counterpart, and return
    ``None`` by construction.
    """
    if not target_id.startswith("G0"):
        return None
    digits = target_id[1:]
    if not digits.isdigit():
        return None
    candidate = "ENSG" + digits.zfill(11)
    return candidate if candidate in cat_gene_ids else None


def evidence_by_ensembl(
    evidence: dict[str, Fantom6Evidence],
) -> dict[str, Fantom6Evidence]:
    """Re-key evidence by unversioned Ensembl gene ID, dropping unmappable targets.

    This is the form MEIVA joins against, since gene models carry Ensembl IDs.
    """
    return {e.ensembl_gene_id: e for e in evidence.values() if e.ensembl_gene_id is not None}


def load_fantom6_evidence(
    deg_matrix: str | Path,
    sample_summary: str | Path,
    *,
    cat_gene_ids: frozenset[str] | None = None,
    min_degs: int = DEFAULT_MIN_DEGS,
) -> dict[str, Fantom6Evidence]:
    """Derive per-lncRNA functional evidence, keyed by FANTOM6 target ID.

    When ``cat_gene_ids`` is supplied (see :func:`load_cat_gene_ids`) each target is
    also resolved to an Ensembl gene ID where possible.

    Raises ``ValueError`` if any ASO column cannot be resolved to a target via the
    sample summary -- a silent mis-join here would quietly corrupt the evidence,
    so it is treated as a hard error rather than dropped.
    """
    pid_to_target, target_to_symbol, target_to_cells = _load_crosswalk(sample_summary)
    aso_cols, counts = _tally_degs(deg_matrix)

    unresolved = [c for c in aso_cols if c not in pid_to_target]
    if unresolved:
        raise ValueError(
            f"{len(unresolved)} ASO column(s) absent from the sample summary, "
            f"e.g. {unresolved[:3]}; the two files are mismatched releases"
        )

    by_target: dict[str, list[str]] = defaultdict(list)
    for col in aso_cols:
        by_target[pid_to_target[col]].append(col)

    evidence: dict[str, Fantom6Evidence] = {}
    for target, cols in by_target.items():
        totals = [counts[c][0] + counts[c][1] for c in cols]
        responding = sum(1 for t in totals if t >= min_degs)
        if responding >= 2:
            tier = EvidenceTier.CONCORDANT
        elif responding == 1:
            tier = EvidenceTier.SINGLE_ASO
        else:
            tier = EvidenceTier.NO_RESPONSE

        evidence[target] = Fantom6Evidence(
            target_id=target,
            target_symbol=target_to_symbol.get(target, ""),
            tier=tier,
            n_aso=len(cols),
            n_aso_responding=responding,
            max_degs=max(totals),
            total_degs=sum(totals),
            n_up=sum(counts[c][0] for c in cols),
            n_down=sum(counts[c][1] for c in cols),
            cell_types=tuple(sorted(target_to_cells.get(target, ()))),
            ensembl_gene_id=(
                target_to_ensembl(target, cat_gene_ids) if cat_gene_ids is not None else None
            ),
        )

    return evidence
