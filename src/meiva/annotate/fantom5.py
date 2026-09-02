"""FANTOM5 regulatory context: transcribed enhancers and CAGE-defined promoters.

Where FANTOM6 asks whether a lncRNA *does* anything, FANTOM5 asks *where the
regulatory elements are*. It is a coordinate resource, so this is an interval
overlap -- but with one wrinkle that shapes the whole design.

**CAGE peaks are tiny.** Across the 209,911 human peaks the median width is 14 bp
(p90 43 bp, max 283 bp), totalling about 2.9 Mb, or 0.09% of the genome. A strict
"the insertion falls inside the peak" test would therefore almost never fire, and
would miss insertions sitting 200 bp from a transcription start site that are
plainly capable of disrupting it. Enhancers are wider (~300 bp, ~19 Mb) but the
same argument applies more weakly.

So every site gets **both** an exact-overlap call and a distance to the nearest
element. The module commits to no window: it reports the evidence and leaves the
threshold to the analysis, the same principle used for the FANTOM6 biotype
disagreement. ``DEFAULT_MAX_DISTANCE`` only bounds how far it bothers to look.

Inputs, all from the ``hg38_latest`` reprocessed release:

``F5.hg38.enhancers.bed.gz``
    BED12 transcribed enhancers. Unstranded, so distances are unsigned.

``hg38_fair+new_CAGE_peaks_phase1and2.bed.gz``
    BED9 CAGE peaks. Stranded, and ``thickStart`` holds the *representative TSS*,
    which is the biologically meaningful point -- distances are measured to it,
    not to the peak edge, and signed by the peak's strand so that negative means
    upstream of the TSS.

``human_phase1and2_CAGE_Peak_name.txt.gz`` (optional)
    Maps a peak to a ``p1@GENE`` label. The number is the promoter's rank for that
    gene, so ``p1`` is its dominant promoter -- an insertion there is a stronger
    claim than one in a minor alternative promoter.

Note the peak identifiers encode the original hg19 coordinates, which are shared
by the ``fair`` and ``liftover`` remappings. The name file can therefore be joined
to the recommended ``fair+new`` BED even though it ships beside the ``liftover``
one: all 209,911 identifiers match.

This layer is purely additive. It does not feed the Layer 2 consequence model,
because an insertion in an enhancer that also sits in an intron is still
intronic; folding regulatory overlap into the region hierarchy would need a
precedence rule we have not yet argued for.
"""

from __future__ import annotations

import bisect
import gzip
import re
from bisect import bisect_left, bisect_right
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

__all__ = [
    "DEFAULT_MAX_DISTANCE",
    "Fantom5Model",
    "RegulatoryContext",
    "load_fantom5",
]

# Beyond this the nearest element is not worth reporting; distance becomes None.
DEFAULT_MAX_DISTANCE = 100_000

# "p1@MTND1P23" -> rank 1, gene MTND1P23. Some entries carry no rank ("p@GENE").
_PEAK_NAME = re.compile(r"^p(\d*)@(.+)$")


@dataclass(frozen=True, slots=True)
class RegulatoryContext:
    """FANTOM5 evidence for one insertion site."""

    enhancer_id: str | None = None
    enhancer_distance: int | None = None
    """Unsigned distance to the nearest enhancer; 0 when inside one."""

    promoter_id: str | None = None
    """CAGE peak identifier, set only when the insertion is inside the peak span."""
    promoter_gene: str | None = None
    promoter_rank: int | None = None
    """Rank of this promoter for its gene; 1 is the gene's dominant promoter."""
    promoter_distance: int | None = None
    """Signed distance to the nearest peak's representative TSS. Negative is
    upstream of the TSS in the peak's own orientation."""

    @property
    def in_enhancer(self) -> bool:
        return self.enhancer_id is not None

    @property
    def in_promoter(self) -> bool:
        return self.promoter_id is not None


@contextmanager
def _open_text(path: str | Path) -> Iterator[TextIO]:
    p = Path(path)
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt", encoding="utf-8") as fh:
        yield fh


@dataclass(frozen=True, slots=True)
class _Track:
    """Per-contig sorted intervals, plus an anchor point used for distances."""

    starts: list[int]
    ends: list[int]
    anchors: list[int]
    names: list[str]
    strands: list[str]
    max_width: int


def _build(raw: dict[str, list[tuple[int, int, int, str, str]]]) -> dict[str, _Track]:
    out: dict[str, _Track] = {}
    for contig, items in raw.items():
        items.sort(key=lambda t: t[0])
        out[contig] = _Track(
            starts=[i[0] for i in items],
            ends=[i[1] for i in items],
            anchors=[i[2] for i in items],
            names=[i[3] for i in items],
            strands=[i[4] for i in items],
            # how far back to scan for a containing interval
            max_width=max((i[1] - i[0] for i in items), default=0),
        )
    return out


def _parse_bed(
    path: str | Path, *, anchor_col: int | None
) -> dict[str, list[tuple[int, int, int, str, str]]]:
    """Read a BED file. ``anchor_col`` picks thickStart as the distance anchor."""
    raw: dict[str, list[tuple[int, int, int, str, str]]] = {}
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                raise ValueError(f"{path}: expected at least 4 BED columns, got {len(f)}")
            start, end = int(f[1]), int(f[2])
            name = f[3]
            strand = f[5] if len(f) > 5 else "."
            if anchor_col is not None and len(f) > anchor_col:
                anchor = int(f[anchor_col])
            else:
                anchor = (start + end) // 2
            raw.setdefault(f[0], []).append((start, end, anchor, name, strand))
    if not raw:
        raise ValueError(f"{path}: no BED records found")
    return raw


def _load_peak_names(path: str | Path) -> dict[str, tuple[str, int | None]]:
    """Map a CAGE peak id to (gene symbol, promoter rank)."""
    out: dict[str, tuple[str, int | None]] = {}
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 or not f[1]:
                continue
            m = _PEAK_NAME.match(f[1])
            if m:
                out[f[0]] = (m.group(2), int(m.group(1)) if m.group(1) else None)
    return out


def _containing(track: _Track, pos0: int) -> int | None:
    """Index of an interval containing this 0-based position, if any."""
    i = bisect_right(track.starts, pos0) - 1
    floor = bisect_left(track.starts, pos0 - track.max_width) - 1
    while i > floor and i >= 0:
        if track.starts[i] <= pos0 < track.ends[i]:
            return i
        i -= 1
    return None


def _nearest_anchor(track: _Track, pos0: int) -> int | None:
    """Index of the interval whose anchor is closest to this position."""
    if not track.anchors:
        return None
    # anchors follow start order closely enough for a local search
    i = bisect.bisect_left(track.starts, pos0)
    best, best_d = None, None
    for j in range(max(0, i - 3), min(len(track.starts), i + 3)):
        d = abs(track.anchors[j] - pos0)
        if best_d is None or d < best_d:
            best, best_d = j, d
    return best


class Fantom5Model:
    """Queryable FANTOM5 enhancer and promoter tracks."""

    def __init__(
        self,
        enhancers: dict[str, _Track] | None = None,
        peaks: dict[str, _Track] | None = None,
        peak_names: dict[str, tuple[str, int | None]] | None = None,
        *,
        max_distance: int = DEFAULT_MAX_DISTANCE,
    ) -> None:
        self._enhancers = enhancers or {}
        self._peaks = peaks or {}
        self._peak_names = peak_names or {}
        self._max_distance = max_distance

    @property
    def n_enhancers(self) -> int:
        return sum(len(t.starts) for t in self._enhancers.values())

    @property
    def n_peaks(self) -> int:
        return sum(len(t.starts) for t in self._peaks.values())

    def annotate(self, chrom: str, pos: int) -> RegulatoryContext:
        """Annotate a 1-based position with enhancer and promoter evidence."""
        pos0 = pos - 1

        enh_id: str | None = None
        enh_dist: int | None = None
        track = self._enhancers.get(chrom)
        if track is not None:
            hit = _containing(track, pos0)
            if hit is not None:
                enh_id, enh_dist = track.names[hit], 0
            else:
                j = _nearest_anchor(track, pos0)
                if j is not None:
                    d = abs(track.anchors[j] - pos0)
                    enh_dist = d if d <= self._max_distance else None

        prom_id: str | None = None
        gene: str | None = None
        rank: int | None = None
        prom_dist: int | None = None
        track = self._peaks.get(chrom)
        if track is not None:
            hit = _containing(track, pos0)
            j = hit if hit is not None else _nearest_anchor(track, pos0)
            if j is not None:
                signed = pos0 - track.anchors[j]
                if track.strands[j] == "-":
                    signed = -signed  # measure in the peak's own orientation
                if abs(signed) <= self._max_distance:
                    prom_dist = signed
                if hit is not None:
                    prom_id = track.names[hit]
                    named = self._peak_names.get(prom_id)
                    if named:
                        gene, rank = named

        return RegulatoryContext(
            enhancer_id=enh_id,
            enhancer_distance=enh_dist,
            promoter_id=prom_id,
            promoter_gene=gene,
            promoter_rank=rank,
            promoter_distance=prom_dist,
        )


def load_fantom5(
    *,
    enhancers: str | Path | None = None,
    cage_peaks: str | Path | None = None,
    peak_names: str | Path | None = None,
    max_distance: int = DEFAULT_MAX_DISTANCE,
) -> Fantom5Model:
    """Build a :class:`Fantom5Model` from whichever tracks are supplied.

    ``peak_names`` is only meaningful alongside ``cage_peaks``.
    """
    if enhancers is None and cage_peaks is None:
        raise ValueError("supply at least one of enhancers or cage_peaks")
    if peak_names is not None and cage_peaks is None:
        raise ValueError("peak_names requires cage_peaks")

    enh = _build(_parse_bed(enhancers, anchor_col=None)) if enhancers is not None else None
    # column 6 (0-based) is thickStart: the representative TSS
    pk = _build(_parse_bed(cage_peaks, anchor_col=6)) if cage_peaks is not None else None
    names = _load_peak_names(peak_names) if peak_names is not None else None
    return Fantom5Model(enh, pk, names, max_distance=max_distance)
