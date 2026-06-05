"""Ingest layer: caller-specific VCF parsers -> normalised MEISite objects.

Parsers OWN normalisation (contig naming, family vocabulary, length sanity);
the model only validates. Use :func:`detect_parser` to auto-select a parser
from a file's header, or instantiate a parser directly.
"""

from __future__ import annotations

from pathlib import Path

from meiva.io.base import Cyvcf2Parser, MEIParser
from meiva.io.xtea import XteaParser

#: registry of available parsers, tried in order during detection
PARSERS: tuple[type[MEIParser], ...] = (XteaParser,)

__all__ = ["PARSERS", "Cyvcf2Parser", "MEIParser", "XteaParser", "detect_parser"]


def detect_parser(path: str | Path) -> MEIParser:
    """Return an instance of the first parser that recognises ``path``.

    Raises :class:`ValueError` if no registered parser claims the file, so an
    unknown format fails loudly instead of being silently mis-parsed.
    """
    for parser_cls in PARSERS:
        if parser_cls.sniff(path):
            return parser_cls()
    raise ValueError(f"no registered MEI parser recognises {path!r}")
