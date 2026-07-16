"""Command-line interface for MEIVA.

Two subcommands:

* ``meiva annotate`` runs the full pipeline (parse -> merge -> annotate -> TSV)
  over a cohort of caller VCFs against a GENCODE GTF.
* ``meiva parse`` normalizes a single caller VCF into a TSV of canonical
  ``MEISite`` records, without annotation -- useful for inspection and for
  demonstrating that different callers reduce to the same model.

Deliberately built on stdlib :mod:`argparse` so the only runtime dependency
stays ``cyvcf2``. Output goes to ``-o/--output`` or, if omitted, to stdout so it
composes with shell pipelines; progress lines go to stderr so they never
contaminate a piped TSV.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import TextIO

from meiva import __version__
from meiva.annotate.fantom6 import Fantom6Evidence
from meiva.model import MEISite

# Columns emitted by ``meiva parse`` (normalized site view, pre-annotation).
PARSE_HEADER = [
    "chrom",
    "pos",
    "ci_lower",
    "ci_upper",
    "family",
    "subfamily",
    "strand",
    "length",
    "tsd",
    "allele_freq",
    "n_carriers",
    "source_caller",
]


def _open_out(output: str | None) -> tuple[TextIO, bool]:
    """Resolve an output target to (stream, should_close). None/'-' -> stdout."""
    if output is None or output == "-":
        return sys.stdout, False
    return open(output, "w", newline="", encoding="utf-8"), True


def _parse_row(site: MEISite) -> list[str]:
    n_carriers = sum(1 for gt in site.genotypes.values() if gt.is_carrier)
    return [
        site.chrom,
        str(site.pos),
        str(site.ci_lower),
        str(site.ci_upper),
        site.family.value,
        site.subfamily or "",
        site.strand.value,
        "" if site.length is None else str(site.length),
        site.tsd or "",
        "" if site.allele_freq is None else f"{site.allele_freq:.4f}",
        str(n_carriers),
        site.source_caller or "",
    ]


def _load_fantom6(args: argparse.Namespace) -> dict[str, Fantom6Evidence] | None:
    """Build the Ensembl-keyed FANTOM6 evidence map, or None when not requested."""
    from meiva.annotate.fantom6 import (
        evidence_by_ensembl,
        load_cat_gene_ids,
        load_fantom6_evidence,
    )

    supplied = [args.fantom6_degs, args.fantom6_samples, args.fantom6_cat]
    if not any(supplied):
        return None
    if not all(supplied):
        raise ValueError(
            "--fantom6-degs, --fantom6-samples and --fantom6-cat must be given together"
        )
    cat = load_cat_gene_ids(args.fantom6_cat)
    evidence = load_fantom6_evidence(args.fantom6_degs, args.fantom6_samples, cat_gene_ids=cat)
    by_gene = evidence_by_ensembl(evidence)
    print(
        f"meiva: FANTOM6: {len(by_gene)} of {len(evidence)} tested lncRNAs mapped to Ensembl",
        file=sys.stderr,
    )
    return by_gene


def _cmd_annotate(args: argparse.Namespace) -> int:
    from meiva.pipeline import annotate_vcfs, write_tsv

    fantom6 = _load_fantom6(args)
    annotated = annotate_vcfs(args.vcf, args.gencode, window=args.merge_window, fantom6=fantom6)
    stream, close = _open_out(args.output)
    try:
        write_tsv(annotated, stream)
    finally:
        if close:
            stream.close()
    print(f"meiva: wrote {len(annotated)} annotated sites", file=sys.stderr)
    return 0


def _cmd_parse(args: argparse.Namespace) -> int:
    from meiva.io import detect_parser

    parser = detect_parser(args.vcf)  # raises ValueError if unrecognized
    sites = list(parser.parse(args.vcf))
    stream, close = _open_out(args.output)
    try:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(PARSE_HEADER)
        for site in sites:
            writer.writerow(_parse_row(site))
    finally:
        if close:
            stream.close()
    print(f"meiva: parsed {len(sites)} sites from {Path(args.vcf).name}", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser (separated out so tests can introspect it)."""
    parser = argparse.ArgumentParser(
        prog="meiva",
        description="MEIVA: caller-agnostic functional annotation of mobile element insertions.",
    )
    parser.add_argument("--version", action="version", version=f"meiva {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="<command>")

    ann = sub.add_parser("annotate", help="annotate a cohort of caller VCFs against GENCODE")
    ann.add_argument(
        "--vcf", nargs="+", required=True, metavar="VCF", help="one or more per-sample caller VCFs"
    )
    ann.add_argument(
        "--gencode", required=True, metavar="GTF", help="GENCODE annotation GTF (.gtf or .gtf.gz)"
    )
    ann.add_argument(
        "-o", "--output", metavar="TSV", default=None, help="output TSV path (default: stdout)"
    )
    ann.add_argument(
        "--merge-window",
        type=int,
        default=None,
        metavar="BP",
        help="cohort merge window in bp (default: 50)",
    )
    ann.add_argument(
        "--fantom6-degs",
        metavar="TSV",
        default=None,
        help="FANTOM6 DESeq2_genes_ASO_signif.tsv[.bz2] (requires the two flags below)",
    )
    ann.add_argument(
        "--fantom6-samples",
        metavar="TSV",
        default=None,
        help="FANTOM6 Published_sample_summary.tsv[.bz2]",
    )
    ann.add_argument(
        "--fantom6-cat",
        metavar="TSV",
        default=None,
        help="FANTOM_CAT.lv3_robust.info_table.ID_mapping.tsv[.gz], used to validate gene IDs",
    )

    par = sub.add_parser("parse", help="normalize a single caller VCF to a TSV (no annotation)")
    par.add_argument("--vcf", required=True, metavar="VCF", help="a single caller VCF")
    par.add_argument(
        "-o", "--output", metavar="TSV", default=None, help="output TSV path (default: stdout)"
    )

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point. Returns a process exit code (0 ok, 1 runtime error, 2 usage)."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help(sys.stderr)
        return 1

    try:
        if args.command == "annotate":
            return _cmd_annotate(args)
        if args.command == "parse":
            return _cmd_parse(args)
    except BrokenPipeError:
        # A downstream consumer closed the pipe (e.g. `meiva parse ... | head`).
        # Redirect stdout to devnull so the interpreter's final flush doesn't
        # raise again, then exit quietly like any well-behaved Unix filter.
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        return 0
    except FileNotFoundError as exc:
        print(f"meiva: file not found: {exc.filename}", file=sys.stderr)
        return 1
    except (ValueError, OSError) as exc:
        print(f"meiva: error: {exc}", file=sys.stderr)
        return 1

    return 1  # unreachable: argparse constrains command to the handled set


if __name__ == "__main__":
    sys.exit(main())
