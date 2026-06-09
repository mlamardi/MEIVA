"""Minimal example: parse an MEI VCF and print a summary.

Auto-detects the caller from the file header, parses every record into a
`MEISite`, and prints a short summary including any records whose reported
length was flagged as unreliable (e.g. xTEA transductions).

Usage:
    python examples/parse_xtea.py path/to/sample.vcf
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from meiva.io import detect_parser


def main() -> None:
    ap = argparse.ArgumentParser(description="Parse an MEI VCF and print a summary.")
    ap.add_argument("vcf", type=Path, help="path to an xTEA (or other supported) VCF")
    args = ap.parse_args()

    parser = detect_parser(args.vcf)
    sites = list(parser.parse(args.vcf))

    families = Counter(site.family.value for site in sites)
    flagged = [s for s in sites if s.raw_info.get("MEIVA_LENGTH_UNRELIABLE")]

    print(f"caller:          {parser.caller}")
    print(f"sites parsed:    {len(sites)}")
    print(f"by family:       {dict(families)}")
    print(f"length-flagged:  {len(flagged)}")
    for site in flagged:
        print(f"  {site}  (raw SVLEN={site.raw_info.get('SVLEN')})")


if __name__ == "__main__":
    main()
