"""Run the MEIVA Layer-1 pipeline: merge per-sample MEI VCFs, annotate against
GENCODE, and write an annotated cohort TSV.

Usage:
    python examples/annotate_cohort.py OUT.tsv GENCODE.gtf.gz SAMPLE1.vcf [SAMPLE2.vcf ...]

Notes:
    * GENCODE.gtf.gz is the GRCh38 GTF (e.g. gencode.v46.annotation.gtf.gz).
    * The output is a Layer-1 preview: cohort genotype summaries + genic context,
      not yet consequence tiers or population frequencies.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from meiva.pipeline import run


def main() -> None:
    ap = argparse.ArgumentParser(description="Merge + genic-annotate MEI VCFs into a TSV.")
    ap.add_argument("out_tsv", type=Path, help="output TSV path")
    ap.add_argument("gencode_gtf", type=Path, help="GRCh38 GENCODE GTF (plain or .gz)")
    ap.add_argument("vcfs", type=Path, nargs="+", help="per-sample MEI VCFs (xTEA, ...)")
    ap.add_argument(
        "--window", type=int, default=None, help="breakpoint merge window in bp (default 50)"
    )
    args = ap.parse_args()

    n = run(args.vcfs, args.gencode_gtf, args.out_tsv, window=args.window)
    print(f"wrote {n} annotated cohort sites to {args.out_tsv}")


if __name__ == "__main__":
    main()
