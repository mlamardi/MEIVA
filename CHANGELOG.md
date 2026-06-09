# Changelog

All notable changes to MEIVA are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- GENCODE GTF loader (`meiva.annotate.gencode`): parses a GRCh38 GENCODE GTF
  (plain or gzipped) into the gene model, recording per-transcript CDS span and
  the MANE Select tag; `IndexedGeneModel` provides a binned index for
  genome-scale lookups. Genic classification now records whether the reported
  transcript is MANE Select.
- Layer 1 genic context (`meiva.annotate.genic`): classifies a site as CDS / UTR /
  non-coding exon / splice-region / intron / promoter / upstream / downstream /
  intergenic against a pluggable `GeneModel`, strand-aware for UTR sides and for
  the promoter/upstream side, reporting the insertion's orientation
  (sense/antisense) relative to the gene and the distance to the nearest gene.
  Includes `InMemoryGeneModel`.
- Cohort merger (`meiva.cohort`): collapses per-sample discovery calls into a
  unified cohort site set with a sparse genotype matrix, using seed-anchored
  interval clustering (no single-linkage chaining) and capturing cross-sample
  breakpoint jitter as the merged site's imprecision interval. Reports a
  discovery-based `carrier_frequency` (explicitly not a true allele frequency).
- Canonical data model (`MEISite`, `SampleGenotype`, `MEIFamily`, `Strand`):
  interval-aware breakpoints, site-level (not per-call) identity, and a
  caller-estimated allele-frequency field.
- Ingest layer (`meiva.io`): a parser ABC with shared cyvcf2 plumbing, the
  xTEA parser (including a guard against the transduction `SVLEN` that encodes
  a genomic coordinate rather than an element length), and `detect_parser`
  header sniffing.
- Project scaffold: hatchling packaging, ruff (lint + format), strict mypy on
  `src`, pytest, GitHub Actions CI (Python 3.10–3.13), and pre-commit hooks.

[Unreleased]: https://github.com/TODO/meiva/commits/main
