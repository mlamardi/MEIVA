# Changelog

All notable changes to MEIVA are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- FANTOM6 evidence wired into the pipeline. `annotate_cohort`, `annotate_vcfs` and
  `run` take an optional `fantom6` mapping of unversioned Ensembl gene ID to
  knockdown evidence; it is joined on the reported gene's ID (version-stripped via
  the new `base_gene_id`), never on symbol, and left absent when the gene was not
  tested. Two new TSV columns: `fantom6_evidence` (the tier) and
  `fantom6_cell_types` -- the cell type is carried because fibroblast and iPSC
  phenotypes agree poorly, so a bare "functional" would overclaim. The CLI grows
  `--fantom6-degs`, `--fantom6-samples` and `--fantom6-cat`, which must be supplied
  together, and reports the Ensembl mapping rate on stderr.
- FANTOM6 to Ensembl crosswalk (`load_cat_gene_ids`, `target_to_ensembl`,
  `evidence_by_ensembl`). FANTOM6 numbers targets `G0<digits>` where the digits are
  the Ensembl accession (`G0214548` is MEG3, `ENSG00000214548`), but the convention is
  inferred and not universal: `G0277925` is TERC (really `ENSG00000270141`),
  `G0278144` is a FANTOM-specific `NEAT1_1` model, and `G0223811` is CAT-novel.
  A derived ID is therefore accepted only if it exists in the FANTOM CAT gene
  universe, read from `FANTOM_CAT.lv3_robust.info_table.ID_mapping.tsv.gz`; the traps
  fall outside it and resolve to `None`. On the real release this maps 146 of 154
  targets (94.8%), leaving 53 concordant, Ensembl-mappable lncRNAs as the usable
  functional-evidence set. Avoids both gene-symbol drift and the hg19-to-hg38 liftover
  FANTOM CAT would otherwise require.
- FANTOM6 lncRNA functional-evidence loader (`meiva.annotate.fantom6`). Reads the
  FANTOM6 significant DEG sign matrix plus the published sample summary and derives,
  per knocked-down lncRNA, an evidence tier (`concordant` = two or more independent
  ASOs produced a molecular response, `single_aso`, `no_response`) alongside DEG
  counts, up/down splits, symbol, and cell types. Joins ASO columns to targets via
  `perturb_id` rather than by parsing column names, because an ASO's identifier does
  not always match its target (`ASO_C013368_02` targets `G0253161`); an unresolvable
  column is a hard error rather than a silent drop. Evidence is keyed by FANTOM6
  target ID; see the crosswalk entry above for resolution to Ensembl gene IDs.
- Command-line interface (`meiva`, stdlib argparse; also `python -m meiva`).
  `meiva annotate --vcf ... --gencode ... [-o out.tsv]` runs the full pipeline and
  writes the annotated table (stdout if `-o` is omitted, so it composes with
  shell pipelines). `meiva parse --vcf ...` normalizes a single caller VCF into a
  TSV of canonical MEISite records for inspection. Registered as a console entry
  point in `pyproject.toml`. Broken-pipe (e.g. `| head`) is handled quietly like a
  well-behaved Unix filter. Extracted `annotate_vcfs` as the I/O-free core of
  `run` so the CLI can stream results without duplicating orchestration.
- Layer 2: MEI-aware consequence model (`meiva.annotate.consequence`). Maps each
  site to a mechanistic consequence term (`coding_disruption`, `splice_disruption`,
  `exonization_candidate`, `polyA_interference`, `promoter_insertion`,
  `utr5_insertion`/`utr3_insertion`, `noncoding_exon_insertion`, `intronic`,
  `upstream_gene`/`downstream_gene`/`intergenic`) and a derived ordinal impact
  tier (HIGH/MODERATE/LOW/MODIFIER), VEP-style. Region sets the ceiling; family,
  orientation, length and biotype modify it. Orientation- and length-dependent
  calls (antisense-Alu exonization, sense full-length-L1 polyA interference)
  degrade gracefully to the generic term with an `orientation_unknown` /
  `length_unknown` flag when xTEA leaves the strand or length undetermined; SVA in
  a regulatory context is flagged `sva_regulatory`. Emitted as new `consequence`,
  `impact`, and `consequence_flags` columns in the output TSV.
- Gene biotype capture: the GENCODE loader now reads each gene's `gene_type`
  (e.g. `protein_coding`, `lncRNA`), carries it on the `Gene` model and the
  `GenicContext`, and emits it as a new `gene_biotype` column in the output TSV.
  This lets lncRNA insertions be distinguished from protein-coding ones — a
  prerequisite for the lncRNA-focused analyses and the Layer 2 consequence model.
- Committed public test fixtures: two openly consented HGDP Tuscan (TSI) xTEA
  VCFs (`HGDP01162`, `HGDP01167`), trimmed to chr1, under `tests/data/`. The
  parser, cohort-merge, and end-to-end pipeline tests now run against these in
  CI instead of skipping on the gitignored VDA cohort data, so the suite is
  green on a fresh clone. The orphan-transduction length guard is now covered by
  a deterministic synthetic-VCF test rather than relying on a fixture record.
- End-to-end Layer-1 pipeline (`meiva.pipeline`): `run()` merges per-sample VCFs,
  annotates each cohort site against a GENCODE model, and writes an annotated TSV
  (cohort genotype summaries + genic context). `annotate_cohort`/`write_tsv` are
  pure and unit-tested; `examples/annotate_cohort.py` exposes it on the command line.
- Reference-data cache (`meiva.cache`): `CacheManager` downloads, checksums,
  and versions external reference files behind a JSON manifest, with atomic
  installs and an injectable fetcher (stdlib `urllib` by default). `gencode_resource`
  pins a GENCODE release (default v47, GRCh38). Cache location honours
  `MEIVA_CACHE_DIR` / `XDG_CACHE_HOME`.
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
