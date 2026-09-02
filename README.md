# MEIVA — Mobile Element Insertion Variant Annotator

[![PyPI](https://img.shields.io/pypi/v/meiva.svg)](https://pypi.org/project/meiva/)
[![Python](https://img.shields.io/pypi/pyversions/meiva.svg)](https://pypi.org/project/meiva/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**VEP for mobile elements.** A caller-agnostic engine that turns raw MEI calls
(Alu, LINE-1, SVA) from tools like **xTEA** and **MELT** into *interpreted*,
*prioritised* variants, with element and orientation-aware functional
consequences, population frequencies, regulatory context, and disease
knowledge.

> Status: **alpha, under active development.** The pipeline runs end to end for
> xTEA cohorts: parse, merge across samples, annotate against GENCODE, and emit a
> per-locus table carrying genic context, MEI-aware consequence calls, FANTOM5
> regulatory overlap and FANTOM6 lncRNA functional evidence. A second caller
> (MELT) and true population frequencies by force-genotyping are the next layers.

## Install

```bash
pip install meiva
```

Or from a checkout, for development:

```bash
git clone https://github.com/mlamardi/MEIVA.git
cd MEIVA
pip install -e ".[dev]"
```

Requires Python 3.10 or newer. The only runtime dependency is `cyvcf2`.
A bioconda recipe is planned.

## Usage

MEIVA installs a `meiva` command with two subcommands (both also work as
`python -m meiva ...`).

Annotate a cohort of caller VCFs against a GENCODE GTF, writing an annotated
per-locus table:

```bash
meiva annotate --vcf sample1.vcf sample2.vcf \
  --gencode gencode.v47.annotation.gtf.gz \
  -o cohort.annotated.tsv
```

Omit `-o` to stream the table to stdout so it composes with shell pipelines.

Optionally add FANTOM6 lncRNA functional evidence (all three files together):

```bash
meiva annotate --vcf sample1.vcf sample2.vcf \
  --gencode gencode.v47.annotation.gtf.gz \
  --fantom6-degs DESeq2_genes_ASO_signif.tsv.bz2 \
  --fantom6-samples Published_sample_summary.tsv.bz2 \
  --fantom6-cat FANTOM_CAT.lv3_robust.info_table.ID_mapping.tsv.gz \
  -o cohort.annotated.tsv
```

Add FANTOM5 regulatory context (either track may be given on its own):

```bash
meiva annotate --vcf sample1.vcf sample2.vcf \
  --gencode gencode.v47.annotation.gtf.gz \
  --fantom5-enhancers F5.hg38.enhancers.bed.gz \
  --fantom5-peaks hg38_fair+new_CAGE_peaks_phase1and2.bed.gz \
  --fantom5-peak-names human_phase1and2_CAGE_Peak_name.txt.gz \
  -o cohort.annotated.tsv
```

### Large cohorts

`meiva annotate` holds every parsed call in memory while merging, so a cohort of
a few thousand genomes can need tens of gigabytes. `scripts/run_meiva_cohort.sh`
splits the work by chromosome, which is exactly equivalent (MEIVA only clusters
calls sharing a contig and family, so no locus can span chromosomes) and brings
peak memory under 2 GB per worker:

```bash
./scripts/run_meiva_cohort.sh \
  --vcf-dir /path/to/vcfs \
  --gencode gencode.v47.annotation.gtf.gz \
  --outdir results \
  --workers 4
```

It is resumable, takes the same optional `--fantom5-*` and `--fantom6-*` inputs,
and `--help` lists every option.

Inspect a single caller's output as normalized `MEISite` records, without
annotation (handy for checking that a new caller parses correctly):

```bash
meiva parse --vcf sample1.vcf | head
```

## Why this exists

MEI *calling* is a crowded, mature space (MELT, xTEA, Mobster, TEBreak, …).
What's missing is everything *downstream*: once you have calls across a cohort,
there is no standard, maintained tool to annotate and prioritise them. People
hand-roll scripts. MEIVA is that missing layer.

The intellectual core is **not** the annotation plumbing — it is the
**MEI-aware consequence model**. General annotators (VEP, SnpEff) treat a
variant as a point. An MEI is not a point: its functional impact depends on the
element family, its orientation, its length, its poly-A tail, and the genic
context it lands in. A sense-oriented full-length L1 in an intron does
something mechanistically different from an AluY in a 3′UTR. MEIVA models that.

## Design principles

- **The parser normalises; the model only validates.** A `MEISite` that exists
  is, by construction, trustworthy. Malformed input fails loudly at ingest, not
  silently three layers later.
- **Sites, not points.** Breakpoints are imprecise; everything downstream
  matches on intervals, never bare coordinates.
- **Caller-agnostic.** One canonical record; per-caller parsers behind it.
- **Reproducible references.** Reference data is versioned and built via a
  cache manager, not bundled and forgotten.

## Status by layer

- **Ingest.** xTEA parser behind a caller-agnostic parser interface. (MELT next.)
- **Cohort merge.** Interval-based merge across samples with breakpoint-jitter
  handling and discovery-based carrier frequencies.
- **Layer 1, genic context.** CDS, UTR, non-coding exon, splice, intron,
  promoter, up/downstream, intergenic; strand-aware, with sense/antisense
  orientation, MANE-Select flag, gene biotype, and nearest-gene distance,
  annotated against a pinned GENCODE release (v47).
- **Layer 2, consequence model.** VEP-style mechanistic consequence term plus an
  ordinal impact tier, derived from element family, orientation, region, length,
  and biotype.

- **Layer 4, lncRNA functional evidence.** FANTOM6 knockdown phenotypes joined onto
  the host gene by Ensembl ID, tiered by how many independent ASOs responded.
- **Layer 4, regulatory context.** FANTOM5 transcribed enhancers and CAGE-defined
  promoters, reported as both exact overlap and distance to the nearest element.

Planned:
- **MELT parser,** to exercise and lock the caller-agnostic design.
- **Layer 3, population frequency** by force-genotyping discovered sites back
  against the samples, for true allele frequencies rather than discovery counts.

## License

MIT, see [LICENSE](LICENSE).
