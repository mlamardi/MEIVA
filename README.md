# MEIVA
Mobile Element Insertion Variant Annotator 
# MEIVA — Mobile Element Insertion Variant Annotator

**VEP for mobile elements.** A caller-agnostic engine that turns raw MEI calls
(Alu, LINE-1, SVA) from tools like **xTEA** and **MELT** into *interpreted*,
*prioritised* variants — with element- and orientation-aware functional
consequences, population frequencies, regulatory context, and disease
knowledge.

> Status: **pre-alpha, under active construction.** The data model and ingest
> layer are landing first. Not yet usable end-to-end.

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

## Roadmap

- **Phase 1 (MVP):** ingest (xTEA + MELT) → genic context → consequence tiers →
  population frequency (gnomAD-SV v4 + 1000G) → annotated VCF/TSV.
- **Phase 2:** regulatory overlap, curated disease-MEI database, gene-constraint
  prioritisation, HTML reports.
- **Phase 3:** more callers, long-read inputs, clinical mode, web front end.

## License

MIT — see [LICENSE](LICENSE).
