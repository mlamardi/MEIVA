# Contributing to MEIVA

Thanks for your interest. This guide covers local setup, the checks your
changes must pass, and the few architectural rules that keep the codebase
coherent.

## Development setup

```bash
git clone https://github.com/TODO/meiva.git
cd meiva
python -m pip install -e ".[dev]"
pre-commit install
```

`cyvcf2` is a compiled dependency; if a wheel isn't available for your
platform it builds against htslib. On conda, `conda install -c bioconda
cyvcf2` first is the easy path.

## The checks (all must pass)

CI runs these on Python 3.10–3.13; run them locally before opening a PR:

```bash
ruff check .            # lint
ruff format --check .   # formatting
mypy                    # strict type-checking on src/
pytest                  # tests
```

`pre-commit` runs ruff and mypy automatically on commit.

## Architectural rules

These are not style preferences — they prevent whole classes of bugs:

1. **Parsers normalise; the model only validates.** All coercion (contig
   naming, family vocabulary, length sanity) lives in `meiva.io`.
   `MEISite.__post_init__` raises on bad input but never silently fixes it. A
   `MEISite` that exists is, by construction, trustworthy.
2. **An MEI is an interval, never a point.** Use `MEISite.search_interval()`
   for any overlap or cross-callset match — never compare bare `pos` values.
   Breakpoints jitter by a few bp between samples and callers.
3. **Identity is the locus, not the payload.** Equality/hashing key on
   chrom/pos/family/etc.; genotypes, qual, and `raw_info` are excluded. Two
   records for the same locus must compare equal.
4. **Never drop information silently.** Unparseable or out-of-scope fields go
   into `raw_info` (lossless) and, where relevant, get a `MEIVA_*` flag —
   they are not discarded.

## Adding a new caller parser

1. Subclass `Cyvcf2Parser` in a new `meiva/io/<caller>.py`.
2. Implement `_build_site(variant, samples) -> MEISite | None` (the
   caller-specific INFO mapping) and `sniff(path) -> bool` (header detection
   on an INFO ID unique to that caller).
3. Register the class in `PARSERS` in `meiva/io/__init__.py`.
4. Add tests under `tests/` that run against a real (or realistic,
   public-data-derived) VCF fixture.

## Tests & data

Test fixtures must be committable — i.e. derived from openly consented public
data (e.g. 1000G/HGDP), never real participant data. Logic that doesn't need a
VCF (model invariants, merge logic) should be unit-tested by constructing
`MEISite` objects directly.

## Commits & PRs

- Keep commits focused; write a clear message explaining the *why*.
- Update `CHANGELOG.md` under `[Unreleased]` for user-visible changes.
- A PR should leave `ruff`, `mypy`, and `pytest` green.
