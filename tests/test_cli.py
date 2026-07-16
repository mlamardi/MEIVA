"""Tests for the ``meiva`` command-line interface.

Driven in-process via ``cli.main([...])`` (fast, and captures the return code and
streams directly) rather than through a subprocess.
"""

from pathlib import Path

import pytest

from meiva import __version__
from meiva.cli import PARSE_HEADER, build_parser, main
from meiva.pipeline import TSV_HEADER

DATA = Path(__file__).parent / "data"
TSI = sorted(str(p) for p in DATA.glob("HGDP*.tsi.vcf"))

# A minimal GTF placing a lncRNA over a real TSI ALU (chr1:5,874,228) so the
# annotate path exercises a non-intergenic consequence end to end.
_GTF = (
    'chr1\tT\tgene\t5870001\t5880000\t.\t+\t.\tgene_id "G1"; gene_name "DEMO"; '
    'gene_type "lncRNA";\n'
    'chr1\tT\ttranscript\t5870001\t5880000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; '
    'gene_type "lncRNA";\n'
    'chr1\tT\texon\t5870001\t5880000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
)


@pytest.fixture
def gtf(tmp_path):
    p = tmp_path / "mini.gtf"
    p.write_text(_GTF)
    return str(p)


# --------------------------------------------------------------------------- #
# Top-level behaviour                                                          #
# --------------------------------------------------------------------------- #
def test_version_prints_and_exits_zero(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_no_command_prints_help_and_returns_one(capsys):
    assert main([]) == 1
    assert "usage" in capsys.readouterr().err.lower()


def test_build_parser_parses_annotate_args():
    args = build_parser().parse_args(["annotate", "--vcf", "a.vcf", "b.vcf", "--gencode", "g.gtf"])
    assert args.command == "annotate"
    assert args.vcf == ["a.vcf", "b.vcf"]
    assert args.gencode == "g.gtf"


# --------------------------------------------------------------------------- #
# annotate                                                                    #
# --------------------------------------------------------------------------- #
def test_annotate_to_file(gtf, tmp_path, capsys):
    out = tmp_path / "out.tsv"
    rc = main(["annotate", "--vcf", *TSI, "--gencode", gtf, "-o", str(out)])
    assert rc == 0
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == TSV_HEADER
    assert len(lines) - 1 == 54  # merged cohort sites
    assert "wrote 54 annotated sites" in capsys.readouterr().err


def test_annotate_to_stdout(gtf, capsys):
    rc = main(["annotate", "--vcf", *TSI, "--gencode", gtf])
    assert rc == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0].split("\t") == TSV_HEADER
    assert len(out) - 1 == 54
    # the lncRNA host over chr1:5,874,228 should yield a real consequence column
    assert any("noncoding_exon_insertion" in ln for ln in out)


def test_annotate_missing_gencode_is_usage_error():
    with pytest.raises(SystemExit) as exc:
        main(["annotate", "--vcf", *TSI])  # no --gencode
    assert exc.value.code == 2  # argparse usage error


def test_annotate_missing_vcf_file_returns_one(gtf, capsys):
    rc = main(["annotate", "--vcf", "/no/such/file.vcf", "--gencode", gtf])
    assert rc == 1
    assert "meiva:" in capsys.readouterr().err


# --------------------------------------------------------------------------- #
# parse                                                                       #
# --------------------------------------------------------------------------- #
def test_parse_to_stdout(capsys):
    rc = main(["parse", "--vcf", TSI[0]])  # HGDP01162 -> 27 records
    assert rc == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0].split("\t") == PARSE_HEADER
    assert len(out) - 1 == 27


def test_parse_to_file(tmp_path):
    out = tmp_path / "sites.tsv"
    rc = main(["parse", "--vcf", TSI[0], "-o", str(out)])
    assert rc == 0
    rows = out.read_text().splitlines()
    assert rows[0].split("\t") == PARSE_HEADER
    # every data row has one field per header column
    assert all(len(r.split("\t")) == len(PARSE_HEADER) for r in rows[1:])


def test_parse_missing_file_returns_one(capsys):
    rc = main(["parse", "--vcf", "/no/such/file.vcf"])
    assert rc == 1
    assert "meiva:" in capsys.readouterr().err


# --------------------------------------------------------------------------- #
# FANTOM6 options                                                              #
# --------------------------------------------------------------------------- #
def test_fantom6_flags_must_be_given_together(gtf, capsys):
    rc = main(["annotate", "--vcf", *TSI, "--gencode", gtf, "--fantom6-degs", "x.tsv"])
    assert rc == 1
    assert "must be given together" in capsys.readouterr().err


def test_annotate_without_fantom6_leaves_columns_blank(gtf, capsys):
    rc = main(["annotate", "--vcf", *TSI, "--gencode", gtf])
    assert rc == 0
    out = capsys.readouterr().out.splitlines()
    hdr = out[0].split("\t")
    assert "fantom6_evidence" in hdr
    idx = hdr.index("fantom6_evidence")
    assert all(ln.split("\t")[idx] == "" for ln in out[1:])
