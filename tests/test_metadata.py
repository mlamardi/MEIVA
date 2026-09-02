"""Guards against version drift between the package and its metadata.

``pyproject.toml`` reads the version from ``meiva.__version__`` via hatchling, so
those two cannot disagree. ``CITATION.cff`` is maintained by hand, so it can --
and a stale citation version is the kind of thing nobody notices until someone
cites the wrong release.
"""

import re
from pathlib import Path

import meiva

ROOT = Path(__file__).resolve().parents[1]


def test_citation_version_matches_package():
    text = (ROOT / "CITATION.cff").read_text()
    match = re.search(r"^version:\s*(\S+)\s*$", text, re.MULTILINE)
    assert match, "CITATION.cff has no version field"
    assert match.group(1) == meiva.__version__, (
        f"CITATION.cff says {match.group(1)}, package says {meiva.__version__}"
    )


def test_pyproject_takes_version_from_the_package():
    text = (ROOT / "pyproject.toml").read_text()
    assert 'dynamic = ["version"]' in text
    assert 'path = "src/meiva/__init__.py"' in text
    assert not re.search(r"^version\s*=", text, re.MULTILINE), (
        "pyproject.toml declares a static version; it should be dynamic"
    )


def test_no_placeholder_metadata_remains():
    for name in ("pyproject.toml", "CITATION.cff"):
        text = (ROOT / name).read_text()
        assert "TODO" not in text, f"{name} still contains a TODO placeholder"
