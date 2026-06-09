"""Reference-data cache.

A small, reproducible manager for the external reference files MEIVA depends
on — starting with the GENCODE GTF, later gnomAD-SV and FANTOM5. It downloads a
pinned version, verifies and records its checksum, and keeps a JSON manifest so
a run can state *exactly* which reference data it used. The annotation layers
then read paths from the cache rather than hardcoding them.

Design notes:

* **The download is injectable.** ``build`` takes a ``fetcher`` callable, so the
  manifest/versioning/integrity logic is testable without network, and the real
  :func:`urllib_fetch` (stdlib only — no ``requests`` dependency) is used by
  default.
* **Atomic installs.** Downloads land in a ``.tmp`` file that is checksum-checked
  and only then renamed into place, so an interrupted download never looks
  installed.
* **Pinned but configurable.** :func:`gencode_resource` defaults to v46 but takes
  any release, so the same machinery serves v33/v45/v46.
* **Relocatable cache.** Honours ``MEIVA_CACHE_DIR`` / ``XDG_CACHE_HOME``; HPC
  users should point it at project or scratch storage, not a quota'd home.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

__all__ = [
    "CacheManager",
    "Resource",
    "default_cache_dir",
    "gencode_resource",
    "urllib_fetch",
]

#: a fetcher downloads ``url`` to the local ``dest`` path
Fetcher = Callable[[str, Path], None]


@dataclass(frozen=True, slots=True)
class Resource:
    """A pinned, downloadable reference file."""

    name: str  # cache key, e.g. "gencode"
    version: str  # e.g. "v46"
    url: str
    filename: str  # local filename within the cache dir


def gencode_resource(release: int = 46, *, basic: bool = False) -> Resource:
    """A GENCODE (human, GRCh38) GTF resource for the given release.

    ``basic=True`` selects the smaller 'basic' annotation; the default is the
    comprehensive set.
    """
    flavour = "basic." if basic else ""
    filename = f"gencode.v{release}.{flavour}annotation.gtf.gz"
    url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{release}/{filename}"
    return Resource(name="gencode", version=f"v{release}", url=url, filename=filename)


def default_cache_dir() -> Path:
    """Resolve the cache directory: ``MEIVA_CACHE_DIR`` > ``XDG_CACHE_HOME`` > ~/.cache."""
    env = os.environ.get("MEIVA_CACHE_DIR")
    if env:
        return Path(env)
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg) if xdg else Path.home() / ".cache"
    return base / "meiva"


def urllib_fetch(url: str, dest: Path) -> None:
    """Default fetcher: stream ``url`` to ``dest`` using the standard library."""
    request: urllib.request.Request | str
    if url.startswith(("http://", "https://")):
        request = urllib.request.Request(url, headers={"User-Agent": "meiva-cache"})
    else:
        request = url  # e.g. file:// — used in tests
    with urllib.request.urlopen(request) as resp, open(dest, "wb") as out:
        shutil.copyfileobj(resp, out)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class CacheManager:
    """Manages a directory of versioned reference files plus a JSON manifest."""

    def __init__(self, cache_dir: str | Path | None = None) -> None:
        self.cache_dir = Path(cache_dir) if cache_dir is not None else default_cache_dir()

    @property
    def _manifest_path(self) -> Path:
        return self.cache_dir / "manifest.json"

    def _load_manifest(self) -> dict[str, Any]:
        if self._manifest_path.exists():
            data: dict[str, Any] = json.loads(self._manifest_path.read_text())
        else:
            data = {}
        data.setdefault("resources", {})
        return data

    def _save_manifest(self, data: dict[str, Any]) -> None:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._manifest_path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")

    def build(
        self,
        resource: Resource,
        *,
        fetcher: Fetcher = urllib_fetch,
        force: bool = False,
        expected_sha256: str | None = None,
    ) -> Path:
        """Download ``resource`` into the cache if needed; return its local path.

        Idempotent: a present file of the matching version is reused unless
        ``force``. If ``expected_sha256`` is given and the download doesn't match,
        nothing is installed and a :class:`ValueError` is raised.
        """
        manifest = self._load_manifest()
        entry = manifest["resources"].get(resource.name)
        dest = self.cache_dir / resource.filename

        if not force and entry and entry.get("version") == resource.version and dest.exists():
            return dest

        self.cache_dir.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_name(dest.name + ".tmp")
        fetcher(resource.url, tmp)
        digest = _sha256(tmp)
        if expected_sha256 is not None and digest.lower() != expected_sha256.lower():
            tmp.unlink(missing_ok=True)
            raise ValueError(
                f"checksum mismatch for {resource.name}: got {digest}, expected {expected_sha256}"
            )
        tmp.replace(dest)

        manifest["resources"][resource.name] = {
            "version": resource.version,
            "filename": resource.filename,
            "url": resource.url,
            "sha256": digest,
            "bytes": dest.stat().st_size,
            "installed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        self._save_manifest(manifest)
        return dest

    def path_for(self, name: str) -> Path:
        """Local path of an installed resource; raises if it isn't cached."""
        entry = self._load_manifest()["resources"].get(name)
        if entry is None:
            raise KeyError(f"resource {name!r} not in cache at {self.cache_dir}; build it first")
        path = self.cache_dir / str(entry["filename"])
        if not path.exists():
            raise FileNotFoundError(f"manifest lists {name!r} but file is missing: {path}")
        return path

    def is_installed(self, name: str) -> bool:
        try:
            self.path_for(name)
        except (KeyError, FileNotFoundError):
            return False
        return True

    def installed(self) -> dict[str, Any]:
        """A copy of the manifest's resource records."""
        return dict(self._load_manifest()["resources"])
