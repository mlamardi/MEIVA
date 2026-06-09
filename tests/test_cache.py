"""Tests for the reference-data cache manager.

The download is exercised two ways: a fake in-memory fetcher (for manifest /
versioning / integrity logic) and the real urllib_fetch against a file:// URL
(so the default fetch path has coverage without external network).
"""

from pathlib import Path

import pytest

from meiva.cache import CacheManager, Resource, default_cache_dir, gencode_resource, urllib_fetch


def _res(version: str = "v46") -> Resource:
    return Resource(
        name="gencode", version=version, url="http://example/x.gtf.gz", filename="x.gtf.gz"
    )


class CountingFetcher:
    """Writes fixed bytes and counts how many times it was invoked."""

    def __init__(self, payload: bytes = b"GENCODE\n") -> None:
        self.payload = payload
        self.calls = 0

    def __call__(self, url: str, dest: Path) -> None:
        self.calls += 1
        dest.write_bytes(self.payload)


# --------------------------------------------------------------------------- #
# Resource construction                                                       #
# --------------------------------------------------------------------------- #
def test_gencode_resource_default():
    r = gencode_resource()
    assert r.version == "v46"
    assert r.filename == "gencode.v46.annotation.gtf.gz"
    assert r.url.endswith("release_46/gencode.v46.annotation.gtf.gz")


def test_gencode_resource_custom_release_and_basic():
    r = gencode_resource(45, basic=True)
    assert r.version == "v45"
    assert r.filename == "gencode.v45.basic.annotation.gtf.gz"
    assert "release_45/" in r.url


# --------------------------------------------------------------------------- #
# build / install lifecycle                                                   #
# --------------------------------------------------------------------------- #
def test_build_installs_and_records(tmp_path):
    cm = CacheManager(tmp_path)
    fetcher = CountingFetcher()
    path = cm.build(_res(), fetcher=fetcher)
    assert path.exists() and path.read_bytes() == b"GENCODE\n"
    assert cm.is_installed("gencode")
    rec = cm.installed()["gencode"]
    assert rec["version"] == "v46"
    assert rec["bytes"] == len(b"GENCODE\n")
    assert len(rec["sha256"]) == 64


def test_build_is_idempotent(tmp_path):
    cm = CacheManager(tmp_path)
    fetcher = CountingFetcher()
    cm.build(_res(), fetcher=fetcher)
    cm.build(_res(), fetcher=fetcher)  # same version, already present
    assert fetcher.calls == 1


def test_force_redownloads(tmp_path):
    cm = CacheManager(tmp_path)
    fetcher = CountingFetcher()
    cm.build(_res(), fetcher=fetcher)
    cm.build(_res(), fetcher=fetcher, force=True)
    assert fetcher.calls == 2


def test_version_change_refetches(tmp_path):
    cm = CacheManager(tmp_path)
    fetcher = CountingFetcher()
    cm.build(_res("v45"), fetcher=fetcher)
    cm.build(_res("v46"), fetcher=fetcher)
    assert fetcher.calls == 2


def test_checksum_mismatch_aborts(tmp_path):
    cm = CacheManager(tmp_path)
    with pytest.raises(ValueError, match="checksum mismatch"):
        cm.build(_res(), fetcher=CountingFetcher(), expected_sha256="deadbeef")
    assert not cm.is_installed("gencode")
    assert not (tmp_path / "x.gtf.gz").exists()


def test_path_for_unknown_resource(tmp_path):
    with pytest.raises(KeyError):
        CacheManager(tmp_path).path_for("gencode")


def test_manifest_persists_across_instances(tmp_path):
    CacheManager(tmp_path).build(_res(), fetcher=CountingFetcher())
    assert CacheManager(tmp_path).is_installed("gencode")


# --------------------------------------------------------------------------- #
# Cache dir resolution                                                        #
# --------------------------------------------------------------------------- #
def test_default_cache_dir_env(monkeypatch, tmp_path):
    monkeypatch.setenv("MEIVA_CACHE_DIR", str(tmp_path / "mycache"))
    assert default_cache_dir() == tmp_path / "mycache"


# --------------------------------------------------------------------------- #
# Real fetcher against a file:// URL                                          #
# --------------------------------------------------------------------------- #
def test_urllib_fetch_file_url(tmp_path):
    src = tmp_path / "source.bin"
    src.write_bytes(b"hello reference data")
    dest = tmp_path / "out.bin"
    urllib_fetch(src.as_uri(), dest)
    assert dest.read_bytes() == b"hello reference data"


def test_build_with_real_fetcher_via_file_url(tmp_path):
    src = tmp_path / "gencode.v46.annotation.gtf.gz"
    src.write_bytes(b"\x1f\x8b\x08fake-gzip")
    resource = Resource(
        name="gencode", version="v46", url=src.as_uri(), filename="gencode.v46.annotation.gtf.gz"
    )
    cm = CacheManager(tmp_path / "cache")
    path = cm.build(resource)  # default urllib_fetch, real code path
    assert path.read_bytes() == b"\x1f\x8b\x08fake-gzip"
