"""S3 / UPath integration tests against a local MinIO instance.

Run via ``test/run_s3_test.sh``, which starts MinIO and sets:

  BGEN_TEST_MINIO_ENDPOINT  e.g. http://127.0.0.1:19000
  BGEN_TEST_MINIO_ACCESS    access key
  BGEN_TEST_MINIO_SECRET    secret key
  BGEN_TEST_S3_METAFILE     path to the pre-built local .metafile

Tests are skipped when those env vars are absent so they never break a
normal ``pytest`` run without MinIO.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

MINIO_ENDPOINT = os.environ.get("BGEN_TEST_MINIO_ENDPOINT", "")
MINIO_ACCESS   = os.environ.get("BGEN_TEST_MINIO_ACCESS", "minio_access")
MINIO_SECRET   = os.environ.get("BGEN_TEST_MINIO_SECRET", "minio_secret")
METAFILE       = os.environ.get("BGEN_TEST_S3_METAFILE", "")
BUCKET_PRIVATE = "bgen-test"
BUCKET_PUBLIC  = "bgen-public"
KEY            = "example.14bits.bgen"

pytestmark = pytest.mark.skipif(
    not MINIO_ENDPOINT or not METAFILE,
    reason="BGEN_TEST_MINIO_ENDPOINT / BGEN_TEST_S3_METAFILE not set",
)

try:
    from upath import UPath
    _has_upath = True
except ImportError:
    _has_upath = False

needs_upath = pytest.mark.skipif(not _has_upath, reason="universal_pathlib not installed")

# Credentials that should NOT be present in the ambient env during UPath tests
_CRED_VARS = (
    "AWS_ACCESS_KEY_ID",
    "AWS_SECRET_ACCESS_KEY",
    "AWS_SESSION_TOKEN",
    "AWS_NO_SIGN_REQUEST",
    "AWS_ENDPOINT_URL",
    "AWS_DEFAULT_REGION",
)


def _check_bgen(bgen_path, metafile_path: str) -> None:
    """Open bgen_path, read a partition from the local metafile, check basic stats."""
    from cbgen import BgenFile, BgenMetafile

    with BgenFile(bgen_path) as bgen:
        assert bgen.nsamples == 500
        assert bgen.nvariants == 199

        with BgenMetafile(metafile_path) as mf:
            part = mf.read_partition(0)
            gt = bgen.read_genotype(part.variants[0].offset)
            assert gt.probability.shape == (500, 3)
            # First sample may be marked missing (NaN) – check non-missing rows
            non_missing = ~np.isnan(gt.probability).any(axis=1)
            assert non_missing.sum() >= 499


# ── Test 1: UPath with explicit key/secret ────────────────────────────────────

@needs_upath
def test_upath_signed_credentials(monkeypatch):
    """UPath carries key+secret; no ambient AWS env vars should be needed."""
    for var in _CRED_VARS:
        monkeypatch.delenv(var, raising=False)

    path = UPath(
        f"s3://{BUCKET_PRIVATE}/{KEY}",
        key=MINIO_ACCESS,
        secret=MINIO_SECRET,
        # MinIO requires path-style; s3fs reads endpoint_url at top-level
        endpoint_url=MINIO_ENDPOINT,
        client_kwargs={"endpoint_url": MINIO_ENDPOINT},
    )
    _check_bgen(path, METAFILE)


# ── Test 2: UPath with anon=True for a public bucket ─────────────────────────

@needs_upath
def test_upath_anon_public_bucket(monkeypatch):
    """UPath anon=True maps to AWS_NO_SIGN_REQUEST; public object accessible."""
    for var in _CRED_VARS:
        monkeypatch.delenv(var, raising=False)

    path = UPath(
        f"s3://{BUCKET_PUBLIC}/{KEY}",
        anon=True,
        endpoint_url=MINIO_ENDPOINT,
        client_kwargs={"endpoint_url": MINIO_ENDPOINT},
    )
    _check_bgen(path, METAFILE)


# ── Test 3: plain s3:// string (credentials via env vars) ────────────────────

def test_plain_s3_string_env_credentials(monkeypatch):
    """Baseline: plain string URI works when AWS env vars are set."""
    monkeypatch.setenv("AWS_ACCESS_KEY_ID",     MINIO_ACCESS)
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY",  MINIO_SECRET)
    monkeypatch.setenv("AWS_ENDPOINT_URL",       MINIO_ENDPOINT)
    monkeypatch.setenv("AWS_DEFAULT_REGION",     "us-east-1")
    monkeypatch.delenv("AWS_NO_SIGN_REQUEST", raising=False)

    _check_bgen(f"s3://{BUCKET_PRIVATE}/{KEY}", METAFILE)


# ── Test 4: UPath str() round-trip preserves the s3:// URI ───────────────────

@needs_upath
def test_upath_str_preserves_uri():
    """Ensure UPath.__str__ keeps the full s3:// URI (pathlib would mangle it)."""
    path = UPath(f"s3://{BUCKET_PRIVATE}/{KEY}")
    assert str(path) == f"s3://{BUCKET_PRIVATE}/{KEY}"
