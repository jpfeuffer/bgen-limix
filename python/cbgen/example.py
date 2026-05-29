"""Download and cache example BGEN files for testing."""

from __future__ import annotations

import logging
import os
from pathlib import Path

import pooch

__all__ = ["get"]

BGEN_CACHE_HOME = Path(
    os.environ.get("BGEN_CACHE_HOME", Path.home() / ".cache" / "bgen")
)
(BGEN_CACHE_HOME / "test_data").mkdir(parents=True, exist_ok=True)

pooch.get_logger().setLevel(logging.ERROR)

_registry = pooch.create(
    path=BGEN_CACHE_HOME / "test_data",
    base_url="https://github.com/fastlmm/bgen-sample-files/raw/main",
    registry={
        "complex.23bits.no.samples.bgen": "25d30a4e489da1aeb05f9893af98e8bf3b09d74db2982bf1828f8c8565886fc6",
        "haplotypes.bgen": "84e0b59efcc83c7c305cf5446e5dc26b49b15aeb4157a9eb36451376ce3efe4c",
        "haplotypes.bgen.metadata.corrupted": "8f55628770c1ae8155c1ced2463f15df80d32bc272a470bb1d6b68225e1604b1",
        "wrong.metadata": "f746345605150076f3234fbeea7c52e86bf95c9329b2f08e1e3e92a7918b98fb",
    },
)


def get(filename: str) -> Path:
    """Get file path to an example BGEN file (downloaded on first use).

    Parameters
    ----------
    filename
        Name of the example file.

    Returns
    -------
    Path to the cached file.
    """
    return Path(_registry.fetch(filename))
