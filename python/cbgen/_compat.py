"""Backward-compatible API shims matching the old cbgen package interface."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Union

import numpy as np
import numpy.typing as npt

from cbgen._core import BgenFile as _BgenFile
from cbgen._core import BgenMetafile as _BgenMetafile
from cbgen._s3_utils import s3_env_from_storage_options
from cbgen.typing import Genotype, Partition, Variants

__all__ = ["bgen_file", "bgen_metafile"]


class bgen_file:
    """
    BGEN file handler (legacy API compatible with the old cbgen package).

    Parameters
    ----------
    filepath
        BGEN file path.
    """

    def __init__(self, filepath: Union[str, "os.PathLike[str]"]):
        s = str(filepath)
        if s.startswith("s3://"):
            self._filepath = filepath
        else:
            self._filepath = Path(filepath)
        with s3_env_from_storage_options(filepath):
            self._impl = _BgenFile(str(self._filepath))

    @property
    def filepath(self) -> "Union[Path, os.PathLike[str], str]":
        """File path."""
        return self._filepath

    @property
    def nvariants(self) -> int:
        """Number of variants."""
        return self._impl.nvariants

    @property
    def nsamples(self) -> int:
        """Number of samples."""
        return self._impl.nsamples

    @property
    def contain_samples(self) -> bool:
        """Check if file contains samples."""
        return self._impl.contain_samples

    def read_samples(self) -> npt.NDArray:
        """
        Read samples as a numpy byte-string array.

        Returns
        -------
        Numpy array of byte-string sample names.
        """
        samples = self._impl.read_samples()
        return np.array([s.encode() for s in samples], dtype="S")

    def create_metafile(self, filepath: Union[str, Path], verbose: bool = False):
        """
        Create a metafile.

        Parameters
        ----------
        filepath
            Metafile path.
        verbose
            Show progress.
        """
        self._impl.create_metafile(str(filepath), verbose)

    def read_genotype(self, offset: int, precision: int = 64) -> Genotype:
        """
        Read genotype at the given variant offset.

        Parameters
        ----------
        offset
            Variant offset.
        precision
            Probability precision in bits: 64 (default) or 32.

        Returns
        -------
        Genotype result.
        """
        if precision not in (64, 32):
            raise ValueError("Precision must be either 64 or 32.")

        if precision == 64:
            d = self._impl.read_genotype(offset)
        else:
            d = self._impl.read_genotype32(offset)

        return Genotype(
            probability=d["probability"],
            phased=d["phased"],
            ploidy=d["ploidy"],
            missing=d["missing"],
        )

    def read_probability(self, offset: int, precision: int = 64) -> npt.NDArray:
        """
        Read genotype probability array only.

        Parameters
        ----------
        offset
            Variant offset.
        precision
            Probability precision in bits: 64 (default) or 32.

        Returns
        -------
        Probability array of shape (nsamples, ncombs).
        """
        if precision not in (64, 32):
            raise ValueError("Precision must be either 64 or 32.")

        if precision == 64:
            d = self._impl.read_genotype(offset)
        else:
            d = self._impl.read_genotype32(offset)

        return d["probability"]

    def close(self):
        """Close file stream."""
        self._impl.close()

    def __del__(self):
        self.close()

    def __enter__(self) -> bgen_file:
        return self

    def __exit__(self, *_):
        self.close()


class bgen_metafile:
    """
    BGEN metafile handler (legacy API compatible with the old cbgen package).

    Parameters
    ----------
    filepath
        BGEN metafile path.
    """

    def __init__(self, filepath: Union[str, Path]):
        self._filepath = Path(filepath)
        self._impl = _BgenMetafile(str(self._filepath))

    @property
    def filepath(self) -> Path:
        """File path."""
        return self._filepath

    @property
    def npartitions(self) -> int:
        """Number of partitions."""
        return self._impl.npartitions

    @property
    def nvariants(self) -> int:
        """Number of variants."""
        return self._impl.nvariants

    @property
    def partition_size(self) -> int:
        """Number of variants per partition."""
        return self._impl.partition_size

    def read_partition(self, index: int) -> Partition:
        """
        Read partition.

        Parameters
        ----------
        index
            Partition index.

        Returns
        -------
        Partition with array-style variants access.
        """
        part = self._impl.read_partition(index)
        variants = Variants(part.variants)
        return Partition(offset=part.offset, variants=variants)

    def close(self):
        """Close file stream."""
        self._impl.close()

    def __del__(self):
        self.close()

    def __enter__(self) -> bgen_metafile:
        return self

    def __exit__(self, *_):
        self.close()
