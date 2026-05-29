"""cbgen - Python wrapper around the BGEN C library."""

from __future__ import annotations

from cbgen import example, typing
from cbgen._compat import bgen_file, bgen_metafile
from cbgen._core import BgenFile as _BgenFile
from cbgen._core import BgenMetafile, PartitionResult, VariantInfo
from cbgen.example import BGEN_CACHE_HOME
from cbgen.typing import Genotype, GenotypeResult, Partition, Variants

__all__ = [
    "BGEN_CACHE_HOME",
    "BgenFile",
    "BgenMetafile",
    "Genotype",
    "GenotypeResult",
    "Partition",
    "PartitionResult",
    "Variants",
    "VariantInfo",
    "__version__",
    "bgen_file",
    "bgen_metafile",
    "example",
    "typing",
]


class BgenFile:
    """BGEN file handler with GenotypeResult support.

    Parameters
    ----------
    filepath
        Path to the BGEN file.
    """

    def __init__(self, filepath: str):
        self._impl = _BgenFile(filepath)

    @property
    def filepath(self) -> str:
        return self._impl.filepath

    @property
    def nvariants(self) -> int:
        return self._impl.nvariants

    @property
    def nsamples(self) -> int:
        return self._impl.nsamples

    @property
    def contain_samples(self) -> bool:
        return self._impl.contain_samples

    def read_samples(self) -> list[str]:
        return self._impl.read_samples()

    def create_metafile(self, filepath: str, verbose: bool = False) -> None:
        self._impl.create_metafile(filepath, verbose)

    def read_genotype(self, offset: int, precision: int = 64) -> GenotypeResult:
        """Read genotype at the given offset.

        Parameters
        ----------
        offset
            Variant offset (from VariantInfo.offset).
        precision
            Probability precision in bits: 64 (default) or 32.

        Returns
        -------
        GenotypeResult with probability, phased, ploidy, and missing arrays.
        """
        if precision not in (64, 32):
            raise ValueError("Precision must be either 64 or 32.")

        if precision == 64:
            d = self._impl.read_genotype(offset)
        else:
            d = self._impl.read_genotype32(offset)

        return GenotypeResult(
            probability=d["probability"],
            phased=d["phased"],
            ploidy=d["ploidy"],
            missing=d["missing"],
        )

    def read_probability(self, offset: int, precision: int = 64):
        """Read genotype probability array only.

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

    def close(self) -> None:
        self._impl.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


try:
    from importlib.metadata import version

    __version__ = version("cbgen")
except Exception:
    __version__ = "0.0.0"
