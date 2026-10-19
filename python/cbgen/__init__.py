"""cbgen - Python wrapper around the BGEN C library."""

from __future__ import annotations

import os
from typing import Union

from cbgen import example, typing
from cbgen._compat import bgen_file, bgen_metafile
from cbgen._core import BgenFile as _BgenFile
from cbgen._core import BgenMetafile, PartitionResult, VariantInfo
from cbgen._s3_utils import s3_env_from_storage_options
from cbgen.example import BGEN_CACHE_HOME
from cbgen.typing import (
    Genotype,
    GenotypeBatchMetadata,
    GenotypeBatchResult,
    GenotypeDimensions,
    GenotypeResult,
    Partition,
    Variants,
)

__all__ = [
    "BGEN_CACHE_HOME",
    "BgenFile",
    "BgenMetafile",
    "Genotype",
    "GenotypeBatchMetadata",
    "GenotypeBatchResult",
    "GenotypeDimensions",
    "GenotypeResult",
    "Partition",
    "PartitionResult",
    "VariantInfo",
    "Variants",
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

    def __init__(self, filepath: Union[str, "os.PathLike[str]"]):
        with s3_env_from_storage_options(filepath):
            self._impl = _BgenFile(str(filepath))

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
    def layout(self) -> int:
        """BGEN layout version (1 or 2), from the native parsed header."""
        return self._impl.layout

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

    def read_ncombs(self, offset: int) -> int:
        """Return the number of genotype combinations for the variant at offset.

        Unlike :meth:`read_probability`, this does **not** allocate or fill a
        probability array — it only opens the genotype header and reads the
        ``ncombs`` field.

        Parameters
        ----------
        offset
            Variant offset (from VariantInfo.offset).

        Returns
        -------
        Number of genotype combinations (int).
        """
        return int(self._impl.read_ncombs(offset))

    def get_ncombs(self, offsets: list[int]):
        """Read uint32 row widths (V,) without decoding probabilities.

        Supports all layouts and compression methods. Compressed layout-2
        blocks are still decompressed. Ordering and duplicates are preserved.
        """
        return self._impl.get_ncombs(offsets)

    def read_genotype_metadata(self, offsets: list[int]) -> GenotypeBatchMetadata:
        """Read metadata for requested offsets without decoding probabilities.

        Supports both layouts and all supported compression methods.
        Compressed layout 2 still decompresses the block; layout 1 reads its
        triplets to determine missingness. Does not allocate probability arrays.
        Validates framing, headers and encoded payload lengths, not numerical
        probability sums or application-specific probability tolerances.
        Output shapes are (V,) and (V, nsamples), including for empty input.
        All arrays remain valid after later reads or closing the file.
        """
        return GenotypeBatchMetadata(**self._impl.read_genotype_metadata(offsets))

    def read_dimensions(
        self,
        offsets: list[int],
        *,
        ploidy_source: str = "actual",
        probability_model: str = "auto",
        max_ploidy: int | None = None,
    ) -> GenotypeDimensions:
        """Discover GP/HP/decode maxima without probability arrays.

        ``actual`` scans sample ploidy; ``declared`` uses header max ploidy.
        ``auto`` reads original phase; ``unphased``/``phased`` assume a model,
        and ``mixed`` reserves both widths for every record. Assumptions are
        not validated here and must be checked during full decode.

        Layout 1 uses fixed K=2/P=2 without genotype payload reads (explicit
        phased reserves HP=4; mixed GP=3/HP=4). Raw/zlib layout 2 read needed
        prefixes: declared+explicit model needs eight uncompressed bytes;
        auto must reach phase after sample ploidies. Zstd loads a full block.
        Prefixes do not validate skipped metadata, trailing data/checksums or
        numerical probabilities. Full reads validate framing/payload integrity;
        consumers must validate numerical probabilities separately.
        hp_width>=1 retains phase presence when all selected ploidies are zero.
        Empty selections return all zeros. Width overflow raises an error.

        max_ploidy in [0,63] overrides the capacity calculation, not sample
        validation. This avoids unused declared-width overflow and, on raw/
        zlib, skips scanning sample ploidy even with source='actual'. AUTO
        still reaches phase; explicit models need only eight header bytes.
        Validate this assumed bound against ploidy during later full decode.
        """
        if ploidy_source not in ("actual", "declared"):
            raise ValueError("ploidy_source must be 'actual' or 'declared'.")
        models = {"auto": 0, "unphased": 1, "phased": 2, "mixed": 3}
        if probability_model not in models:
            raise ValueError("Invalid probability_model.")
        if max_ploidy is not None and not 0 <= max_ploidy <= 63:
            raise ValueError("max_ploidy must be in [0, 63].")
        return GenotypeDimensions(
            **self._impl.read_dimensions(
                offsets, ploidy_source == "actual", models[probability_model],
                -1 if max_ploidy is None else max_ploidy,
            )
        )

    def read_genotypes_batch(
        self,
        offsets: list[int],
        max_ncomb: int | None = None,
        precision: int = 64,
        *,
        max_actual_ncomb: int | None = None,
    ) -> GenotypeBatchResult:
        """Decode a native structured batch across all supported encodings.

        ``offsets`` select original genotype blocks, preserving order and
        duplicates. ``max_ncomb`` is the padded width W; None scans requested
        widths first, then reads each block again to decode. An explicit
        uint32 width avoids that scan and reads/decompresses each block once.
        A width smaller than any requested single-reader width raises an error.
        Probability dtype is float64 (default) or float32.

        Returns probabilities (V, nsamples, W) and per-variant/sample metadata;
        see GenotypeBatchResult for field shapes and padding semantics.
        Empty input returns (0, nsamples, W), with W=0 if auto-sized. A wholly
        zero-ploidy phased batch may also have W=0. Allocation is limited to
        the requested batch and checked for dimension/byte-count overflow.
        Arrays are C-contiguous, independently owned and stable after close.
        Output storage is V*nsamples*W*itemsize bytes plus 6*V*nsamples+8*V
        metadata bytes. Native scratch holds at most one genotype block and
        one variant's decoded probabilities; no second full batch is built.
        Request float64 when deriving dosage before a later float32 cast.

        Alternatively, ``max_actual_ncomb`` is a compact capacity bounding
        actual structural sample widths, including missing samples. It is
        mutually exclusive with max_ncomb and never pre-scans. Full records
        validate before each variant's output is written; bounds are checked,
        not used to truncate. Declared width may exceed capacity. No declared-
        width probability scratch is allocated, and ncombs reports effective
        actual widths (zero for zero samples) in this mode. Earlier variants
        may have been decoded when a later variant raises an exception; no
        partial batch is returned. Default behavior is unchanged.
        """
        if precision not in (64, 32):
            raise ValueError("Precision must be either 64 or 32.")
        if max_ncomb is not None and not 0 <= max_ncomb <= 2**32 - 1:
            raise ValueError("max_ncomb must be an unsigned 32-bit width.")
        if max_actual_ncomb is not None:
            if max_ncomb is not None:
                raise ValueError("max_ncomb and max_actual_ncomb are mutually exclusive.")
            if not 0 <= max_actual_ncomb <= 2**32 - 1:
                raise ValueError("max_actual_ncomb must be an unsigned 32-bit width.")
        return GenotypeBatchResult(
            **self._impl.read_genotypes_batch_result(
                offsets, -1 if max_ncomb is None else max_ncomb, precision,
                -1 if max_actual_ncomb is None else max_actual_ncomb,
            )
        )

    def close(self) -> None:
        self._impl.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


try:
    from importlib.metadata import version

    __version__ = version("bgen-limix")
except Exception:
    __version__ = "0.0.0"
