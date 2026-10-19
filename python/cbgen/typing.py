from __future__ import annotations

from dataclasses import dataclass
from typing import Union

import numpy as np
import numpy.typing as npt

__all__ = [
    "Genotype",
    "GenotypeBatchMetadata",
    "GenotypeBatchResult",
    "GenotypeDimensions",
    "GenotypeResult",
    "Partition",
    "Variants",
]


@dataclass
class GenotypeResult:
    """
    Genotype result.

    Attributes
    ----------
    probability
        Probability array of shape (nsamples, ncombs).
    phased
        Phasedness.
    ploidy
        Ploidy array of shape (nsamples,).
    missing
        Missingness array of shape (nsamples,).
    """

    probability: Union[npt.NDArray[np.float64], npt.NDArray[np.float32]]
    phased: bool
    ploidy: npt.NDArray[np.uint8]
    missing: npt.NDArray[np.bool_]


# Legacy alias for API compatibility with the old cbgen package.
Genotype = GenotypeResult


@dataclass(frozen=True)
class GenotypeDimensions:
    """Probability capacities, with no per-sample metadata arrays.

    gp_width/hp_width are maxima under the selected model assumption.
    Zero means that model is absent; hp_width is at least one if selected,
    even for zero-ploidy phased records. decode_width is their maximum.
    Empty selections have all widths zero. These dimensions do not certify
    probability payload integrity or explicit model assumptions.
    """

    gp_width: int
    hp_width: int
    decode_width: int


@dataclass
class GenotypeBatchMetadata:
    """Metadata for V requested variants and S samples, in requested order.

    ``phased``, ``ncombs``, ``nalleles`` and ``nbits`` have shape (V,).
    ``ploidy``, ``missing`` and ``probability_count`` have shape (V, S).
    ``ncombs`` is the single-reader row width for each variant, computed from
    declared maximum ploidy (which may exceed all actual sample ploidies).
    ``probability_count`` is each sample's structural width, independent of
    missingness: ploidy * nalleles if phased, otherwise the number of unordered
    genotypes (one for zero-ploidy unphased samples). Missingness is never
    inferred from NaNs or padding. ``nbits`` is fixed at 16 for layout 1.
    All arrays are C-contiguous and own their data.
    """

    phased: npt.NDArray[np.bool_]
    ploidy: npt.NDArray[np.uint8]
    missing: npt.NDArray[np.bool_]
    ncombs: npt.NDArray[np.uint32]
    probability_count: npt.NDArray[np.uint32]
    nalleles: npt.NDArray[np.uint16]
    nbits: npt.NDArray[np.uint8]


@dataclass
class GenotypeBatchResult(GenotypeBatchMetadata):
    """Decoded batch with probabilities of shape (V, S, W).

    The first ncombs[v] columns match a single-record read exactly. Remaining
    batch-width slots and unused sample-width slots are NaN. Missing samples
    are also NaN: use missing/probability_count to distinguish missingness and
    structural padding. Phased rows flatten haplotype-major allele
    probabilities, just like GenotypeResult.
    In compact-capacity mode, ncombs is the maximum actual sample width,
    not the declared maximum-ploidy width (zero for zero samples).
    """

    probability: npt.NDArray[np.float64] | npt.NDArray[np.float32]


class Variants:
    """
    Array-style access to variant metadata within a partition.

    Supports both old-style indexed access (``variants.id[i]``) returning numpy
    byte-string arrays, and new-style list access (``variants[i]``) returning
    ``VariantInfo`` objects.

    Attributes
    ----------
    id
        Variant IDs as numpy byte-string array.
    rsid
        Reference SNP cluster IDs as numpy byte-string array.
    chromosome
        Chromosomes as numpy byte-string array.
    chrom
        Alias for chromosome (new-style name).
    position
        Positions as numpy uint32 array.
    nalleles
        Number of alleles as numpy uint16 array.
    allele_ids
        Allele IDs as numpy byte-string array (comma-joined).
    offset
        Variant byte offsets as numpy uint64 array.
    """

    def __init__(self, variant_list: list):
        self._variants = variant_list
        self._id: npt.NDArray | None = None
        self._rsid: npt.NDArray | None = None
        self._chrom: npt.NDArray | None = None
        self._position: npt.NDArray | None = None
        self._nalleles: npt.NDArray | None = None
        self._allele_ids: npt.NDArray | None = None
        self._offset: npt.NDArray | None = None

    def __len__(self) -> int:
        return len(self._variants)

    def __getitem__(self, index):
        return self._variants[index]

    @property
    def size(self) -> int:
        """Number of variants."""
        return len(self._variants)

    @property
    def id(self) -> npt.NDArray:
        """Variant IDs as numpy byte-string array."""
        if self._id is None:
            self._id = np.array([v.id_bytes for v in self._variants], dtype="S")
        return self._id

    @property
    def rsid(self) -> npt.NDArray:
        """RSIDs as numpy byte-string array."""
        if self._rsid is None:
            self._rsid = np.array([v.rsid_bytes for v in self._variants], dtype="S")
        return self._rsid

    @property
    def chromosome(self) -> npt.NDArray:
        """Chromosomes as numpy byte-string array."""
        if self._chrom is None:
            self._chrom = np.array([v.chrom_bytes for v in self._variants], dtype="S")
        return self._chrom

    @property
    def chrom(self) -> npt.NDArray:
        """Alias for chromosome."""
        return self.chromosome

    @property
    def position(self) -> npt.NDArray[np.uint32]:
        """Positions as numpy uint32 array."""
        if self._position is None:
            self._position = np.array(
                [v.position for v in self._variants], dtype=np.uint32
            )
        return self._position

    @property
    def nalleles(self) -> npt.NDArray[np.uint16]:
        """Number of alleles as numpy uint16 array."""
        if self._nalleles is None:
            self._nalleles = np.array(
                [v.nalleles for v in self._variants], dtype=np.uint16
            )
        return self._nalleles

    @property
    def allele_ids(self) -> npt.NDArray:
        """Allele IDs as numpy byte-string array (comma-joined)."""
        if self._allele_ids is None:
            self._allele_ids = np.array(
                [v.allele_ids_bytes for v in self._variants], dtype="S"
            )
        return self._allele_ids

    @property
    def offset(self) -> npt.NDArray[np.uint64]:
        """Variant byte offsets as numpy uint64 array."""
        if self._offset is None:
            self._offset = np.array(
                [v.offset for v in self._variants], dtype=np.uint64
            )
        return self._offset


@dataclass
class Partition:
    """
    Partition of variants.

    Supports both old-style access (``part.variants.offset[i]``) and new-style
    access (``part.variants[i].offset``).

    Attributes
    ----------
    offset
        Partition offset (partition_size * index).
    variants
        Variants accessor.
    """

    offset: int
    variants: Variants
