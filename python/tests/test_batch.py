"""Native batched decoding must match the public single-record API."""

import math
import struct

import numpy as np
import pytest
from cbgen import (
    BgenFile,
    BgenMetafile,
    GenotypeBatchMetadata,
    GenotypeBatchResult,
    bgen_file,
)

from .test_native_layouts import _layout2, _write_bgen

METADATA_DTYPES = {
    "phased": np.bool_,
    "ploidy": np.uint8,
    "missing": np.bool_,
    "ncombs": np.uint32,
    "probability_count": np.uint32,
    "nalleles": np.uint16,
    "nbits": np.uint8,
}


def _check_arrays(result):
    assert isinstance(result, GenotypeBatchMetadata)
    for field, dtype in METADATA_DTYPES.items():
        array = getattr(result, field)
        assert array.dtype == dtype
        assert array.flags.c_contiguous
        assert array.flags.owndata
    if isinstance(result, GenotypeBatchResult):
        assert result.probability.flags.c_contiguous
        assert result.probability.flags.owndata
    else:
        assert not hasattr(result, "probability")


@pytest.mark.parametrize("layout,compression", [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2)])
@pytest.mark.parametrize("precision", [32, 64])
def test_structured_batch(tmp_path, layout, compression, precision):
    nsamples = 4
    if layout == 1:
        rows = np.array(
            [[6553, 10922, 15292], [0, 32768, 0], [0, 0, 0], [32768, 0, 0]],
            dtype="<u2",
        )
        blocks = [rows.tobytes(), bytes(nsamples * 6), rows.tobytes()]
    else:
        blocks = [
            _layout2(32, False, [0, 1, 3, 4], 3, [False, True, False, False])[0],
            _layout2(7, True, [2, 0, 1, 3], 4, [False, False, True, False])[0],
            _layout2(8, True, [0] * 4, 2, [True, False, True, False])[0],
        ]
    path = tmp_path / "batch.bgen"
    offsets = _write_bgen(path, layout, compression, blocks, nsamples)
    selected = [offsets[2], offsets[0], offsets[1], offsets[0]]
    retained = []
    with BgenFile(path) as bgen:
        bgen.create_metafile(str(tmp_path / "index"))
        with BgenMetafile(str(tmp_path / "index")) as index:
            variants = [
                variant
                for p in range(index.npartitions)
                for variant in index.read_partition(p).variants
            ]
        indexed_alleles = {v.offset: v.nalleles for v in variants}
        singles = [bgen.read_genotype(offset, precision) for offset in selected]
        metadata = bgen.read_genotype_metadata(selected)
        _check_arrays(metadata)
        width = max(single.probability.shape[1] for single in singles)
        for explicit_width in (None, width, width + 2):
            batch = bgen.read_genotypes_batch(selected, explicit_width, precision)
            _check_arrays(batch)
            assert batch.probability.shape == (4, nsamples, explicit_width or width)
            assert batch.probability.dtype == np.dtype(f"float{precision}")
            for i, single in enumerate(singles):
                nc = single.probability.shape[1]
                np.testing.assert_array_equal(
                    batch.probability[i, :, :nc], single.probability
                )
                assert np.isnan(batch.probability[i, :, nc:]).all()
                assert batch.phased[i] == single.phased
                assert batch.ncombs[i] == nc
                assert batch.nalleles[i] == indexed_alleles[selected[i]]
                assert batch.nbits[i] == (16 if layout == 1 else [8, 32, 7, 32][i])
                np.testing.assert_array_equal(batch.ploidy[i], single.ploidy)
                np.testing.assert_array_equal(batch.missing[i], single.missing)
                for s, ploidy in enumerate(single.ploidy):
                    k = int(batch.nalleles[i])
                    count = (
                        int(ploidy) * k
                        if single.phased
                        else math.comb(int(ploidy) + k - 1, k - 1)
                    )
                    assert batch.probability_count[i, s] == count
            for field in (
                "phased",
                "ploidy",
                "missing",
                "ncombs",
                "probability_count",
                "nalleles",
                "nbits",
            ):
                np.testing.assert_array_equal(
                    getattr(batch, field), getattr(metadata, field)
                )
            retained.append(batch)
        with pytest.raises(RuntimeError, match="max_ncomb"):
            bgen.read_genotypes_batch(selected, width - 1)
        np.testing.assert_array_equal(bgen.get_ncombs(selected), metadata.ncombs)
        expected64 = bgen.read_genotypes_batch(selected).probability
        snapshots = {
            field: getattr(metadata, field).copy() for field in METADATA_DTYPES
        }
        probability_snapshots = [batch.probability.copy() for batch in retained]
        bgen.read_genotype(offsets[-1])
    assert retained[0].probability.shape == (4, nsamples, width)
    for batch, snapshot in zip(retained, probability_snapshots):
        np.testing.assert_array_equal(batch.probability, snapshot)
    for field, snapshot in snapshots.items():
        np.testing.assert_array_equal(getattr(metadata, field), snapshot)
        for batch in retained:
            np.testing.assert_array_equal(getattr(batch, field), snapshot)
    with bgen_file(path) as legacy:
        np.testing.assert_array_equal(legacy.read_genotypes_batch(selected), expected64)
        np.testing.assert_array_equal(legacy.get_ncombs(selected), metadata.ncombs)
        padded = legacy.read_genotypes_batch(selected, width + 2)
        np.testing.assert_array_equal(padded[:, :, :width], expected64)
        assert np.isnan(padded[:, :, width:]).all()


@pytest.mark.parametrize("compression", [0, 1, 2])
@pytest.mark.parametrize("precision", [32, 64])
@pytest.mark.parametrize("bits", range(1, 33))
def test_batch_all_bit_depths(tmp_path, compression, precision, bits):
    models = [
        (False, [0, 1, 3, 7], 4),
        (True, [2, 0, 1, 3], 3),
        (False, [0, 0, 0, 0], 2),
        (True, [0, 0, 0, 0], 2),
    ]
    blocks = [
        _layout2(bits, phased, ploidies, k, [False, True, False, True])[0]
        for phased, ploidies, k in models
    ]
    path = tmp_path / "bits.bgen"
    offsets = _write_bgen(path, 2, compression, blocks, 4)
    with BgenFile(path) as bgen:
        for selected in (offsets, offsets[1:3], offsets[-1:]):
            singles = [bgen.read_genotype(o, precision) for o in selected]
            width = max(s.probability.shape[1] for s in singles)
            batch = bgen.read_genotypes_batch(
                np.array(selected, dtype=np.uint64), width, precision
            )
            _check_arrays(batch)
            for i, single in enumerate(singles):
                n = single.probability.shape[1]
                np.testing.assert_array_equal(
                    batch.probability[i, :, :n], single.probability
                )
                assert np.isnan(batch.probability[i, :, n:]).all()
                np.testing.assert_array_equal(batch.ploidy[i], single.ploidy)
                np.testing.assert_array_equal(batch.missing[i], single.missing)
                assert batch.phased[i] == single.phased
                assert batch.nbits[i] == bits


@pytest.mark.parametrize("compression", [0, 1, 2])
@pytest.mark.parametrize(
    "phased,ploidies,k",
    [
        (False, [0, 63, 1], 3),
        (True, [63, 0, 1], 4),
        (False, [12, 1, 0], 4),
        (False, [63, 0, 1], 1),
        (True, [63, 0, 1], 1),
    ],
)
def test_batch_large_and_single_allele_models(
    tmp_path, compression, phased, ploidies, k
):
    block, expected = _layout2(32, phased, ploidies, k, [False, True, False])
    path = tmp_path / "large.bgen"
    offsets = _write_bgen(path, 2, compression, [block, block], 3)
    with BgenFile(path) as bgen:
        for precision in (32, 64):
            batch = bgen.read_genotypes_batch(offsets, precision=precision)
            np.testing.assert_allclose(batch.probability[0], expected, atol=1e-7)
            np.testing.assert_array_equal(
                batch.probability[1],
                bgen.read_genotype(offsets[1], precision).probability,
            )
            counts = [
                p * k if phased else math.comb(p + k - 1, k - 1) for p in ploidies
            ]
            np.testing.assert_array_equal(batch.probability_count, [counts, counts])


@pytest.mark.parametrize("compression", [0, 1, 2])
def test_declared_width_and_empty_samples(tmp_path, compression):
    block = bytearray(_layout2(8, False, [0, 1], 3)[0])
    block[7] = 4  # Declared max ploidy may exceed every sample's actual ploidy.
    path = tmp_path / "declared.bgen"
    offsets = _write_bgen(path, 2, compression, [block], 2)
    with BgenFile(path) as bgen:
        batch = bgen.read_genotypes_batch(offsets)
        np.testing.assert_array_equal(batch.ncombs, [15])
        np.testing.assert_array_equal(batch.probability_count, [[1, 3]])
        np.testing.assert_array_equal(
            batch.probability[0], bgen.read_genotype(offsets[0]).probability
        )
    for phased in (False, True):
        block = struct.pack("<IHBBBB", 0, 2, 0, 0, phased, 8)
        offsets = _write_bgen(path, 2, compression, [block], 0)
        with BgenFile(path) as bgen:
            batch = bgen.read_genotypes_batch(offsets)
            assert batch.probability.shape == (1, 0, 0 if phased else 1)
            assert batch.ploidy.shape == (1, 0)
            _check_arrays(batch)


def test_empty_zero_and_errors(tmp_path):
    block, _ = _layout2(8, True, [0, 0], missing=[True, False])
    path = tmp_path / "zero.bgen"
    offsets = _write_bgen(path, 2, 2, [block], 2)
    bgen = BgenFile(path)
    for width in (None, 0, 4):
        batch = bgen.read_genotypes_batch([], width)
        _check_arrays(batch)
        assert batch.probability.shape == (0, 2, width or 0)
        assert (
            batch.ploidy.shape
            == batch.missing.shape
            == batch.probability_count.shape
            == (0, 2)
        )
        assert batch.phased.shape == batch.ncombs.shape == (0,)
        zero = bgen.read_genotypes_batch(offsets, width)
        assert zero.probability.shape == (1, 2, width or 0)
        assert np.isnan(zero.probability).all()
        np.testing.assert_array_equal(zero.missing, [[True, False]])
    assert bgen.read_genotype_metadata([]).ploidy.shape == (0, 2)
    for width in (-1, 2**32, 2**80):
        with pytest.raises((ValueError, OverflowError)):
            bgen.read_genotypes_batch(offsets, width)
    with pytest.raises(ValueError, match="Precision"):
        bgen.read_genotypes_batch(offsets, precision=16)
    for offset in (-1, 0, len(path.read_bytes()) + 1, 2**63, 2**64):
        with pytest.raises((RuntimeError, TypeError, ValueError, OverflowError)):
            bgen.read_genotypes_batch([offset], 0)
        with pytest.raises((RuntimeError, TypeError, ValueError, OverflowError)):
            bgen.read_genotype_metadata([offset])
    bgen.close()
    bgen.close()
    with pytest.raises(RuntimeError, match="closed"):
        bgen.read_genotypes_batch(offsets)
    with pytest.raises(RuntimeError, match="closed"):
        bgen.read_genotype_metadata(offsets)


@pytest.mark.parametrize("layout,compression", [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2)])
def test_batch_truncated_record_and_recovery(tmp_path, layout, compression):
    block = (
        np.array([[6553, 10922, 15292], [0, 0, 0]], dtype="<u2").tobytes()
        if layout == 1
        else _layout2(8, False, [1, 2])[0]
    )
    path = tmp_path / "truncated.bgen"
    offsets = _write_bgen(path, layout, compression, [block, block], 2)
    path.write_bytes(path.read_bytes()[:-1])
    with BgenFile(path) as bgen:
        for method in (
            bgen.read_genotype_metadata,
            bgen.get_ncombs,
            bgen.read_genotypes_batch,
        ):
            with pytest.raises(RuntimeError):
                method(offsets)
        with pytest.raises(RuntimeError):
            bgen.read_genotypes_batch(offsets, 3)
        valid = bgen.read_genotypes_batch(offsets[:1], 3)
        np.testing.assert_array_equal(
            valid.probability[0], bgen.read_genotype(offsets[0]).probability
        )


@pytest.mark.parametrize("compression", [0, 1, 2])
def test_metadata_does_not_validate_probability_values(tmp_path, compression):
    block = bytearray(_layout2(8, False, [2])[0])
    block[-2:] = b"\xff\xff"  # Well-framed but impossible sum of stored probabilities.
    path = tmp_path / "values.bgen"
    offsets = _write_bgen(path, 2, compression, [block], 1)
    with BgenFile(path) as bgen:
        metadata = bgen.read_genotype_metadata(offsets)
        _check_arrays(metadata)
        np.testing.assert_array_equal(metadata.ncombs, [3])
        np.testing.assert_array_equal(metadata.probability_count, [[3]])


@pytest.mark.parametrize("compression", [0, 1, 2])
@pytest.mark.parametrize("phased", [False, True])
def test_all_missing_retains_structural_counts(tmp_path, compression, phased):
    block, _ = _layout2(11, phased, [0, 1, 2], 3, [True] * 3)
    path = tmp_path / "missing.bgen"
    offsets = _write_bgen(path, 2, compression, [block], 3)
    with BgenFile(path) as bgen:
        batch = bgen.read_genotypes_batch(offsets)
        assert np.isnan(batch.probability).all()
        np.testing.assert_array_equal(batch.missing, [[True] * 3])
        np.testing.assert_array_equal(batch.ploidy, [[0, 1, 2]])
        np.testing.assert_array_equal(
            batch.probability_count, [[0, 3, 6] if phased else [1, 3, 6]]
        )
        zero, _ = _layout2(8, False, [0] * 3)
    offsets = _write_bgen(path, 2, compression, [zero], 3)
    with BgenFile(path) as bgen:
        with pytest.raises(RuntimeError, match="max_ncomb"):
            bgen.read_genotypes_batch(offsets, 0)
        batch = bgen.read_genotypes_batch(offsets, 1)
        np.testing.assert_array_equal(batch.probability, np.ones((1, 3, 1)))
        empty = bgen.read_genotypes_batch([], 2**32 - 1)
        assert empty.probability.shape == (0, 3, 2**32 - 1)
