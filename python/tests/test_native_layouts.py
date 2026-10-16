"""Generated BGEN records exercising the compiled public reader and indexer."""

import math
import struct
import zlib

import numpy as np
import pytest
from cbgen import BgenFile, BgenMetafile, bgen_file


def _string(value, width):
    data = value.encode()
    return len(data).to_bytes(width, "little") + data


def _compressed(data, compression):
    if compression == 1:
        return zlib.compress(data)
    if compression == 2:
        # Single-segment zstd frame with one final raw block, no external encoder.
        assert len(data) <= 128 * 1024
        return (
            b"\x28\xb5\x2f\xfd\xa0"
            + struct.pack("<I", len(data))
            + (len(data) << 3 | 1).to_bytes(3, "little")
            + data
        )
    return data


def _write_bgen(path, layout, compression, blocks, nsamples):
    data = bytearray(
        struct.pack(
            "<IIII4sI",
            20,
            20,
            len(blocks),
            nsamples,
            b"bgen",
            layout << 2 | compression,
        )
    )
    offsets = []
    for index, block in enumerate(blocks):
        if layout == 1:
            data += struct.pack("<I", nsamples)
        data += _string(f"snp_{index}", 2) + _string(f"rs{index}", 2)
        data += _string("chr1", 2) + struct.pack("<I", index + 1)
        nalleles = 2 if layout == 1 else struct.unpack_from("<H", block, 4)[0]
        if layout == 2:
            data += struct.pack("<H", nalleles)
        for allele in range(nalleles):
            data += _string("ACGT"[allele % 4], 4)
        offsets.append(len(data))
        payload = _compressed(block, compression)
        if layout == 2 and compression:
            payload = struct.pack("<I", len(block)) + payload
        if layout == 2 or compression:
            data += struct.pack("<I", len(payload))
        data += payload
    path.write_bytes(data)
    return offsets


def _layout2(bits, phased, ploidies, nalleles=2, missing=None):
    missing = missing if missing is not None else [False] * len(ploidies)
    block = struct.pack("<IHBB", len(ploidies), nalleles, min(ploidies), max(ploidies))
    block += bytes(p | (128 if m else 0) for p, m in zip(ploidies, missing))
    block += bytes([phased, bits])
    denominator = (1 << bits) - 1
    packed = 0
    shift = 0
    expected = []
    width = (
        max(ploidies) * nalleles
        if phased
        else math.comb(max(ploidies) + nalleles - 1, nalleles - 1)
    )
    for ploidy, is_missing in zip(ploidies, missing):
        groups = ploidy if phased else 1
        count = nalleles if phased else math.comb(ploidy + nalleles - 1, ploidy)
        row = []
        for _ in range(groups):
            values = [denominator // count] * (count - 1)
            for value in values:
                packed |= value << shift
                shift += bits
            row += [value / denominator for value in values]
            row += [(denominator - sum(values)) / denominator]
        row += [np.nan] * (width - len(row))
        expected.append([np.nan] * width if is_missing else row)
    block += packed.to_bytes((shift + 7) // 8, "little")
    return block, np.array(expected).reshape(len(ploidies), width)


def _index(bgen, path, offsets):
    bgen.create_metafile(str(path))
    with BgenMetafile(str(path)) as metafile:
        variants = [
            v
            for i in range(metafile.npartitions)
            for v in metafile.read_partition(i).variants
        ]
    assert [v.offset for v in variants] == offsets
    for i, variant in enumerate(variants):
        assert (variant.id, variant.rsid, variant.chrom, variant.position) == (
            f"snp_{i}",
            f"rs{i}",
            "chr1",
            i + 1,
        )
    return variants


@pytest.mark.parametrize("compression", [0, 1])
def test_layout1_index_and_probabilities(tmp_path, compression):
    rows = np.array(
        [[32768, 0, 0], [0, 32768, 0], [6553, 10922, 15292], [0, 0, 0]], dtype="<u2"
    )
    expected = rows.astype(float) / 32768
    expected[-1] = np.nan
    path = tmp_path / "layout1.bgen"
    offsets = _write_bgen(
        path, 1, compression, [rows.tobytes(), bytes(24), rows.tobytes()], 4
    )
    retained = []
    with BgenFile(str(path)) as bgen:
        assert bgen.layout == 1
        variants = _index(bgen, tmp_path / "index", offsets)
        for variant in variants:
            assert variant.nalleles == 2
            assert variant.allele_ids == ["A", "C"]
        for offset in (offsets[0], offsets[1], offsets[2], offsets[0]):
            gt = bgen.read_genotype(offset)
            want = np.full((4, 3), np.nan) if offset == offsets[1] else expected
            np.testing.assert_array_equal(gt.probability, want)
            np.testing.assert_array_equal(gt.ploidy, [2] * 4)
            np.testing.assert_array_equal(gt.missing, np.isnan(want[:, 0]))
            assert gt.phased is False
            retained.append(gt)
    with bgen_file(path) as legacy:
        assert legacy.layout == 1
        for precision in (32, 64):
            for offset in offsets:
                want = np.full((4, 3), np.nan) if offset == offsets[1] else expected
                gt = legacy.read_genotype(offset, precision)
                np.testing.assert_array_equal(gt.probability, want)
                np.testing.assert_array_equal(
                    legacy.read_probability(offset, precision), want
                )
                np.testing.assert_array_equal(gt.ploidy, [2] * 4)
                np.testing.assert_array_equal(gt.missing, np.isnan(want[:, 0]))
                assert gt.probability.dtype == np.dtype(f"float{precision}")
    np.testing.assert_array_equal(retained[0].probability, expected)
    np.testing.assert_array_equal(retained[0].missing, [False, False, False, True])


@pytest.mark.parametrize("compression", [0, 1, 2])
@pytest.mark.parametrize("phased", [False, True])
@pytest.mark.parametrize("bits", range(1, 33))
@pytest.mark.parametrize("nalleles", [2, 3])
def test_layout2_variable_blocks(tmp_path, compression, phased, bits, nalleles):
    ploidies = [0, 1, 3, 4]
    missing = [False, True, False, False]
    block, expected = _layout2(bits, phased, ploidies, nalleles, missing)
    path = tmp_path / "layout2.bgen"
    offsets = _write_bgen(path, 2, compression, [block, block], 4)
    with BgenFile(str(path)) as bgen:
        assert bgen.layout == 2
        _index(bgen, tmp_path / "index", offsets)
        for offset in reversed(offsets):
            gt = bgen.read_genotype(offset)
            np.testing.assert_allclose(gt.probability, expected, rtol=0, atol=1e-15)
            np.testing.assert_array_equal(gt.ploidy, ploidies)
            np.testing.assert_array_equal(gt.missing, missing)
            assert gt.phased == phased
            gt32 = bgen.read_genotype(offset, precision=32)
            np.testing.assert_allclose(gt32.probability, expected, atol=1e-7)
            np.testing.assert_array_equal(gt32.ploidy, ploidies)
            np.testing.assert_array_equal(gt32.missing, missing)


@pytest.mark.parametrize("compression", [0, 1, 2])
def test_layout2_zero_ploidy_phased(tmp_path, compression):
    empty, expected = _layout2(8, True, [0, 0], missing=[False, True])
    normal, normal_expected = _layout2(8, True, [2, 2])
    path = tmp_path / "empty.bgen"
    offsets = _write_bgen(path, 2, compression, [empty, normal, empty], 2)
    with BgenFile(str(path)) as bgen:
        _index(bgen, tmp_path / "index", offsets)
        for offset in offsets:
            gt = bgen.read_genotype(offset)
            want = normal_expected if offset == offsets[1] else expected
            np.testing.assert_array_equal(gt.probability, want)
            assert gt.phased is True
            if offset != offsets[1]:
                np.testing.assert_array_equal(gt.ploidy, [0, 0])
                np.testing.assert_array_equal(gt.missing, [False, True])
                assert gt.probability.shape == (2, 0)
                assert bgen.read_probability(offset, precision=32).shape == (2, 0)
                assert bgen.read_ncombs(offset) == 0


def test_layout2_many_combinations(tmp_path):
    block, expected = _layout2(32, False, [12], 4)
    path = tmp_path / "many.bgen"
    offsets = _write_bgen(path, 2, 1, [block], 1)
    assert expected.shape == (1, 455)
    with BgenFile(str(path)) as bgen:
        np.testing.assert_allclose(bgen.read_genotype(offsets[0]).probability, expected)


@pytest.mark.parametrize(
    "flags", [0, 3 | (2 << 2), 3 << 2, (1 << 2) | 2, (2 << 2) | (1 << 6)]
)
def test_unsupported_header_flags(tmp_path, flags):
    path = tmp_path / "bad.bgen"
    path.write_bytes(struct.pack("<IIII4sI", 20, 20, 0, 0, b"bgen", flags))
    with pytest.raises(RuntimeError, match="Failed to open"):
        BgenFile(path)


@pytest.mark.parametrize("length", range(24))
def test_truncated_header(tmp_path, length):
    path = tmp_path / "truncated.bgen"
    path.write_bytes(struct.pack("<IIII4sI", 20, 20, 0, 0, b"bgen", 8)[:length])
    with pytest.raises(RuntimeError, match="Failed to open"):
        BgenFile(path)


@pytest.mark.parametrize("count", [2**31, 2**32 - 1])
def test_unsigned_variant_count(tmp_path, count):
    path = tmp_path / "count.bgen"
    path.write_bytes(struct.pack("<IIII4sI", 20, 20, count, 0, b"bgen", 8))
    with BgenFile(path) as bgen:
        assert bgen.nvariants == count


@pytest.mark.parametrize("compression", [0, 1, 2])
@pytest.mark.parametrize(
    "mutation",
    [
        "nsamples",
        "phase",
        "bits",
        "ploidy",
        "short_header",
        "short_payload",
        "long_payload",
    ],
)
def test_malformed_layout2(tmp_path, compression, mutation):
    block, _ = _layout2(16, False, [2, 2])
    block = bytearray(block)
    if mutation == "nsamples":
        struct.pack_into("<I", block, 0, 3)
    elif mutation == "phase":
        block[10] = 2
    elif mutation == "bits":
        block[11] = 33
    elif mutation == "ploidy":
        block[8] = 63
    elif mutation == "short_header":
        block = block[:9]
    elif mutation == "short_payload":
        block = block[:-1]
    elif mutation == "long_payload":
        block += b"\0"
    path = tmp_path / "bad.bgen"
    offsets = _write_bgen(path, 2, compression, [block], 2)
    with (
        BgenFile(path) as bgen,
        pytest.raises(RuntimeError, match="Could not open genotype"),
    ):
        bgen.read_genotype(offsets[0])


@pytest.mark.parametrize("layout,compression", [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2)])
def test_truncated_variant_index(tmp_path, layout, compression):
    block = bytes(12) if layout == 1 else _layout2(8, False, [2, 2])[0]
    path = tmp_path / "truncated.bgen"
    _write_bgen(path, layout, compression, [block], 2)
    path.write_bytes(path.read_bytes()[:-1])
    with (
        BgenFile(path) as bgen,
        pytest.raises(RuntimeError, match="Error creating metafile"),
    ):
        bgen.create_metafile(str(tmp_path / "index"))


def test_layout1_variant_sample_count(tmp_path):
    path = tmp_path / "wrong_nsamples.bgen"
    _write_bgen(path, 1, 0, [bytes(12)], 2)
    data = bytearray(path.read_bytes())
    struct.pack_into("<I", data, 24, 3)
    path.write_bytes(data)
    with (
        BgenFile(path) as bgen,
        pytest.raises(RuntimeError, match="Error creating metafile"),
    ):
        bgen.create_metafile(str(tmp_path / "index"))


@pytest.mark.parametrize("layout,compression", [(1, 1), (2, 1), (2, 2)])
@pytest.mark.parametrize("change", [-1, 1])
def test_decompressed_length_mismatch(tmp_path, layout, compression, change):
    block = bytes(12 + change) if layout == 1 else _layout2(8, False, [2, 2])[0]
    path = tmp_path / "length.bgen"
    offsets = _write_bgen(path, layout, compression, [block], 2)
    if layout == 2:
        data = bytearray(path.read_bytes())
        struct.pack_into("<I", data, offsets[0] + 4, len(block) + change)
        path.write_bytes(data)
    with (
        BgenFile(path) as bgen,
        pytest.raises(RuntimeError, match="Could not open genotype"),
    ):
        bgen.read_genotype(offsets[0])


@pytest.mark.parametrize("compression", [0, 1, 2])
@pytest.mark.parametrize("phased", [False, True])
@pytest.mark.parametrize("ploidies,nalleles", [([0, 0], 2), ([0, 63], 2), ([0, 63], 1)])
def test_layout2_ploidy_boundaries(tmp_path, compression, phased, ploidies, nalleles):
    block, expected = _layout2(32, phased, ploidies, nalleles)
    path = tmp_path / "boundary.bgen"
    offsets = _write_bgen(path, 2, compression, [block], 2)
    with BgenFile(path) as bgen:
        np.testing.assert_array_equal(
            bgen.read_genotype(offsets[0]).probability, expected
        )
