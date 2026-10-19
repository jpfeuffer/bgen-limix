"""Public dimension discovery and single-pass compact capacity contracts."""

import math
import random
import struct

import numpy as np
import pytest
from cbgen import BgenFile

from .test_native_layouts import _layout2, _write_bgen


@pytest.mark.parametrize("compression", [0, 1, 2])
@pytest.mark.parametrize("source", ["actual", "declared"])
@pytest.mark.parametrize("model", ["auto", "unphased", "phased", "mixed"])
def test_dimensions_strategies(tmp_path, compression, source, model):
    models = [(False, [0, 1, 2], 3), (True, [0, 1, 2], 4)]
    blocks = []
    for phase, ploidies, k in models:
        block = bytearray(_layout2(11, phase, ploidies, k)[0])
        block[7] = 4
        blocks.append(block)
    path = tmp_path / "dimensions.bgen"
    offsets = _write_bgen(path, 2, compression, blocks, 3)
    p = 2 if source == "actual" else 4
    gp = hp = 0
    for phase, _, k in models:
        if model in ("unphased", "mixed") or (model == "auto" and not phase):
            gp = max(gp, math.comb(p + k - 1, k - 1))
        if model in ("phased", "mixed") or (model == "auto" and phase):
            hp = max(hp, 1, p * k)
    with BgenFile(path) as bgen:
        result = bgen.read_dimensions(
            offsets, ploidy_source=source, probability_model=model
        )
        assert (result.gp_width, result.hp_width, result.decode_width) == (
            gp,
            hp,
            max(gp, hp),
        )
        result = bgen.read_dimensions([], ploidy_source=source, probability_model=model)
        assert (result.gp_width, result.hp_width, result.decode_width) == (0, 0, 0)


@pytest.mark.parametrize("compression", [0, 1, 2])
@pytest.mark.parametrize("precision", [32, 64])
def test_compact_capacity(tmp_path, compression, precision):
    blocks = []
    for phase, ploidies, k in [
        (False, [0, 1, 2], 3),
        (True, [2, 0, 1], 4),
        (True, [0, 0, 0], 2),
    ]:
        block = bytearray(_layout2(13, phase, ploidies, k, [False, True, False])[0])
        block[7] = 4
        blocks.append(block)
    path = tmp_path / "compact.bgen"
    offsets = _write_bgen(path, 2, compression, blocks, 3)
    with BgenFile(path) as bgen:
        dimensions = bgen.read_dimensions(offsets)
        assert dimensions.decode_width == 8
        batch = bgen.read_genotypes_batch(
            offsets, precision=precision, max_actual_ncomb=dimensions.decode_width
        )
        np.testing.assert_array_equal(batch.ncombs, [6, 8, 0])
        assert batch.probability.shape == (3, 3, 8)
        for i, offset in enumerate(offsets):
            single = bgen.read_genotype(offset, precision)
            n = int(batch.ncombs[i])
            np.testing.assert_array_equal(
                batch.probability[i, :, :n], single.probability[:, :n]
            )
            assert np.isnan(batch.probability[i, :, n:]).all()
            np.testing.assert_array_equal(batch.ploidy[i], single.ploidy)
            np.testing.assert_array_equal(batch.missing[i], single.missing)
            assert batch.phased[i] == single.phased
        with pytest.raises(RuntimeError):
            bgen.read_genotypes_batch(offsets, max_actual_ncomb=7)
        with pytest.raises(ValueError):
            bgen.read_genotypes_batch(offsets, 8, max_actual_ncomb=8)
        with pytest.raises(RuntimeError):
            bgen.read_genotypes_batch(offsets, 8)


@pytest.mark.parametrize("compression", [0, 1, 2])
@pytest.mark.parametrize("phase", [False, True])
@pytest.mark.parametrize("bits", range(1, 33))
def test_compact_all_bits_and_models(tmp_path, compression, phase, bits):
    ploidies = [0, 1, 3, 2]
    block, _ = _layout2(bits, phase, ploidies, 3, [False, True, False, False])
    block = bytearray(block)
    block[7] = 63
    path = tmp_path / "bits.bgen"
    offsets = _write_bgen(path, 2, compression, [block] * 3, 4)
    selected = [offsets[2], offsets[0], offsets[2]]
    with BgenFile(path) as bgen:
        width = 9 if phase else 10
        dimensions = bgen.read_dimensions(selected)
        assert dimensions.decode_width == width
        for precision in (32, 64):
            batch = bgen.read_genotypes_batch(
                selected, precision=precision, max_actual_ncomb=width + 2
            )
            assert batch.probability.dtype == np.dtype(f"float{precision}")
            assert batch.probability.flags.c_contiguous
            assert batch.probability.flags.owndata
            for i, offset in enumerate(selected):
                single = bgen.read_genotype(offset, precision)
                np.testing.assert_array_equal(
                    batch.probability[i, :, :width], single.probability[:, :width]
                )
                assert np.isnan(batch.probability[i, :, width:]).all()
                assert batch.ncombs[i] == width
            snapshot = batch.probability.copy()
            mask = batch.missing.copy()
            bgen.read_genotype_metadata(offsets)
        bgen.close()
        np.testing.assert_array_equal(batch.probability, snapshot)
        np.testing.assert_array_equal(batch.missing, mask)


@pytest.mark.parametrize("compression", [0, 1, 2])
def test_unused_declared_overflow_and_override(tmp_path, compression):
    block = bytearray(_layout2(8, False, [0, 0], 16, [False, True])[0])
    block[7] = 63
    path = tmp_path / "overflow.bgen"
    offsets = _write_bgen(path, 2, compression, [block], 2)
    with BgenFile(path) as bgen:
        dimensions = bgen.read_dimensions(offsets)
        assert dimensions.gp_width == dimensions.decode_width == 1
        with pytest.raises(RuntimeError):
            bgen.read_dimensions(offsets, ploidy_source="declared")
        for source in ("actual", "declared"):
            bound = bgen.read_dimensions(offsets, ploidy_source=source, max_ploidy=0)
            assert bound == dimensions
            mixed = bgen.read_dimensions(
                offsets,
                ploidy_source=source,
                max_ploidy=1,
                probability_model="mixed",
            )
            assert (mixed.gp_width, mixed.hp_width, mixed.decode_width) == (16, 16, 16)
        for precision in (32, 64):
            batch = bgen.read_genotypes_batch(
                offsets, precision=precision, max_actual_ncomb=1
            )
            np.testing.assert_array_equal(batch.probability, [[[1], [np.nan]]])
            np.testing.assert_array_equal(batch.ncombs, [1])
            np.testing.assert_array_equal(batch.probability_count, [[1, 1]])
        with pytest.raises(RuntimeError):
            bgen.read_genotypes_batch(offsets)
        with pytest.raises(RuntimeError):
            bgen.read_genotypes_batch(offsets, max_actual_ncomb=0)


@pytest.mark.parametrize("compression", [0, 1])
def test_layout1_fixed_dimensions_without_payload(tmp_path, compression):
    rows = np.array([[6553, 10922, 15292], [0, 0, 0]], dtype="<u2")
    path = tmp_path / "layout1.bgen"
    offsets = _write_bgen(path, 1, compression, [rows.tobytes()] * 2, 2)
    with BgenFile(path) as bgen:
        for source in ("actual", "declared"):
            for model, widths in [
                ("auto", (3, 0)),
                ("unphased", (3, 0)),
                ("phased", (0, 4)),
                ("mixed", (3, 4)),
            ]:
                dims = bgen.read_dimensions(
                    offsets, ploidy_source=source, probability_model=model
                )
                assert (dims.gp_width, dims.hp_width) == widths
        batch = bgen.read_genotypes_batch(offsets, max_actual_ncomb=3)
        np.testing.assert_array_equal(batch.probability[0, 0], rows[0] / 32768)
        assert np.isnan(batch.probability[:, 1]).all()
    # No genotype bytes are necessary to discover the intrinsic layout-1 model.
    path.write_bytes(path.read_bytes()[: offsets[0]])
    with BgenFile(path) as bgen:
        assert bgen.read_dimensions(offsets[:1]).gp_width == 3
        with pytest.raises(RuntimeError):
            bgen.read_genotypes_batch(offsets[:1], max_actual_ncomb=3)


@pytest.mark.parametrize("compression", [0, 1, 2])
@pytest.mark.parametrize("phase", [False, True])
def test_zero_ploidy_and_empty_sample_semantics(tmp_path, compression, phase):
    path = tmp_path / "zero.bgen"
    block, _ = _layout2(32, phase, [0, 0], 3, [False, True])
    offsets = _write_bgen(path, 2, compression, [block], 2)
    with BgenFile(path) as bgen:
        dims = bgen.read_dimensions(offsets)
        assert dims.gp_width == (0 if phase else 1)
        assert dims.hp_width == (1 if phase else 0)
        assert dims.decode_width == 1
        batch = bgen.read_genotypes_batch(offsets, max_actual_ncomb=0 if phase else 1)
        np.testing.assert_array_equal(batch.ncombs, [0 if phase else 1])
        np.testing.assert_array_equal(batch.missing, [[False, True]])
        if phase:
            assert batch.probability.shape == (1, 2, 0)
        else:
            np.testing.assert_array_equal(batch.probability, [[[1], [np.nan]]])
    block = struct.pack("<IHBBBB", 0, 3, 0, 0, phase, 8)
    offsets = _write_bgen(path, 2, compression, [block], 0)
    with BgenFile(path) as bgen:
        batch = bgen.read_genotypes_batch(offsets, max_actual_ncomb=0)
        assert batch.probability.shape == (1, 0, 0)
        np.testing.assert_array_equal(batch.ncombs, [0])
        np.testing.assert_array_equal(batch.probability_count, np.empty((1, 0)))


@pytest.mark.parametrize("compression", [0, 1, 2])
def test_errors_and_assumptions(tmp_path, compression):
    path = tmp_path / "errors.bgen"
    block, _ = _layout2(8, True, [1, 2], 3, [False, True])
    offsets = _write_bgen(path, 2, compression, [block], 2)
    with BgenFile(path) as bgen:
        for source in ("actual", "declared"):
            # Wrong explicit models/bounds are assumptions, not validations.
            dims = bgen.read_dimensions(
                offsets,
                ploidy_source=source,
                probability_model="unphased",
                max_ploidy=0,
            )
            assert dims.gp_width == 1 and dims.hp_width == 0
        # Missing samples still count towards the capacity bound.
        with pytest.raises(RuntimeError):
            bgen.read_genotypes_batch(offsets, max_actual_ncomb=3)
        for value in (-1, 2**32, 2**80):
            with pytest.raises((ValueError, OverflowError)):
                bgen.read_genotypes_batch(offsets, max_actual_ncomb=value)
        for kwargs in [
            {"ploidy_source": "wrong"},
            {"probability_model": "wrong"},
            {"max_ploidy": -1},
            {"max_ploidy": 64},
        ]:
            with pytest.raises(ValueError):
                bgen.read_dimensions(offsets, **kwargs)
        for offset in (0, -1, 2**63, len(path.read_bytes()) + 1):
            with pytest.raises((RuntimeError, TypeError)):
                bgen.read_dimensions([offset])
        batch = bgen.read_genotypes_batch([], max_actual_ncomb=4)
        assert batch.probability.shape == (0, 2, 4)
    with pytest.raises(RuntimeError, match="closed"):
        bgen.read_dimensions(offsets)
    with pytest.raises(RuntimeError, match="closed"):
        bgen.read_genotypes_batch(offsets, max_actual_ncomb=6)


@pytest.mark.parametrize("compression", [0, 1, 2])
def test_prefix_validation_contract(tmp_path, compression):
    block = bytearray(_layout2(8, False, [1, 2], 3)[0])
    block[9] = 64  # Invalid second sample ploidy.
    path = tmp_path / "bad-ploidy.bgen"
    offsets = _write_bgen(path, 2, compression, [block], 2)
    with BgenFile(path) as bgen:
        if compression != 2:
            declared = bgen.read_dimensions(
                offsets, ploidy_source="declared", probability_model="unphased"
            )
            assert declared.gp_width == 6  # Skipped sample metadata, by contract.
        with pytest.raises(RuntimeError):
            bgen.read_dimensions(offsets)
        with pytest.raises(RuntimeError):
            bgen.read_genotypes_batch(offsets, max_actual_ncomb=6)


@pytest.mark.parametrize("compression", [0, 1, 2])
def test_phase_assumptions_and_malformed_headers(tmp_path, compression):
    good = bytearray(_layout2(8, False, [1, 2], 3)[0])
    path = tmp_path / "malformed.bgen"
    for position, value in [(0, 3), (4, 0), (6, 3), (7, 64), (10, 2), (11, 0)]:
        block = good.copy()
        block[position] = value
        offsets = _write_bgen(path, 2, compression, [block], 2)
        with BgenFile(path) as bgen:
            with pytest.raises(RuntimeError):
                bgen.read_dimensions(offsets)
            with pytest.raises(RuntimeError):
                bgen.read_genotypes_batch(offsets, max_actual_ncomb=6)
            if compression != 2 and position in (10, 11):
                # Explicit model intentionally does not read trailing phase/bits.
                assert (
                    bgen.read_dimensions(
                        offsets, ploidy_source="declared", probability_model="unphased"
                    ).gp_width
                    == 6
                )
    # Encoded block length is too small even for required sample header metadata.
    offsets = _write_bgen(path, 2, compression, [good], 2)
    data = bytearray(path.read_bytes())
    struct.pack_into("<I", data, offsets[0], 3)
    path.write_bytes(data)
    with BgenFile(path) as bgen, pytest.raises(RuntimeError):
        bgen.read_dimensions(offsets, ploidy_source="declared")


@pytest.mark.parametrize("compression", [0, 1, 2])
def test_prefix_is_not_full_payload_validation(tmp_path, compression):
    n = 10000
    rng = random.Random(2314)
    payload = bytearray()
    for _ in range(n):
        a = rng.randrange(256)
        payload += bytes((a, rng.randrange(256 - a)))
    block = struct.pack("<IHBB", n, 2, 2, 2) + bytes([2]) * n + bytes([0, 8]) + payload
    path = tmp_path / "trailing.bgen"
    offsets = _write_bgen(path, 2, compression, [block], n)
    path.write_bytes(path.read_bytes()[:-32])
    with BgenFile(path) as bgen:
        if compression != 2:
            assert bgen.read_dimensions(offsets).gp_width == 3
            assert (
                bgen.read_dimensions(
                    offsets, ploidy_source="declared", probability_model="unphased"
                ).gp_width
                == 3
            )
        else:
            with pytest.raises(RuntimeError):
                bgen.read_dimensions(offsets)
        with pytest.raises(RuntimeError):
            bgen.read_genotypes_batch(offsets, max_actual_ncomb=3)


@pytest.mark.parametrize("compression", [0, 1, 2])
def test_actual_max63_and_single_allele(tmp_path, compression):
    path = tmp_path / "max.bgen"
    for phase, k in [(False, 3), (True, 4), (False, 1), (True, 1)]:
        block, _ = _layout2(8, phase, [63, 0], k, [True, False])
        offsets = _write_bgen(path, 2, compression, [block], 2)
        width = 63 * k if phase else math.comb(63 + k - 1, k - 1)
        with BgenFile(path) as bgen:
            assert bgen.read_dimensions(offsets).decode_width == width
            batch = bgen.read_genotypes_batch(offsets, max_actual_ncomb=width)
            np.testing.assert_array_equal(
                batch.probability[0], bgen.read_genotype(offsets[0]).probability
            )
            with pytest.raises(RuntimeError):
                bgen.read_genotypes_batch(offsets, max_actual_ncomb=width - 1)
