"""Independent wire oracles for specialized hard-call and fused encoding."""

import concurrent.futures
import itertools
import struct
import sys
import threading
import zlib

import numpy as np
import pytest
from cbgen import (
    BgenFile,
    BufferTooSmall,
    encode_layout2_fixed,
    encode_layout2_hardcalls,
    encode_layout2_hardcalls_legacy,
    layout2_fixed_size,
)


def _block(gt, phase, legacy=True):
    ns = len(gt)
    ploidies, flags, values = [], [], []
    for a, b in gt:
        p = 0 if a == -2 else 1 if b == -2 else 2
        missing = a == -1 or b == -1
        ploidies.append(p)
        flags.append(p | (128 if missing else 0))
        if p == 2:
            values.extend(
                [int(a == 0 and not missing) * 255, int(b == 0 and not missing) * 255]
                if phase
                else [
                    int(a == b == 0 and not missing) * 255,
                    int(a != b and not missing) * 255,
                ]
            )
        elif p == 1:
            values.append(int(a == 0) * 255)
    default = 2 if legacy else 0
    return (
        struct.pack(
            "<IHBB",
            ns,
            2,
            min(ploidies, default=default),
            max(ploidies, default=default),
        )
        + bytes(flags)
        + bytes([phase, 8])
        + bytes(values)
    )


def _stored(raw):
    chunks = [raw[i : i + 65535] for i in range(0, len(raw), 65535)] or [b""]
    return (
        b"\x78\x01"
        + b"".join(
            bytes([i == len(chunks) - 1])
            + struct.pack("<HH", len(x), 65535 - len(x))
            + x
            for i, x in enumerate(chunks)
        )
        + struct.pack(">I", zlib.adler32(raw))
    )


def _field(data, width=0):
    return np.frombuffer(
        data.ljust(max(len(data), width), b"\0"), dtype=np.uint8
    ).reshape(1, -1)


def _fixed_inputs(ns=3, nv=2, ploidy=2):
    rng = np.random.default_rng(132)
    gt = rng.integers(0, 2, (nv, ns, 2), dtype=np.int8)
    if ploidy == 1:
        gt[:, :, 1] = -2
    if ns:
        gt[:, 0, 0] = -1
    fields = [
        np.repeat(_field(x, 7), nv, axis=0) for x in (b"id", b"rs", b"chr", b"A", b"T")
    ]
    position = np.arange(nv, dtype=np.int32)
    phase = np.arange(nv) % 2 == 1
    size = layout2_fixed_size(ns, ploidy, 9)
    out = np.full((nv, size), 0xCC, dtype=np.uint8)
    return [*fields, position, gt, phase, out, ploidy, 9]


def _variant(args, v):
    fields = [bytes(a[v]).split(b"\0")[0] for a in args[:5]]
    raw = _block(args[6][v], bool(args[7][v]))
    compressed = _stored(raw)
    record = b"".join(struct.pack("<H", len(s)) + s for s in fields[:3])
    record += struct.pack("<IH", int(args[5][v]), 2)
    record += b"".join(struct.pack("<I", len(s)) + s for s in fields[3:])
    record += struct.pack("<II", 4 + len(compressed), len(raw)) + compressed
    return record


@pytest.mark.parametrize("phase", [False, True])
def test_all_valid_pairs_and_general_selection(phase):
    pairs = list(itertools.product((-1, 0, 1), (-2, -1, 0, 1)))
    gt = np.array([pairs, pairs[::-1]], dtype=np.int8)
    phases = np.full(2, phase)
    padded, lengths = encode_layout2_hardcalls_legacy(gt, phases)
    raw, offsets = encode_layout2_hardcalls(
        genotypes=gt, phased=phases, allele_counts=np.array([2, 2], dtype=np.uint16)
    )
    for v in range(2):
        expected = _block(gt[v], phase)
        assert bytes(padded[v, : lengths[v]]) == expected
        assert raw[int(offsets[v]) : int(offsets[v + 1])] == expected
    assert padded.dtype == np.uint8 and lengths.dtype == np.uint32


@pytest.mark.parametrize("phase", [False, True])
def test_general_fastpath_preserves_p0_and_haploid(phase):
    for width in (1, 2):
        gt = np.array([[[0, -2], [-1, -2], [-2, -2], [1, -2]]], dtype=np.int8)
        gt = gt[:, :, :width].copy()
        raw, _ = encode_layout2_hardcalls(
            genotypes=gt,
            phased=np.array([phase]),
            allele_counts=np.array([2], dtype=np.uint16),
        )
        padded = np.pad(gt, ((0, 0), (0, 0), (0, 2 - width)), constant_values=-2)
        assert raw == _block(padded[0], phase, legacy=False)
    with pytest.raises(ValueError, match="zero-ploidy"):
        encode_layout2_hardcalls_legacy(padded, np.array([phase]))


@pytest.mark.parametrize("slot", [0, 1])
def test_every_invalid_int8_value(slot):
    for invalid in range(-128, 128):
        if invalid in (-2, -1, 0, 1):
            continue
        gt = np.zeros((1, 3, 2), dtype=np.int8)
        gt[0, 2, slot] = invalid
        with pytest.raises(ValueError, match="variant 0, sample 2.*out of range"):
            encode_layout2_hardcalls_legacy(gt, np.array([False]))


@pytest.mark.parametrize("nv,ns", [(0, 3), (3, 0), (0, 0)])
def test_empty_legacy(nv, ns):
    gt = np.empty((nv, ns, 2), dtype=np.int8)
    data, lengths = encode_layout2_hardcalls_legacy(gt, np.zeros(nv, dtype=bool))
    assert data.shape == (nv, 10 + 3 * ns)
    assert lengths.shape == (nv,)
    for v in range(nv):
        assert bytes(data[v, : lengths[v]]) == _block(gt[v], False)


@pytest.mark.parametrize(
    "ns,ploidy",
    [
        (0, 1),
        (0, 2),
        (1, 1),
        (4, 2),
        (21841, 2),
        (21842, 2),
        (32762, 1),
        (32763, 1),
        (43687, 2),
        (43688, 2),
    ],
)
def test_fixed_bytes_stored_boundaries_and_reader(tmp_path, ns, ploidy):
    args = _fixed_inputs(ns=ns, ploidy=ploidy)
    assert encode_layout2_fixed(*args) is None
    expected = b"".join(_variant(args, v) for v in range(2))
    assert args[8].tobytes() == expected
    assert args[8].shape[1] == layout2_fixed_size(ns, ploidy, 9)
    path = tmp_path / "fixed.bgen"
    header = struct.pack("<IIII4sI", 20, 20, 2, ns, b"bgen", 2 << 2 | 1)
    path.write_bytes(header + expected)
    with BgenFile(path) as reader:
        for v in range(2):
            # Position after variable identifying fields, just before C.
            offset = 24 + v * args[8].shape[1] + 20 + 9
            result = reader.read_genotype(offset)
            np.testing.assert_array_equal(result.ploidy, np.full(ns, ploidy))
            if ns:
                assert result.missing[0]
                np.testing.assert_array_equal(result.missing[1:], False)
            assert result.phased == bool(v)
    # 1D caller capacity may exceed the encoded region; suffix remains untouched.
    args[8] = np.full(len(expected) + 11, 0xAB, dtype=np.uint8)
    encode_layout2_fixed(*args)
    assert bytes(args[8][:-11]) == expected
    np.testing.assert_array_equal(args[8][-11:], 0xAB)


def test_fixed_empty_strings_nuls_and_positions():
    args = _fixed_inputs(ns=1, nv=1)
    args[:5] = [
        _field(b""),
        _field("é".encode()),
        _field(b"x\0ignored"),
        _field(b"A"),
        _field(b"T"),
    ]
    args[10] = 5
    args[5][:] = np.iinfo(np.int32).max
    args[8] = np.empty((1, layout2_fixed_size(1, 2, 5)), dtype=np.uint8)
    encode_layout2_fixed(*args)
    assert args[8].tobytes() == _variant(args, 0)
    args[5][:] = -1
    with pytest.raises(ValueError, match="position must be non-negative"):
        encode_layout2_fixed(*args)
    args = _fixed_inputs(ns=1, nv=1)
    args[0] = np.full((1, 65536), ord("a"), dtype=np.uint8)
    with pytest.raises(ValueError, match="string exceeds BGEN length"):
        encode_layout2_fixed(*args)


def test_fixed_validation_capacity_alias_and_missing():
    args = _fixed_inputs(ns=3, nv=2)
    args[8] = np.zeros(1, dtype=np.uint8)
    with pytest.raises(BufferTooSmall):
        encode_layout2_fixed(*args)
    args = _fixed_inputs()
    args[8].flags.writeable = False
    with pytest.raises(ValueError, match="writeable"):
        encode_layout2_fixed(*args)
    args = _fixed_inputs()
    args[6][1, 2, 1] = -2
    with pytest.raises(NotImplementedError, match="variant 1, sample 2.*mixed ploidy"):
        encode_layout2_fixed(*args)
    args = _fixed_inputs()
    args[0][:] = 0
    with pytest.raises(ValueError, match="string byte sum"):
        encode_layout2_fixed(*args)
    args = _fixed_inputs(ns=100)
    args[8] = np.empty(1000, dtype=np.uint8)
    args[6] = args[8][:400].view(np.int8).reshape(2, 100, 2)
    with pytest.raises(ValueError, match="overlap"):
        encode_layout2_fixed(*args)
    args = _fixed_inputs(nv=0)
    args[8] = np.full(20, 7, dtype=np.uint8)
    encode_layout2_fixed(*args)
    np.testing.assert_array_equal(args[8], 7)


@pytest.mark.parametrize(
    "ns,p,strings",
    [
        (-1, 2, 1),
        (1, 0, 0),
        (1, 3, 0),
        (1, 2, -1),
        (2**32, 2, 0),
        (2**32 - 1, 2, 0),
    ],
)
def test_invalid_fixed_sizes(ns, p, strings):
    with pytest.raises(ValueError):
        layout2_fixed_size(ns, p, strings)


def test_legacy_allocation_limit():
    gt = np.zeros((2, 3, 2), dtype=np.int8)
    with pytest.raises(ValueError, match="max_output_bytes"):
        encode_layout2_hardcalls_legacy(
            gt, np.zeros(2, dtype=bool), max_output_bytes=37
        )


@pytest.mark.parametrize("index", range(9))
def test_fixed_strict_arrays(index):
    args = _fixed_inputs(ns=3, nv=2)
    args[index] = args[index].tolist()
    with pytest.raises(TypeError):
        encode_layout2_fixed(*args)
    args = _fixed_inputs(ns=3, nv=2)
    args[index] = args[index].astype(np.float64)
    with pytest.raises(ValueError, match="Wrong dtype"):
        encode_layout2_fixed(*args)
    args = _fixed_inputs(ns=3, nv=2)
    args[index] = args[index][::-1]
    with pytest.raises(ValueError, match="C-contiguous"):
        encode_layout2_fixed(*args)
    args = _fixed_inputs(ns=3, nv=2)
    args[index] = args[index][None]
    with pytest.raises(ValueError, match="dimension"):
        encode_layout2_fixed(*args)


def test_specialized_bad_phase_and_legacy_shapes():
    args = _fixed_inputs()
    args[7].view(np.uint8)[0] = 2
    with pytest.raises(ValueError, match="variant 0.*phased must be"):
        encode_layout2_fixed(*args)
    with pytest.raises(ValueError, match="variant 0.*phased must be"):
        encode_layout2_hardcalls_legacy(args[6], args[7])
    for gt in (np.zeros((1, 1), dtype=np.int8), np.zeros((1, 1, 1), dtype=np.int8)):
        with pytest.raises(ValueError):
            encode_layout2_hardcalls_legacy(gt, np.zeros(1, dtype=bool))
    for gt, phase in [
        (np.zeros((1, 1, 2), dtype=np.int16), np.zeros(1, dtype=bool)),
        (np.zeros((1, 1, 2), dtype=np.int8), np.zeros(1, dtype=np.uint8)),
    ]:
        with pytest.raises(ValueError, match="Wrong dtype"):
            encode_layout2_hardcalls_legacy(gt, phase)
    with pytest.raises(ValueError, match="phased.shape"):
        encode_layout2_hardcalls_legacy(
            np.zeros((1, 1, 2), dtype=np.int8), np.zeros(2, dtype=bool)
        )
    gt = np.zeros((1, 2, 2), dtype=np.int8)
    with pytest.raises(ValueError, match="C-contiguous"):
        encode_layout2_hardcalls_legacy(gt[:, :, ::-1], np.zeros(1, dtype=bool))
    gt.flags.writeable = False
    encode_layout2_hardcalls_legacy(gt, np.zeros(1, dtype=bool))


@pytest.mark.parametrize("which", ["fixed", "legacy"])
def test_specialized_threads_and_gil(which):
    args = _fixed_inputs(ns=8192, nv=16)

    def run():
        if which == "legacy":
            out, lengths = encode_layout2_hardcalls_legacy(args[6], args[7])
            return b"".join(bytes(row[:n]) for row, n in zip(out, lengths))
        current = list(args)
        current[8] = np.empty_like(args[8])
        encode_layout2_fixed(*current)
        return current[8].tobytes()

    expected = run()
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        assert all(result == expected for result in pool.map(lambda _: run(), range(8)))
    ready, request = threading.Event(), threading.Event()
    progress = []

    def worker():
        ready.set()
        request.wait()
        progress.append(True)

    thread = threading.Thread(target=worker)
    thread.start()
    ready.wait()
    interval = sys.getswitchinterval()
    try:
        sys.setswitchinterval(100)
        request.set()
        run()
        assert progress
    finally:
        sys.setswitchinterval(interval)
        thread.join()
