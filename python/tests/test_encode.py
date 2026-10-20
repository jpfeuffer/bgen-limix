"""Public encoder tests with independently assembled layout-2 byte oracles."""

import concurrent.futures
import itertools
import math
import struct
import sys
import threading

import numpy as np
import pytest
from cbgen import BgenFile, encode_layout2, encode_layout2_hardcalls

from .test_native_layouts import _write_bgen


def _simplex(p, scale):
    values = [float(x) for x in p]
    total = sum(values)
    scaled = [x / total * scale for x in values]
    rounded = [math.floor(x) for x in scaled]
    order = sorted(range(len(values)), key=lambda j: (rounded[j] - scaled[j], j))
    for j in order[: scale - sum(rounded)]:
        rounded[j] += 1
    return rounded


def _oracle(k, ploidy, missing, phased, probabilities, bits):
    scale = (1 << bits) - 1
    stored, decoded = [], []
    for s, raw_p in enumerate(ploidy):
        p = int(raw_p)
        count = k if phased else math.comb(p + k - 1, p)
        groups = p if phased else 1
        decoded_row = []
        for h in range(groups):
            values = probabilities[s, h * count : (h + 1) * count]
            q = [0] * count if missing[s] else _simplex(values, scale)
            stored.extend(q[:-1])
            decoded_row.extend([np.nan] * count if missing[s] else np.array(q) / scale)
        decoded.append(decoded_row)
    # A single arbitrary-width integer is independent of native streaming packing.
    packed = sum(x << (bits * j) for j, x in enumerate(stored))
    prefix = struct.pack(
        "<IHBB", len(ploidy), k, min(ploidy, default=0), max(ploidy, default=0)
    )
    prefix += bytes(int(p) | (128 if m else 0) for p, m in zip(ploidy, missing))
    prefix += bytes([phased, bits])
    return prefix + packed.to_bytes((len(stored) * bits + 7) // 8, "little"), decoded


def _inputs(k=3, ploidies=(0, 1, 2, 4), dtype=np.float64):
    p = np.array([ploidies, ploidies], dtype=np.uint8)
    missing = np.zeros(p.shape, dtype=bool)
    if len(ploidies) > 1:
        missing[:, 1] = True
    width = max(math.comb(int(x) + k - 1, int(x)) for x in ploidies)
    hp_width = max(ploidies) * k
    gp = np.full((2, len(ploidies), width), np.nan, dtype=dtype)
    hp = np.full((2, len(ploidies), hp_width), np.nan, dtype=dtype)
    rng = np.random.default_rng(42)
    for v, s in itertools.product(range(2), range(len(ploidies))):
        count = math.comb(int(p[v, s]) + k - 1, int(p[v, s]))
        values = rng.random(count)
        gp[v, s, :count] = values / values.sum()
        for h in range(int(p[v, s])):
            values = rng.random(k)
            hp[v, s, h * k : (h + 1) * k] = values / values.sum()
        if missing[v, s]:
            gp[v, s] = np.nan
            hp[v, s] = np.inf
    return {
        "allele_counts": np.full(2, k, dtype=np.uint16),
        "ploidy": p,
        "missing": missing,
        "phased": np.array([False, True]),
        "gp": gp,
        "hp": hp,
    }


@pytest.mark.parametrize("bits", range(1, 33))
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_probability_bytes_and_native_read(tmp_path, bits, dtype):
    args = _inputs(dtype=dtype)
    payload, offsets = encode_layout2(**args, bits=bits)
    assert isinstance(payload, bytes)
    assert offsets.dtype == np.uint64 and offsets.flags.c_contiguous
    assert offsets[0] == 0 and offsets[-1] == len(payload)
    blocks, expected = [], []
    for v in range(2):
        block, decoded = _oracle(
            3,
            args["ploidy"][v],
            args["missing"][v],
            bool(v),
            args["hp" if v else "gp"][v],
            bits,
        )
        assert payload[int(offsets[v]) : int(offsets[v + 1])] == block
        blocks.append(block)
        expected.append(decoded)
    path = tmp_path / "encoded.bgen"
    file_offsets = _write_bgen(path, 2, bits % 3, blocks, 4)
    with BgenFile(path) as reader:
        for v, offset in enumerate(file_offsets):
            result = reader.read_genotype(offset)
            for s, values in enumerate(expected[v]):
                np.testing.assert_allclose(
                    result.probability[s, : len(values)], values, rtol=0, atol=1e-15
                )
            np.testing.assert_array_equal(result.ploidy, args["ploidy"][v])
            np.testing.assert_array_equal(result.missing, args["missing"][v])
            assert result.phased == bool(v)


@pytest.mark.parametrize(
    "k,ploidies", [(1, (0, 1, 63)), (2, (0, 1, 63)), (8, (0, 1, 4)), (65535, (0, 1))]
)
def test_extreme_dimensions(k, ploidies):
    args = _inputs(k=k, ploidies=ploidies)
    payload, offsets = encode_layout2(**args, bits=23)
    for v in range(2):
        block, _ = _oracle(
            k,
            args["ploidy"][v],
            args["missing"][v],
            bool(v),
            args["hp" if v else "gp"][v],
            23,
        )
        assert payload[int(offsets[v]) : int(offsets[v + 1])] == block


def test_stable_ties_normalization_and_missing():
    args = _inputs(k=3, ploidies=(1, 1))
    args["gp"][:] = [1 / 3, 1 / 3, 1 / 3]
    args["hp"][:] = [1 / 3, 1 / 3, 1 / 3]
    args["missing"][:] = False
    payload, offsets = encode_layout2(**args, bits=1)
    # Scale=1: the tied first state wins, second and implicit third are zero.
    assert payload[int(offsets[1]) - 1] == 0b0101
    args["gp"][:] = [0.5, 0.25, 0.2500005]
    encode_layout2(**args)  # normalize a permitted small sum error
    args["gp"][:] = [0.5, 0.25, 0.250002]
    with pytest.raises(ValueError, match="variant 0, sample 0.*sum to 1"):
        encode_layout2(**args)


@pytest.mark.parametrize("dtype", [np.int8, np.int16, np.int32, np.int64])
@pytest.mark.parametrize("phased", [False, True])
@pytest.mark.parametrize("bits", [1, 7, 8, 15, 16, 23, 31, 32])
def test_hardcalls_colex_and_missing(dtype, phased, bits):
    gt = np.array(
        [
            [
                [-2, -2, -2, -2],
                [2, -2, -2, -2],
                [1, 0, -2, -2],
                [2, 0, 1, 2],
                [-1, 1, -2, -2],
            ]
        ],
        dtype=dtype,
    )
    p = np.array([[0, 1, 2, 4, 2]], dtype=np.uint8)
    missing = np.array([[False, False, False, False, True]])
    width = 12 if phased else 15
    probs = np.zeros((1, 5, width))
    for s, ploidy in enumerate(p[0]):
        alleles = gt[0, s, :ploidy]
        if missing[0, s]:
            continue
        if phased:
            for h, a in enumerate(alleles):
                probs[0, s, h * 3 + a] = 1
        else:
            states = sorted(
                itertools.combinations_with_replacement(range(3), int(ploidy)),
                key=lambda a: tuple(reversed(a)),
            )
            probs[0, s, states.index(tuple(sorted(alleles)))] = 1
    metadata = {
        "allele_counts": np.array([3], dtype=np.uint16),
        "phased": np.array([phased]),
    }
    expected = encode_layout2(
        **metadata,
        ploidy=p,
        missing=missing,
        bits=bits,
        **{"hp" if phased else "gp": probs},
    )
    actual = encode_layout2_hardcalls(**metadata, genotypes=gt, bits=bits)
    assert actual[0] == expected[0]
    np.testing.assert_array_equal(actual[1], expected[1])
    explicit = encode_layout2_hardcalls(
        **metadata, genotypes=gt, ploidy=p, missing=missing, bits=bits
    )
    assert explicit[0] == expected[0]


@pytest.mark.parametrize(
    "k,dtype", [(129, np.int16), (65535, np.int32), (65535, np.int64)]
)
def test_hardcall_large_allele_indexes(k, dtype):
    gt = np.array([[[k - 1], [k - 2]]], dtype=dtype)
    payload, _ = encode_layout2_hardcalls(
        genotypes=gt,
        allele_counts=np.array([k], dtype=np.uint16),
        phased=np.array([False]),
        bits=8,
    )
    assert len(payload) == 12 + 2 * (k - 1)
    assert payload[12 : 12 + k - 1] == bytes(k - 1)
    assert payload[-1] == 255
    assert payload[12 + k - 1 : -1] == bytes(k - 2)


@pytest.mark.parametrize("k", range(1, 6))
@pytest.mark.parametrize("ploidy", range(6))
def test_exhaustive_small_colex_states(k, ploidy):
    states = sorted(
        itertools.combinations_with_replacement(range(k), ploidy),
        key=lambda a: tuple(reversed(a)),
    )
    gt = np.array(states, dtype=np.int16).reshape(1, len(states), ploidy)
    payload, _ = encode_layout2_hardcalls(
        genotypes=gt,
        allele_counts=np.array([k], dtype=np.uint16),
        phased=np.array([False]),
        bits=8,
    )
    stored = payload[10 + len(states) :]
    matrix = np.frombuffer(stored, dtype=np.uint8).reshape(len(states), len(states) - 1)
    np.testing.assert_array_equal(
        matrix, np.eye(len(states), dtype=np.uint8)[:, :-1] * 255
    )


@pytest.mark.parametrize("value", [np.nan, np.inf, -0.01, 1.01])
def test_invalid_numeric(value):
    args = _inputs()
    args["gp"][0, 2, 0] = value
    with pytest.raises(ValueError, match=r"variant 0, sample 2.*finite.*\[0, 1\]"):
        encode_layout2(**args)


@pytest.mark.parametrize(
    "name", ["allele_counts", "ploidy", "missing", "phased", "gp", "hp"]
)
def test_strict_dtype_and_shape(name):
    args = _inputs()
    args[name] = np.zeros(args[name].shape, dtype=np.int32)
    with pytest.raises(TypeError):
        encode_layout2(**args)
    args = _inputs()
    args[name] = args[name].reshape(-1)
    if name in ("allele_counts", "phased"):
        args[name] = args[name][:-1]
    with pytest.raises(ValueError):
        encode_layout2(**args)


@pytest.mark.parametrize(
    "name", ["gp", "hp", "ploidy", "missing", "phased", "allele_counts"]
)
def test_strides_and_buffer_objects(name):
    args = _inputs()
    args[name] = args[name][::-1]
    with pytest.raises(ValueError, match="C-contiguous"):
        encode_layout2(**args)
    args = _inputs()
    args[name] = args[name].tolist()
    with pytest.raises(TypeError, match="NumPy array"):
        encode_layout2(**args)


def test_endianness_unaligned_and_readonly():
    args = _inputs()
    args["gp"] = args["gp"].astype(">f8" if sys.byteorder == "little" else "<f8")
    with pytest.raises(ValueError, match="native-endian"):
        encode_layout2(**args)
    args = _inputs()
    original = args["gp"]
    args["gp"] = np.ndarray(
        original.shape,
        dtype=np.float64,
        buffer=bytearray(original.nbytes + 1),
        offset=1,
    )
    args["gp"][:] = original
    with pytest.raises(ValueError, match="aligned"):
        encode_layout2(**args)
    args = _inputs()
    for array in args.values():
        array.flags.writeable = False
    encode_layout2(**args)


def test_contiguous_singleton_axis_with_zero_stride():
    args = _inputs()
    args = {name: value[:1] for name, value in args.items()}
    original = encode_layout2(**args)
    args["gp"] = np.lib.stride_tricks.as_strided(
        args["gp"], shape=args["gp"].shape, strides=(0, *args["gp"].strides[1:])
    )
    assert args["gp"].flags.c_contiguous
    actual = encode_layout2(**args)
    assert actual[0] == original[0]
    np.testing.assert_array_equal(actual[1], original[1])


@pytest.mark.parametrize("bits", [0, 33, -1, 256])
def test_invalid_bits(bits):
    with pytest.raises(ValueError, match="bits"):
        encode_layout2(**_inputs(), bits=bits)


def test_width_missing_arrays_and_limits():
    args = _inputs()
    args["gp"] = args["gp"][:, :, :1].copy()
    with pytest.raises(ValueError, match="variant 0, sample 1.*insufficient"):
        encode_layout2(**args)
    for name in ("gp", "hp"):
        args = _inputs()
        args.pop(name)
        with pytest.raises(ValueError, match=f"requires {name}"):
            encode_layout2(**args)
    args = _inputs()
    args["ploidy"][1, 0] = 64
    with pytest.raises(ValueError, match="variant 1, sample 0.*ploidy"):
        encode_layout2(**args)
    args = _inputs()
    args["allele_counts"][0] = 0
    with pytest.raises(ValueError, match="variant 0.*allele count"):
        encode_layout2(**args)
    args = _inputs()
    payload, _ = encode_layout2(**args)
    assert encode_layout2(**args, max_output_bytes=len(payload))[0] == payload
    with pytest.raises(ValueError, match="max_output_bytes"):
        encode_layout2(**args, max_output_bytes=len(payload) - 1)
    # This count overflows combinatorial arithmetic long before allocation.
    args["allele_counts"][0] = 65535
    args["ploidy"][0, 0] = 63
    with pytest.raises(ValueError, match="variant 0, sample 0.*overflow"):
        encode_layout2(**args)


@pytest.mark.parametrize(
    "gt,ploidy,missing,match",
    [
        ([0, -2, 1], None, None, "fill must be trailing"),
        ([0, -3, -2], None, None, "allele outside"),
        ([0, 3, -2], None, None, "allele outside"),
        ([0, -1, -2], None, False, "nonmissing"),
        ([0, 1, -2], 1, None, "ploidy disagrees"),
        ([0] * 64, None, None, "ploidy must"),
    ],
)
def test_invalid_hardcalls(gt, ploidy, missing, match):
    with pytest.raises(ValueError, match=f"variant 0, sample 0.*{match}"):
        encode_layout2_hardcalls(
            genotypes=np.array([[gt]], dtype=np.int32),
            allele_counts=np.array([3], dtype=np.uint16),
            phased=np.array([False]),
            ploidy=None if ploidy is None else np.array([[ploidy]], dtype=np.uint8),
            missing=None if missing is None else np.array([[missing]]),
        )


@pytest.mark.parametrize("dtype", [np.uint8, np.uint16, np.float64, bool])
def test_hardcall_dtype(dtype):
    with pytest.raises(TypeError, match="signed"):
        encode_layout2_hardcalls(
            genotypes=np.zeros((1, 1, 2), dtype=dtype),
            allele_counts=np.array([2], dtype=np.uint16),
            phased=np.array([False]),
        )


@pytest.mark.parametrize("phased", [False, True])
def test_empty_batches_samples_and_p0(phased):
    for nv, ns in [(0, 2), (2, 0), (2, 3)]:
        args = {
            "allele_counts": np.ones(nv, dtype=np.uint16),
            "phased": np.full(nv, phased),
            "ploidy": np.zeros((nv, ns), dtype=np.uint8),
            "missing": np.zeros((nv, ns), dtype=bool),
        }
        payload, offsets = encode_layout2(
            **args,
            **{"hp" if phased else "gp": np.ones((nv, ns, 0 if phased else 1))},
        )
        assert len(offsets) == nv + 1 and offsets[-1] == len(payload)
        hard, hard_offsets = encode_layout2_hardcalls(
            **args,
            genotypes=np.empty((nv, ns, 0), dtype=np.int8),
        )
        assert hard == payload
        np.testing.assert_array_equal(hard_offsets, offsets)


@pytest.mark.parametrize("phased", [False, True])
def test_explicit_missing_at_p0(phased):
    args = {
        "allele_counts": np.array([1], dtype=np.uint16),
        "phased": np.array([phased]),
        "ploidy": np.zeros((1, 2), dtype=np.uint8),
        "missing": np.array([[False, True]]),
    }
    data = np.array([[[1.0], [np.nan]]]) if not phased else np.empty((1, 2, 0))
    payload, _ = encode_layout2(**args, **{"hp" if phased else "gp": data})
    hard, _ = encode_layout2_hardcalls(
        **args, genotypes=np.empty((1, 2, 0), dtype=np.int8)
    )
    assert hard == payload
    assert payload[8:10] == b"\x00\x80"
    if not phased:
        data[0, 0, 0] = 0
        with pytest.raises(ValueError, match="variant 0, sample 0.*sum to 1"):
            encode_layout2(**args, gp=data)


def test_threads_and_gil_release():
    args = _inputs(k=8, ploidies=(4,) * 128)
    baseline = encode_layout2(**args, bits=32)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        outputs = list(pool.map(lambda _: encode_layout2(**args, bits=32), range(16)))
    for payload, offsets in outputs:
        assert payload == baseline[0]
        np.testing.assert_array_equal(offsets, baseline[1])

    # A worker must execute while the native call runs, not just around its edges.
    started = threading.Event()
    requested = threading.Event()
    counter = [0]

    def worker():
        started.set()
        requested.wait()
        counter[0] += 1

    old_interval = sys.getswitchinterval()
    thread = threading.Thread(target=worker)
    thread.start()
    started.wait()
    try:
        # No interpreter thread switch between requesting work and entering C.
        sys.setswitchinterval(100)
        requested.set()
        encode_layout2(**args, bits=32)
        assert counter[0] == 1
    finally:
        sys.setswitchinterval(old_interval)
        thread.join()
