"""Compare two Release builds through their exported native public APIs.

Generate deterministic synthetic files with --generate DIR --zstd LIBZSTD.
Then compare --old LIBBGEN --new LIBBGEN --manifest DIR/fixtures.json.
Fixture generation, process startup, output-buffer setup and warm-up are not timed.
"""

import argparse
import ctypes as ct
import json
import math
import os
import random
import statistics
import struct
import subprocess
import sys
import sysconfig
import time
import zlib
from pathlib import Path


def string(value, width):
    value = value.encode()
    return len(value).to_bytes(width, "little") + value


def fixtures(root, zstd_library, nsamples, nvariants):
    zstd = ct.CDLL(zstd_library)
    zstd.ZSTD_compressBound.argtypes = [ct.c_size_t]
    zstd.ZSTD_compressBound.restype = ct.c_size_t
    zstd.ZSTD_compress.argtypes = [
        ct.c_void_p,
        ct.c_size_t,
        ct.c_void_p,
        ct.c_size_t,
        ct.c_int,
    ]
    zstd.ZSTD_compress.restype = ct.c_size_t
    rng = random.Random(93275)
    cases = []
    for bits, phased, compression in [
        (8, False, 1),
        (8, False, 2),
        (16, False, 2),
        (32, False, 2),
        (8, True, 2),
    ]:
        name = f"{bits}bit-{'phased' if phased else 'unphased'}-{'zlib' if compression == 1 else 'zstd'}"
        path = root / f"{name}.bgen"
        denominator = (1 << bits) - 1
        pairs = [
            (denominator, 0),
            (0, denominator),
            (0, 0),
            (denominator * 3 // 4, denominator // 8),
            (denominator // 8, denominator * 3 // 4),
        ]
        encoded = [
            struct.pack({8: "<BB", 16: "<HH", 32: "<II"}[bits], *p) for p in pairs
        ]
        offsets = []
        with path.open("wb") as output:
            output.write(
                struct.pack(
                    "<IIII4sI", 20, 20, nvariants, nsamples, b"bgen", 8 | compression
                )
            )
            for variant in range(nvariants):
                output.write(string(f"snp{variant}", 2) + string(f"rs{variant}", 2))
                output.write(string("1", 2) + struct.pack("<IH", variant + 1, 2))
                output.write(string("A", 4) + string("T", 4))
                offsets.append(output.tell())
                data = struct.pack("<IHBB", nsamples, 2, 2, 2)
                data += bytes([2]) * nsamples + bytes([phased, bits])
                data += b"".join(
                    rng.choices(encoded, weights=[70, 15, 8, 5, 2], k=nsamples)
                )
                if compression == 1:
                    compressed = zlib.compress(data)
                else:
                    bound = zstd.ZSTD_compressBound(len(data))
                    dest = ct.create_string_buffer(bound)
                    size = zstd.ZSTD_compress(dest, bound, data, len(data), 3)
                    assert size <= bound
                    compressed = dest.raw[:size]
                output.write(struct.pack("<II", len(compressed) + 4, len(data)))
                output.write(compressed)
        cases.append(
            {
                "name": name,
                "path": str(path),
                "offsets": offsets,
                "nsamples": nsamples,
                "nvariants": nvariants,
                "ncombs": 4 if phased else 3,
            }
        )
    manifest = root / "fixtures.json"
    manifest.write_text(json.dumps(cases))
    return manifest


def worker(library, case, mode, loops):
    lib = ct.CDLL(library)
    signatures = {
        "bgen_file_open": ([ct.c_char_p], ct.c_void_p),
        "bgen_file_close": ([ct.c_void_p], None),
        "bgen_file_open_genotype": ([ct.c_void_p, ct.c_uint64], ct.c_void_p),
        "bgen_genotype_close": ([ct.c_void_p], None),
        "bgen_genotype_read64": ([ct.c_void_p, ct.c_void_p], ct.c_int),
        "bgen_genotype_read32": ([ct.c_void_p, ct.c_void_p], ct.c_int),
        "bgen_file_read_genotypes_batch": (
            [ct.c_void_p, ct.c_void_p, ct.c_uint32, ct.c_uint32, ct.c_void_p],
            ct.c_int,
        ),
        "bgen_file_read_ncombs_batch": (
            [ct.c_void_p, ct.c_void_p, ct.c_uint32, ct.c_void_p],
            ct.c_int,
        ),
        "bgen_file_read_genotypes_batch_padded": (
            [ct.c_void_p, ct.c_void_p, ct.c_uint32, ct.c_uint32, ct.c_void_p],
            ct.c_int,
        ),
        "bgen_metafile_create": (
            [ct.c_void_p, ct.c_char_p, ct.c_uint32, ct.c_int],
            ct.c_void_p,
        ),
        "bgen_metafile_close": ([ct.c_void_p], ct.c_int),
    }
    for name, (args, result) in signatures.items():
        function = getattr(lib, name)
        function.argtypes, function.restype = args, result
    handle = lib.bgen_file_open(case["path"].encode())
    assert handle
    offsets = (ct.c_uint64 * case["nvariants"])(*case["offsets"])
    width = case["ncombs"] + (2 if mode == "padded" else 0)
    n = case["nsamples"] * width
    scalar = ct.c_float if mode == "single32" else ct.c_double
    result = (
        scalar * (n * (case["nvariants"] if mode in ("batch", "padded") else 1))
    )()
    counts = (ct.c_uint32 * case["nvariants"])()
    index = (case["path"] + ".metafile").encode()

    def run():
        if mode in ("single", "single32"):
            for offset in offsets:
                genotype = lib.bgen_file_open_genotype(handle, offset)
                assert genotype
                decode = (
                    lib.bgen_genotype_read32
                    if mode == "single32"
                    else lib.bgen_genotype_read64
                )
                assert decode(genotype, result) == 0
                lib.bgen_genotype_close(genotype)
        elif mode == "batch":
            assert (
                lib.bgen_file_read_genotypes_batch(
                    handle, offsets, case["nvariants"], case["ncombs"], result
                )
                == 0
            )
        elif mode == "padded":
            assert (
                lib.bgen_file_read_genotypes_batch_padded(
                    handle, offsets, case["nvariants"], width, result
                )
                == 0
            )
        elif mode == "scan":
            assert (
                lib.bgen_file_read_ncombs_batch(
                    handle, offsets, case["nvariants"], counts
                )
                == 0
            )
            assert all(count == case["ncombs"] for count in counts)
        elif mode == "index":
            metafile = lib.bgen_metafile_create(handle, index, 1, 0)
            assert metafile
            assert lib.bgen_metafile_close(metafile) == 0

    run()
    start = time.perf_counter()
    for _ in range(loops):
        run()
    elapsed = (time.perf_counter() - start) / loops
    lib.bgen_file_close(handle)
    print(
        json.dumps(
            {
                "seconds": elapsed,
                "checksum": sum(x for x in result[:16] if not math.isnan(x)),
            }
        )
    )


def python_worker(package, case, mode, loops):
    """Run with -S and an explicit installed package, bypassing editable hooks."""
    sys.path.insert(0, str(package))
    from cbgen import BgenFile, bgen_file

    with BgenFile(case["path"]) as modern, bgen_file(case["path"]) as legacy:

        def run():
            if mode in ("single", "single32"):
                for offset in case["offsets"]:
                    result = modern.read_genotype(
                        offset, precision=32 if mode == "single32" else 64
                    ).probability
                return float(result[0].sum())
            if mode in ("batch", "padded"):
                width = case["ncombs"] + (2 if mode == "padded" else 0)
                result = legacy.read_genotypes_batch(case["offsets"], max_ncomb=width)
                return float(result[0, 0, : case["ncombs"]].sum())
            if mode == "structured":
                result = modern.read_genotypes_batch(
                    case["offsets"], max_ncomb=case["ncombs"], precision=64
                )
                return float(result.probability[0, 0].sum())
            if mode == "metadata":
                result = modern.read_genotype_metadata(case["offsets"])
                return int(result.ncombs.sum())
            if mode == "scan":
                result = legacy.get_ncombs(case["offsets"])
                return int(result.sum())
            raise ValueError(f"Unsupported Python benchmark mode: {mode}")

        checksum = run()
        start = time.perf_counter()
        for _ in range(loops):
            run()
        print(
            json.dumps(
                {
                    "seconds": (time.perf_counter() - start) / loops,
                    "checksum": checksum,
                }
            )
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--generate", type=Path)
    parser.add_argument("--zstd")
    parser.add_argument("--nsamples", type=int, default=20000)
    parser.add_argument("--nvariants", type=int, default=128)
    parser.add_argument("--worker")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--case", type=int)
    parser.add_argument("--mode")
    parser.add_argument("--loops", type=int, default=3)
    parser.add_argument("--old")
    parser.add_argument("--new")
    parser.add_argument(
        "--python-old", type=Path, help="Installed baseline package directory"
    )
    parser.add_argument(
        "--python-new", type=Path, help="Installed current package directory"
    )
    parser.add_argument("--python-worker", type=Path)
    parser.add_argument(
        "--modes",
        nargs="+",
        default=["single", "single32", "batch", "padded", "scan", "index"],
    )
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if min(args.nsamples, args.nvariants, args.loops, args.repeats) < 1:
        parser.error("sample, variant, loop and repetition counts must be positive")
    if args.generate:
        if not args.zstd:
            parser.error("--generate requires --zstd")
        args.generate.mkdir(parents=True, exist_ok=True)
        print(fixtures(args.generate, args.zstd, args.nsamples, args.nvariants))
        return
    if not args.manifest:
        parser.error("--manifest is required")
    python_mode = args.python_old is not None or args.python_new is not None
    if python_mode and not (args.python_old and args.python_new):
        parser.error("Python comparison requires both installed package directories")
    if (
        not args.worker
        and not args.python_worker
        and not python_mode
        and (not args.old or not args.new)
    ):
        parser.error("comparison requires --old and --new")
    worker_modes = {"single", "single32", "batch", "padded", "scan", "index"}
    if args.python_worker:
        worker_modes = worker_modes - {"index"} | {"structured", "metadata"}
    if (args.worker or args.python_worker) and (
        args.case is None or args.mode not in worker_modes
    ):
        parser.error("worker requires --case and a valid --mode")
    cases = json.loads(args.manifest.read_text())
    if args.worker:
        worker(args.worker, cases[args.case], args.mode, args.loops)
        return
    if args.python_worker:
        python_worker(args.python_worker, cases[args.case], args.mode, args.loops)
        return
    valid_modes = {"single", "single32", "batch", "padded", "scan", "index"}
    if not set(args.modes) <= valid_modes or (python_mode and "index" in args.modes):
        parser.error("Invalid modes (Python supports all modes except index)")
    report = []
    for i, case in enumerate(cases):
        for mode in args.modes:
            times = {"old": [], "new": []}
            checksums = set()
            for repeat in range(args.repeats):
                order = ["old", "new"] if repeat % 2 == 0 else ["new", "old"]
                for label in order:
                    command = [
                        sys.executable,
                        *(["-S"] if python_mode else []),
                        __file__,
                        "--python-worker" if python_mode else "--worker",
                        str(getattr(args, f"python_{label}" if python_mode else label)),
                        "--manifest",
                        str(args.manifest),
                        "--case",
                        str(i),
                        "--mode",
                        mode,
                        "--loops",
                        str(args.loops),
                    ]
                    env = None
                    if python_mode:
                        env = dict(
                            os.environ,
                            PYTHONPATH=os.pathsep.join(
                                {
                                    sysconfig.get_path("platlib"),
                                    sysconfig.get_path("purelib"),
                                }
                            ),
                        )
                    measurement = json.loads(
                        subprocess.check_output(command, text=True, env=env)
                    )
                    times[label].append(measurement["seconds"])
                    checksums.add(measurement["checksum"])
            assert len(checksums) == 1, (case["name"], mode, checksums)
            old, new = (statistics.median(times[label]) for label in ["old", "new"])
            row = {
                "case": case["name"],
                "mode": mode,
                "old_ms": old * 1000,
                "new_ms": new * 1000,
                "ratio": new / old,
                "samples": times,
            }
            report.append(row)
            print(
                f"{case['name']:24} {mode:6} {old * 1000:8.3f} -> {new * 1000:8.3f} ms "
                f"({(new / old - 1) * 100:+.1f}%)",
                flush=True,
            )
    if args.output:
        args.output.write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
