"""Compare two Release builds through their exported native public APIs.

Generate deterministic synthetic files with --generate DIR --zstd LIBZSTD.
Then compare --old LIBBGEN --new LIBBGEN --manifest DIR/fixtures.json.
Fixture generation, process startup, output-buffer setup and warm-up are not timed.
"""

import argparse
import ctypes as ct
import json
import random
import statistics
import struct
import subprocess
import sys
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
        "bgen_file_read_genotypes_batch": (
            [ct.c_void_p, ct.c_void_p, ct.c_uint32, ct.c_uint32, ct.c_void_p],
            ct.c_int,
        ),
        "bgen_file_read_ncombs_batch": (
            [ct.c_void_p, ct.c_void_p, ct.c_uint32, ct.c_void_p],
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
    n = case["nsamples"] * case["ncombs"]
    result = (ct.c_double * (n * (case["nvariants"] if mode == "batch" else 1)))()
    counts = (ct.c_uint32 * case["nvariants"])()
    index = (case["path"] + ".metafile").encode()

    def run():
        if mode == "single":
            for offset in offsets:
                genotype = lib.bgen_file_open_genotype(handle, offset)
                assert genotype
                assert lib.bgen_genotype_read64(genotype, result) == 0
                lib.bgen_genotype_close(genotype)
        elif mode == "batch":
            assert (
                lib.bgen_file_read_genotypes_batch(
                    handle, offsets, case["nvariants"], case["ncombs"], result
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
    print(json.dumps({"seconds": elapsed, "checksum": sum(result[:16])}))


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
    if not args.worker and (not args.old or not args.new):
        parser.error("comparison requires --old and --new")
    if args.worker and (
        args.case is None or args.mode not in ("single", "batch", "scan", "index")
    ):
        parser.error("worker requires --case and a valid --mode")
    cases = json.loads(args.manifest.read_text())
    if args.worker:
        worker(args.worker, cases[args.case], args.mode, args.loops)
        return
    report = []
    for i, case in enumerate(cases):
        for mode in ["single", "batch", "scan", "index"]:
            times = {"old": [], "new": []}
            checksums = set()
            for repeat in range(args.repeats):
                order = ["old", "new"] if repeat % 2 == 0 else ["new", "old"]
                for label in order:
                    command = [
                        sys.executable,
                        __file__,
                        "--worker",
                        getattr(args, label),
                        "--manifest",
                        str(args.manifest),
                        "--case",
                        str(i),
                        "--mode",
                        mode,
                        "--loops",
                        str(args.loops),
                    ]
                    measurement = json.loads(
                        subprocess.check_output(command, text=True)
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
