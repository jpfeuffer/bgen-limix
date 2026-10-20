"""Isolate the installed specialized kernels, including allocation and binding costs."""

import argparse
import hashlib
import importlib.metadata
import json
import statistics
import time

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=["native", "vcz"], default="native")
    parser.add_argument("--variants", type=int, default=128)
    parser.add_argument("--samples", type=int, default=8192)
    parser.add_argument("--repeats", type=int, default=31)
    parser.add_argument("--ploidy", type=int, choices=[1, 2], default=2)
    parser.add_argument("--phased", action="store_true")
    args = parser.parse_args()
    if args.variants < 1 or args.samples < 1 or args.repeats < 1:
        parser.error("dimensions and repeats must be positive")
    if args.backend == "native":
        from cbgen import (
            encode_layout2_fixed as fixed,
        )
        from cbgen import (
            encode_layout2_hardcalls_legacy as variable,
        )
        from cbgen import (
            layout2_fixed_size as fixed_size,
        )

        version = importlib.metadata.version("bgen-limix")
    else:
        from vcztools._vcztools import (
            bgen_variant_block_size as fixed_size,
        )
        from vcztools._vcztools import (
            encode_bgen_chunk_slice_level0 as fixed,
        )
        from vcztools._vcztools import (
            encode_bgen_geno_blocks as variable,
        )

        version = importlib.metadata.version("vcztools")
    rng = np.random.default_rng(201026)
    gt = rng.integers(0, 2, (args.variants, args.samples, 2), dtype=np.int8)
    if args.ploidy == 1:
        gt[:, :, 1] = -2
    gt[:, ::17, 0] = -1
    phase = np.full(args.variants, args.phased)
    fields = [np.full((args.variants, 1), ord("."), dtype=np.uint8) for _ in range(5)]
    position = np.arange(args.variants, dtype=np.int32)
    row = fixed_size(args.samples, args.ploidy, 5)
    out = np.empty(args.variants * row, dtype=np.uint8)

    def run_variable():
        return variable(gt, phase)

    def run_fixed():
        fixed(*fields, position, gt, phase, out, args.ploidy, 5)

    variable_result = run_variable()
    run_fixed()
    timings = {"variable": [], "fixed": []}
    for j in range(args.repeats):
        calls = [("variable", run_variable), ("fixed", run_fixed)]
        if j % 2:
            calls.reverse()
        for name, function in calls:
            start = time.perf_counter()
            function()
            timings[name].append(time.perf_counter() - start)
    rows, lengths = variable_result
    variable_bytes = b"".join(bytes(r[:n]) for r, n in zip(rows, lengths))
    print(
        json.dumps(
            {
                "version": version,
                "parameters": vars(args),
                "results": {
                    name: {
                        "median_seconds": statistics.median(times),
                        "durations": times,
                        "sha256": hashlib.sha256(
                            variable_bytes if name == "variable" else out.tobytes()
                        ).hexdigest(),
                    }
                    for name, times in timings.items()
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
