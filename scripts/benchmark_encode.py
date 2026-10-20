"""Time the installed public encoder; run with `pixi run -e wheel python`."""

import argparse
import hashlib
import json
import math
import statistics
import time

import cbgen
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=32)
    parser.add_argument("--samples", type=int, default=1024)
    parser.add_argument("--alleles", type=int, default=2)
    parser.add_argument("--ploidy", type=int, default=2)
    parser.add_argument("--bits", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--dtype", choices=["float32", "float64"], default="float64")
    args = parser.parse_args()
    if (
        args.variants < 1
        or args.samples < 1
        or args.repeats < 1
        or not 1 <= args.alleles <= 65535
        or not 0 <= args.ploidy <= 63
        or not 1 <= args.bits <= 32
    ):
        parser.error("invalid dimensions, repeats or bit depth")
    rng = np.random.default_rng(42)
    metadata = {
        "allele_counts": np.full(args.variants, args.alleles, dtype=np.uint16),
        "ploidy": np.full((args.variants, args.samples), args.ploidy, dtype=np.uint8),
        "missing": np.zeros((args.variants, args.samples), dtype=bool),
    }
    results = {}
    for model in ("gp", "hp"):
        phased = model == "hp"
        groups = args.ploidy if phased else 1
        states = (
            args.alleles
            if phased
            else math.comb(args.ploidy + args.alleles - 1, args.ploidy)
        )
        data = rng.random((args.variants, args.samples, groups, states))
        data /= data.sum(axis=-1, keepdims=True)
        data = data.reshape(args.variants, args.samples, groups * states).astype(
            args.dtype
        )
        inputs = dict(
            **metadata,
            phased=np.full(args.variants, phased),
            bits=args.bits,
            **{model: data},
        )
        payload, offsets = cbgen.encode_layout2(**inputs)
        timings = []
        for _ in range(args.repeats):
            start = time.perf_counter()
            current, bounds = cbgen.encode_layout2(**inputs)
            timings.append(time.perf_counter() - start)
            if current != payload or not np.array_equal(bounds, offsets):
                raise RuntimeError("encoder returned nondeterministic output")
        results[model] = {
            "median_seconds": statistics.median(timings),
            "payload_bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
        }
    print(
        json.dumps(
            {
                "version": cbgen.__version__,
                "parameters": vars(args),
                "results": results,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
