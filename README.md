# bgen

C library and Python bindings for reading [BGEN files](https://www.well.ox.ac.uk/~gav/bgen_format/) (format specifications 1.2 and 1.3).

This repository is a maintained fork of the original limix/bgen and limix/cbgen projects, merged into a single codebase with a modern build system, nanobind-based Python bindings, and S3 support.

## Features

- Fast C library for parsing BGEN genotype files
- Python bindings via [nanobind](https://github.com/wjakob/nanobind)
- Supports BGEN 1.2 and 1.3 format specifications
- zlib and zstd decompression support
- Metafile indexing for efficient random access

## Quick start (pixi)

```bash
pixi install
pixi run test
```

## Building the C library

```bash
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

System dependencies: `zlib`, `zstd` (installed via your system package manager or conda).

## Installing the Python package

```bash
pip install .
```

Or in development mode:

```bash
pip install -e ".[test]"
pytest python/tests
```

## Python usage

```python
from cbgen import BgenFile, BgenMetafile

# Open a BGEN file
with BgenFile("example.bgen") as bgen:
    print(f"Variants: {bgen.nvariants}, Samples: {bgen.nsamples}")

    if bgen.contain_samples:
        samples = bgen.read_samples()

    # Create a metafile index for efficient access
    bgen.create_metafile("example.bgen.metafile")

# Read variant data via metafile
with BgenFile("example.bgen") as bgen:
    with BgenMetafile("example.bgen.metafile") as mf:
        partition = mf.read_partition(0)
        for variant in partition.variants:
            gt = bgen.read_genotype(variant.offset)
            print(f"{variant.rsid}: {gt.probability.shape}")
```

`bgen.layout` reports the native header's layout version (`1` or `2`), also
available on the legacy `bgen_file` API. Sample and variant counts are unsigned
32-bit values exposed as Python integers. Unsupported header flags and truncated
headers are rejected by the native reader.

Layout 1 supports raw and zlib-compressed, biallelic, diploid, unphased records,
including metafile indexing. Its three independently rounded probabilities are
decoded as `q / 32768` without renormalization; all-zero triplets are marked
missing and returned as NaNs. Layout 2 supports raw, zlib and zstd blocks with
1--32-bit probabilities, variable ploidy, and multiple alleles. Phased records
whose samples all have zero ploidy return shape `(nsamples, 0)`, retaining
their ploidy and missingness arrays. Both layouts support `precision=32` and
`precision=64`; returned Python arrays remain valid after subsequent reads and
after closing the file.

### Performance checks

`scripts/benchmark_decode.py` compares two Release shared-library builds through
their exported native public APIs. It generates deterministic, compressed,
biallelic diploid fixtures and measures single reads, batch reads, combination
scans and indexing. Timings exclude fixture generation, initial output-buffer
allocation and warm-up; old/new process order alternates across repetitions.
For example, on macOS:

```bash
pixi run python scripts/benchmark_decode.py --generate build-benchmark/fixtures \
  --zstd .pixi/envs/default/lib/libzstd.dylib
pixi run python scripts/benchmark_decode.py \
  --manifest build-benchmark/fixtures/fixtures.json \
  --old /path/to/before/libbgen.dylib --new /path/to/after/libbgen.dylib \
  --repeats 11 --loops 10 --output build-benchmark/results.json
```

Use matching compilers, optimization flags and compression libraries, and run
without competing builds/tests. Results are local, warm-cache measurements, not
a guarantee for other hardware, remote storage or arbitrary datasets.

## C API

Link against the `bgen` library:

```cmake
find_package(bgen REQUIRED)
target_link_libraries(mytarget PRIVATE BGEN::bgen)
```

See `include/bgen/bgen.h` for the full C API.

## Acknowledgments

- [Marc Jan Bonder](https://github.com/Bonder-MJ) for bug-reporting and improvement suggestions.
- [Yan Wong](https://github.com/hyanwong) for bug-reporting and improvement suggestions.
- [Carl Kadie](https://github.com/CarlKCarlK) for bug-reporting and improvement suggestions.

## Authors

- [Danilo Horta](https://github.com/horta)
- [Julianus Pfeuffer](https://github.com/jpfeuffer)

## License

This project is licensed under the [MIT License](https://raw.githubusercontent.com/jpfeuffer/bgen-limix/main/LICENSE.md).
