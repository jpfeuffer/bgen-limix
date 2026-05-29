# bgen

C library and Python bindings for reading [BGEN files](https://www.well.ox.ac.uk/~gav/bgen_format/) (format specifications 1.2 and 1.3).

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
