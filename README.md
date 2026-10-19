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

### Structured batch reads

The modern `BgenFile` supports native batches for both layouts and every
supported compression method, including mixed phase, allele counts and ploidy.
One batch can supply multiple conversion fields without rereading each variant:

```python
with BgenFile("example.bgen") as bgen:
    with BgenMetafile("example.bgen.metafile") as mf:
        offsets = [v.offset for v in mf.read_partition(0).variants]
    # Use bounded chunks, not the entire file, for large datasets.
    metadata = bgen.read_genotype_metadata(offsets)
    width = int(metadata.ncombs.max()) if len(offsets) else 0
    batch = bgen.read_genotypes_batch(offsets, max_ncomb=width, precision=64)
    probabilities = batch.probability  # (variants, samples, padded width)
    ploidy, missing = batch.ploidy, batch.missing
    for v in range(len(offsets)):
        row = probabilities[v, :, :batch.ncombs[v]]
        if batch.phased[v]:
            haplotypes = row.reshape(bgen.nsamples, -1, int(batch.nalleles[v]))
```

`GenotypeBatchResult` exposes columnar, C-contiguous, NumPy-owned arrays:

| Field | Shape | Dtype |
| --- | --- | --- |
| `probability` | `(V, S, W)` | float64 (default), or float32 |
| `phased` | `(V,)` | bool |
| `ploidy`, `missing` | `(V, S)` | uint8, bool |
| `ncombs` | `(V,)` | uint32 |
| `probability_count` | `(V, S)` | uint32 |
| `nalleles`, `nbits` | `(V,)` | uint16, uint8 |

`ncombs` is each variant's single-reader width, based on its declared maximum
ploidy, not the caller's padded width. `probability_count` is each sample's
structural width even when missing: `ploidy * nalleles` if phased, otherwise
`choose(ploidy + nalleles - 1, nalleles - 1)`. Zero-ploidy unphased samples have
one state; zero-ploidy phased samples have no probabilities. All-missing
variants retain their structural widths. Missing samples and unused columns
are NaN; distinguish them using `missing` and `probability_count`, not NaNs.
Phased probabilities are haplotype-major; unphased probabilities retain BGEN's
colex order. Layout-1 `nbits` is 16 and its exact independent rounding is
preserved.

Offsets preserve order and duplicates. `max_ncomb` is capacity, never
truncation: too-small or invalid widths raise an exception. An explicit width
reads/decompresses and decodes each selected block once, returning metadata from
that same read. `None` first scans widths, then rereads blocks to decode;
compressed layout-2 blocks therefore decompress twice, but probability values
are decoded only once. `read_genotype_metadata` returns `GenotypeBatchMetadata`
without allocating or decoding probabilities. It still decompresses compressed
layout-2 blocks and reads layout-1 triplets for missingness, and validates
framing/header/payload lengths, not numerical probability sums. `get_ncombs`
provides a narrower width-only scan.

Empty batches return `(0, S, W)`, with `W=0` for automatic sizing. Explicit
`W=0` is valid for wholly zero-ploidy phased batches. All arrays survive later
reads and reader closure. Compute dosage from float64 probability views before
casting output to float32 if rounding accuracy matters.

Output storage is `V*S*W*itemsize + 6*V*S + 8*V` bytes, including probability
padding and metadata, plus offsets and native reusable one-variant encoded/
decoded scratch. There is no second full probability tensor or Python object
per sample/variant. Choose `V` from a byte budget when `S` or `W` is large;
sample slicing is not supported. File handles are not thread-safe.
Legacy `bgen_file.read_genotypes_batch` remains probability-only float64 and
does not allocate the new metadata arrays.

### Lightweight dimensions and compact capacity (4.14.0)

```python
with BgenFile("example.bgen") as bgen:
    dimensions = bgen.read_dimensions(
        offsets, ploidy_source="actual", probability_model="auto"
    )
    batch = bgen.read_genotypes_batch(
        offsets, max_actual_ncomb=dimensions.decode_width, precision=64
    )
    # All output fields can reuse this one batch.
```

`GenotypeDimensions` contains Python integer maxima `gp_width`, `hp_width` and
`decode_width = max(gp_width, hp_width)`; native counts must fit uint32.
It does not return per-sample arrays. GP width zero means no unphased model
was selected; HP width zero means no phased model was selected. **HP width is
at least one if phased is selected**, even if all haplotypes have ploidy zero.
This preserves model presence for output schema creation; actual decoded phased
rows may still have zero probability entries. Empty selections return all zeros.

| Discovery strategy | Raw layout 2 | Zlib layout 2 | Zstd layout 2 |
| --- | --- | --- | --- |
| Actual + auto | Read header, sample ploidies and phase | Inflate prefix and scan ploidies | Load full validated block and scan ploidies |
| Declared + auto | Read header, seek over ploidies, read phase | Inflate/discard ploidies to reach phase | Load full validated block |
| Declared + explicit model | Read first eight bytes | Inflate first eight bytes | Load full validated block |

`ploidy_source="actual"` bounds dimensions using all sample ploidies, including
missing samples. `"declared"` uses the encoded maximum ploidy and may reserve
more space than any sample needs. `probability_model="auto"` reads original
phase. `"unphased"` or `"phased"` assumes that model, without validating it
during prefix discovery; `"mixed"` reserves both GP and HP bounds for every
record. The phase byte occurs **after** the sample-ploidy bytes, so it is not
available from the first eight bytes. Layout 1 always has intrinsic K=2/P=2,
so no genotype payload is needed; explicit phased/mixed assumptions still
reserve HP=4 (mixed also GP=3), and the full reader reports original unphasedness.

Optional `max_ploidy=bound` (0--63) substitutes a caller-provided capacity
assumption instead of either source of ploidy. Raw/zlib then skip scanning
sample ploidies: auto still reads phase, while explicit models only need the
header. This avoids overflow from an unused, loose declared maximum. The bound
is **not validated against samples in discovery**; consumers must check decoded
ploidy and model against their assumptions. With known allele counts, a supplied
ploidy bound and explicit model, consumers can calculate capacity and call
compact batch directly without invoking discovery at all.

Discovery is **not an integrity-validation pass**. Raw/zlib prefix strategies
may accept a truncated tail, skipped invalid ploidies/phase or a missing checksum
that full reads reject. Zstd deliberately uses
the full-block native header path; actual/override strategies validate compact
actual counts, while declared discovery validates declared counts too.
Layout-1 fixed dimensions do not validate payloads or authenticate index offsets.
Numerical probability sums are not validated by discovery or native decoding;
consumers requiring numerical validation must check decoded probabilities.

`max_actual_ncomb=W` is mutually exclusive with `max_ncomb`. It reads each block
once, validates the full native header/payload lengths and bounds **actual**
structural sample widths before writing that variant. Missing samples count
towards capacity. Too-small capacity raises, never truncates. Even overflowing
unused declared combination counts do not force a huge intermediate allocation.
In this mode only, `batch.ncombs[v]` is the effective maximum actual sample width,
not the declared maximum-ploidy width; it is zero for zero-sample variants.
Probabilities, precision, metadata ownership, ordering, duplicate offsets and
NaN padding otherwise retain their usual semantics. There is no width pre-scan
and no declared-width probability scratch; output and optional one-variant
decoded scratch are bounded by W. All default APIs keep their original behavior.

### Performance checks

`scripts/benchmark_decode.py` compares two Release shared-library builds through
their exported native public APIs. It generates deterministic, compressed,
biallelic diploid fixtures and measures single reads, batch reads, combination
scans and indexing, including float32 single reads and wider padded batches.
Timings exclude fixture generation, initial native output-buffer
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
To compare actual installed Python wheels with the same interpreter and
dependencies, use `--python-old /path/to/old-install --python-new
/path/to/new-install --modes single single32 batch padded scan` instead of
`--old`/`--new`. These measurements include NumPy output allocation and bypass
editable-install hooks.

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
