# bgen

C library and Python bindings for reading [BGEN files](https://www.well.ox.ac.uk/~gav/bgen_format/) and encoding layout-2 genotype payloads (format specifications 1.2 and 1.3).

This repository is a maintained fork of the original limix/bgen and limix/cbgen projects, merged into a single codebase with a modern build system, nanobind-based Python bindings, and S3 support.

## Features

- Fast C library for parsing BGEN genotype files
- Python bindings via [nanobind](https://github.com/wjakob/nanobind)
- Supports BGEN 1.2 and 1.3 format specifications
- zlib and zstd decompression support
- Metafile indexing for efficient random access
- Native bounded-batch layout-2 probability and hard-call encoding

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

For a self-contained local CPython 3.12 wheel, including the private nanobind
backend and bundled native libraries:

```bash
pixi run -e wheel python -m pip wheel --no-deps --no-build-isolation -w dist \
  -Cbuild-dir=build-wheel \
  -Ccmake.define.CMAKE_DISABLE_FIND_PACKAGE_bgen=ON \
  -Ccmake.define.BGEN_NANOBIND_BACKEND=private -Cwheel.py-api= .
```

The empty `wheel.py-api` override is intentional: a private backend is specific
to its build interpreter, so its wheel must be tagged `cp312-cp312`, not `abi3`.
Bundled dependencies' license notices are installed with the wheel.

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

## Encoding layout-2 genotype blocks (4.15 API)

The encoder returns **uncompressed genotype payloads**, not complete BGEN
files. File headers, variant identifiers, length prefixes and compression
remain the caller's responsibility. It does not depend on VCZ or interpret
format-specific NaN sentinel bit patterns.

```python
import numpy as np
from cbgen import encode_layout2, encode_layout2_hardcalls

payload, offsets = encode_layout2(
    allele_counts=np.array([2], dtype=np.uint16),
    ploidy=np.array([[2, 1]], dtype=np.uint8),
    missing=np.array([[False, True]], dtype=bool),
    phased=np.array([False], dtype=bool),
    gp=np.array([[[0.5, 0.25, 0.25], [np.nan, np.nan, np.nan]]]),
    bits=16,
)
first_block = payload[int(offsets[0]):int(offsets[1])]

hard_payload, hard_offsets = encode_layout2_hardcalls(
    genotypes=np.array([[[0, 1], [-1, -2]]], dtype=np.int8),
    allele_counts=np.array([2], dtype=np.uint16),
    phased=np.array([False], dtype=bool),
)
```

Both interfaces are keyword-only:

```text
encode_layout2(*, allele_counts, ploidy, missing, phased,
               gp=None, hp=None, bits=8, max_output_bytes=268435456)
encode_layout2_hardcalls(*, genotypes, allele_counts, phased,
                         ploidy=None, missing=None, bits=8,
                         max_output_bytes=268435456)
    -> (bytes, uint64[nvariants + 1] offsets)
```

All inputs must be aligned, native-endian, C-contiguous NumPy arrays with
matching variant/sample dimensions. No implicit dtype or layout conversions
are made. Read-only inputs are accepted; do not mutate inputs concurrently.
Offsets start at zero and end at the payload length; every variant starts
on a byte boundary. Empty batches, zero samples and zero ploidy are supported.
The native measurement and encoding passes release the GIL.

| Input | Dtype | Shape / meaning |
| --- | --- | --- |
| `allele_counts` | `uint16` | `(variants,)`, K in 1..65535 |
| `ploidy` | `uint8` | `(variants, samples)`, P in 0..63 |
| `missing` | `bool` | `(variants, samples)`, explicit whole-call missingness |
| `phased` | `bool` | `(variants,)`, one phase flag per variant |
| `gp` | `float32` or `float64` | `(variants, samples, gp_width)`, unphased states in BGEN colex order |
| `hp` | `float32` or `float64` | `(variants, samples, hp_width)`, flattened haplotype-major simplexes |
| `genotypes` | `int8`, `int16`, `int32` or `int64` | `(variants, samples, max_ploidy)`, allele indexes with -1 missing and -2 trailing fill |

For each unphased sample, the first `C(P+K-1,P)` GP entries are active.
For phased samples, the first `P*K` HP entries are active, arranged as P
simplexes of K alleles. Width must cover all samples, including missing ones.
Padding beyond active width is **ignored**, not validated. The caller must
validate any external format's fill/missing conventions. Only the selected
GP or HP field is read for a variant; an unused field may be omitted.

Active nonmissing simplexes must be finite, in [0,1], and sum to one within
absolute tolerance 1e-6 (no relative tolerance). They are normalized and
quantized with stable largest-remainder rounding; exact ties favor earlier
state indexes. The last state of each simplex is omitted. All bit depths
1..32 use little-endian least-significant-bit-first packing. Explicitly
missing probability samples skip numeric validation and store only zeros.
At P=0, unphased GP still has one state (which must be one if nonmissing);
phased HP has no states. Missingness remains explicit even at P=0.

For hard calls, ploidy is derived from trailing -2 fill. Optional explicit
ploidy must agree. Without an explicit `missing` array, any -1 allele marks
the whole sample missing. With one, a nonmissing call containing -1 is an
error; missing calls may contain valid alleles or -1, but all allele indexes
and trailing fill are still checked. Multiallelic unphased hard calls map to
the BGEN colex genotype rank. Hard calls are packed directly without
materializing a dense probability array.

The default allocation limit is 256 MiB of payload per batch, overridable by
`max_output_bytes`. Each individual block must also fit a BGEN uint32 length
field. Offsets and one variant descriptor per input variant are allocated
separately; quantization scratch is reused and scales with the largest active
simplex (approximately 24 bytes per state on 64-bit platforms). Chunk inputs
to bound both input and scratch memory. Invalid values, insufficient widths,
overflow and allocation limits raise errors; per-call data errors identify
the zero-based variant and sample. No partial payload is returned.

## C API

Link against the `bgen` library:

```cmake
find_package(bgen REQUIRED)
target_link_libraries(mytarget PRIVATE BGEN::bgen)
```

See `include/bgen/bgen.h` for the full C API.
The public `bgen/encode.h` declares `bgen_layout2_size` and
`bgen_layout2_encode`, supporting the same numeric and hard-call models.
Each `bgen_layout2_variant` contains explicit input buffer byte lengths and
padding width. Measure first, allocate the payload and `nvariants+1` uint64
offsets, then encode. Returns are 0 on success and -1 on error; an optional
per-call `bgen_encode_error` gives the location and message (`SIZE_MAX` means
not applicable). The size pass validates structure and hard calls; numeric
probabilities are validated by the encode pass. On encode failure, C output
buffers may be partially written and must be discarded. Output/input aliasing
and concurrent input mutation are unsupported. No global encoder state is used.

### Specialized biallelic compatibility and fixed-stride encoding

The following public APIs keep the specialized legacy biallelic hard-call
paths in the native library, rather than in downstream format adapters:

```text
encode_layout2_hardcalls_legacy(genotypes, phased, *, max_output_bytes=268435456)
    -> (uint8[variants, 10 + 3*samples] padded_blocks, uint32[variants] lengths)
layout2_fixed_size(num_samples, uniform_ploidy, total_string_length) -> int
encode_layout2_fixed(varid, rsid, chrom, allele1, allele2, position,
                     genotypes, phased, out, uniform_ploidy, total_string_length)
    -> None
```

`genotypes` must be contiguous `int8[variants,samples,2]`, `phased` is
`bool[variants]`. Haploid calls use -2 in slot 1; -1 in either allele makes
the whole call missing. These **compatibility** APIs reject -2 in slot 0,
including all-fill calls. Use the general hard-call API for P=0. Empty sample
axes retain the legacy declared ploidy range 2/2. Only `padded_blocks[v,
:lengths[v]]` is encoded; trailing bytes are unspecified. No Python payload
copy is required.

`encode_layout2_fixed` fuses genotype encoding, identifying-field framing and
zlib level-0 stored-block serialization into caller-owned output, with one
reusable native genotype scratch buffer. Every sample must have the declared
`uniform_ploidy` (1 or 2); mixed ploidy raises `NotImplementedError`.
String fields are contiguous `uint8[variants,width]` padded byte rows, ending
at the first NUL or row width. The sum of their actual lengths must equal
`total_string_length` **for every variant**; callers can pad an identifier
field to achieve a fixed total. Identifier lengths must fit uint16 and allele
lengths uint32. `position` is contiguous non-negative `int32[variants]`.
Invalid alleles, lengths, metadata and string budgets raise `ValueError`
with variant/sample context where applicable.

`out` must be a writable contiguous uint8 array, either 1D with at least
`variants * layout2_fixed_size(...)` capacity, or exactly
`(variants, layout2_fixed_size(...))`. Input/output overlap is rejected.
Insufficient 1D capacity raises `cbgen.BufferTooSmall`, a `ValueError`
subclass; malformed shapes/dtypes and read-only buffers raise `ValueError`.
On failure, caller-owned output may be partially written and must be
discarded. On success, any excess 1D capacity remains untouched. The output
contains complete identifying fields and C/D length framing, but **no file
header**; it requires layout-2/zlib flags in the caller's BGEN header.
Stored DEFLATE chunks are deterministically at most 65535 bytes, with the
RFC1950 Adler-32 supplied by the optimized zlib-ng dependency.

Both specialized calls release the GIL and validate and encode hard calls
in one native sample pass, without general combinatorial calculations or
dense probability arrays. The general `encode_layout2_hardcalls` also selects
the direct biallelic kernel for int8 inputs with K=2, bits=8 and width 1/2
when explicit ploidy/missing arrays are absent; all other inputs use the
general checked encoder. Its P=0 and missingness semantics do not change.
The public C counterparts and bounded buffer descriptors are declared in
`bgen/encode.h`.

Reproduce installed-kernel timings with
`pixi run -e wheel python scripts/benchmark_encode_hardcalls.py`.
The optional `--backend vcz` runs an installed older VCZ implementation in its
own environment for byte-hash and throughput comparisons; it is not a library
dependency. The specialized implementation is independently written against
the BGEN and RFC1950/1951 formats, preserving legacy observable behavior
without incorporating Apache-licensed VCZ implementation source.

## Acknowledgments

- [Marc Jan Bonder](https://github.com/Bonder-MJ) for bug-reporting and improvement suggestions.
- [Yan Wong](https://github.com/hyanwong) for bug-reporting and improvement suggestions.
- [Carl Kadie](https://github.com/CarlKCarlK) for bug-reporting and improvement suggestions.

## Authors

- [Danilo Horta](https://github.com/horta)
- [Julianus Pfeuffer](https://github.com/jpfeuffer)

## License

This project is licensed under the [MIT License](https://raw.githubusercontent.com/jpfeuffer/bgen-limix/main/LICENSE.md).
