# Changelog

## 4.13.0 - 2026-10-19

### Added

- Modern `BgenFile.read_genotypes_batch(offsets, max_ncomb=None, precision=64)`
  returning `GenotypeBatchResult`: owned contiguous probabilities and columnar
  phase, ploidy, missingness, allele/bit counts, declared variant widths and
  structural per-sample probability lengths. Supports both layouts, every
  supported compression method, float32/64 and zero-width phased records.
- Public native structured float32/64 batch and metadata APIs sharing the
  single-record header/decoding core, with checked dimensions, bounded
  one-variant scratch, explicit padding capacity and NaN padding.
- `BgenFile.read_genotype_metadata` and `get_ncombs` for discovery without
  probability decoding. Compressed layout-2 blocks still decompress; automatic
  batch sizing scans widths before the decode pass.

### Fixed

- Support layout-1 raw/zlib and raw layout-2 records in existing native and legacy
  Python batch/count APIs, preserving probability-only legacy return values and
  avoiding metadata work on paths that do not request it.
- Reject invalid batch offsets, insufficient/overflowing widths and use after
  reader closure explicitly. Preserve requested ordering and duplicates, exact
  layout-1 rounding, original missingness and result ownership after close.

### Performance

- Keep probability-only and width-only paths free of new metadata allocation
  and copying. Reuse native decode/decompression scratch and initialize padded
  output with lazy contiguous NaN fills rather than slower per-row tail loops.
- Extend the reproducible comparison tool to float32 single reads, wider padded
  batches and actual installed Python packages, using pre-change 4.12.0 as the
  baseline for this API work.

## 4.12.0 - 2026-10-16

### Fixed

- Read uncompressed layout-2 genotype blocks using their encoded byte lengths,
  rather than a fixed six-bytes-per-sample assumption. Validate sample counts,
  ploidy ranges, bit depths and probability payload boundaries.
- Support layout-1 metafile creation and partition reads for raw and
  zlib-compressed files, preserving original genotype offsets and variant
  metadata and validating each record's sample count.
- Populate layout-1 native ploidy and missingness metadata, preventing NULL
  buffer access from Python genotype/probability reads. Layout-1 samples are
  diploid and unphased; all-zero probability triplets are missing. Preserve
  independently rounded probabilities as `q / 32768` without renormalization.
- Accept valid phased layout-2 records with zero ploidy for every sample,
  returning probability arrays of shape `(nsamples, 0)` while retaining
  missingness and ploidy metadata.
- Remove eight-bit truncation of unphased genotype combination loop counts,
  widen probability bit-offset arithmetic, and compute implied probabilities
  using integer residuals before floating-point conversion.
- Reject unsupported header flags, truncated records and decompressed-length
  mismatches explicitly.
- Expose native unsigned 32-bit sample and variant counts without signed
  narrowing in Python, and report the installed `bgen-limix` distribution
  version through `cbgen.__version__`.

### Added

- Public `layout` properties on `BgenFile` and legacy `bgen_file`, backed by
  the native parsed header.
- Generated native and Python regression tests covering raw/zlib/zstd blocks,
  1--32-bit probabilities, phase, multiple alleles, variable and zero ploidy,
  layout-1 indexing, both floating-point precisions, missingness and malformed
  inputs.
- A reproducible native performance comparison tool,
  `scripts/benchmark_decode.py`, with deterministic synthetic fixtures.

### Performance

- Fuse fixed-ploidy validation with metadata copying and precompute genotype
  combination counts outside per-sample decode loops, retaining bounds checks
  without the initial correctness patch's normal-use slowdown.
- On an Apple M2 Max, matched Release builds with 20,000 samples and 128
  variants measured end-to-end unphased Python reads approximately 4--6% faster
  than 4.11.0 for the tested 8/16/32-bit compressed workloads. Phased 8-bit
  reads were essentially unchanged; combination-only scans retained about
  0--3% validation overhead. These are local warm-cache measurements, not
  guarantees for other systems or datasets.
