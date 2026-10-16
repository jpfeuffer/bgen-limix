# Changelog

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
