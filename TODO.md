# TODO

## Expose compressed-block access so consumers can parallelise

`bgen_file_read_genotypes_batch` is serial. Its wins are real but orthogonal to
threading: seek elision on sorted offsets, reusable decompression scratch
buffers ("Plan A") and persistent decompressor contexts ("Plan B"). The header
already states one `bgen_file` handle per thread.

regenie does not use it. It reads raw compressed genotype blocks itself in one
serial pass, then inflates and bit-unpacks them OpenMP-parallel across a block
of variants. That two-phase split is where its 1.88x scaling from 1 to 4
threads comes from; a serial batch call cannot provide it.

Making the batch API internally parallel is the wrong fix: it would impose
OpenMP on every consumer, including the Python bindings, and force the library
to own per-thread decompressor contexts.

Better: expose the pieces, keep threading with the caller.

- an accessor giving the compressed block for an offset, plus its compressed
  and uncompressed lengths
- a pure decode entry point taking that buffer and a caller-supplied scratch
  or decompressor context, so N threads can decode concurrently

Then regenie's two-phase pattern is library-supported, and it can delete its
hand-written zlib/zstd inflate and probability bit-unpacker.

The stronger motivation is coverage rather than speed: regenie's own decoder
only implements layout-2 with 8-bit probabilities, so v1.1 files and 16/32-bit
probabilities silently fall back to the serial path. A library-provided block
API would let it use one code path for every encoding.

Measure against the existing baseline before and after: 109.2k variants/sec at
4 threads on 1M variants, run-to-run spread 0.15%.
