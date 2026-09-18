# TODO

## Drop the athr (almosthere) dependency

`athr` provides one thing: the progress bar shown while `bgen_metafile_create`
writes the index. It is used in a single file, `src/metafile_write.h`, at three
call sites:

```c
at = athr_create((long)nvariants, "Writing variants", ATHR_BAR | ATHR_ETA);
athr_consume(at, 1);
athr_finish(at);
```

For that, it costs:

- a **PUBLIC** dependency of `libbgen` (`target_link_libraries(bgen PUBLIC ... ATHR::athr)`),
  so it lands on every consumer's link line
- transitive curses + form + pthread
- the conda-forge `athr` package hardcodes its build machine's Xcode SDK path
  (`/Applications/Xcode_12.4.app/.../libcurses.tbd`) into `ATHR::athr`, which
  breaks downstream builds. We repair it in two places now — `CMakeLists.txt`
  for our own build, and `bgen-config.cmake.in` for consumers, because
  `find_dependency()` re-imports the target in the consumer's scope.
- a macOS `find_library(CURSES_LIBRARY curses)` patch
- `C11THREADS_NO_TIMED_MUTEX`, because almosthere v2.0.2 compiles with strict
  POSIX flags that hide `PTHREAD_MUTEX_TIMED_NP` on some manylinux/aarch64
  toolchains
- a FetchContent fallback path and an extra install target

Metafile creation takes about 2 s for 1M variants, so the bar earns very little.

**Plan:** replace with a few lines of `fprintf(stderr, ...)` gated on the
existing `verbose` parameter, then delete `BGEN_USE_SYSTEM_ATHR`, the
`find_package(athr)` block and all its workarounds, `ATHR::athr` from the link
line, `find_dependency(athr)` from `bgen-config.cmake.in`, and the install
target entry.

Two things to settle first:

1. It changes the progress output for anyone calling
   `bgen_metafile_create(..., verbose=1)` — including the Python bindings, where
   `create_metafile` is user-facing. Check whether any cbgen test asserts on it.
2. A TTY check needs a Windows shim (`_isatty`/`_fileno`), or just always print
   when `verbose` is set.

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
