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
