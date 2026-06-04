C library + Python bindings for BGEN format (specs 1.2 & 1.3). C11/C++17, CMake + Ninja, pixi for env management.

Build & test: `pixi run test` (runs both C and Python test suites).
C only: `cmake --preset pixi -DCMAKE_BUILD_TYPE=Release && cmake --build build && ctest --test-dir build`.
Python: `pixi run build-python && pixi run test-python`. The Python build (`python/CMakeLists.txt`, driven by scikit-build via `cmake.source-dir = "python"`) reuses a pre-installed bgen via `find_package(bgen)` when available (the `install-c` task installs it into `$CONDA_PREFIX`), and otherwise builds the C library from source for a self-contained wheel.

Key CMake options: `BGEN_ENABLE_S3` (needs libcurl + OpenSSL), `BGEN_BUILD_TESTS`.
System deps resolved via conda-forge (pixi): zlib-ng, zstd, almosthere, nanobind.
Pre-built build dirs: `build` (pixi preset), `build-s3` (S3 enabled), `build-no-s3`.
