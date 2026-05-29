C library + Python bindings for BGEN format (specs 1.2 & 1.3). C11/C++17, CMake + Ninja, pixi for env management.

Build & test: `pixi run test` (runs both C and Python test suites).
C only: `cmake --preset pixi -DCMAKE_BUILD_TYPE=Release && cmake --build build && ctest --test-dir build`.
Python: `pip install --no-build-isolation -e . && pytest python/tests`.

Key CMake options: `BGEN_ENABLE_S3` (needs libcurl + OpenSSL), `BGEN_BUILD_PYTHON`, `BGEN_BUILD_TESTS`.
System deps resolved via conda-forge (pixi): zlib-ng, zstd, almosthere, nanobind.
Pre-built build dirs: `build` (pixi preset), `build-s3` (S3 enabled), `build-no-s3`.
