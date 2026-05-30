# Stable ABI (abi3) analysis for `bgen-limix`

## Scope
This report analyzes likely runtime impact of enabling nanobind stable ABI (`Py_LIMITED_API`) for this repository, based on the current implementation and PR #3 (`feat: nanobind stable ABI (abi3) wheels`).

## 1) Current nanobind / NumPy usage in this codebase

### What is used now
The Python extension is concentrated in `python/bindings.cpp` and is built via `nanobind_add_module` in `python/CMakeLists.txt`.

Key patterns:
- **No `nb::ndarray` objects are constructed/returned in the current code path**.
- Genotype output arrays are allocated with the **NumPy C API** (`PyArray_EMPTY`) and raw pointers are obtained via `PyArray_DATA`.
- Arrays are wrapped into nanobind objects with `nb::steal<nb::object>`.
- Return values are Python containers (`nb::dict`) plus scalars.
- Binding surface uses `NB_MODULE`, `nb::class_`, `.def*`, `nb::bytes`, `nb::object`, `nb::args`.

Relevant locations:
- NumPy allocation/data access: `python/bindings.cpp:26-55`
- Genotype hot path (`read_genotype`, `read_genotype32`): `python/bindings.cpp:163-239`
- Module/class bindings: `python/bindings.cpp:320-385`
- Current module build invocation: `python/CMakeLists.txt:4`

### Performance-relevant shape
In the hot methods, most work is native-side:
- `bgen_genotype_read64/read32` fills large probability buffers in C (`python/bindings.cpp:179,218`)
- per-sample ploidy/missing loops run in C++ (`python/bindings.cpp:187-190,226-229`)

So Python API overhead is mostly at array/object construction boundaries, not inner numeric loops.

## 2) APIs in this repo that are affected by `Py_LIMITED_API`

If `STABLE_ABI` is enabled in nanobind (as in PR #3), the following used APIs are in the affected boundary layer:

- nanobind binding/object APIs:
  - `NB_MODULE`, `nb::class_`, `.def*`, `nb::init`, `nb::dict`, `nb::bytes`, `nb::object`, `nb::steal`, `nb::args`
- NumPy C-API entry points used directly:
  - `_import_array()`, `PyArray_EMPTY`, `PyArray_DATA`

Important nuance for this repo:
- The common nanobind stable-ABI performance concern is often around **`nb::ndarray` metadata/data access and conversion paths**.
- This code currently avoids that path and uses NumPy C API allocation directly.
- Therefore, the most sensitive `nb::ndarray`-specific overhead is largely **not present** in current hot paths.

## 3) Expected runtime cost for this specific library

### Which operations are slower under `Py_LIMITED_API` here?

Likely slower (small absolute cost):
1. **Python object/binding boundary operations**
   - constructing return dicts and class-bound objects
   - attribute/property plumbing done by nanobind
2. **NumPy API call indirection at creation time**
   - `PyArray_EMPTY` and `PyArray_DATA` are called through API entry points
3. **Module import/type registration one-time cost**
   - `NB_MODULE` setup and class registration

Likely unchanged in dominant runtime:
- `bgen_genotype_read64/read32` decode and buffer fill
- C++ loops writing ploidy/missing arrays
- overall memory traffic for large genotype matrices

### Estimated impact ranges (for this codebase)

Given the implementation shape above, expected impact is:

- **`read_genotype` / `read_genotype32` on realistic variant chunks (moderate/large arrays):**
  - **~0% to 3%** typical
  - worst case usually still **<5%**
- **Very small calls dominated by Python object creation (tiny arrays, very high call counts):**
  - could reach **~5% to 15%** for those micro-cases
- **One-time module import cost:**
  - measurable but not user-relevant in steady state

Why these ranges are small here:
- heavy work is done in C/C++ over contiguous buffers
- array metadata/pointer extraction happens once per returned array, not per element
- no current `nb::ndarray`-heavy per-access path in tight loops

## 4) PR #3 change review

PR #3 (`feat: nanobind stable ABI (abi3) wheels`) changes:

- `python/CMakeLists.txt`
  - `nanobind_add_module(_core STABLE_ABI NB_STATIC ...)`
- `python/bindings.cpp`
  - adds `#define NPY_TARGET_VERSION NPY_2_0_API_VERSION`
- `pyproject.toml`
  - adds `wheel.py-api = "cp311"`
- `.github/workflows/wheels.yml`
  - builds only `cp311-*` wheels (abi3 wheel intended for 3.11+)

Net effect:
- broader wheel compatibility (one wheel per platform for 3.11+), less CI build matrix pressure
- expected runtime penalty in this repo is low because hot paths are native and not `nb::ndarray`-centric

## Recommendation: is the tradeoff worth it?

For this library, **yes, likely worth it**.

Reasoning:
- The project is primarily an I/O/decode + native-array fill workload.
- Current Python boundary usage is relatively thin.
- Estimated steady-state runtime impact is small for realistic usage.
- Packaging/ops gains are material: fewer wheels to build, test, publish, and maintain.

Caveat:
- If future changes move to heavy `nb::ndarray` metadata/data access in tight binding-layer loops, re-benchmark. That is where stable-ABI overhead can become more visible.

## Baseline validation notes from this investigation

- `pytest python/tests` passes locally (`6 passed, 4 skipped`).
- Full C test path via `cmake -B build ...` failed in this environment due to missing system `libcurl` dependency (environment issue, not due to report changes).
- GitHub workflow view for PR #3 branch showed:
  - `CI`: success
  - `Build wheels`: in progress (at inspection time)
