cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(cbgen LANGUAGES C CXX)

# Reuse a pre-built bgen C library if one is discoverable (e.g. installed into
# CMAKE_PREFIX_PATH / CONDA_PREFIX by a prior C-library build, or provided by a
# conda `libbgen` package). Otherwise fall back to building it from source from
# the parent repository so that PyPI sdist / cibuildwheel produce a
# self-contained wheel.
find_package(bgen 4.15 QUIET)

if(TARGET BGEN::bgen)
  message(STATUS "cbgen: using pre-built bgen (BGEN::bgen target found)")
else()
  message(STATUS "cbgen: no pre-built bgen found, building C library from source")
  # Self-contained wheel: vendor all C dependencies.
  set(BGEN_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(BGEN_USE_SYSTEM_ZLIB_NG OFF CACHE BOOL "" FORCE)
  set(BGEN_USE_SYSTEM_ZSTD OFF CACHE BOOL "" FORCE)
  set(BGEN_USE_SYSTEM_CURL OFF CACHE BOOL "" FORCE)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/.." "${CMAKE_BINARY_DIR}/_bgen")
endif()

# Development.SABIModule is required by nanobind split mode on Windows.
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy
  OPTIONAL_COMPONENTS Development.SABIModule)
find_package(nanobind CONFIG REQUIRED)

# Binary distributions must retain the binding runtime's license notices too.
if(EXISTS "${NB_DIR}/LICENSE")
  set(_nanobind_license "${NB_DIR}/LICENSE")
else()
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c
      "from importlib.metadata import distribution; d = distribution('nanobind'); print(next(str(d.locate_file(f)) for f in d.files if f.name == 'LICENSE'))"
    OUTPUT_VARIABLE _nanobind_license
    OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY)
endif()
file(TO_CMAKE_PATH "${_nanobind_license}" _nanobind_license)
install(FILES "${_nanobind_license}" DESTINATION cbgen/licenses RENAME nanobind.txt)
file(READ "${NB_DIR}/ext/robin_map/include/tsl/robin_map.h" _robin_header)
string(FIND "${_robin_header}" "*/" _robin_license_end)
if(_robin_license_end LESS 0)
  message(FATAL_ERROR "Cannot locate robin-map license notice")
endif()
math(EXPR _robin_license_length "${_robin_license_end} + 2")
string(SUBSTRING "${_robin_header}" 0 ${_robin_license_length} _robin_license)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/robin-map-license.txt" "${_robin_license}\n")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/robin-map-license.txt"
  DESTINATION cbgen/licenses RENAME robin-map.txt)

# "shared" links against the externally hosted nanobind-backend wheel/package,
# a separately versioned package that pip resolves per-interpreter, giving a
# single abi3 wheel that covers Python 3.10+ on PyPI. "private" builds the
# backend module into our own binary instead: this is NOT itself abi3 (the
# backend shim touches parts of the C API outside Py_LIMITED_API), so a
# "private" build only works for the exact Python version it was compiled
# against -- one conda package per Python version, same as any other
# extension module. conda-forge has no nanobind-backend package to depend on,
# so it uses "private" and pays for that with the usual per-version matrix.
set(BGEN_NANOBIND_BACKEND "shared" CACHE STRING "nanobind split-mode backend: shared (external nanobind-backend) or private (self-contained)")
set_property(CACHE BGEN_NANOBIND_BACKEND PROPERTY STRINGS shared private)

if(BGEN_NANOBIND_BACKEND STREQUAL "private")
  nanobind_add_backend(_bgen_nb_backend)
  # BACKEND_MODULE must match the module's actual import path: it's installed
  # inside the cbgen package (DESTINATION cbgen below), not at the top level.
  nanobind_add_module(_core STABLE_ABI BACKEND_MODULE cbgen._bgen_nb_backend ${CMAKE_CURRENT_SOURCE_DIR}/bindings.cpp)
  install(TARGETS _bgen_nb_backend LIBRARY DESTINATION cbgen)
else()
  nanobind_add_module(_core STABLE_ABI BACKEND_MODULE nanobind_backend ${CMAKE_CURRENT_SOURCE_DIR}/bindings.cpp)
endif()
target_link_libraries(_core PRIVATE BGEN::bgen)
target_sources(_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/encoding.cpp)
target_include_directories(_core PRIVATE ${Python_NumPy_INCLUDE_DIRS})

install(TARGETS _core LIBRARY DESTINATION cbgen)
