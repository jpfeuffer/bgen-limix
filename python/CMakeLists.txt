cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(cbgen LANGUAGES C CXX)

# Reuse a pre-built bgen C library if one is discoverable (e.g. installed into
# CMAKE_PREFIX_PATH / CONDA_PREFIX by a prior C-library build, or provided by a
# conda `libbgen` package). Otherwise fall back to building it from source from
# the parent repository so that PyPI sdist / cibuildwheel produce a
# self-contained wheel.
find_package(bgen QUIET)

if(TARGET BGEN::bgen)
  message(STATUS "cbgen: using pre-built bgen (BGEN::bgen target found)")
  # The conda-forge athr package hardcodes Xcode SDK paths in its link
  # interface (e.g. .../MacOSX11.0.sdk/usr/lib/libcurses.tbd). When we reuse a
  # pre-built bgen, ATHR::athr is pulled in transitively, so apply the same
  # fix the C build uses.
  if(TARGET ATHR::athr)
    get_target_property(_athr_libs ATHR::athr INTERFACE_LINK_LIBRARIES)
    if(_athr_libs)
      list(FILTER _athr_libs EXCLUDE REGEX "Xcode.*\\.tbd$")
      set_target_properties(ATHR::athr PROPERTIES INTERFACE_LINK_LIBRARIES "${_athr_libs}")
    endif()
    if(APPLE)
      find_library(CURSES_LIBRARY curses)
      if(CURSES_LIBRARY)
        target_link_libraries(ATHR::athr INTERFACE ${CURSES_LIBRARY})
      endif()
    endif()
  endif()
else()
  message(STATUS "cbgen: no pre-built bgen found, building C library from source")
  # Self-contained wheel: vendor all C dependencies.
  set(BGEN_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(BGEN_USE_SYSTEM_ZLIB_NG OFF CACHE BOOL "" FORCE)
  set(BGEN_USE_SYSTEM_ZSTD OFF CACHE BOOL "" FORCE)
  set(BGEN_USE_SYSTEM_ATHR OFF CACHE BOOL "" FORCE)
  set(BGEN_USE_SYSTEM_CURL OFF CACHE BOOL "" FORCE)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/.." "${CMAKE_BINARY_DIR}/_bgen")
endif()

# Development.SABIModule is required by nanobind split mode on Windows.
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy
  OPTIONAL_COMPONENTS Development.SABIModule)
find_package(nanobind CONFIG REQUIRED)

# "shared" links against the externally hosted nanobind-backend wheel/package
# so a single abi3 wheel covers Python 3.10+ on PyPI. "private" builds our own
# backend module instead: conda-forge builds one package per Python version
# regardless of abi3 (no shared-backend infra there), so a "shared" build
# would just add an unmet nanobind-backend dependency for no benefit; conda
# packaging should use "private".
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
target_include_directories(_core PRIVATE ${Python_NumPy_INCLUDE_DIRS})

install(TARGETS _core LIBRARY DESTINATION cbgen)
