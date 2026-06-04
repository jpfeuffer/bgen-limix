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

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_core STABLE_ABI NB_STATIC ${CMAKE_CURRENT_SOURCE_DIR}/bindings.cpp)
target_link_libraries(_core PRIVATE BGEN::bgen)
target_include_directories(_core PRIVATE ${Python_NumPy_INCLUDE_DIRS})

install(TARGETS _core LIBRARY DESTINATION cbgen)
