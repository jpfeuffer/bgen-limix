cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(cbgen LANGUAGES C CXX)

# Reuse a pre-built bgen C library if one is discoverable (e.g. installed into
# CMAKE_PREFIX_PATH / CONDA_PREFIX by a prior C-library build, or provided by a
# conda `libbgen` package). Otherwise fall back to building it from source from
# the parent repository so that PyPI sdist / cibuildwheel produce a
# self-contained wheel.
find_package(bgen QUIET)

# Verify compiler compatibility on Windows. If the pre-built library was built
# with GCC (e.g. MinGW) but the Python bindings are being built with MSVC (or vice
# versa), we CANNOT link them due to ABI / object format mismatch (LNK1143).
# In such cases, we discard the pre-built target and fall back to from-source.
set(_bgen_compat_found FALSE)
if(TARGET BGEN::bgen)
  set(_bgen_compat_found TRUE)
  if(WIN32)
    # Check imported location file extension or config properties
    get_target_property(_bgen_loc BGEN::bgen IMPORTED_LOCATION_RELEASE)
    if(NOT _bgen_loc)
      get_target_property(_bgen_loc BGEN::bgen IMPORTED_LOCATION)
    endif()
    
    if(_bgen_loc)
      # GCC uses .a extension on Windows, MSVC uses .lib
      if(MSVC AND _bgen_loc MATCHES "\\.a$")
        message(STATUS "cbgen: Discarding pre-built bgen because of compiler mismatch (Pre-built is GCC, but this build is MSVC).")
        set(_bgen_compat_found FALSE)
      elseif(NOT MSVC AND _bgen_loc MATCHES "\\.lib$")
        message(STATUS "cbgen: Discarding pre-built bgen because of compiler mismatch (Pre-built is MSVC, but this build is GCC).")
        set(_bgen_compat_found FALSE)
      endif()
    endif()
  endif()
endif()

if(_bgen_compat_found)
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
