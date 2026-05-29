cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(
  bgen
  VERSION 4.2.0
  LANGUAGES C CXX
  DESCRIPTION "C library and Python bindings for the BGEN format (specifications 1.2 and 1.3)")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(BGEN_BUILD_TESTS "Build C tests" ON)
option(BGEN_BUILD_PYTHON "Build Python bindings" OFF)
option(BGEN_ENABLE_SANITIZERS "Enable sanitizers (Debug builds)" OFF)
option(BGEN_USE_SYSTEM_ZLIB "Use system-provided zlib" ON)
option(BGEN_USE_SYSTEM_ZSTD "Use system-provided zstd" ON)
option(BGEN_USE_SYSTEM_ATHR "Use system-provided athr (almosthere)" ON)

include(FetchContent)

# --- ZLIB ---
if(BGEN_USE_SYSTEM_ZLIB)
  find_package(ZLIB REQUIRED)
else()
  message(STATUS "Fetching zlib-ng via FetchContent ...")
  set(ZLIB_COMPAT ON CACHE BOOL "" FORCE)
  set(ZLIB_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
  set(ZLIBNG_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
  set(WITH_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(zlib-ng
    GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git
    GIT_TAG 2.3.3
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(zlib-ng)
  get_target_property(_zlib_real zlibstatic ALIASED_TARGET)
  if(_zlib_real)
    add_library(ZLIB::ZLIB ALIAS ${_zlib_real})
  else()
    add_library(ZLIB::ZLIB ALIAS zlibstatic)
  endif()
endif()

# --- zstd ---
if(BGEN_USE_SYSTEM_ZSTD)
  find_package(zstd REQUIRED)
else()
  message(STATUS "Fetching zstd via FetchContent ...")
  set(ZSTD_BUILD_PROGRAMS OFF CACHE BOOL "" FORCE)
  set(ZSTD_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
  set(ZSTD_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(ZSTD_BUILD_SHARED OFF CACHE BOOL "" FORCE)
  set(ZSTD_BUILD_STATIC ON CACHE BOOL "" FORCE)
  FetchContent_Declare(zstd
    GIT_REPOSITORY https://github.com/facebook/zstd.git
    GIT_TAG v1.5.7
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR build/cmake)
  FetchContent_MakeAvailable(zstd)
  if(NOT TARGET zstd::libzstd_static)
    add_library(zstd::libzstd_static ALIAS libzstd_static)
  endif()
endif()

# --- athr (almosthere) ---
if(BGEN_USE_SYSTEM_ATHR)
  find_package(athr REQUIRED)
  # Fix broken conda-forge athr package that hardcodes Xcode SDK paths
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
  message(STATUS "Fetching almosthere via FetchContent ...")
  FetchContent_Declare(almosthere
    GIT_REPOSITORY https://github.com/horta/almosthere.git
    GIT_TAG v2.0.2
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(almosthere)
  # ATHR::athr alias is defined by almosthere's CMakeLists.txt
endif()

# --- Compiler warnings ---
if(MSVC)
  set(BGEN_WARNING_FLAGS /W4)
else()
  set(BGEN_WARNING_FLAGS -Wall -Wextra -Wpedantic -Wstrict-prototypes
                         -Wshadow -Wconversion -Wmissing-prototypes)
endif()

# --- Sanitizers ---
if(BGEN_ENABLE_SANITIZERS AND NOT MSVC)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()

# --- bgen C library ---
add_library(
  bgen
  src/file.c
  src/genotype.c
  src/io.c
  src/layout1.c
  src/layout2.c
  src/metafile.c
  src/report.c
  src/samples.c
  src/variant.c
  src/partition.c
  src/bstring.c
  src/unzlib.c
  src/unzstd.c
  src/strdup.c)
add_library(BGEN::bgen ALIAS bgen)

include(GenerateExportHeader)
generate_export_header(
  bgen
  BASE_NAME BGEN
  INCLUDE_GUARD_NAME BGEN_EXPORT_H
  EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/bgen/export.h)

target_include_directories(
  bgen
  PUBLIC $<INSTALL_INTERFACE:include>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

set_target_properties(
  bgen PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})

target_compile_definitions(
  bgen PRIVATE
  $<$<BOOL:${WIN32}>:_CRT_SECURE_NO_WARNINGS>
  $<$<BOOL:${WIN32}>:_CRT_NONSTDC_NO_DEPRECATE>)

target_link_libraries(bgen PUBLIC ZLIB::ZLIB ATHR::athr)
if(TARGET zstd::libzstd_shared)
  target_link_libraries(bgen PUBLIC zstd::libzstd_shared)
elseif(TARGET zstd::libzstd_static)
  target_link_libraries(bgen PUBLIC zstd::libzstd_static)
endif()

target_compile_options(bgen PRIVATE ${BGEN_WARNING_FLAGS})

# --- Feature checks ---
if(NOT c_restrict IN_LIST CMAKE_C_COMPILE_FEATURES)
  message(WARNING "restrict feature is not supported")
  target_compile_definitions(bgen PUBLIC restrict=)
endif()

include(CheckCSourceCompiles)
check_c_source_compiles(
  "
  static void f(char const* err, ...) __attribute__((format(printf, 1, 2)));
  int main(void) {return 0;}
  "
  HAVE_ATTR_FORMAT)
target_compile_definitions(
  bgen PRIVATE $<$<BOOL:${HAVE_ATTR_FORMAT}>:HAVE_ATTR_FORMAT>)

include(TestBigEndian)
test_big_endian(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
  message(FATAL_ERROR "Sorry, big-endian platforms are not supported yet.")
endif()

# --- Install ---
include(GNUInstallDirs)

set(_bgen_install_targets bgen)
if(NOT BGEN_USE_SYSTEM_ZLIB)
  get_target_property(_zlib_real zlibstatic ALIASED_TARGET)
  if(_zlib_real)
    list(APPEND _bgen_install_targets ${_zlib_real})
  else()
    list(APPEND _bgen_install_targets zlibstatic)
  endif()
endif()
if(NOT BGEN_USE_SYSTEM_ZSTD)
  list(APPEND _bgen_install_targets libzstd_static)
endif()
if(NOT BGEN_USE_SYSTEM_ATHR)
  list(APPEND _bgen_install_targets athr)
endif()

install(
  TARGETS ${_bgen_install_targets}
  EXPORT bgen-targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/bgen/export.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/bgen)

install(
  EXPORT bgen-targets
  FILE bgen-targets.cmake
  NAMESPACE BGEN::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bgen)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  bgen-config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/bgen-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bgen
  NO_SET_AND_CHECK_MACRO)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/bgen-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/bgen-config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/bgen-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bgen)

# --- Tests ---
if(BGEN_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

# --- Python bindings ---
if(BGEN_BUILD_PYTHON)
  add_subdirectory(python)
endif()
