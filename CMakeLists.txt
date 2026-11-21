cmake_minimum_required(VERSION 3.24)
project(terra LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(TERRA_BUILD_TESTS "Build unit tests" ON)
# e.g. -DTERRA_SANITIZERS=address,undefined
set(TERRA_SANITIZERS "" CACHE STRING "Comma-separated -fsanitize= list")

if(TERRA_SANITIZERS)
  add_compile_options(-fsanitize=${TERRA_SANITIZERS} -fno-omit-frame-pointer -fno-sanitize-recover=all)
  add_link_options(-fsanitize=${TERRA_SANITIZERS})
endif()

set(TERRA_WARNINGS -Wall -Wextra -Wpedantic -Werror)

add_library(terra_core src/ellipsoid.cpp src/tiles.cpp src/terrain.cpp src/camera.cpp src/select.cpp src/mesh.cpp src/render.cpp)
target_include_directories(terra_core PUBLIC include)
target_compile_options(terra_core PRIVATE ${TERRA_WARNINGS})

add_executable(terra tools/terra.cpp)
target_link_libraries(terra PRIVATE terra_core)
target_compile_options(terra PRIVATE ${TERRA_WARNINGS})

if(TERRA_BUILD_TESTS)
  include(FetchContent)
  # Vendored tarball when present (offline machines), else GitHub; same hash.
  set(_gtest_tarball ${CMAKE_CURRENT_SOURCE_DIR}/third_party/googletest-1.15.2.tar.gz)
  if(EXISTS ${_gtest_tarball})
    set(_gtest_url ${_gtest_tarball})
  else()
    set(_gtest_url https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz)
  endif()
  FetchContent_Declare(googletest
    URL ${_gtest_url}
    URL_HASH SHA256=7b42b4d6ed48810c5362c265a17faebe90dc2373c885e5216439d37927f02926
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  enable_testing()
  add_executable(terra_tests tests/geo_test.cpp tests/tiles_test.cpp tests/raster_test.cpp)
  target_link_libraries(terra_tests PRIVATE terra_core GTest::gtest_main)
  target_compile_options(terra_tests PRIVATE ${TERRA_WARNINGS})
  include(GoogleTest)
  gtest_discover_tests(terra_tests DISCOVERY_MODE PRE_TEST DISCOVERY_TIMEOUT 60)
  add_test(NAME engine_checks COMMAND terra --check)
endif()
