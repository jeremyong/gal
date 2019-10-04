cmake_minimum_required(VERSION 3.15)
project(gal LANGUAGES C CXX)

option(GAL_TESTS_ENABLED "Enable GAL test compilation" ON)
option(GAL_SAMPLES_ENABLED "Enable GAL samples compilation" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Requires libfmt
add_subdirectory(src)

if (GAL_TESTS_ENABLED)
  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup()

  enable_testing()
  add_subdirectory(test)
endif()

if (GAL_SAMPLES_ENABLED)
  add_subdirectory(samples)
endif()

