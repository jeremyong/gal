# NOTE: A recent cmake is required due to the dependency on the FetchContent module
# used to fetch doctest. If you are not building the tests, you can (probably) relax
# this requirement.
cmake_minimum_required(VERSION 3.15)

project(gal LANGUAGES C CXX)

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(GAL_STANDALONE ON)
else()
  set(GAL_STANDALONE OFF)
endif()

# NOTE: When GAL is included as a transitive dependency, these switches are disabled
option(GAL_TESTS_ENABLED "Enable GAL test compilation" ON)
option(GAL_SAMPLES_ENABLED "Enable GAL samples compilation" ON)
option(GAL_FORMATTERS_ENABLED "Enable formatters for use with fmtlib" ON)
option(GAL_PROFILE_COMPILATION_ENABLED "Enable use of the compiler time trace facilities if available" OFF)
option(GAL_TEST_IK_ENABLED "Enable benchmark ik test compilation" ON)

# NEVER mutate global cmake state unless we are building as a standalone project
if (GAL_STANDALONE)
  # Export compile commands for IDE/language-server/editor support
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

  # Output to a few canonical directories for convenience
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif()

add_subdirectory(src)

if (GAL_TESTS_ENABLED AND GAL_STANDALONE)
  enable_testing()
  add_subdirectory(test)
endif()

if (GAL_SAMPLES_ENABLED AND GAL_STANDALONE)
  # add_subdirectory(samples)
endif()

