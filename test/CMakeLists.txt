include(FetchContent)

FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/onqtam/doctest.git
    GIT_TAG 2.3.5
)
FetchContent_MakeAvailable(doctest)

list(APPEND CMAKE_MODULE_PATH ${doctest_SOURCE_DIR}/scripts/cmake)

add_executable(gal_test
    test.cpp
    test_algebra.cpp
    test_algorithm.cpp
    test_cga.cpp
    test_ega.cpp
    test_dfa.cpp
    test_pga.cpp)

if (GAL_TEST_IK_ENABLED)
    # target_sources(gal_test PUBLIC test_ik.cpp)
endif()


target_link_libraries(gal_test PRIVATE gal doctest)
target_compile_definitions(gal_test PRIVATE
    GAL_DEBUG
    DOCTEST_CONFIG_SUPER_FAST_ASSERTS # uses a function call for asserts to speed up compilation
    DOCTEST_CONFIG_USE_STD_HEADERS # prevent non-standard overloading of std declarations
    DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS # enable doctest::Approx() to take any argument explicitly convertible to a double
    )
# Uncomment for profiling with clang
if (GAL_PROFILE_COMPILATION_ENABLED)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_compile_options(gal_test PRIVATE -ftime-trace)
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(gal_test PRIVATE -ftime-report -ftime-report-details -Q)
    endif()
endif()

include(doctest)

doctest_discover_tests(gal_test)
