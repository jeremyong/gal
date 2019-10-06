add_executable(gal_test
    test.cpp
    test_algebra.cpp
    test_engine.cpp
    test_ega.cpp
    test_ga.cpp
    test_cga.cpp
    test_pga2.cpp
    test_pga.cpp)
target_link_libraries(gal_test PRIVATE gal)
target_include_directories(gal_test PRIVATE ${CONAN_INCLUDE_DIRS_DOCTEST})
target_compile_definitions(gal_test PRIVATE
    DOCTEST_CONFIG_SUPER_FAST_ASSERTS # uses a function call for asserts to speed up compilation
    DOCTEST_CONFIG_USE_STD_HEADERS # prevent non-standard overloading of std declarations
    DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS # enable doctest::Approx() to take any argument explicitly convertible to a double
    )

list(APPEND CMAKE_MODULE_PATH ${CONAN_LIB_DIRS_DOCTEST}/cmake/doctest)
include(doctest)

doctest_discover_tests(gal_test)