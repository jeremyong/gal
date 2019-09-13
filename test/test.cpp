#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

// Cheat sheet
// TEST_SUITE_BEGIN defines the start of a scope of tests grouped under a suite
// TEST_SUITE_END ends the test suite scope
// TEST_CASE defines a test case and takes a string name
//   After the test name, the * operator can be used to add decorators
//     - skip(bool = true)
//     - may_fail(bool = true)
//     - should_fail(bool = true)
//     - expected_failures(int)
//     - timeout(double)
//     - description("text")
// TEST_CASE_TEMPLATE takes a comma separated list of types after the test name and type parameter name
// SUBCASE defines a sub-test case within the scope of a test case

// CHECK checks a condition
// REQUIRE checks a condition and stops the test if the condition fails
// REQUIRE_FALSE ensures that a condition fails
// For most asserts though, use [REQUIRE|CHECK|WARN]_[EQ|NE|GT|LT|GE|LE|UNARY|UNARY_FALSE] as these compile faster.

// For floating point comparison, there is the handy helper `doctest::Approx(2, 1)` which
// when compared to another float, checks the comparison to within a percentage-based tolerance.
// You can specify the tolerance using `doctest::Approx::epsilon(float percentage)`.

// Command-line
// -ltc list-test-cases
// -lts list-test-suites
// -tc <filters> (comma separated test cases)
// -tce negated -tc
// -ts <filters> (comma separated test suites)
// -tse negated -ts
// -rs <seed> randomizes test order
// -d time each test and print duration