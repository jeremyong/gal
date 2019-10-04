#include "test_util.hpp"

#include <doctest/doctest.h>
#include <ega.hpp>
#include <formatters.hpp>

using namespace gal;
using namespace gal::ega;

TEST_SUITE_BEGIN("euclidean-geometric-algebra");

TEST_CASE("rotors")
{
    SUBCASE("vector-reflection-through-vector")
    {
        auto v = e0{} + e1{};
        auto reflected = conjugate(e0{}, v);
        static_assert(std::is_same<decltype(e0{} - e1{}), decltype(reflected)>::value);

        static_assert(std::is_same<decltype(-e0{} - e1{}), decltype(conjugate(e2{}, v))>::value);
    }

    SUBCASE("blade-reflection-through-vector")
    {
        static_assert(std::is_same<e01, decltype(conjugate(e2{}, e01{}))>::value);
        static_assert(std::is_same<e01, decltype(conjugate(-e2{}, e01{}))>::value);
        static_assert(std::is_same<typename scale<2, e12>::type, decltype(conjugate(e0{} + e2{}, e01{}))>::value);
    }
}

TEST_SUITE_END();