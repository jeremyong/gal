#include "test_util.hpp"

#include <doctest/doctest.h>
#include <ega.hpp>
#include <formatters.hpp>

using namespace gal;
using namespace gal::ega;

TEST_SUITE_BEGIN("euclidean-geometric-algebra");

TEST_CASE("rotors")
{
    multivector<void, term<element<0>, monomial<multiplier<1>, factor<degree<1>, 0>>>> e;
    multivector<void, term<element<0b1>, monomial<multiplier<1>, factor<degree<1>, 1>>>> e0;
    multivector<void, term<element<0b10>, monomial<multiplier<1>, factor<degree<1>, 2>>>> e1;
    multivector<void, term<element<0b100>, monomial<multiplier<1>, factor<degree<1>, 3>>>> e2;
    multivector<void, term<element<0b11>, monomial<multiplier<1>, factor<degree<1>, 4>>>> e01;
    multivector<void, term<element<0b110>, monomial<multiplier<1>, factor<degree<1>, 4>>>> e12;

    multivector<void, term<element<0>, monomial<multiplier<1>>>> b;
    multivector<void, term<element<0b1>, monomial<multiplier<1>>>> b0;
    multivector<void, term<element<0b10>, monomial<multiplier<1>>>> b1;
    multivector<void, term<element<0b100>, monomial<multiplier<1>>>> b2;

    SUBCASE("vector-reflection-through-vector")
    {
        auto v = e0 + e1;
        auto reflected = conjugate(b0, v);
        static_assert(std::is_same<decltype(e0 - e1), decltype(reflected)>::value);

        static_assert(std::is_same<decltype(-e0 - e1), decltype(conjugate(b2, v))>::value);
    }

    SUBCASE("blade-reflection-through-vector")
    {
        static_assert(std::is_same<decltype(e01), decltype(conjugate(b2, e01))>::value);
        static_assert(std::is_same<decltype(e01), decltype(conjugate(-b2, e01))>::value);
        static_assert(std::is_same<typename scale<2, decltype(e12)>::type, decltype(conjugate(b0 + b2, e01))>::value);
        fmt::print("{}\n", conjugate(b0 + b1 + b2, e01));
    }
}

TEST_SUITE_END();