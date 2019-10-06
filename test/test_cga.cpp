#include "test_util.hpp"

#include <doctest/doctest.h>
#include <formatters.hpp>
#include <cga.hpp>

using namespace gal;
using namespace gal::cga;

TEST_SUITE_BEGIN("conformal-geometric-algebra");

TEST_CASE("metric")
{
    SUBCASE("timelike-units")
    {
        multivector<void, term<element<0b10000>, monomial<one>>> m;
        static_assert(std::is_same<decltype(m >> m), multivector<void, term<element<0>, monomial<minus_one>>>>::value);

        // Point at (1, 2, 1)
        point_t<1, 2, 1> p1;
        // In the Conformal Geometric Algebra, the contraction of a point onto itself is exactly 0
        static_assert(std::is_same<decltype(p1 >> p1), multivector<void>>::value);
    }
}

TEST_SUITE_END();