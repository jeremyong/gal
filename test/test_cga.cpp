#include "test_util.hpp"

#include <doctest/doctest.h>
#include <gal/cga.hpp>

using namespace gal;
using namespace gal::cga;

TEST_SUITE_BEGIN("conformal-geometric-algebra");

TEST_CASE("null-basis-conversion")
{
    SUBCASE("to-natural-basis")
    {
        auto p = point<>::ie(0);
        auto pn = gal::detail::to_natural_basis(p);
        // TODO: verify
    }
}

TEST_CASE("point-norm")
{
    point<float> p{3.9f, 1.2f, -29.f};
    auto p_norm = compute([](auto p) { return p >> p; }, p);
    // Notice that we can compute that the norm of a point is exactly zero at compile time.
    static_assert(p_norm.size() == 0);
    CHECK_EQ(p_norm.size(), 0);
}

TEST_SUITE_END();