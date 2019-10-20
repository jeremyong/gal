#include "test_util.hpp"

#include <doctest/doctest.h>

#include <gal/cga.hpp>
#include <gal/engine.hpp>


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

TEST_SUITE_END();