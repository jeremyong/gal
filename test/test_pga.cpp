#include <doctest/doctest.h>
#include <gal/pga.hpp>
#include <gal/format.hpp>
#include <gal/engine.hpp>

#include <iostream>

using namespace gal;
using namespace gal::pga;

TEST_SUITE_BEGIN("projective-geometric-algebra");

TEST_CASE("incidence")
{
    SUBCASE("point-to-line-construction")
    {
        point<> p1{0, 0, 1};
        point<> p2{1, 0, 1};

        auto line = compute([](auto p1, auto p2) { return p1 & p2; }, p1, p2);

        // auto l = evaluate<point<>, point<>>{}.debug([](auto p1, auto p2) { return p1 & p2; });

        std::cout << to_string(line) << std::endl;
    }
}

TEST_SUITE_END();