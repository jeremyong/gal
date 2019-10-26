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

TEST_CASE("motors")
{
    SUBCASE("line-exp")
    {
        // Join two points to create a line
        point<> p1{1, 0, 1};
        point<> p2{1, 4, 1};
        // This will be a horizontal line in the +y direction passing through (1, 0, 1)
        auto l = compute([](auto p1, auto p2) { return p1 & p2; }, p1, p2);
        // motor<> axis{0, l[0], l[1], l[2], l[3], l[4], l[5], 0};
        // motor<> m = exp(axis);
        // Take the exponential to find a motor which should generate a rotation about the axis
        motor<> m{1.09165, 0.362261, -0.360923, -0.210585, -0.190524, 0.284162, 0.00798034, -0.0545484};
        std::cout << "m: " << to_string(m) << std::endl;
        auto logm = log(m);
        // std::cout << "logm: " << to_string(logm) << std::endl;
        auto explogm = exp(logm);
        std::cout << "m?: " << to_string(explogm) << std::endl;

        line<> logm2{0.3214595, -0.3244836, -0.1841874, -0.1691257, 0.2522469, 0.0079841};
        // std::cout << "logm2: " << to_string(m2) << std::endl;
        auto expm2 = exp(logm2);
        std::cout << "m2?: " << to_string(expm2) << std::endl;
    }
}

TEST_SUITE_END();