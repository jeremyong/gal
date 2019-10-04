#include "test_util.hpp"

#include <doctest/doctest.h>
#include <pga.hpp>
#include <ring_generator.hpp>
#include <formatters.hpp>

using namespace gal::pga;

TEST_SUITE_BEGIN("projective-geometric-algebra");

TEST_CASE("primitives")
{
    SUBCASE("line-construction")
    {
        // A line is either a join of two points, or the wedge of two planes

        ipoint<-1, 0, -1> point1;
        ipoint<1, 0, 1> point2;
        // Make a line that cuts from the -xz quadrant to the +xz quadrant

        // Meet the xz plane and the xy plane rotated into the +xz quadrant
        // Extra scaling to account for unnormalized points
        iplane<0, -2, 0, -2> plane1;
        iplane<0, 0, 1, 0> plane2;

        static_assert(std::is_same<decltype(point1 | point2), decltype(plane1 ^ plane2)>::value);
    }

    SUBCASE("plane-construction")
    {
        SUBCASE("at-origin")
        {
            // Construct the xy plane from points
            ipoint<1, 0, 0> point1;
            ipoint<-1, 0, 0> point2;
            // First two points constructs the x-axis
            ipoint<0, -1, 0> point3;
            auto plane = point1 | point2 | point3;

            static_assert(std::is_same<decltype(gal::detail::simplify(iplane<0, 0, 0, 2>{})), decltype(plane)>::value);
        }

        SUBCASE("translated-from-origin")
        {
            // Construct the yz plane from points offset one unit from the origin
            ipoint<1, 1, 0> point1;
            ipoint<1, 0, 0> point2;
            ipoint<1, 1, 1> point3;
            auto plane = point1 | point2 | point3;
            static_assert(std::is_same<decltype(gal::detail::simplify(iplane<-1, 1, 0, 0>{})), decltype(plane)>::value);
        }
    }
}

TEST_SUITE_END();