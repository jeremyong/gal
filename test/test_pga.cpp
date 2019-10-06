#include "test_util.hpp"

#include <doctest/doctest.h>
#include <formatters.hpp>
#include <pga.hpp>
#include <engine.hpp>

using namespace gal::pga;

TEST_SUITE_BEGIN("projective-geometric-algebra");

TEST_CASE("primitives")
{
    SUBCASE("line-construction")
    {
        // A line is either a join of two points, or the wedge of two planes

        point_t<-1, 0, -1> point1;
        point_t<1, 0, 1> point2;
        // Make a line that cuts from the -xz quadrant to the +xz quadrant

        // Meet the xz plane and the xy plane rotated into the +xz quadrant
        // Extra scaling to account for unnormalized points
        plane_t<0, -1, 0, 1> plane1;
        plane_t<0, 0, 1, 0> plane2;

        // The two lines constructed this way are off by a scalar multiple
        static_assert(std::is_same<decltype(point1 & point2), gal::scale<-2, decltype(plane1 ^ plane2)>::type>::value);
    }

    SUBCASE("plane-construction")
    {
        SUBCASE("at-origin")
        {
            // Construct the xy plane from points
            point_t<1, 0, 0> point1;
            point_t<-1, 0, 0> point2;
            // First two points constructs the x-axis
            point_t<0, -1, 0> point3;
            auto plane = point1 & point2 & point3;

            // The plane returned satisfies the equation z = 0 (times a scalar factor)
            static_assert(std::is_same<decltype(gal::simplify(plane_t<0, 0, 0, -2>{})), decltype(plane)>::value);
        }

        SUBCASE("at-origin-2")
        {
            // Construct the xy plane from points
            point_t<1, 0, 0> point1;
            point_t<0, 1, 0> point2;
            point_t<1, 1, 0> point3;
            auto plane = point1 & point2 & point3;

            static_assert(std::is_same<decltype(gal::simplify(plane_t<0, 0, 0, 1>{})), decltype(plane)>::value);
        }

        SUBCASE("at-origin-3")
        {
            // Construct the yz plane from points
            point_t<0, 1, 1> point1; // -e012 + e013 + e123
            point_t<0, 1, 0> point2; // e013 + e123
            point_t<0, -1, 1> point3; // -e012 - e013 + e123
            auto plane = point1 & point2 & point3;
            static_assert(std::is_same<decltype(gal::simplify(plane_t<0, 2, 0, 0>{})), decltype(plane)>::value);
        }

        SUBCASE("translated-from-origin")
        {
            // Construct the yz plane from points offset one unit from the origin
            point_t<1, 1, 0> point1;
            point_t<1, 0, 0> point2;
            point_t<1, 1, 1> point3;
            auto plane = point1 & point2 & point3;
            static_assert(std::is_same<decltype(gal::simplify(plane_t<-1, 1, 0, 0>{})), decltype(plane)>::value);
        }
    }

    SUBCASE("distance-between-points")
    {
        point<> point1{0.0f, 0.0f, 1.0f};
        point<> point2{.70710678118f, .70710678118f, 1.0f};
        gal::engine engine{point1, point2};

        auto distance = engine.compute<gal::scalar<>>([](auto point1, auto point2) {
            auto l = point1 & point2;
            return l >> l;
        });
        CHECK_EQ(std::abs(distance), doctest::Approx(1.0f));
    }
}

TEST_SUITE_END();