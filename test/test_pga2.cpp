#include "test_util.hpp"

#include <doctest/doctest.h>
#include <gal/pga2.hpp>
#include <gal/formatters.hpp>

using namespace gal::pga2;

TEST_SUITE_BEGIN("projective-2D-algebra");

TEST_CASE("incidence")
{
    SUBCASE("points-and-lines")
    {
        {
            // Construct two points
            point_t<0, 1> point1;
            point_t<1, 1> point2;

            // Join them to construct the line y = 1
            static_assert(std::is_same<decltype(point1 & point2), decltype(simplify(line_t<0, 1, -1>{}))>::value);
        }

        {
            point_t<1, 2> point1;
            point_t<4, 0> point2;
            // The slope of the line constructed here is -2/3
            constexpr auto line = point1 & point2;
            static_assert(std::is_same<gal::rational<-2, 3>, decltype(line_slope(line))>::value);
        }

        {
            // Join two lines at a point
            line_t<1, -1, 2> line1;
            line_t<-1, 0, 1> line2;

            constexpr auto point = line1 ^ line2;
            // The lines should meet at (1, 3)
            constexpr auto coords = cartesian_point(point);
            static_assert(coords.first == gal::one{});
            static_assert(coords.second == gal::rational<3>{});
        }
    }
}

TEST_SUITE_END();