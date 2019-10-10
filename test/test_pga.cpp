#include "test_util.hpp"

#include <doctest/doctest.h>
#include <gal/engine.hpp>
#include <gal/formatters.hpp>
#include <gal/pga.hpp>

TEST_SUITE_BEGIN("projective-geometric-algebra");

TEST_CASE("primitives")
{
    using namespace gal::pga;
    SUBCASE("simple-line-construction")
    {
        point_t<0, 0, 0> origin;
        point_t<1, 0, 0> px;
        point_t<0, 1, 0> py;
        point_t<0, 0, 1> pz;

        // x-direction := e23
        // fmt::print("PGA line along x:");
        // print(origin & px);

        // // y-direction := e31
        // fmt::print("PGA line along y:");
        // print(origin & py);

        // // z-direction := e31
        // fmt::print("PGA line along z:");
        // print(origin & pz);
    }

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
        static_assert(std::is_same<decltype(point1 & point2), decltype(gal::rational<-2>{} * plane1 ^ plane2)>::value);
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
            static_assert(std::is_same<decltype(simplify(plane_t<0, 0, 0, -2>{})), decltype(plane)>::value);
        }

        SUBCASE("at-origin-2")
        {
            // Construct the xy plane from points
            point_t<1, 0, 0> point1;
            point_t<0, 1, 0> point2;
            point_t<1, 1, 0> point3;
            auto plane = point1 & point2 & point3;

            static_assert(std::is_same<decltype(simplify(plane_t<0, 0, 0, 1>{})), decltype(plane)>::value);
        }

        SUBCASE("at-origin-3")
        {
            // Construct the yz plane from points
            point_t<0, 1, 1> point1;  // -e012 + e013 + e123
            point_t<0, 1, 0> point2;  // e013 + e123
            point_t<0, -1, 1> point3; // -e012 - e013 + e123
            auto plane = point1 & point2 & point3;
            static_assert(std::is_same<decltype(simplify(plane_t<0, 2, 0, 0>{})), decltype(plane)>::value);
        }

        SUBCASE("translated-from-origin")
        {
            // Construct the yz plane from points offset one unit from the origin
            point_t<1, 1, 0> point1;
            point_t<1, 0, 0> point2;
            point_t<1, 1, 1> point3;
            auto plane = point1 & point2 & point3;
            static_assert(std::is_same<decltype(simplify(plane_t<-1, 1, 0, 0>{})), decltype(plane)>::value);
        }
    }

    SUBCASE("distance-between-points")
    {
        point<> point1{0.0f, 0.0f, 1.0f};
        point<> point2{.70710678118f, .70710678118f, 1.0f};
        gal::engine engine{point1, point2};

        scalar<> distance = engine.compute([](auto point1, auto point2) {
            auto l = point1 & point2;
            return l >> l;
        });
        CHECK_EQ(std::abs(distance), doctest::Approx(1.0f));
    }
}

TEST_CASE("rotors")
{
    using namespace gal::pga;
    using point = point<float>;
    using rotor = rotor<float>;
    using translator = translator<float>;
    SUBCASE("point-rotation-types")
    {
        point_t<0, 0, 1> p;
        rotor_t<0, 1, 0> r;
        // fmt::print("Rotating (0, 0, 1) about the y-axis gives: {}\n", conjugate(r, p));
    }

    SUBCASE("point-rotation")
    {
        // Produce a pi/2 radian rotation of a point (0, 0, 1) about the y-axis
        point p{0.0f, 0.0f, 1.0f};
        rotor r{M_PI * 0.5f, 0.0f, 1.0f, 0.0f};
        gal::engine engine{p, r};

        point p2 = engine.compute([](auto p, auto r) { return conjugate(r, p); });

        CHECK_EQ(p2.x, doctest::Approx(-1.0f));
        CHECK_EQ(p2.y, epsilon);
    }

    SUBCASE("point-translation")
    {
        // Move a point along the xy-line toward (1, 1, 1);
        point p{0.0f, 0.0f, 1.0f};
        translator t{std::sqrt(2.0f), 1.0f, 1.0f, 0.0f};
        t.normalize();

        gal::engine engine{p, t};

        point p2 = engine.compute([](auto p, auto t) { return conjugate(t, p); });

        CHECK_EQ(p2.x, doctest::Approx(1.0f));
        CHECK_EQ(p2.y, doctest::Approx(1.0f));
        CHECK_EQ(p2.z, doctest::Approx(1.0f));
    }
}

TEST_SUITE_END();