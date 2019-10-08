#include "test_util.hpp"

#include <doctest/doctest.h>
#include <formatters.hpp>
#include <cga.hpp>
#include <ega.hpp>
#include <pga.hpp>
#include <pga2.hpp>
#include <engine.hpp>

TEST_SUITE_BEGIN("engine");

TEST_CASE("basic-computation")
{
    SUBCASE("variadic-return")
    {
        using namespace gal::ega;

        point<> p1{2.4f, 3.6f, 4.2};
        point<> p2{-1.1f, 2.7f, 1.6};

        gal::engine engine{p1, p2};
        auto&& [q1, q2] = (std::tuple<point<>, point<>>)engine.compute([](auto p1, auto p2)
        {
            return std::tuple{p1, p2};
        });

        for (size_t i = 0; i != 3; ++i)
        {
            CHECK_EQ(p1[i], q1[i]);
            CHECK_EQ(p2[i], q2[i]);
        }
    }

    SUBCASE("pga2-expression-evaluation")
    {
        using namespace gal::pga2;

        point<> p1{2.4f, 3.6f};
        point<> p2{-1.1f, 2.7f};

        gal::engine engine{p1, p2};
        line<> l = engine.compute([](auto p1, auto p2)
        {
            return p1 & p2;
        });

        CHECK_EQ(p1.x * l.a + p1.y * l.b + l.c, epsilon);
        CHECK_EQ(p2.x * l.a + p2.y * l.b + l.c, epsilon);
    }

    SUBCASE("pga3-expression-evaluation")
    {
        using namespace gal::pga;
        point<double> p1{2.4f, 3.6f, 1.3f};
        point<double> p2{-1.1f, 11.7f, 5.0f};
        point<double> p3{-1.8f, -2.7f, -4.3f};

        gal::engine engine{p1, p2, p3};
        plane<double> p = engine.compute([](auto p1, auto p2, auto p3)
        {
            return p1 & p2 & p3;
        });

        CHECK_EQ(p1.x * p.x + p1.y * p.y + p1.z * p.z + p.d, epsilon);
        CHECK_EQ(p2.x * p.x + p2.y * p.y + p2.z * p.z + p.d, epsilon);
        CHECK_EQ(p3.x * p.x + p3.y * p.y + p3.z * p.z + p.d, epsilon);
    }

    SUBCASE("cga-expression-evaluation")
    {
        using namespace gal::cga;
        point<> p{3.5f, 0.2f, -23.9f};

        gal::engine engine{p};
        // Compute the contraction of a point as itself
        gal::scalar<> n = engine.compute([](auto p)
        {
            return p >> p;
        });
        CHECK_EQ(n, epsilon);
    }
}

TEST_SUITE_END();