#include <doctest/doctest.h>

#include <formatters.hpp>
#include <pga2.hpp>
#include <engine.hpp>

using namespace gal::pga2;

TEST_SUITE_BEGIN("engine");

TEST_CASE("basic-computation")
{
    SUBCASE("expression-evaluation")
    {
        point<> p1{2.4f, 3.6f};
        point<> p2{-1.1f, 2.7f};

        gal::engine engine{p1, p2};
        auto l = engine.compute<line<>>([](auto p1, auto p2)
        {
            return p1 | p2;
        });
        CHECK_EQ(p1.x * l.a + p1.y * l.b + l.c, doctest::Approx(0.f));
        CHECK_EQ(p2.x * l.a + p2.y * l.b + l.c, doctest::Approx(0.f));
    }
}

TEST_SUITE_END();