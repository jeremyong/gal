#include <doctest/doctest.h>
#include <gal/pga.hpp>
#include <gal/format.hpp>
#include <gal/expression_debug.hpp>

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

        std::cout << "line: " <<  to_string(line) << std::endl;
    }

    SUBCASE("plane-construction")
    {
        point<> p1{1, 0, 0};
        point<> p2{0, 1, 0};
        point<> p3{0, 0, 1};
        plane<> p = compute([](auto pl1, auto pl2, auto pl3) { return pl1 & pl2 & pl3; }, p1, p2, p3);
        std::printf("plane: %f + %f*x + %f*y + %fz)\n", p.d, p.x, p.y, p.z);
    }
}

TEST_CASE("bivector-norm")
{
    auto l2 = evaluate<line<>>{}.debug([](auto l){ return ((l | l) + (l ^ l)); });
    CHECK_EQ(l2.size.term, 2);
}

TEST_CASE("motors")
{
    SUBCASE("simple-motor")
    {
        // Join two points to create a line
        point<> p1{0, 0, 0};
        point<> p2{0, 0, M_PI / 4};
        // This will be a vertical line in the +z direction through the origin (pi/2 rotation)
        line<> l = compute([](auto p1, auto p2) { return p1 & p2; }, p1, p2);
        motor<> m = exp(l);

        point<> p3{1, 0, 0};

        point<> p3_motor = compute([](auto p3, auto m) { return p3 % m; }, p3, m);
        CHECK_EQ(p3_motor.x, doctest::Approx(0.0));
        CHECK_EQ(p3_motor.y, doctest::Approx(1.0));
        CHECK_EQ(p3_motor.z, doctest::Approx(0.0));

        line<> l2 = log(m);
        CHECK_EQ(l2.dz, doctest::Approx(M_PI / 4));
    }

    SUBCASE("line-exp")
    {
        // Random line :)
        line<> l{3.234, -12.3, 4.2, 1.293, -3.58, -1.1};
        motor<> m = exp(l);
        line<> l2 = log(m);
        motor<> m2 = exp(l2);
        line<> l3 = log(m2);
        for (size_t i = 0; i != 6; ++i)
        {
            CHECK_EQ(l2[i], doctest::Approx(l3[i]));
        }
        for (size_t i = 0; i != 8; ++i)
        {
            CHECK_EQ(m[i], doctest::Approx(m2[i]));
        }

        auto m_norm = compute([](auto m) { return m * ~m; }, m2);
        printf("m_norm: %f, %f\n", m_norm[0], m_norm[1]);
        CHECK_EQ(m_norm[0], doctest::Approx(1));
        CHECK_EQ(m_norm[1], doctest::Approx(0));
    }

    SUBCASE("normalize-motor")
    {
        motor<> m{{293.2, -39.3, 59.3, -1.04, 434.3, 23.0, 72.874}};
        m.normalize();
        auto m_norm = compute([](auto m) { return m * ~m; }, m);
        printf("m_norm: %f, %f\n", m_norm[0], m_norm[1]);
        CHECK_EQ(m_norm[0], doctest::Approx(1));
        CHECK_EQ(m_norm[1], doctest::Approx(0));
    }
}

TEST_SUITE_END();