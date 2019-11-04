#include <doctest/doctest.h>
#include <gal/format.hpp>
#include <gal/pga.hpp>

#include <iostream>

using namespace gal;
using namespace gal::pga;

template <typename Line>
constexpr auto expp(Line l)
{
    auto s        = (l | l)[0];
    auto p        = (l ^ l)[0b1111];
    auto u        = sqrt(-s);
    auto v        = -p / (2 * u);
    auto inv_norm = 1 / u + v / s * 1_ps;
    auto cos_u    = cos(u);
    auto sin_u    = sin(u);
    auto real     = cos_u - v * sin_u * 1_ps;
    auto ideal    = sin_u + v * cos_u * 1_ps;
    return real + ideal * inv_norm * l;
}

template <typename Line>
constexpr auto normalized_line(Line l)
{
    auto norm = sqrt(l[0b110] * l[0b110] + l[0b1010] * l[0b1010] + l[0b1100] * l[0b1100]);
    return l / norm;
}

template <typename T, typename Line>
constexpr auto rotor(T ang, Line l)
{
    return cos(ang / 2) + sin(ang / 2) * normalized_line(l);
}

template <typename T, typename Line>
constexpr auto translator(T d, Line l)
{
    return 1 + d / 2 * l;
}

template <typename T, typename Line>
constexpr auto circle(T t, T r, Line l)
{
    return ::rotor(2 * PI * t, l) * ::translator(r, -1_e01);
}

template <typename T, typename Line>
constexpr auto torus(T s, T t, T r1, Line l1, T r2, Line l2)
{
    return circle(s, r2, l2) * circle(t, r1, l1);
}

using sc = scalar<pga_algebra, float>;

TEST_SUITE_BEGIN("projective-geometric-algebra");
TEST_CASE("point-translation")
{
    sc d{2};
    point<> p1{0, 0, 0};
    point<> p2 = compute(
        [](auto p, auto d) {
            auto t = ::translator(d, 1_e01);
            return t * p * ~t;
        },
        p1,
        d);
    CHECK_EQ(p2.x, doctest::Approx(2));
    CHECK_EQ(p2.y, doctest::Approx(0));
    CHECK_EQ(p2.z, doctest::Approx(0));
}

TEST_CASE("line-norm")
{
    sc a{M_PI / 2.0};
    auto rpn = evaluate<sc>::rpnf_reshaped([](auto a) {
        auto nl = normalized_line(1_e12 + 1_e13);
        return nl;
    });
    printf("line-norm: %s\n", to_string(rpn).c_str());
    auto ie = evaluate<sc>::ie_reshaped([](auto a) {
        auto nl = normalized_line(1_e12 + 1_e13);
        return nl;
    });
    auto l2 = compute(
        [](auto a) {
            auto nl = normalized_line(1_e12 + 1_e13);
            return nl;
        },
        a);
    CHECK_EQ(l2[0], doctest::Approx(1.0 / std::sqrt(2)));
    CHECK_EQ(l2[1], doctest::Approx(1.0 / std::sqrt(2)));
}

TEST_CASE("point-rotation")
{
    sc a{M_PI / 2.0};
    point<> p1{1, 0, 0};
    auto a1 = compute([](auto a) { return a / 2; }, a);
    printf("half-angle: %f\n", a1[0]);
    auto r = compute(
        [](auto a) {
            auto r = ::rotor(a, 1_e12);
            return r;
        },
        a);
    auto r1 = compute(
        [](auto a) {
            auto r = ::rotor(a, 1_e12);
            return ~r;
        },
        a);
    CHECK_EQ(r[0], r1[0]);
    CHECK_EQ(r[1], doctest::Approx(-r1[1]));
    point<> p2 = compute(
        [](auto p, auto a) {
            auto r = ::rotor(a, 1_e12);
            return r * p * ~r;
        },
        p1,
        a);
}

TEST_CASE("torus-construction")
{
    sc s{.2};
    sc t{.3};
    sc r1{0.25};
    sc r2{0.6};

    auto rpn = evaluate<sc, sc, sc, sc>::rpnf_reshaped([](auto s, auto t, auto r1, auto r2) {
        auto to = torus(s, t, r1, 1_e12, r2, 1_e13);
        return to * 1_e123 * ~to;
    });
    printf("torus exp: %s\n", to_string(rpn).c_str());

    point<> p = compute(
        [](auto s, auto t, auto r1, auto r2) {
            auto to = torus(s, t, r1, 1_e12, r2, 1_e13);
            return to * 1_e123 * ~to;
        },
        s,
        t,
        r1,
        r2);
}

TEST_CASE("incidence")
{
    SUBCASE("point-to-line-construction")
    {
        point<> p1{0, 0, 1};
        point<> p2{1, 0, 1};

        auto line = compute([](auto p1, auto p2) { return p1 & p2; }, p1, p2);

        // auto l = evaluate<point<>, point<>>{}.debug([](auto p1, auto p2) { return p1 & p2; });

        std::cout << "line: " << to_string(line) << std::endl;
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
    auto rpn = evaluate<line<>>::rpnf([](auto l) { return ((l | l) + (l ^ l)); });
    std::printf("bivector norm: %s\n", gal::to_string(rpn).c_str());
    auto l2 = evaluate<line<>>::ie([](auto l) { return ((l | l) + (l ^ l)); });
    CHECK_EQ(l2.size.term, 2);
}

/*
TEST_CASE("motors")
{
    SUBCASE("simple-motor")
    {
        // Join two points to create a line
        point<> p1{0, 0, 0};
        point<> p2{0, 0, M_PI / 4};
        // This will be a vertical line in the +z direction through the origin (pi/2 rotation)
        auto [l, m] = std::tuple<line<>, motor<>>(compute(
            [](auto p1, auto p2) {
                auto l = p1 & p2;
                return gal::make_tuple(l, exp(l));
            },
            p1,
            p2));

        point<> p3{1, 0, 0};

        point<> p3_motor = compute([](auto p3, auto m) { return p3 % m; }, p3, m);
        CHECK_EQ(p3_motor.x, doctest::Approx(0.0));
        CHECK_EQ(p3_motor.y, doctest::Approx(1.0));
        CHECK_EQ(p3_motor.z, doctest::Approx(0.0));

        // line<> l2 = compute([](auto m) { return log(m); }, m);
        // CHECK_EQ(l2.dz, doctest::Approx(M_PI / 4));
    }

    SUBCASE("line-exp")
    {
        // Random line :)
        line<> l{3.234, -12.3, 4.2, 1.293, -3.58, -1.1};

        auto rpn = evaluate<line<>>::rpnf([](auto l) {
            auto m  = exp(l);
            auto l2 = log(m);
            auto m2 = exp(l2);
            // auto l3 = log(m2);
            // return gal::make_tuple(m, l2, m2);
            // return gal::make_tuple(m, l2);
            return gal::make_tuple(l2);
        });
        printf("line-exp: %s\n", gal::to_string(rpn).c_str());
        // TODO: the test below works on clang but not gcc
        auto result = compute(
            [](auto l) {
                auto m  = exp(l);
                auto l2 = log(m);
                auto m2 = exp(l2);
                auto l3 = log(m2);
                // return gal::make_tuple(m, l2, m2);
                return m;
            },
            l);
        // motor<> m  = std::get<0>(result);
        // line<> l2  = std::get<1>(result);
        // motor<> m2 = std::get<2>(result);
        // // line<> l3  = std::get<3>(result);
        // line<> l3 = l2;

        // for (size_t i = 0; i != 6; ++i)
        // {
        //     CHECK_EQ(l2[i], doctest::Approx(l3[i]));
        // }
        // for (size_t i = 0; i != 8; ++i)
        // {
        //     CHECK_EQ(m[i], doctest::Approx(m2[i]));
        // }

        // auto m_norm = compute([](auto m) { return m * ~m; }, m2);
        // printf("m_norm: %f, %f\n", m_norm[0], m_norm[1]);
        // CHECK_EQ(m_norm[0], doctest::Approx(1));
        // CHECK_EQ(m_norm[1], doctest::Approx(0));
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
*/

TEST_SUITE_END();
