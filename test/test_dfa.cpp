#include <doctest/doctest.h>

#include <gal/pga.hpp>

TEST_SUITE_BEGIN("dfa");

using namespace gal::pga;
using namespace gal;

TEST_CASE("shift")
{
    auto rpn = gal::evaluate<plane<>>::rpnf([](auto p) { return 1 + p; });
    printf("%s\n", ::gal::to_string(rpn).c_str());
    CHECK_EQ(rpn.count, 5);
    CHECK_EQ(rpn.nodes[3].o, detail::c_scalar);
    CHECK_EQ(rpn.nodes[2].q.num, 1);
    CHECK_EQ(rpn.nodes[2].q.den, 1);

    auto rpn2 = gal::evaluate<plane<>>::rpnf([](auto p) { return 1 + p + 2; });
    printf("%s\n", ::gal::to_string(rpn2).c_str());
    CHECK_EQ(rpn2.nodes[3].o, detail::c_scalar);
    CHECK_EQ(rpn2.nodes[2].q.num, 3);
    CHECK_EQ(rpn2.nodes[2].q.den, 1);
    auto rpn3 = gal::evaluate<plane<>>::rpnf([](auto p) { return 1 + 2 * p - 1; });
    printf("%s\n", ::gal::to_string(rpn3).c_str());
    CHECK_EQ(rpn3.count, 1);
    CHECK_EQ(rpn3.q.num, 2);
}

TEST_CASE("scale")
{
    auto rpn = gal::evaluate<plane<>>::rpnf([](auto p) { return 2 * p / 3; });
    printf("%s\n", ::gal::to_string(rpn).c_str());
    CHECK_EQ(rpn.q.num, 2);
    CHECK_EQ(rpn.q.den, 3);
}

TEST_CASE("sum-dependent")
{
    auto rpn = gal::evaluate<plane<>>::rpnf([](auto p) { return p / 3 + 2 * p / 3; });
    printf("%s\n", ::gal::to_string(rpn).c_str());
    CHECK_EQ(rpn.count, 1);
    CHECK_EQ(rpn.q.num, 1);
    CHECK_EQ(rpn.q.den, 1);
    CHECK_EQ(rpn.nodes[0].o, detail::op_id);
    CHECK_EQ(rpn.nodes[0].checksum, 0);
}

TEST_CASE("sum-independent")
{
    auto rpn = gal::evaluate<plane<>, plane<>>::rpnf([](auto p1, auto p2) { return p1 / 3 + p2; });
    printf("%s\n", ::gal::to_string(rpn).c_str());
    auto rpn2 = gal::evaluate<plane<>, plane<>, plane<>, plane<>>{}.rpnf(
        [](auto p1, auto p2, auto p3, auto p4) { return p1 + p2 + p3 + p4; });
    printf("%s\n", ::gal::to_string(rpn2).c_str());
    CHECK_EQ(rpn2.count, 9);
}

TEST_CASE("sum-cancellation")
{
    auto rpn = gal::evaluate<plane<>>::rpnf([](auto p) { return p - p; });
    printf("%s\n", ::gal::to_string(rpn).c_str());
}

TEST_CASE("unary-ops")
{
    auto rpn = gal::evaluate<plane<>>::rpnf([](auto p) { return !p; });
    printf("%s\n", ::gal::to_string(rpn).c_str());
    CHECK_EQ(rpn.count, 2);
    rpn = gal::evaluate<plane<>>::rpnf([](auto p) { return ~p; });
    printf("%s\n", ::gal::to_string(rpn).c_str());
    auto rpn2 = gal::evaluate<plane<>>::rpnf([](auto p) { return ~~p; });
    printf("%s\n", ::gal::to_string(rpn2).c_str());
    CHECK_EQ(rpn2.count, 1);
}

TEST_CASE("geometric-product")
{
    auto rpn  = gal::evaluate<plane<>, plane<>>::rpnf([](auto p1, auto p2) { return p1 * p2; });
    auto rpn1 = gal::evaluate<plane<>, plane<>>::rpnf([](auto p1, auto p2) { return p2 * p1; });
    printf("%s\n", ::gal::to_string(rpn).c_str());
    printf("%s\n", ::gal::to_string(rpn1).c_str());
    // // The geometric product does not commute in general
    CHECK_NE(rpn.nodes[1].checksum, rpn1.nodes[1].checksum);
}

TEST_CASE("exterior-product")
{
    auto rpn = gal::evaluate<plane<>>::rpnf([](auto p) { return p ^ p; });
    printf("%s\n", ::gal::to_string(rpn).c_str());
    auto rpn1 = gal::evaluate<plane<>, line<>>::rpnf([](auto p, auto l) { return p ^ l; });
    printf("%s\n", ::gal::to_string(rpn1).c_str());
}

TEST_CASE("regressive-product")
{
    auto rpn = gal::evaluate<gal::pga::point<>, gal::pga::line<>>::rpnf(
        [](auto p, auto l) { return p & l; });
    printf("%s\n", ::gal::to_string(rpn).c_str());
}

TEST_CASE("left-contraction")
{
    auto rpn = gal::evaluate<gal::pga::point<>, gal::pga::line<>>::rpnf(
        [](auto p, auto l) { return p >> l; });
    printf("%s\n", ::gal::to_string(rpn).c_str());
    auto rpn1 = gal::evaluate<gal::pga::point<>, gal::pga::line<>>::rpnf(
        [](auto p, auto l) { return l >> p; });
    printf("%s\n", ::gal::to_string(rpn1).c_str());
    // // The left contraction doesn't commute
    CHECK_NE(rpn.back().checksum, rpn1.back().checksum);
}

TEST_CASE("various-cancellations")
{
    auto rpn = gal::evaluate<gal::pga::point<>, gal::pga::plane<>>::rpnf([](auto point, auto plane) {
        return (-1 + 2 * point + 1)
               ^ (2 * plane + 3 * plane * point + point + 2 * plane - plane
                  - 3 * (plane * point + plane));
    });
    printf("%s\n", ::gal::to_string(rpn).c_str());
}

TEST_CASE("multivector-construction")
{
    auto rpn = gal::evaluate<gal::pga::point<>>::rpnf([](auto p) { return 1 + 2_e1; });
    printf("%s\n", ::gal::to_string(rpn).c_str());

    auto mv = gal::evaluate<gal::pga::point<>>::ie([](auto p) { return 1 + 2_e1; });
    CHECK_EQ(mv.size.ind, 0);
    CHECK_EQ(mv.size.mon, 2);
    CHECK_EQ(mv.size.term, 2);
    CHECK_EQ(mv.terms[0].element, 0);
    CHECK_EQ(mv.terms[1].element, 2);
    CHECK_EQ(mv.mons[0].q.num, 1);
    CHECK_EQ(mv.mons[1].q.num, 2);

    auto mv2 = gal::evaluate<gal::pga::point<>>::ie([](auto p) { return 2_e1; });
    // The scaling factor here is left on the outer RPN
    CHECK_EQ(mv2.mons[0].q.num, 1);

    auto mv3 = gal::evaluate<gal::pga::point<>>::ie([](auto p) { return 2_e1 + 3_e0; });
    CHECK_EQ(mv3.mons[0].q.num, 3);
    CHECK_EQ(mv3.mons[1].q.num, 2);
}

TEST_CASE("compound-sum")
{
    auto mv
        = gal::evaluate<gal::pga::point<>, gal::pga::point<>, gal::pga::point<>, gal::pga::point<>>::ie(
            [](auto p1, auto p2, auto p3, auto p4) { return p1 + p2 * p3 * p4; });
}

TEST_CASE("compute-simple-expression")
{
    gal::pga::plane<> p{1, 2, 3, 4};
    gal::pga::plane<> p1 = gal::compute([](auto p) { return 2 * p; }, p);

    CHECK_EQ(2 * p.d, p1.d);
    CHECK_EQ(2 * p.x, p1.x);
    CHECK_EQ(2 * p.y, p1.y);
    CHECK_EQ(2 * p.z, p1.z);
}

/*
TEST_CASE("compute-with-temporary")
{
    gal::pga::line<> l{1, 2, 3, 4, 5, 6};
    gal::pga::motor<> m = gal::compute([](auto l) { return gal::exp(l); }, l);

    auto l2 = compute([](auto l) { return ((l | l) + (l ^ l)); }, l);
    auto s  = l2.template select<0>();
    auto p  = l2.template select<0b1111>();
    auto u  = std::sqrt(-s);
    auto v  = -p / (2 * u);

    gal::entity<gal::pga::pga_algebra, float, 0, 0b1111> inv_norm{1.0f / u, v / s};
    auto cos_u = std::cos(u);
    auto sin_u = std::sin(u);
    using in_t = decltype(inv_norm);
    decltype(inv_norm) real{cos_u, -v * sin_u};
    decltype(inv_norm) ideal{sin_u, v * cos_u};

    gal::pga::motor<> m1 = compute(
        [](auto real, auto ideal, auto inv_norm, auto l) { return real + ideal * inv_norm * l; },
        real,
        ideal,
        inv_norm,
        l);

    CHECK_EQ(m[0], m1[0]);
    CHECK_EQ(m[1], m1[1]);
    CHECK_EQ(m[2], m1[2]);
    CHECK_EQ(m[3], m1[3]);
    CHECK_EQ(m[4], m1[4]);
    CHECK_EQ(m[5], m1[5]);
    CHECK_EQ(m[6], m1[6]);
    CHECK_EQ(m[7], m1[7]);
}
*/

TEST_CASE("rpne-literals")
{
    gal::pga::motor<> ps = compute([](auto l) { return 1_ps; }, gal::pga::line<>{1, 2, 3, 4, 5, 6});
    for (size_t i = 0; i != 6; ++i)
    {
        CHECK_EQ(ps[i], 0);
    }
    CHECK_EQ(ps[7], 1);
}

TEST_CASE("extract-elements")
{
    auto rpn = evaluate<gal::pga::line<>>::rpnf([](auto l) { return l[0b1100]; });
    printf("%s\n", ::gal::to_string(rpn).c_str());
    auto ie = evaluate<gal::pga::line<>>::ie([](auto l) { return l[0b1100]; });

    gal::pga::motor<> m{1, 2, 3, 4, 5, 6, 7, 8};
    auto result = compute([](auto m) { return m[0] + m[0b1111] * 1_ps; }, m);
    CHECK_EQ(decltype(result)::size(), 2);
    CHECK_EQ(result[0], 1);
    CHECK_EQ(result[1], 8);
}

TEST_CASE("scalar-transcendentals")
{
    gal::pga::motor<> m{1, 2, -M_PI / 2.0f, 2.0f * M_PI, 5, 6, 7, 8};
    auto result = compute([](auto m) { return sqrt(m[0b11]); }, m);
    CHECK_EQ(result[0], doctest::Approx(std::sqrt(2.0f)));

    auto rpn = evaluate<motor<>>::rpnf_reshaped(
        [](auto m) { return sin(m[0b101]) + cos(m[0b110]) * 1_ps; });
    printf("scalar transcendentals: %s\n", gal::to_string(rpn).c_str());

    auto ie
        = evaluate<motor<>>::ie_reshaped([](auto m) { return sin(m[0b101]) + cos(m[0b110]) * 1_ps; });

    auto s = compute([](auto m) { return sin(m[0b101]) + cos(m[0b110]) * 1_ps; }, m);
    CHECK_EQ(s[0], doctest::Approx(-1.0f));
    CHECK_EQ(s[1], doctest::Approx(1.0f));
}

TEST_CASE("exterior-product")
{
    SUBCASE("vector")
    {
        auto rpn = evaluate<gal::pga::line<>>::rpnf([](auto dummy) {
            auto l = 1_e0 + 1_e1;
            return l ^ l;
        });
        printf("%s\n", gal::to_string(rpn).c_str());
        auto ie = evaluate<gal::pga::line<>>::ie([](auto dummy) {
            auto l = 1_e0 + 1_e1;
            return l ^ l;
        });
        CHECK_EQ(ie.size.term, 0);
    }

    SUBCASE("bivector")
    {
        auto rpn = evaluate<gal::pga::line<>>::rpnf_reshaped([](auto dummy) {
            auto l = 1_e01 + 1_e23;
            return l ^ l;
        });
        printf("%s\n", gal::to_string(rpn).c_str());
        auto ie = evaluate<gal::pga::line<>>::ie([](auto dummy) {
            auto l = 1_e01 + 1_e23;
            return l ^ l;
        });
        CHECK_EQ(ie.terms[0].count, 1);
        CHECK_EQ(ie.mons[0].q.num, 2);
    }
}

/*
TEST_CASE("exp-and-reverse")
{
    gal::pga::line<double> l{1, 2, 3, 4, 5, 6};
    // After exponentiation, logarithm, and exponentiation again, we expect the same motor
    gal::pga::motor<double> m  = compute([](auto l) { return gal::exp(l); }, l);
    gal::pga::line<double> l1  = compute([](auto m) { return gal::log(m); }, m);
    gal::pga::motor<double> m1 = compute([](auto l) { return gal::exp(l); }, l1);
    // The two normalized lines should also be equal
    gal::pga::line<double> l2  = compute([](auto m) { return gal::log(m); }, m1);
    gal::pga::motor<double> m2 = compute([](auto l) { return gal::exp(l); }, l2);

    CHECK_EQ(m[0], doctest::Approx(m1[0]));
    CHECK_EQ(m[1], doctest::Approx(m1[1]));
    CHECK_EQ(m[2], doctest::Approx(m1[2]));
    CHECK_EQ(m[3], doctest::Approx(m1[3]));
    CHECK_EQ(m[4], doctest::Approx(m1[4]));
    CHECK_EQ(m[5], doctest::Approx(m1[5]));
    CHECK_EQ(m[6], doctest::Approx(m1[6]));
    CHECK_EQ(m[7], doctest::Approx(m1[7]));

    CHECK_EQ(l1[0], doctest::Approx(l2[0]));
    CHECK_EQ(l1[1], doctest::Approx(l2[1]));
    CHECK_EQ(l1[2], doctest::Approx(l2[2]));
    CHECK_EQ(l1[3], doctest::Approx(l2[3]));
    CHECK_EQ(l1[4], doctest::Approx(l2[4]));
    CHECK_EQ(l1[5], doctest::Approx(l2[5]));
}
*/

TEST_CASE("variadic-return")
{
    gal::pga::plane<> p{1, 2, 3, 4};
    auto&& [p1, p2] = compute([](auto p) { return gal::tuple{p, 2 * p}; }, p);
    CHECK_EQ(p1[0], 1);
    CHECK_EQ(p1[1], 2);
    CHECK_EQ(p1[2], 3);
    CHECK_EQ(p1[3], 4);
    CHECK_EQ(p2[0], 2);
    CHECK_EQ(p2[1], 4);
    CHECK_EQ(p2[2], 6);
    CHECK_EQ(p2[3], 8);
}

TEST_CASE("scalar-quantities")
{
    gal::pga::plane<> p{1, 2, 3, 4};
    float d  = 2.3f;
    auto rpn = evaluate<gal::pga::plane<>, float>::rpnf([](auto p, auto d) { return p * d; });
    std::printf("scaled: %s\n", gal::to_string(rpn).c_str());

    gal::pga::plane<> p2 = compute([](auto p, auto d) { return p * d; }, p, d);
    CHECK_EQ(p2.d, d);
}

TEST_CASE("cse")
{
    auto rpn = evaluate<gal::pga::plane<>, gal::pga::plane<>>::rpnf_reshaped(
        [](auto p1, auto p2) { return (p1 + p2) * (p1 + p2); });
    std::printf("cse: %s\n", gal::to_string(rpn).c_str());
}

TEST_CASE("scalar-division")
{
    auto rpn = evaluate<gal::pga::plane<>>::rpnf_reshaped([](auto p) {
        auto l = 1_e12;
        return p / l[0b110];
    });
    std::printf("scalar-division: %s\n", gal::to_string(rpn).c_str());
    auto ie = evaluate<gal::pga::plane<>>::ie_reshaped([](auto p) {
        auto l = 1_e12 / 2;
        return p / l[0b110];
    });

    gal::pga::plane<> p{1, 2, 3, 4};
    gal::pga::plane<> scale_2 = compute(
        [](auto p) {
            auto l = 1_e12 / 2;
            return p / l[0b110];
        },
        p);
    CHECK_EQ(scale_2.d, doctest::Approx(2.0f));

    using sc = gal::scalar<gal::pga::pga_algebra, float>;
    sc s{2.0f};
    auto rpn1 = evaluate<sc>::rpnf_reshaped([](auto s) { return s / 2; });
    std::printf("scalar-division2: %s\n", gal::to_string(rpn1).c_str());

    auto ie1 = evaluate<sc>::ie_reshaped([](auto s) { return s / 2; });

    auto half_s = compute([](auto s) { return s / 2; }, s);
    CHECK_EQ(half_s[0], doctest::Approx(1.0f));

    SUBCASE("in-transcendental")
    {
        sc p{M_PI};
        auto c_rpn = evaluate<sc>::rpnf_reshaped([](auto p) { return sin(p / 2); });
        std::printf("sin pi/2: %s\n", gal::to_string(c_rpn).c_str());
        auto c_ie = evaluate<sc>::ie([](auto p) { return sin(p / 2); });
        auto c    = compute([](auto p) { return sin(p / 2); }, p);
        CHECK_EQ(c[0], doctest::Approx(1.0f));
    }
}

TEST_CASE("special-constant-evaluation")
{
    using sc = gal::scalar<gal::pga::pga_algebra, float>;
    sc s{2.0f};
    auto rpn = evaluate<sc>::rpnf_reshaped([](auto s) { return s * gal::PI; });
    std::printf("2 * pi: %s\n", gal::to_string(rpn).c_str());

    auto ie = evaluate<sc>::ie_reshaped([](auto s) { return s * gal::PI; });

    auto two_pi = compute([](auto s) { return gal::PI * s; }, s);
    CHECK_EQ(two_pi[0], doctest::Approx(2.0 * M_PI));

    auto sin_34_pi = compute([](auto s) { return sin(3 * gal::PI * s / 4); }, s);
    CHECK_EQ(sin_34_pi[0], doctest::Approx(-1));

    auto foo = compute([](auto s) { return cos(3 * gal::PI * s / 4) + sin(3 * gal::PI * s / 4); }, s);
    CHECK_EQ(foo[0], doctest::Approx(-1));
}

TEST_SUITE_END();
