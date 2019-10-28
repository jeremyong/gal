#include "test_util.hpp"

#include <doctest/doctest.h>
#include <fmt/core.h>
#include <gal/engine.hpp>
#include <gal/cga.hpp>

using namespace gal;
using algebra_t = gal::cga::cga_algebra;

TEST_SUITE_BEGIN("finite-algebra");

TEST_CASE("multivector-arithmetic")
{
    using S = gal::scalar<algebra_t, float>;

    SUBCASE("scalar-sum")
    {
        S s1{1};
        S s2{2};

        auto sum = compute([](auto const& s1, auto const& s2) { return s1 + s2; }, s1, s2);
        CHECK_EQ(sum[0], doctest::Approx(3.0));
    }

    SUBCASE("scalar-negation")
    {
        S s{2};
        auto negation = compute([](auto const& s) { return -s; }, s);
        CHECK_EQ(negation[0], doctest::Approx(-2));
    }

    SUBCASE("scalar-sum-cancellation")
    {
        S s{1};

        auto sum = compute([](auto const& s) { return s + -s; }, s);
        CHECK_EQ(sum.size(), 0);
    }

    SUBCASE("scalar-subtraction")
    {
        S s1{1};
        S s2{2};

        auto diff = compute([](auto const& s1, auto const& s2) { return s1 - s2; }, s1, s2);
        // auto diff = compute_rt<([](auto const& s1, auto const& s2) { return s1 - s2; })>(s1, s2);
        CHECK_EQ(diff[0], doctest::Approx(-1.0));

        auto diff2 = compute([](auto const& s1, auto const& s2) { return s2 - s1 + s1; }, s1, s2);
        CHECK_EQ(diff2[0], doctest::Approx(2.0));
    }

    SUBCASE("independent-term-sum")
    {
        using e1_t = entity<algebra_t, float, 1>;
        e1_t e1{1.0};
        S s1{1.0};
        auto sum = compute([](auto s1, auto e1) { return s1 + e1; }, s1, e1);
        CHECK_EQ(sum[0], doctest::Approx(1.0));
        CHECK_EQ(sum[1], doctest::Approx(1.0));
    }

    SUBCASE("coincident-term-sum")
    {
        using e1_t = entity<algebra_t, float, 0b0, 0b1>;
        auto sum   = evaluate<e1_t, e1_t>{}.debug([](auto v1, auto v2) { return v1 + v2; });
        CHECK_EQ(sum.inds[0].id, 0);
        CHECK_EQ(sum.inds[1].id, 2);
        CHECK_EQ(sum.inds[2].id, 1);
        CHECK_EQ(sum.inds[3].id, 3);
        CHECK_EQ(sum.mons[0].count, 1);
        CHECK_EQ(sum.mons[1].count, 1);
        CHECK_EQ(sum.mons[2].count, 1);
        CHECK_EQ(sum.mons[3].count, 1);
        CHECK_EQ(sum.terms[0].element, 0);
        CHECK_EQ(sum.terms[0].count, 2);
        CHECK_EQ(sum.terms[0].mon_offset, 0);
        CHECK_EQ(sum.terms[1].element, 1);
        CHECK_EQ(sum.terms[1].count, 2);
        CHECK_EQ(sum.terms[1].mon_offset, 2);
    }

    SUBCASE("coincident-term-sum-with-scalar")
    {
        using e1_t = entity<algebra_t, float, 0b0, 0b1>;
        auto sum   = evaluate<e1_t, e1_t>{}.debug([](auto v1, auto v2) { return frac<1> + v1 + v2; });
        CHECK_EQ(sum.inds[0].id, 0);
        CHECK_EQ(sum.inds[1].id, 2);
        CHECK_EQ(sum.inds[2].id, 1);
        CHECK_EQ(sum.inds[3].id, 3);
        CHECK_EQ(sum.mons[0].q.num, 1);
        CHECK_EQ(sum.mons[0].q.den, 1);
        CHECK_EQ(sum.mons[0].count, 0);
        CHECK_EQ(sum.mons[1].count, 1);
        CHECK_EQ(sum.mons[1].ind_offset, 0);
        CHECK_EQ(sum.mons[2].count, 1);
        CHECK_EQ(sum.mons[3].count, 1);
        CHECK_EQ(sum.mons[4].count, 1);
        CHECK_EQ(sum.terms[0].element, 0);
        CHECK_EQ(sum.terms[0].count, 3);
        CHECK_EQ(sum.terms[0].mon_offset, 0);
        CHECK_EQ(sum.terms[1].element, 1);
        CHECK_EQ(sum.terms[1].count, 2);
        CHECK_EQ(sum.terms[1].mon_offset, 3);
    }

    SUBCASE("noncoincident-term-sum-with-scalar")
    {
        using e1_t = entity<algebra_t, float, 0b1, 0b10>;
        auto sum   = evaluate<e1_t>{}([](auto v) { return frac<2> + v; });
        CHECK_EQ(sum.inds[0].id, 0);
        CHECK_EQ(sum.inds[1].id, 1);
        CHECK_EQ(sum.mons[0].q.num, 2);
        CHECK_EQ(sum.mons[0].q.den, 1);
        CHECK_EQ(sum.mons[0].count, 0);
        CHECK_EQ(sum.mons[1].count, 1);
        CHECK_EQ(sum.mons[1].ind_offset, 0);
        CHECK_EQ(sum.mons[2].count, 1);
        CHECK_EQ(sum.mons[2].ind_offset, 1);
        CHECK_EQ(sum.terms[0].element, 0);
        CHECK_EQ(sum.terms[0].count, 1);
        CHECK_EQ(sum.terms[0].mon_offset, 0);
        CHECK_EQ(sum.terms[1].element, 1);
        CHECK_EQ(sum.terms[1].count, 1);
        CHECK_EQ(sum.terms[1].mon_offset, 1);
        CHECK_EQ(sum.terms[2].element, 2);
        CHECK_EQ(sum.terms[2].count, 1);
        CHECK_EQ(sum.terms[2].mon_offset, 2);
    }
}

struct sm
{
    constexpr static uint8_t dimension = 1;
};

// Simple algebra that operates only on scalar quantities
struct sa
{
    using metric_t = sm;

    [[nodiscard]] constexpr static std::pair<uint8_t, int> product(uint8_t g1, uint8_t g2)
    {
        return {0, 1};
    }
};

using gal::detail::product;

TEST_CASE("polynomial-expansion")
{
    SUBCASE("binomial-expansion")
    {
        mv<sa, 2, 2, 1> b1{
            mv_size{2, 2, 1},
            {ind{0, 1}, ind{1, 2}},
            {mon{one, one, 1, 0}, mon{one, rat{2}, 1, 1}},
            {term{2, 0, 0}}
        };
        mv<sa, 2, 2, 1> b2{
            mv_size{2, 2, 1},
            {ind{2, 1}, ind{3, 1}},
            {mon{one, one, 1, 0}, mon{one, one, 1, 1}},
            {term{2, 0, 0}}
        };
        auto b12 = product(sa{}, b1, b2);
        CHECK_EQ(b12.size.term, 1);
        CHECK_EQ(b12.size.mon, 4);
        CHECK_EQ(b12.size.ind, 8);
        CHECK_EQ(b12.inds[0].id, 0);
        CHECK_EQ(b12.inds[1].id, 2);
        CHECK_EQ(b12.inds[2].id, 0);
        CHECK_EQ(b12.inds[3].id, 3);
        CHECK_EQ(b12.inds[4].id, 1);
        CHECK_EQ(b12.inds[4].degree.num, 2);
        CHECK_EQ(b12.inds[5].id, 2);
        CHECK_EQ(b12.inds[6].id, 1);
        CHECK_EQ(b12.inds[6].degree.num, 2);
        CHECK_EQ(b12.inds[7].id, 3);
        CHECK_EQ(b12.mons[1].ind_offset, 2);
        CHECK_EQ(b12.mons[1].degree.num, 2);
        CHECK_EQ(b12.mons[2].count, 2);
        CHECK_EQ(b12.mons[2].degree.num, 3);
        CHECK_EQ(b12.mons[3].count, 2);
        CHECK_EQ(b12.mons[3].degree.num, 3);
        CHECK_EQ(b12.terms[0].count, 4);
        CHECK_EQ(b12.terms[0].mon_offset, 0);
    }

    SUBCASE("difference-of-squares")
    {
        mv<sa, 2, 2, 1> b1{
            mv_size{2, 2, 1},
            {ind{0, 1}, ind{1, 1}},
            {mon{one, one, 1, 0}, mon{one, one, 1, 1}},
            {term{2, 0, 0}}
        };
        mv<sa, 2, 2, 1> b2{
            mv_size{2, 2, 1},
            {ind{0, 1}, ind{1, 1}},
            {mon{one, one, 1, 0}, mon{minus_one, one, 1, 1}},
            {term{2, 0, 0}}
        };
        auto b12 = product(sa{}, b1, b2);
        CHECK_EQ(b12.size.term, 1);
        CHECK_EQ(b12.size.mon, 2);
        CHECK_EQ(b12.size.ind, 2);
        CHECK_EQ(b12.inds[0].id, 0);
        CHECK_EQ(b12.inds[0].degree.num, 2);
        CHECK_EQ(b12.inds[1].id, 1);
        CHECK_EQ(b12.inds[1].degree.num, 2);
        CHECK_EQ(b12.mons[0].count, 1);
        CHECK_EQ(b12.mons[0].ind_offset, 0);
        CHECK_EQ(b12.mons[0].degree.num, 2);
        CHECK_EQ(b12.mons[1].count, 1);
        CHECK_EQ(b12.mons[1].ind_offset, 1);
        CHECK_EQ(b12.mons[1].degree.num, 2);
    }

    SUBCASE("monomial-term-cancellation")
    {
        mv<sa, 1, 1, 1> m1{
            mv_size{1, 1, 1},
            {ind{0, 1}},
            {mon{one, one, 1, 0}},
            {term{1, 0, 0}}
        };
        mv<sa, 1, 1, 1> m2{
            mv_size{1, 1, 1},
            {ind{0, -1}},
            {mon{one, minus_one, 1, 0}},
            {term{1, 0, 0}}
        };
        auto m12 = product(sa{}, m1, m2);
        CHECK_EQ(m12.size.term, 1);
        CHECK_EQ(m12.size.mon, 1);
        CHECK_EQ(m12.size.ind, 0);
        CHECK_EQ(m12.mons[0].count, 0);
        CHECK_EQ(m12.mons[0].degree.num, 0);
        CHECK_EQ(m12.mons[0].q.num, 1);
        CHECK_EQ(m12.mons[0].q.den, 1);
        CHECK_EQ(m12.terms[0].count, 1);
        CHECK_EQ(m12.terms[0].mon_offset, 0);
    }

    SUBCASE("monomial-term-multiplication")
    {
        mv<sa, 1, 1, 1> m1{
            mv_size{1, 1, 1},
            {ind{1, 1}},
            {mon{one, one, 1, 0}},
            {term{1, 0, 0}}
        };
        mv<sa, 1, 1, 1> m2{
            mv_size{1, 1, 1},
            {ind{1, 2}},
            {mon{one, rat{2}, 1, 0}},
            {term{1, 0, 0}}
        };
        auto m12 = product(sa{}, m1, m2);
        CHECK_EQ(m12.size.term, 1);
        CHECK_EQ(m12.size.mon, 1);
        CHECK_EQ(m12.size.ind, 1);
        CHECK_EQ(m12.inds[0].id, 1);
        CHECK_EQ(m12.inds[0].degree.num, 3);
        CHECK_EQ(m12.mons[0].count, 1);
        CHECK_EQ(m12.mons[0].degree.num, 3);
        CHECK_EQ(m12.mons[0].q.num, 1);
        CHECK_EQ(m12.mons[0].q.den, 1);
        CHECK_EQ(m12.terms[0].count, 1);
        CHECK_EQ(m12.terms[0].mon_offset, 0);
    }
}

TEST_SUITE_END();