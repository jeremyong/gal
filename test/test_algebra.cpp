#include <doctest/doctest.h>

#include <finite_algebra.hpp>
#include <formatters.hpp>

using namespace gal;

TEST_SUITE_BEGIN("finite-algebra");

// The majority of these tests are actually compile-time tests (compilation and correctness are equivalent)

TEST_CASE("term-arithmetic")
{
    SUBCASE("different-degree-and-factor")
    {
        term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>> t1;
        term<element<1>, monomial<multiplier<2>, factor<degree<2>, 2>>> t2;
        constexpr auto t12 = t1 + t2;
        static_assert(decltype(t12)::size == 2);
        static_assert(decltype(t12)::first_t::multiplier == 1);
        static_assert(decltype(t12)::subsequent_t::first_t::multiplier == 2);
    }

    SUBCASE("term-cancellation")
    {
        term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>> t1;
        term<element<1>, monomial<multiplier<-1>, factor<degree<1>, 1>>> t2;
        constexpr auto t12 = t1 + t2;
        static_assert(decltype(t12)::size == 0);
    }

    {
        term<element<1>, monomial<multiplier<-2>, factor<degree<1>, 4>>, monomial<multiplier<1>, factor<degree<1>, 1>, factor<degree<2>, 3>>> t1;
        term<element<1>, monomial<multiplier<2>, factor<degree<1>, 4>>, monomial<multiplier<-1>, factor<degree<1>, 1>, factor<degree<2>, 3>>> t2;
        constexpr auto t12 = t1 + t2;
        static_assert(decltype(t12)::size == 0);
    }

    SUBCASE("successive-addition")
    {
        term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>> t1;
        constexpr auto t111 = t1 + t1 + t1;
        static_assert(decltype(t111)::size == 1);
        static_assert(decltype(t111)::first_t::multiplier == 3);
    }

    SUBCASE("different-monomial-addition")
    {
        term<element<1>, monomial<multiplier<-2>, factor<degree<1>, 4>>, monomial<multiplier<1>, factor<degree<1>, 1>, factor<degree<2>, 3>>> t1;
        term<element<1>, monomial<multiplier<2>, factor<degree<1>, 4>>, monomial<multiplier<1>, factor<degree<1>, 1>, factor<degree<2>, 3>>> t2;
        term<element<1>, monomial<multiplier<2>, factor<degree<2>, 4>>, monomial<multiplier<1>, factor<degree<1>, 1>, factor<degree<2>, 3>>> t3;
        constexpr auto t12 = t1 + t2;
        static_assert(decltype(t12)::size == 1);
        static_assert(
            std::is_same<monomial<multiplier<2>, factor<degree<1>, 1>, factor<degree<2>, 3>>, typename decltype(t12)::first_t>::value);

        constexpr auto t13 = t1 + t3;
        static_assert(t13.size == 3);
        static_assert(std::is_same<monomial<multiplier<-2>, factor<degree<1>, 4>>, typename decltype(t13)::first_t>::value);
        static_assert(std::is_same<monomial<multiplier<2>, factor<degree<2>, 4>>, typename decltype(t13)::subsequent_t::first_t>::value);
        static_assert(std::is_same<monomial<multiplier<2>, factor<degree<1>, 1>, factor<degree<2>, 3>>,
                                   typename decltype(t13)::subsequent_t::subsequent_t::first_t>::value);
    }
}

TEST_CASE("multivector-addition")
{
    SUBCASE("single-term-addition")
    {
        multivector<void, term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>>> v1;
        multivector<void, term<element<1>, monomial<multiplier<2>, factor<degree<2>, 2>>>> v2;
        auto v12 = v1 + v2;
        static_assert(v12.size == 1);
        static_assert(
            std::is_same<term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>, monomial<multiplier<2>, factor<degree<2>, 2>>>,
                         decltype(v12)::first_t>::value);

        multivector<void, term<element<1>, monomial<multiplier<2>, factor<degree<1>, 1>>>> v3;
        auto v13 = v1 + v3;
        static_assert(v13.size == 1);
        static_assert(std::is_same<term<element<1>, monomial<multiplier<3>, factor<degree<1>, 1>>>, decltype(v13)::first_t>::value);
    }

    SUBCASE("multiple-term-addition")
    {
        multivector<void, term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>>> v1;
        multivector<void, term<element<2>, monomial<multiplier<1>, factor<degree<1>, 1>>>> v2;
        auto v12 = v1 + v2;
        static_assert(v12.size == 2);
        static_assert(std::is_same<term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>>, decltype(v12)::first_t>::value);
        static_assert(
            std::is_same<term<element<2>, monomial<multiplier<1>, factor<degree<1>, 1>>>, decltype(v12)::subsequent_t::first_t>::value);

        multivector<void, term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>>, term<element<2>, monomial<multiplier<1>, factor<degree<1>, 1>>>> v3;
        auto v13 = v1 + v3;
        static_assert(v13.size == 2);
        static_assert(std::is_same<term<element<1>, monomial<multiplier<2>, factor<degree<1>, 1>>>, decltype(v13)::first_t>::value);
        static_assert(
            std::is_same<term<element<2>, monomial<multiplier<1>, factor<degree<1>, 1>>>, decltype(v13)::subsequent_t::first_t>::value);

        multivector<void, term<element<2>, monomial<multiplier<1>, factor<degree<1>, 1>>>, term<element<3>, monomial<multiplier<1>, factor<degree<1>, 1>>>> v4;
        auto v14 = v1 + v4;
        static_assert(v14.size == 3);
        static_assert(std::is_same<term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>>, decltype(v14)::first_t>::value);
        static_assert(
            std::is_same<term<element<2>, monomial<multiplier<1>, factor<degree<1>, 1>>>, decltype(v14)::subsequent_t::first_t>::value);
        static_assert(
            std::is_same<term<element<3>, monomial<multiplier<1>, factor<degree<1>, 1>>>, decltype(v14)::subsequent_t::subsequent_t::first_t>::value);
    }

    SUBCASE("term-cancellation")
    {
        multivector<void, term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>>> v1;
        multivector<void, term<element<1>, monomial<multiplier<-1>, factor<degree<1>, 1>>>> v2;
        auto v12 = v1 + v2;
        static_assert(v12.size == 0);

        multivector<void, term<element<1>, monomial<multiplier<1>, factor<degree<1>, 2>>>> v3;
        multivector<void, term<element<1>, monomial<multiplier<-1>, factor<degree<1>, 1>>, monomial<multiplier<-1>, factor<degree<1>, 2>>>> v4;
        auto v134 = v1 + v3 + v4;
        static_assert(v134.size == 0);
    }
}

TEST_CASE("term-multiplication")
{
    SUBCASE("factor-product")
    {
        factor<degree<1>, 1> f1;
        factor<degree<1>, 2> f2;
        static_assert(std::is_same<monomial<multiplier<1>, factor<degree<1>, 1>, factor<degree<1>, 2>>, decltype(f1 * f2)>::value);

        factor<degree<2>, 1> f3;
        static_assert(std::is_same<monomial<multiplier<1>, factor<degree<3>, 1>>, decltype(f1 * f3)>::value);
    }

    SUBCASE("monomial-product")
    {
        monomial<multiplier<1>, factor<degree<1>, 1>> m1;
        monomial<multiplier<2>, factor<degree<1>, 2>> m2;
        auto m12 = m1 * m2;
        static_assert(std::is_same<monomial<multiplier<2>, factor<degree<1>, 1>, factor<degree<1>, 2>>, decltype(m12)>::value);
        // monomial multiplication commutes (order should be the same)
        static_assert(std::is_same<monomial<multiplier<2>, factor<degree<1>, 1>, factor<degree<1>, 2>>, decltype(m2 * m1)>::value);
        CHECK(1 == 1);

        monomial<multiplier<3>, factor<degree<2>, 1>> m3;
        auto m13 = m1 * m3;
        static_assert(std::is_same<monomial<multiplier<3>, factor<degree<3>, 1>>, decltype(m13)>::value);

        monomial<multiplier<1>, factor<degree<2>, 1>, factor<degree<1>, 3>, factor<degree<3>, 5>> m4;
        auto m124 = m1 * m2 * m4;
        static_assert(std::is_same<monomial<multiplier<2>, factor<degree<3>, 1>, factor<degree<1>, 2>, factor<degree<1>, 3>, factor<degree<3>, 5>>, decltype(m124)>::value);
        static_assert(std::is_same<decltype(m1 * m4 * m2), decltype(m124)>::value);
        static_assert(std::is_same<decltype(m2 * m1 * m4), decltype(m124)>::value);
        static_assert(std::is_same<decltype(m4 * m1 * m2), decltype(m124)>::value);
        static_assert(std::is_same<decltype(m4 * m2 * m1), decltype(m124)>::value);

        monomial<multiplier<2>> m5;
        auto m15 = m1 * m5;
        static_assert(std::is_same<monomial<multiplier<2>, factor<degree<1>, 1>>, decltype(m15)>::value);
        monomial<multiplier<1>> m6;
        auto m16 = m1 * m6;
        static_assert(std::is_same<monomial<multiplier<1>, factor<degree<1>, 1>>, decltype(m16)>::value);
    }

    SUBCASE("term-product")
    {
        // Single term product
        term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>> t1;
        term<element<1>, monomial<multiplier<2>, factor<degree<2>, 2>>> t2;
        auto t12 = t1 * t2;
        static_assert(std::is_same<term<element<1>, monomial<multiplier<2>, factor<degree<1>, 1>, factor<degree<2>, 2>>>, decltype(t12)>::value);

        // One term into two
        term<element<1>, monomial<multiplier<-2>, factor<degree<1>, 4>>, monomial<multiplier<1>, factor<degree<1>, 1>, factor<degree<2>, 3>>> t3;
        auto t13 = t1 * t3;
        static_assert(std::is_same<term<element<1>, monomial<multiplier<-2>, factor<degree<1>, 1>, factor<degree<1>, 4>>,
                                        monomial<multiplier<1>, factor<degree<2>, 1>, factor<degree<2>, 3>>>,
                                   decltype(t13)>::value);

        // Binomial product
        term<element<1>, monomial<multiplier<3>, factor<degree<2>, 2>>, monomial<multiplier<2>, factor<degree<1>, 2>, factor<degree<2>, 4>>> t4;
        auto t34 = t3 * t4;
        static_assert(
            std::is_same<term<element<1>, monomial<multiplier<-6>, factor<degree<2>, 2>, factor<degree<1>, 4>>,
                              monomial<multiplier<-4>, factor<degree<1>, 2>, factor<degree<3>, 4>>,
                              monomial<multiplier<3>, factor<degree<1>, 1>, factor<degree<2>, 2>, factor<degree<2>, 3>>,
                              monomial<multiplier<2>, factor<degree<1>, 1>, factor<degree<1>, 2>, factor<degree<2>, 3>, factor<degree<2>, 4>>>,
                         decltype(t34)>::value);
        static_assert(std::is_same<decltype(t4 * t3), decltype(t34)>::value);

        // Term cancellation
        term<element<2>, monomial<multiplier<1>, factor<degree<1>, 1>>, monomial<multiplier<1>, factor<degree<1>, 2>>> t5;
        term<element<2>, monomial<multiplier<1>, factor<degree<1>, 1>>, monomial<multiplier<-1>, factor<degree<1>, 2>>> t6;
        auto t56 = t5 * t6;
        static_assert(
            std::is_same<term<element<2>, monomial<multiplier<1>, factor<degree<2>, 1>>, monomial<multiplier<-1>, factor<degree<2>, 2>>>,
                         decltype(t56)>::value);
        static_assert(std::is_same<decltype(t6 * t5), decltype(t56)>::value);

        // Permutation test
        static_assert(std::is_same<decltype(t1 * t2 * t3 * t4), decltype(t4 * t3 * t2 * t1)>::value);
        static_assert(std::is_same<decltype(t1 * t3 * t2 * t4), decltype(t4 * t2 * t3 * t1)>::value);
        static_assert(std::is_same<decltype(t3 * t3 * t4), decltype(t4 * t3 * t3)>::value);
        static_assert(std::is_same<decltype(t3 * t3 * t4), decltype(t3 * t4 * t3)>::value);
    }
}

TEST_SUITE_END();