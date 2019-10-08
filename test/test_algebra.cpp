#include <doctest/doctest.h>

#include <finite_algebra.hpp>
#include <formatters.hpp>

#include "test_util.hpp"

using namespace gal;

TEST_SUITE_BEGIN("finite-algebra");

// The majority of these tests are actually compile-time tests (compilation and correctness are equivalent)
TEST_CASE("rational-perturbation")
{
    SUBCASE("small-undisturbed")
    {
        constexpr auto q = detail::overflow_gate(rational<1, 128>{});
        static_assert(q.num == 1);
        static_assert(q.den == 128);
    }

    SUBCASE("perturb-irreducible-rationals")
    {
        // TODO: Ensure accuracy of rational mediant approximation
    }
}

TEST_CASE("term-arithmetic")
{
    SUBCASE("different-degree-and-generator")
    {
        term<element<1>, monomial<one, generator<tag<1>, degree<1>>>> t1;
        term<element<1>, monomial<rational<2>, generator<tag<2>, degree<2>>>> t2;
        constexpr auto t12 = t1 + t2;
        static_assert(decltype(t12)::size == 2);
        static_assert(decltype(t12)::first_t::rational_t::num == 1);
        static_assert(decltype(t12)::subsequent_t::first_t::rational_t::num == 2);
    }

    SUBCASE("term-cancellation")
    {
        term<element<1>, monomial<one, generator<tag<1>, degree<1>>>> t1;
        term<element<1>, monomial<rational<-1>, generator<tag<1>, degree<1>>>> t2;
        constexpr auto t12 = t1 + t2;
        static_assert(decltype(t12)::size == 0);
    }

    SUBCASE("term-addition")
    {
        term<element<1>, monomial<rational<-2>, generator<tag<4>, degree<1>>>,
             monomial<one, generator<tag<1>, degree<1>>, generator<tag<3>, degree<2>>>>
            t1;
        term<element<1>, monomial<rational<2>, generator<tag<4>, degree<1>>>,
             monomial<rational<-1>, generator<tag<1>, degree<1>>, generator<tag<3>, degree<2>>>>
            t2;
        constexpr auto t12 = t1 + t2;
        static_assert(decltype(t12)::size == 0);
    }

    SUBCASE("successive-addition")
    {
        term<element<1>, monomial<one, generator<tag<1>, degree<1>>>> t1;
        constexpr auto t111 = t1 + t1 + t1;
        static_assert(decltype(t111)::size == 1);
        static_assert(decltype(t111)::first_t::rational_t::num == 3);
    }

    SUBCASE("different-monomial-addition")
    {
        term<element<1>, monomial<rational<-2>, generator<tag<4>, degree<1>>>, monomial<one, generator<tag<1>, degree<1>>, generator<tag<3>, degree<2>>>> t1;
        term<element<1>, monomial<rational<2>, generator<tag<4>, degree<1>>>, monomial<one, generator<tag<1>, degree<1>>, generator<tag<3>, degree<2>>>> t2;
        term<element<1>, monomial<rational<2>, generator<tag<4>, degree<2>>>, monomial<one, generator<tag<1>, degree<1>>, generator<tag<3>, degree<2>>>> t3;
        constexpr auto t12 = t1 + t2;
        static_assert(decltype(t12)::size == 1);
        static_assert(
            std::is_same<monomial<rational<2>, generator<tag<1>, degree<1>>, generator<tag<3>, degree<2>>>, typename decltype(t12)::first_t>::value);

        constexpr auto t13 = t1 + t3;
        static_assert(t13.size == 3);
        static_assert(std::is_same<monomial<rational<-2>, generator<tag<4>, degree<1>>>, typename decltype(t13)::first_t>::value);
        static_assert(std::is_same<monomial<rational<2>, generator<tag<4>, degree<2>>>, typename decltype(t13)::subsequent_t::first_t>::value);
        static_assert(std::is_same<monomial<rational<2>, generator<tag<1>, degree<1>>, generator<tag<3>, degree<2>>>,
                                   typename decltype(t13)::subsequent_t::subsequent_t::first_t>::value);
    }
}

TEST_CASE("multivector-addition")
{
    SUBCASE("single-term-addition")
    {
        multivector<void, term<element<1>, monomial<one, generator<tag<1>, degree<1>>>>> v1;
        multivector<void, term<element<1>, monomial<rational<2>, generator<tag<2>, degree<2>>>>> v2;
        auto v12 = v1 + v2;
        static_assert(v12.size == 1);
        static_assert(
            std::is_same<term<element<1>, monomial<one, generator<tag<1>, degree<1>>>, monomial<rational<2>, generator<tag<2>, degree<2>>>>,
                         decltype(v12)::first_t>::value);

        multivector<void, term<element<1>, monomial<rational<2>, generator<tag<1>, degree<1>>>>> v3;
        auto v13 = v1 + v3;
        static_assert(v13.size == 1);
        static_assert(std::is_same<term<element<1>, monomial<rational<3>, generator<tag<1>, degree<1>>>>, decltype(v13)::first_t>::value);
    }

    SUBCASE("multiple-term-addition")
    {
        multivector<void, term<element<1>, monomial<one, generator<tag<1>, degree<1>>>>> v1;
        multivector<void, term<element<2>, monomial<one, generator<tag<1>, degree<1>>>>> v2;
        auto v12 = v1 + v2;
        static_assert(v12.size == 2);
        static_assert(std::is_same<term<element<1>, monomial<one, generator<tag<1>, degree<1>>>>, decltype(v12)::first_t>::value);
        static_assert(
            std::is_same<term<element<2>, monomial<one, generator<tag<1>, degree<1>>>>, decltype(v12)::subsequent_t::first_t>::value);

        multivector<void, term<element<1>, monomial<one, generator<tag<1>, degree<1>>>>, term<element<2>, monomial<one, generator<tag<1>, degree<1>>>>> v3;
        auto v13 = v1 + v3;
        static_assert(v13.size == 2);
        static_assert(std::is_same<term<element<1>, monomial<rational<2>, generator<tag<1>, degree<1>>>>, decltype(v13)::first_t>::value);
        static_assert(
            std::is_same<term<element<2>, monomial<one, generator<tag<1>, degree<1>>>>, decltype(v13)::subsequent_t::first_t>::value);

        multivector<void, term<element<2>, monomial<one, generator<tag<1>, degree<1>>>>, term<element<3>, monomial<one, generator<tag<1>, degree<1>>>>> v4;
        auto v14 = v1 + v4;
        static_assert(v14.size == 3);
        static_assert(std::is_same<term<element<1>, monomial<one, generator<tag<1>, degree<1>>>>, decltype(v14)::first_t>::value);
        static_assert(
            std::is_same<term<element<2>, monomial<one, generator<tag<1>, degree<1>>>>, decltype(v14)::subsequent_t::first_t>::value);
        static_assert(
            std::is_same<term<element<3>, monomial<one, generator<tag<1>, degree<1>>>>, decltype(v14)::subsequent_t::subsequent_t::first_t>::value);
    }

    SUBCASE("term-cancellation")
    {
        multivector<void, term<element<1>, monomial<one, generator<tag<1>, degree<1>>>>> v1;
        multivector<void, term<element<1>, monomial<rational<-1>, generator<tag<1>, degree<1>>>>> v2;
        auto v12 = v1 + v2;
        static_assert(v12.size == 0);

        multivector<void, term<element<1>, monomial<one, generator<tag<2>, degree<1>>>>> v3;
        multivector<void, term<element<1>, monomial<rational<-1>, generator<tag<1>, degree<1>>>, monomial<rational<-1>, generator<tag<2>, degree<1>>>>> v4;
        auto v134 = v1 + v3 + v4;
        static_assert(v134.size == 0);
    }
}

TEST_CASE("term-multiplication")
{
    SUBCASE("generator-product")
    {
        monomial<one, generator<tag<1>, degree<1>>> f1;
        monomial<one, generator<tag<2>, degree<1>>> f2;
        static_assert(std::is_same<monomial<one, generator<tag<1>, degree<1>>, generator<tag<2>, degree<1>>>, decltype(f1 * f2)>::value);

        monomial<one, generator<tag<1>, degree<2>>> f3;
        static_assert(std::is_same<monomial<one, generator<tag<1>, degree<3>>>, decltype(f1 * f3)>::value);
    }

    SUBCASE("monomial-product")
    {
        monomial<one, generator<tag<1>, degree<1>>> m1;
        monomial<rational<2>, generator<tag<2>, degree<1>>> m2;
        auto m12 = m1 * m2;
        static_assert(
            std::is_same<monomial<rational<2>, generator<tag<1>, degree<1>>, generator<tag<2>, degree<1>>>, decltype(m12)>::value);
        // monomial multiplication commutes (order should be the same)
        static_assert(
            std::is_same<monomial<rational<2>, generator<tag<1>, degree<1>>, generator<tag<2>, degree<1>>>, decltype(m2 * m1)>::value);

        monomial<rational<3>, generator<tag<1>, degree<2>>> m3;
        auto m13 = m1 * m3;
        static_assert(std::is_same<monomial<rational<3>, generator<tag<1>, degree<3>>>, decltype(m13)>::value);

        monomial<one, generator<tag<1>, degree<2>>, generator<tag<3>, degree<1>>, generator<tag<5>, degree<3>>> m4;
        auto m124 = m1 * m2 * m4;
        static_assert(std::is_same<monomial<rational<2>, generator<tag<1>, degree<3>>, generator<tag<2>, degree<1>>,
                                            generator<tag<3>, degree<1>>, generator<tag<5>, degree<3>>>,
                                   decltype(m124)>::value);
        static_assert(std::is_same<decltype(m1 * m4 * m2), decltype(m124)>::value);
        static_assert(std::is_same<decltype(m2 * m1 * m4), decltype(m124)>::value);
        static_assert(std::is_same<decltype(m4 * m1 * m2), decltype(m124)>::value);
        static_assert(std::is_same<decltype(m4 * m2 * m1), decltype(m124)>::value);

        monomial<rational<2>> m5;
        auto m15 = m1 * m5;
        static_assert(std::is_same<monomial<rational<2>, generator<tag<1>, degree<1>>>, decltype(m15)>::value);
        monomial<one> m6;
        auto m16 = m1 * m6;
        static_assert(std::is_same<monomial<one, generator<tag<1>, degree<1>>>, decltype(m16)>::value);
    }

    SUBCASE("term-product")
    {
        // Single term product
        term<element<1>, monomial<one, generator<tag<1>, degree<1>>>> t1;
        term<element<1>, monomial<rational<2>, generator<tag<2>, degree<2>>>> t2;
        auto t12 = t1 * t2;
        static_assert(std::is_same<term<element<1>, monomial<rational<2>, generator<tag<1>, degree<1>>, generator<tag<2>, degree<2>>>>, decltype(t12)>::value);

        // One term into two
        term<element<1>, monomial<rational<-2>, generator<tag<4>, degree<1>>>, monomial<one, generator<tag<1>, degree<1>>, generator<tag<3>, degree<2>>>> t3;
        auto t13 = t1 * t3;
        static_assert(
            std::is_same<term<element<1>, monomial<rational<-2>, generator<tag<1>, degree<1>>, generator<tag<4>, degree<1>>>,
                              monomial<one, generator<tag<1>, degree<2>>, generator<tag<3>, degree<2>>>>,
                         decltype(t13)>::value);

        // Binomial product
        term<element<1>, monomial<rational<3>, generator<tag<2>, degree<2>>>, monomial<rational<2>, generator<tag<2>, degree<1>>, generator<tag<4>, degree<2>>>> t4;
        auto t34 = t3 * t4;
        static_assert(
            std::is_same<term<element<1>, monomial<rational<-6>, generator<tag<2>, degree<2>>, generator<tag<4>, degree<1>>>,
                              monomial<rational<-4>, generator<tag<2>, degree<1>>, generator<tag<4>, degree<3>>>,
                              monomial<rational<3>, generator<tag<1>, degree<1>>, generator<tag<2>, degree<2>>, generator<tag<3>, degree<2>>>,
                              monomial<rational<2>, generator<tag<1>, degree<1>>, generator<tag<2>, degree<1>>,
                                       generator<tag<3>, degree<2>>, generator<tag<4>, degree<2>>>>,
                         decltype(t34)>::value);
        static_assert(std::is_same<decltype(t4 * t3), decltype(t34)>::value);

        // Term cancellation
        term<element<2>, monomial<one, generator<tag<1>, degree<1>>>, monomial<one, generator<tag<2>, degree<1>>>> t5;
        term<element<2>, monomial<one, generator<tag<1>, degree<1>>>, monomial<rational<-1>, generator<tag<2>, degree<1>>>> t6;
        auto t56 = t5 * t6;
        static_assert(
            std::is_same<term<element<2>, monomial<one, generator<tag<1>, degree<2>>>, monomial<rational<-1>, generator<tag<2>, degree<2>>>>,
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