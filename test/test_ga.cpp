#include "test_util.hpp"

#include <doctest/doctest.h>
#include <pga.hpp>
#include <formatters.hpp>

using namespace gal;
using namespace gal::pga;

using e     = multivector<void, term<element<0>, monomial<one>>>;
using e0    = multivector<void, term<element<0b1>, monomial<one>>>;
using e1    = multivector<void, term<element<0b10>, monomial<one>>>;
using e2    = multivector<void, term<element<0b100>, monomial<one>>>;
using e3    = multivector<void, term<element<0b1000>, monomial<one>>>;
using e012  = multivector<void, term<element<0b111>, monomial<one>>>;
using e013  = multivector<void, term<element<0b1011>, monomial<one>>>;
using e023  = multivector<void, term<element<0b1101>, monomial<one>>>;
using e123  = multivector<void, term<element<0b1110>, monomial<one>>>;
using e0123 = multivector<void, term<element<0b1111>, monomial<one>>>;

TEST_SUITE_BEGIN("geometric-algebra");

TEST_CASE("symmetric-inner-product")
{
    using inner = ga::algebra<gal::pga::metric>::inner;
    SUBCASE("inner-same")
    {
        auto [e1, p1] = inner::inner_product<2, 2>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 1);
    }

    SUBCASE("inner-orthogonal")
    {
        auto [e1, p1] = inner::inner_product<1, 2>();
        auto [e2, p2] = inner::inner_product<2, 1>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);
        CHECK_EQ(e2, 0);
        CHECK_EQ(p2, 0);
    }

    SUBCASE("inner-mixed-grade")
    {
        auto [e1, p1] = inner::inner_product<0b110, 0b1>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);

        auto [e2, p2] = inner::inner_product<0b1, 0b110>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);

        auto [e3, p3] = inner::inner_product<0b110, 0b10>();
        CHECK_EQ(e3, 0b100);
        CHECK_EQ(p3, -1);

        auto [e4, p4] = inner::inner_product<0b10, 0b110>();
        CHECK_EQ(e4, 0b100);
        CHECK_EQ(p4, 1);
    }
}

TEST_CASE("blade-contraction")
{
    using contract = ga::algebra<gal::pga::metric>::contract;
    SUBCASE("contract-same")
    {
        auto [e1, p1] = contract::contract_product<2, 2>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 1);
    }

    SUBCASE("contract-orthogonal")
    {
        auto [e1, p1] = contract::contract_product<1, 2>();
        auto [e2, p2] = contract::contract_product<2, 1>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);
        CHECK_EQ(e2, 0);
        CHECK_EQ(p2, 0);
    }

    SUBCASE("contract-higher-grade-to-lower-grade")
    {
        auto [e1, p1] = contract::contract_product<0b110, 0b1>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);

        auto [e2, p2] = contract::contract_product<0b11, 0b1>();
        CHECK_EQ(e2, 0);
        CHECK_EQ(p2, 0);
    }

    SUBCASE("contract-lower-grade-to-higher-grade")
    {
        auto [e1, p1] = contract::contract_product<0b1, 0b110>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);

        auto [e2, p2] = contract::contract_product<0b10, 0b110>();
        CHECK_EQ(e2, 0b100);
        CHECK_EQ(p2, 1);

        auto [e3, p3] = contract::contract_product<0b1000, 0b1100>();
        CHECK_EQ(e3, 0b100);
        CHECK_EQ(p3, -1);
    }

    SUBCASE("contract-blade-to-blade")
    {
        auto [e1, p1] = contract::contract_product<0b110, 0b1110>();
        CHECK_EQ(e1, 0b1000);
        CHECK_EQ(p1, -1);

        auto [e2, p2] = contract::contract_product<0b11, 0b1100>();
        CHECK_EQ(e2, 0);
        CHECK_EQ(p2, 0);

        auto [e3, p3] = contract::contract_product<0b110, 0b110>();
        CHECK_EQ(e3, 0);
        CHECK_EQ(p3, -1);

        auto [e4, p4] = contract::contract_product<0b110, 0b111>();
        CHECK_EQ(e4, 1);
        CHECK_EQ(p4, -1);
    }
}

TEST_CASE("blade-wedge")
{
    using exterior = ga::algebra<gal::pga::metric>::exterior;
    SUBCASE("wedge-same")
    {
        auto [e1, p1] = exterior::exterior_product<1, 1>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);
    }

    SUBCASE("wedge-orthogonal")
    {
        auto [e1, p1] = exterior::exterior_product<0b1, 0b10>();
        CHECK_EQ(e1, 0b11);
        CHECK_EQ(p1, 1);
        auto [e2, p2] = exterior::exterior_product<0b10, 0b1>();
        CHECK_EQ(e2, 0b11);
        CHECK_EQ(p2, -1);
        auto [e3, p3] = exterior::exterior_product<0b10, 0b100>();
        CHECK_EQ(e3, 0b110);
        CHECK_EQ(p3, 1);
        auto [e4, p4] = exterior::exterior_product<0b100, 0b10>();
        CHECK_EQ(e4, 0b110);
        CHECK_EQ(p4, -1);
    }

    SUBCASE("wedge-dependent-blades")
    {
        auto [e1, p1] = exterior::exterior_product<0b11, 0b101>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);
        auto [e2, p2] = exterior::exterior_product<0b11, 0b110>();
        CHECK_EQ(e2, 0);
        CHECK_EQ(p2, 0);
        auto [e3, p3] = exterior::exterior_product<0b1100, 0b1100>();
        CHECK_EQ(e3, 0);
        CHECK_EQ(p3, 0);
        auto [e4, p4] = exterior::exterior_product<0b1110, 0b1001>();
        CHECK_EQ(e4, 0);
        CHECK_EQ(p4, 0);
    }

    SUBCASE("wedge-independent-blades")
    {
        auto [e1, p1] = exterior::exterior_product<0b11, 0b1100>();
        CHECK_EQ(e1, 0b1111);
        CHECK_EQ(p1, 1);
        auto [e2, p2] = exterior::exterior_product<0b1100, 0b11>();
        CHECK_EQ(e2, 0b1111);
        CHECK_EQ(p2, 1);
        auto [e3, p3] = exterior::exterior_product<0b1110, 0b1>();
        CHECK_EQ(e3, 0b1111);
        CHECK_EQ(p3, -1);
        auto [e4, p4] = exterior::exterior_product<0b1, 0b1110>();
        CHECK_EQ(e4, 0b1111);
        CHECK_EQ(p4, 1);
    }
}

TEST_CASE("multivector-contraction")
{
    SUBCASE("single-element")
    {
        multivector<void, term<element<2>, monomial<rational<1>, generator<tag<1>, degree<1>>>>> m1{};
        multivector<void, term<element<2>, monomial<rational<2>, generator<tag<1>, degree<1>>, generator<tag<2>, degree<2>>>>>
            m2{};
        auto m12 = m1 >> m2;
        static_assert(
            std::is_same<
                multivector<void,
                            term<element<0>, monomial<rational<2>, generator<tag<1>, degree<2>>, generator<tag<2>, degree<2>>>>>,
                decltype(m12)>::value);
    }

    SUBCASE("single-element-orthogonal")
    {
        multivector<void, term<element<1>, monomial<rational<1>, generator<tag<1>, degree<1>>>>> m1{};
        multivector<void, term<element<2>, monomial<rational<1>, generator<tag<2>, degree<1>>>>> m2{};
        auto m12 = m1 >> m2;
        static_assert(std::is_same<multivector<void>, decltype(m12)>::value);
    }

    SUBCASE("distribute-one-to-many")
    {
        multivector<void, term<element<0b10>, monomial<rational<1>, generator<tag<1>, degree<1>>>>> m1{};
        multivector<void,
                    term<element<0b10>, monomial<rational<1>, generator<tag<2>, degree<1>>>>,
                    term<element<0b100>, monomial<rational<1>, generator<tag<3>, degree<1>>>>>
            m2{};
        auto m12 = m1 >> m2;
        static_assert(
            std::is_same<
                multivector<void,
                            term<element<0>, monomial<rational<1>, generator<tag<1>, degree<1>>, generator<tag<2>, degree<1>>>>>,
                decltype(m12)>::value);
        // The contraction operator does not generally commute but it does in this case
        static_assert(std::is_same<decltype(m1 >> m2), decltype(m2 >> m1)>::value);
    }

    SUBCASE("polynomial")
    {
        auto p1 = plane<>::type<1>{};
        auto p2 = plane<>::type<2>{};
        // Equivalent to the inner vector product
        auto p12 = p1 >> p2;
        static_assert(
            std::is_same<
                multivector<void,
                            term<element<0>,
                                 monomial<rational<1>, generator<tag<1, 1>, degree<1>>, generator<tag<2, 1>, degree<1>>>,
                                 monomial<rational<1>, generator<tag<1, 2>, degree<1>>, generator<tag<2, 2>, degree<1>>>,
                                 monomial<rational<1>, generator<tag<1, 3>, degree<1>>, generator<tag<2, 3>, degree<1>>>>>,
                decltype(p12)>::value);
        static_assert(std::is_same<decltype(p2 >> p1), decltype(p1 >> p2)>::value);
        // Vector norm
        static_assert(std::is_same<multivector<void,
                                               term<element<0>,
                                                    monomial<rational<1>, generator<tag<1, 1>, degree<2>>>,
                                                    monomial<rational<1>, generator<tag<1, 2>, degree<2>>>,
                                                    monomial<rational<1>, generator<tag<1, 3>, degree<2>>>>>,
                                   decltype(p1 >> p1)>::value);
    }
}

TEST_CASE("multivector-wedge")
{
    SUBCASE("single-element")
    {
        multivector<void, term<element<1>, monomial<rational<1>, generator<tag<1>, degree<1>>>>> m1{};
        multivector<void, term<element<1>, monomial<rational<1>, generator<tag<2>, degree<1>>>>> m2{};
        auto m12 = m1 ^ m2;
        static_assert(std::is_same<multivector<void>, decltype(m12)>::value);
    }

    SUBCASE("single-element-orthogonal")
    {
        multivector<void, term<element<1>, monomial<rational<1>, generator<tag<1>, degree<1>>>>> m1{};
        multivector<void, term<element<2>, monomial<rational<1>, generator<tag<2>, degree<1>>>>> m2{};
        auto m12 = m1 ^ m2;
        static_assert(
            std::is_same<multivector<void,
                                     term<element<0b11>,
                                          monomial<rational<1>, generator<tag<1>, degree<1>>, generator<tag<2>, degree<1>>>>>,
                         decltype(m12)>::value);
        static_assert(std::is_same<decltype(-(m2 ^ m1)), decltype(m12)>::value);
        static_assert(std::is_same<decltype(-m2 ^ m1), decltype(m12)>::value);
        static_assert(std::is_same<decltype(m2 ^ -m1), decltype(m12)>::value);
    }

    SUBCASE("binomial")
    {
        multivector<void, term<element<0b1>, monomial<rational<1>>>, term<element<0b10>, monomial<rational<1>>>> m1{};
        multivector<void, term<element<0b10>, monomial<rational<1>>>, term<element<0b100>, monomial<rational<1>>>> m2{};
        auto m12 = m1 ^ m2;
        static_assert(std::is_same<multivector<void,
                                               term<element<0b11>, monomial<rational<1>>>,
                                               term<element<0b101>, monomial<rational<1>>>,
                                               term<element<0b110>, monomial<rational<1>>>>,
                                   decltype(m12)>::value);
        static_assert(std::is_same<decltype(m2 ^ -m1), decltype(m12)>::value);
    }
}

TEST_CASE("reversion")
{
    using namespace gal::pga;
    SUBCASE("single-element")
    {
        multivector<void, term<element<1>, monomial<rational<1>>>> t;
        static_assert(std::is_same<decltype(~t), decltype(t)>::value);
    }

    SUBCASE("2-blade")
    {
        multivector<void, term<element<0b11>, monomial<rational<1>>>> t;
        static_assert(std::is_same<decltype(~t), decltype(-t)>::value);
    }

    SUBCASE("3-blade")
    {
        multivector<void, term<element<0b111>, monomial<rational<1>>>> t;
        static_assert(std::is_same<decltype(~t), decltype(-t)>::value);
    }

    SUBCASE("4-blade")
    {
        multivector<void, term<element<0b1111>, monomial<rational<1>>>> t;
        static_assert(std::is_same<decltype(~t), decltype(t)>::value);
    }

    SUBCASE("multivector")
    {
        multivector<void, term<element<0b1>, monomial<rational<1>>>, term<element<0b11>, monomial<rational<1>>>> m{};
        static_assert(
            std::is_same<
                multivector<void, term<element<0b1>, monomial<rational<1>>>, term<element<0b11>, monomial<rational<-1>>>>,
                decltype(~m)>::value);
    }
}

TEST_CASE("dual")
{
    SUBCASE("pointcare dual parity")
    {
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 0>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 1>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 2>(), -1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 3>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 4>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 5>(), -1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 6>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 7>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 8>(), -1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 9>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 10>(), -1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 11>(), -1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 12>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 13>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 14>(), -1);
        CHECK_EQ(gal::ga::poincare_complement_parity<4, 15>(), 1);

        CHECK_EQ(gal::ga::poincare_complement_parity<3, 0>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<3, 1>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<3, 2>(), -1);
        CHECK_EQ(gal::ga::poincare_complement_parity<3, 3>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<3, 4>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<3, 5>(), -1);
        CHECK_EQ(gal::ga::poincare_complement_parity<3, 6>(), 1);
        CHECK_EQ(gal::ga::poincare_complement_parity<3, 7>(), 1);
    }

    SUBCASE("single-element")
    {
        multivector<void, term<element<1>, monomial<one>>> m;
        static_assert(std::is_same<multivector<void, term<element<0b1110>, monomial<one>>>, decltype(!m)>::value);
    }

    SUBCASE("2-blade")
    {
        multivector<void, term<element<0b101>, monomial<one>>> m;
        static_assert(std::is_same<multivector<void, term<element<0b1010>, monomial<rational<-1>>>>, decltype(!m)>::value);
    }
}

TEST_CASE("geometric-product")
{
    SUBCASE("self-inverse")
    {
        static_assert(std::is_same<decltype(e1{} * e1{}), e>::value);
    }

    SUBCASE("grade-raising")
    {
        static_assert(std::is_same<decltype(e1{} * e2{}), multivector<void, term<element<0b110>, monomial<one>>>>::value);
    }

    SUBCASE("non-euclidean-metric")
    {
        static_assert(std::is_same<decltype(e0{} * e0{}), multivector<void>>::value);
    }

    SUBCASE("vector-inverse")
    {
        auto v = e2{} + e3{};
        static_assert(std::is_same<decltype(v * v), multivector<void, term<element<0>, monomial<rational<2>>>>>::value);
    }
}

TEST_SUITE_END();
