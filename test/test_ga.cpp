#include "test_util.hpp"

#include <doctest/doctest.h>
#include <pga.hpp>
#include <formatters.hpp>

using namespace gal;
using namespace gal::pga;

static multivector<void, term<element<0>, monomial<multiplier<1>>>> e;
static multivector<void, term<element<1>, monomial<multiplier<1>>>> e0;
static multivector<void, term<element<0b10>, monomial<multiplier<1>>>> e1;
static multivector<void, term<element<0b100>, monomial<multiplier<1>>>> e2;
static multivector<void, term<element<0b1000>, monomial<multiplier<1>>>> e3;

TEST_SUITE_BEGIN("geometric-algebra");

TEST_CASE("blade-contraction")
{
    SUBCASE("contract-same")
    {
        auto [e1, p1] = ga::algebra<gal::pga::metric>::contract<1, 1>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 1);
    }

    SUBCASE("contract-orthogonal")
    {
        auto [e1, p1] = ga::algebra<gal::pga::metric>::contract<1, 2>();
        auto [e2, p2] = ga::algebra<gal::pga::metric>::contract<2, 1>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);
        CHECK_EQ(e2, 0);
        CHECK_EQ(p2, 0);
    }

    SUBCASE("contract-higher-grade-to-lower-grade")
    {
        auto [e1, p1] = ga::algebra<gal::pga::metric>::contract<0b110, 0b1>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);

        auto [e2, p2] = gal::ga::algebra<gal::pga::metric>::contract<0b11, 0b1>();
        CHECK_EQ(e2, 0);
        CHECK_EQ(p2, 0);
    }

    SUBCASE("contract-lower-grade-to-higher-grade")
    {
        auto [e1, p1] = ga::algebra<gal::pga::metric>::contract<0b1, 0b110>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);

        auto [e2, p2] = ga::algebra<gal::pga::metric>::contract<0b1, 0b11>();
        CHECK_EQ(e2, 0b10);
        CHECK_EQ(p2, 1);

        auto [e3, p3] = ga::algebra<gal::pga::metric>::contract<0b10, 0b11>();
        CHECK_EQ(e3, 0b01);
        CHECK_EQ(p3, -1);
    }

    SUBCASE("contract-blade-to-blade")
    {
        auto [e1, p1] = ga::algebra<gal::pga::metric>::contract<0b11, 0b111>();
        CHECK_EQ(e1, 0b100);
        CHECK_EQ(p1, -1);

        auto [e2, p2] = ga::algebra<gal::pga::metric>::contract<0b11, 0b1100>();
        CHECK_EQ(e2, 0);
        CHECK_EQ(p2, 0);

        auto [e3, p3] = ga::algebra<gal::pga::metric>::contract<0b11, 0b11>();
        CHECK_EQ(e3, 0);
        CHECK_EQ(p3, -1);

        auto [e4, p4] = ga::algebra<gal::pga::metric>::contract<0b110, 0b111>();
        CHECK_EQ(e4, 1);
        CHECK_EQ(p4, 1);
    }
}

TEST_CASE("blade-wedge")
{
    SUBCASE("wedge-same")
    {
        auto [e1, p1] = ga::algebra<gal::pga::metric>::exterior<1, 1>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);
    }

    SUBCASE("wedge-orthogonal")
    {
        auto [e1, p1] = ga::algebra<gal::pga::metric>::exterior<0b1, 0b10>();
        CHECK_EQ(e1, 0b11);
        CHECK_EQ(p1, 1);
        auto [e2, p2] = ga::algebra<gal::pga::metric>::exterior<0b10, 0b1>();
        CHECK_EQ(e2, 0b11);
        CHECK_EQ(p2, -1);
        auto [e3, p3] = ga::algebra<gal::pga::metric>::exterior<0b10, 0b100>();
        CHECK_EQ(e3, 0b110);
        CHECK_EQ(p3, 1);
        auto [e4, p4] = ga::algebra<gal::pga::metric>::exterior<0b100, 0b10>();
        CHECK_EQ(e4, 0b110);
        CHECK_EQ(p4, -1);
    }

    SUBCASE("wedge-dependent-blades")
    {
        auto [e1, p1] = ga::algebra<gal::pga::metric>::exterior<0b11, 0b101>();
        CHECK_EQ(e1, 0);
        CHECK_EQ(p1, 0);
        auto [e2, p2] = ga::algebra<gal::pga::metric>::exterior<0b11, 0b110>();
        CHECK_EQ(e2, 0);
        CHECK_EQ(p2, 0);
        auto [e3, p3] = ga::algebra<gal::pga::metric>::exterior<0b1100, 0b1100>();
        CHECK_EQ(e3, 0);
        CHECK_EQ(p3, 0);
        auto [e4, p4] = ga::algebra<gal::pga::metric>::exterior<0b1110, 0b1001>();
        CHECK_EQ(e4, 0);
        CHECK_EQ(p4, 0);
    }

    SUBCASE("wedge-independent-blades")
    {
        auto [e1, p1] = ga::algebra<gal::pga::metric>::exterior<0b11, 0b1100>();
        CHECK_EQ(e1, 0b1111);
        CHECK_EQ(p1, 1);
        auto [e2, p2] = ga::algebra<gal::pga::metric>::exterior<0b1100, 0b11>();
        CHECK_EQ(e2, 0b1111);
        CHECK_EQ(p2, 1);
        auto [e3, p3] = ga::algebra<gal::pga::metric>::exterior<0b1110, 0b1>();
        CHECK_EQ(e3, 0b1111);
        CHECK_EQ(p3, -1);
        auto [e4, p4] = ga::algebra<gal::pga::metric>::exterior<0b1, 0b1110>();
        CHECK_EQ(e4, 0b1111);
        CHECK_EQ(p4, 1);
    }
}

TEST_CASE("multivector-contraction")
{
    SUBCASE("single-element")
    {
        multivector<void, term<element<2>, monomial<multiplier<1>, factor<degree<1>, 1>>>> m1{};
        multivector<void, term<element<2>, monomial<multiplier<2>, factor<degree<1>, 1>, factor<degree<2>, 2>>>> m2{};
        auto m12 = m1 >> m2;
        static_assert(
            std::is_same<multivector<void, term<element<0>, monomial<multiplier<2>, factor<degree<2>, 1>, factor<degree<2>, 2>>>>, decltype(m12)>::value);
    }

    SUBCASE("single-element-orthogonal")
    {
        multivector<void, term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>>> m1{};
        multivector<void, term<element<2>, monomial<multiplier<1>, factor<degree<1>, 2>>>> m2{};
        auto m12 = m1 >> m2;
        static_assert(std::is_same<multivector<void>, decltype(m12)>::value);
    }

    SUBCASE("distribute-one-to-many")
    {
        multivector<void, term<element<0b10>, monomial<multiplier<1>, factor<degree<1>, 1>>>> m1{};
        multivector<void, term<element<0b10>, monomial<multiplier<1>, factor<degree<1>, 2>>>,
                    term<element<0b100>, monomial<multiplier<1>, factor<degree<1>, 3>>>>
            m2{};
        auto m12 = m1 >> m2;
        static_assert(
            std::is_same<multivector<void, term<element<0>, monomial<multiplier<1>, factor<degree<1>, 1>, factor<degree<1>, 2>>>>, decltype(m12)>::value);
        // The contraction operator does not generally commute but it does in this case
        static_assert(std::is_same<decltype(m1 >> m2), decltype(m2 >> m1)>::value);
    }

    SUBCASE("trinomial")
    {
        auto m1 = vector3<1>::multivector_t{};
        auto m2 = vector3<4>::multivector_t{};
        // Equivalent to the inner vector product
        auto m12 = m1 >> m2;
        static_assert(
            std::is_same<multivector<void, term<element<0>, monomial<multiplier<1>, factor<degree<1>, 1>, factor<degree<1>, 4>>,
                                                monomial<multiplier<1>, factor<degree<1>, 2>, factor<degree<1>, 5>>,
                                                monomial<multiplier<1>, factor<degree<1>, 3>, factor<degree<1>, 6>>>>,
                         decltype(m12)>::value);
        static_assert(std::is_same<decltype(m2 >> m1), decltype(m1 >> m2)>::value);
        // Vector norm
        static_assert(
            std::is_same<multivector<void, term<element<0>, monomial<multiplier<1>, factor<degree<2>, 1>>,
                                                monomial<multiplier<1>, factor<degree<2>, 2>>, monomial<multiplier<1>, factor<degree<2>, 3>>>>,
                         decltype(m1 >> m1)>::value);
    }
}

TEST_CASE("multivector-wedge")
{
    SUBCASE("single-element")
    {
        multivector<void, term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>>> m1{};
        multivector<void, term<element<1>, monomial<multiplier<1>, factor<degree<1>, 2>>>> m2{};
        auto m12 = m1 ^ m2;
        static_assert(std::is_same<multivector<void>, decltype(m12)>::value);
    }

    SUBCASE("single-element-orthogonal")
    {
        multivector<void, term<element<1>, monomial<multiplier<1>, factor<degree<1>, 1>>>> m1{};
        multivector<void, term<element<2>, monomial<multiplier<1>, factor<degree<1>, 2>>>> m2{};
        auto m12 = m1 ^ m2;
        static_assert(
            std::is_same<multivector<void, term<element<0b11>, monomial<multiplier<1>, factor<degree<1>, 1>, factor<degree<1>, 2>>>>, decltype(m12)>::value);
        static_assert(std::is_same<decltype(-(m2 ^ m1)), decltype(m12)>::value);
        static_assert(std::is_same<decltype(-m2 ^ m1), decltype(m12)>::value);
        static_assert(std::is_same<decltype(m2 ^ -m1), decltype(m12)>::value);
    }

    SUBCASE("binomial")
    {
        multivector<void, term<element<0b1>, monomial<multiplier<1>>>, term<element<0b10>, monomial<multiplier<1>>>> m1{};
        multivector<void, term<element<0b10>, monomial<multiplier<1>>>, term<element<0b100>, monomial<multiplier<1>>>> m2{};
        auto m12 = m1 ^ m2;
        static_assert(
            std::is_same<multivector<void, term<element<0b11>, monomial<multiplier<1>>>,
                                     term<element<0b101>, monomial<multiplier<1>>>, term<element<0b110>, monomial<multiplier<1>>>>,
                         decltype(m12)>::value);
        static_assert(std::is_same<decltype(m2 ^ -m1), decltype(m12)>::value);
    }
}

TEST_CASE("reversion")
{
    SUBCASE("single-element")
    {
        term<element<1>, monomial<multiplier<1>>> t;
        static_assert(std::is_same<decltype(~t), decltype(t)>::value);
    }

    SUBCASE("2-blade")
    {
        term<element<0b11>, monomial<multiplier<1>>> t;
        static_assert(std::is_same<decltype(~t), decltype(-t)>::value);
    }

    SUBCASE("3-blade")
    {
        term<element<0b111>, monomial<multiplier<1>>> t;
        static_assert(std::is_same<decltype(~t), decltype(-t)>::value);
    }

    SUBCASE("4-blade")
    {
        term<element<0b1111>, monomial<multiplier<1>>> t;
        static_assert(std::is_same<decltype(~t), decltype(t)>::value);
    }

    SUBCASE("multivector")
    {
        multivector<void, term<element<0b1>, monomial<multiplier<1>>>, term<element<0b11>, monomial<multiplier<1>>>> m{};
        static_assert(
            std::is_same<multivector<void, term<element<0b1>, monomial<multiplier<1>>>, term<element<0b11>, monomial<multiplier<-1>>>>,
                         decltype(~m)>::value);
    }
}

TEST_CASE("dual")
{
    SUBCASE("single-element")
    {
        multivector<void, term<element<1>, monomial<multiplier<1>>>> m;
        static_assert(std::is_same<multivector<void, term<element<0b1110>, monomial<multiplier<1>>>>, decltype(dual(m))>::value);
    }

    SUBCASE("2-blade")
    {
        multivector<void, term<element<0b11>, monomial<multiplier<1>>>> m;
        static_assert(std::is_same<multivector<void, term<element<0b1100>, monomial<multiplier<-1>>>>, decltype(dual(m))>::value);
    }
}

TEST_CASE("geometric-product")
{
    multivector<void, term<element<0>, monomial<multiplier<1>>>> s;
    multivector<void, term<element<1>, monomial<multiplier<1>>>> e1;
    multivector<void, term<element<0b10>, monomial<multiplier<1>>>> e2;
    multivector<void, term<element<0b100>, monomial<multiplier<1>>>> e3;
    multivector<void, term<element<0b1000>, monomial<multiplier<1>>>> e4;

    SUBCASE("self-inverse")
    {
        static_assert(std::is_same<decltype(e0 * e0), decltype(s)>::value);
    }

    SUBCASE("grade-raising")
    {
        static_assert(std::is_same<decltype(e1 * e2), multivector<void, term<element<0b11>, monomial<multiplier<1>>>>>::value);
    }

    SUBCASE("non-euclidean-metric")
    {
        static_assert(std::is_same<decltype(e4 * e4), multivector<void>>::value);
    }

    SUBCASE("vector-inverse")
    {
        auto v = e1 + e2;
        static_assert(std::is_same<decltype(v * v), multivector<void, term<element<0>, monomial<multiplier<2>>>>>::value);
    }
}

TEST_SUITE_END();
