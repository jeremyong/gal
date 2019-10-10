#include "test_util.hpp"

#include <doctest/doctest.h>
#include <gal/formatters.hpp>
#include <gal/cga.hpp>


TEST_SUITE_BEGIN("conformal-geometric-algebra");

TEST_CASE("metric")
{
    using namespace gal;
    using namespace gal::cga;
    using contract = cga_algebra::contract;
    SUBCASE("cga-contraction")
    {
        static_assert(std::is_same<decltype(e_inf >> e_inf), multivector<void>>::value);

        static_assert(
            std::is_same<decltype(e_o >> e_inf), multivector<void, term<element<0>, monomial<minus_one>>>>::value);

        // Point at (1, 2, 1)
        point_t<1, 2, 1> p1;
        // In the Conformal Geometric Algebra, the contraction of a point onto itself is exactly 0
        static_assert(std::is_same<decltype(p1 >> p1), multivector<void>>::value);
    }

    SUBCASE("blade-reverse")
    {
        static_assert(!cga_metric::is_diagonal);
        multivector<void, term<element<0b11>, monomial<one>>> mv;
        static_assert(std::is_same<decltype(mv * ~mv), multivector<void, term<element<0>, monomial<one>>>>::value);
        static_assert(std::is_same<decltype(-cga::pseudoscalar::inverse), std::decay<decltype(cga::pseudoscalar::value)>::type>::value);
    }

    SUBCASE("rotors")
    {
        auto r = e + e_x * e_y;
        auto line = e_x * e_y + e_x * e_inf + e_y * e_inf;
        // fmt::print("{}, {}\n", r, r * ~r);
        // fmt::print("{} -> {}\n", line, conjugate(r, line));
    }
}

TEST_SUITE_END();