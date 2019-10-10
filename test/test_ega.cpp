#include "test_util.hpp"

#include <doctest/doctest.h>
#include <gal/ega.hpp>
#include <gal/engine.hpp>
#include <gal/formatters.hpp>

using namespace gal::ega;

TEST_SUITE_BEGIN("euclidean-geometric-algebra");

TEST_CASE("rotors")
{
    using namespace gal;
    using e    = multivector<void, term<element<0>, monomial<one>>>;
    using e0   = multivector<void, term<element<0b1>, monomial<one>>>;
    using e1   = multivector<void, term<element<0b10>, monomial<one>>>;
    using e2   = multivector<void, term<element<0b100>, monomial<one>>>;
    using e01  = multivector<void, term<element<0b11>, monomial<one>>>;
    using e12  = multivector<void, term<element<0b110>, monomial<one>>>;
    using e012 = multivector<void, term<element<0b111>, monomial<one>>>;

    // TODO: modernize these tests to use the newer syntax
    SUBCASE("vector-reflection-through-vector")
    {
        auto v = e0{} + e1{};
        auto reflected = conjugate(e0{}, v);
        static_assert(std::is_same<decltype(e0{} - e1{}), decltype(reflected)>::value);

        static_assert(std::is_same<decltype(-e0{} - e1{}), decltype(conjugate(e2{}, v))>::value);
    }

    SUBCASE("blade-reflection-through-vector")
    {
        static_assert(std::is_same<e01, decltype(conjugate(e2{}, e01{}))>::value);
        static_assert(std::is_same<e01, decltype(conjugate(-e2{}, e01{}))>::value);
        static_assert(std::is_same<decltype(rational<2>{} * e12{}), decltype(conjugate(e0{} + e2{}, e01{}))>::value);
    }
}

TEST_SUITE_END();