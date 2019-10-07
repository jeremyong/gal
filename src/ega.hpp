#pragma once

#include "finite_algebra.hpp"
#include "ga.hpp"

// The 3D Euclidean Geometric Algebra

namespace gal
{
namespace ega
{
    using metric = ::gal::metric<3, 0, 0>;

    using algebra = ga::algebra<metric>;

    GAL_OPERATORS(algebra);

    template <size_t ID>
    using t = term<element<0>, monomial<one, generator<tag<ID>>>>;
    template <size_t ID>
    using t0 = term<element<0b1>, monomial<one, generator<tag<ID>>>>;
    template <size_t ID>
    using t1 = term<element<0b10>, monomial<one, generator<tag<ID>>>>;
    template <size_t ID>
    using t2 = term<element<0b100>, monomial<one, generator<tag<ID>>>>;
    template <size_t ID>
    using t012 = term<element<0b111>, monomial<one, generator<tag<ID>>>>;

    using e    = multivector<void, term<element<0>, monomial<one>>>;
    using e0   = multivector<void, term<element<0b1>, monomial<one>>>;
    using e1   = multivector<void, term<element<0b10>, monomial<one>>>;
    using e2   = multivector<void, term<element<0b100>, monomial<one>>>;
    using e01  = multivector<void, term<element<0b11>, monomial<one>>>;
    using e12  = multivector<void, term<element<0b110>, monomial<one>>>;
    using e012 = multivector<void, term<element<0b111>, monomial<one>>>;

    template <size_t ID, typename F = float>
    struct vector3
    {
        constexpr static size_t id = ID;
        using type                 = ::gal::multivector<void, t0<ID>, t1<ID + 1>, t2<ID + 2>>;
        F x;
        F y;
        F z;
    };

    template <size_t ID, typename F = float>
    struct point3
    {
        constexpr static size_t id = ID;
        using type                 = ::gal::multivector<void, t0<ID>, t1<ID + 1>, t2<ID + 2>>;
        F x;
        F y;
        F z;
    };
    // TODO: this file isn't as standardized as the other geometries due to the restrictions in R3
    // and lack of use
} // namespace ega
} // namespace gal