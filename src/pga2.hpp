#pragma once

#include "finite_algebra.hpp"
#include "ga.hpp"

// Projective Geometric Algebra for representing Euclidean 2-space

namespace gal
{
namespace pga2
{
// NOTE: the inner product of e0 can be set to +1 or -1 without any change in the algebra's geometric
// interpretation. Here, we opt to define e0^2 := 1 by convention
using metric = ::gal::metric<0, 2, 1>;

// The PGA is a graded algebra with 16 basis elements
struct algebra : public ::gal::ga::algebra<metric>
{
};

// Left contraction
// NOTE: We choose this operator because like the left contraction, >> is right associative
template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator>>(::gal::multivector<void, I...> lhs, ::gal::multivector<void, J...> rhs) noexcept
{
    return ::gal::impl::product<::gal::ga::module<algebra>::contract>(lhs, rhs);
}

template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator^(::gal::multivector<void, I...> lhs, ::gal::multivector<void, J...> rhs) noexcept
{
    return ::gal::impl::product<::gal::ga::module<algebra>::exterior>(lhs, rhs);
}

template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator*(::gal::multivector<void, I...> lhs, ::gal::multivector<void, J...> rhs) noexcept
{
    return ::gal::impl::product<::gal::ga::module<algebra>::geometric>(lhs, rhs);
}

template <typename M>
[[nodiscard]] constexpr auto dual(M input) noexcept
{
    return ::gal::dual<metric>(input);
}

template <size_t ID, typename F = float>
struct vector2
{
    constexpr static size_t id = ID;
    using multivector_t = ::gal::multivector<void,
                                             term<element<0b10>, monomial<multiplier<1>, factor<degree<1>, ID>>>,
                                             term<element<0b100>, monomial<multiplier<1>, factor<degree<1>, ID + 1>>>>;
    F x;
    F y;
};

template <typename F>
struct vector2<0, F>
{
    constexpr static size_t id = ID;
    using multivector_t = ::gal::multivector<void,
                                             term<element<0b10>, monomial<multiplier<1>, factor<degree<1>, 0>>>,
                                             term<element<0b100>, monomial<multiplier<1>, factor<degree<1>, 0>>>>;
    F x;
    F y;
};

template <size_t ID = 0, typename F = float>
struct point2
{
    constexpr static size_t id = ID;
    using multivector_t = ::gal::multivector<void,
                                             term<element<0b1>, monomial<multiplier<1>>>,
                                             term<element<0b10>, monomial<multiplier<1>, factor<degree<1>, ID>>>,
                                             term<element<0b100>, monomial<multiplier<1>, factor<degree<1>, ID + 1>>>>;
    F x;
    F y;
};

template <typename F>
struct point2<0, F>
{
    constexpr static size_t id = ID;
    using multivector_t = ::gal::multivector<void,
                                             term<element<0b1>, monomial<multiplier<1>>>,
                                             term<element<0b10>, monomial<multiplier<1>, factor<degree<1>, 0>>>,
                                             term<element<0b100>, monomial<multiplier<1>, factor<degree<1>, 0>>>>;
    F x;
    F y;
};
} // namespace pga2
} // namespace gal