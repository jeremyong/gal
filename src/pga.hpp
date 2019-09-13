#pragma once

#include "ga.hpp"

// The Projective Geometric Algebra for representing Euclidean 3-space
// NOTE: In comments, we often write "the PGA" to mean literally "the projective geometric algebra."

namespace gal
{
namespace pga
{
// NOTE: the inner product of e0 can be set to +1 or -1 without any change in the algebra's geometric
// interpretation. Here, we opt to define e0^2 := 1 by convention
using metric = ::gal::metric<0, 3, 1>;

// The PGA is a graded algebra with 16 basis elements
struct algebra : public ga::algebra<metric>
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

template <typename V, typename T>
[[nodiscard]] constexpr auto conjugate(V versor, T value) noexcept
{
    return versor * value * ~versor;
}

template <typename M>
[[nodiscard]] constexpr auto dual(M input) noexcept
{
    return ::gal::dual<metric>(input);
}

template <size_t ID, typename F = float>
struct vector3
{
    constexpr static size_t id = ID;
    using multivector_t = ::gal::multivector<void,
                                             term<element<0b1>, monomial<multiplier<1>, factor<degree<1>, ID>>>,
                                             term<element<0b10>, monomial<multiplier<1>, factor<degree<1>, ID + 1>>>,
                                             term<element<0b100>, monomial<multiplier<1>, factor<degree<1>, ID + 2>>>>;
    F x;
    F y;
    F z;
};

template <typename F>
struct vector3<0, F>
{
    constexpr static size_t id = 0;
    using multivector_t = ::gal::multivector<void,
                                             term<element<0b1>, monomial<multiplier<1>, factor<degree<1>, 0>>>,
                                             term<element<0b10>, monomial<multiplier<1>, factor<degree<1>, 0>>>,
                                             term<element<0b100>, monomial<multiplier<1>, factor<degree<1>, 0>>>>;
    F x;
    F y;
    F z;
};

template <size_t ID = 0, typename F = float>
struct point3
{
    constexpr static size_t id = ID;
    using multivector_t = ::gal::multivector<void,
                                             term<element<0b1>, monomial<multiplier<1>, factor<degree<1>, ID>>>,
                                             term<element<0b10>, monomial<multiplier<1>, factor<degree<1>, ID + 1>>>,
                                             term<element<0b100>, monomial<multiplier<1>, factor<degree<1>, ID + 2>>>,
                                             term<element<0b1000>, monomial<multiplier<1>>>>;
    F x;
    F y;
    F z;
};

template <typename F>
struct point3<0, F>
{
    constexpr static size_t id = 0;
    using multivector_t = ::gal::multivector<void,
                                             term<element<0b1>, monomial<multiplier<1>, factor<degree<1>, 0>>>,
                                             term<element<0b10>, monomial<multiplier<1>, factor<degree<1>, 0>>>,
                                             term<element<0b100>, monomial<multiplier<1>, factor<degree<1>, 0>>>,
                                             term<element<0b1000>, monomial<multiplier<1>>>>;
    F x;
    F y;
    F z;
};
} // namespace pga
} // namespace gal