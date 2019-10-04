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
using metric = ::gal::metric<3, 0, 1>;

// The PGA is a graded algebra with 16 basis elements
struct algebra : public ga::algebra<metric>
{
};

// Left contraction
// NOTE: We choose this operator because like the left contraction, >> is right associative
template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator>>(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept
{
    return detail::product<ga::module<algebra>::contract>(lhs, rhs);
}

template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator^(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept
{
    return detail::product<ga::module<algebra>::exterior>(lhs, rhs);
}

template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator*(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept
{
    return detail::product<ga::module<algebra>::geometric>(lhs, rhs);
}

template <typename V, typename T>
[[nodiscard]] constexpr auto conjugate(V versor, T value) noexcept
{
    return versor * value * ~versor;
}

template <typename... I>
[[nodiscard]] constexpr auto operator!(multivector<void, I...> input) noexcept
{
    return dual<metric>(input);
}

template <typename M1, typename M2>
[[nodiscard]] constexpr auto operator|(M1 lhs, M2 rhs) noexcept
{
    return !(!lhs ^ !rhs);
}
template <size_t ID>
using t = term<element<0>, monomial<identity, generator<tag<ID>>>>;
template <size_t ID>
using t0 = term<element<0b1>, monomial<identity, generator<tag<ID>>>>;
template <size_t ID>
using t1 = term<element<0b10>, monomial<identity, generator<tag<ID>>>>;
template <size_t ID>
using t2 = term<element<0b100>, monomial<identity, generator<tag<ID>>>>;
template <size_t ID>
using t3 = term<element<0b1000>, monomial<identity, generator<tag<ID>>>>;
template <size_t ID>
using t012 = term<element<0b111>, monomial<identity, generator<tag<ID>>>>;
template <size_t ID>
using t013 = term<element<0b1011>, monomial<identity, generator<tag<ID>>>>;
template <size_t ID>
using t023 = term<element<0b1101>, monomial<identity, generator<tag<ID>>>>;
template <size_t ID>
using t123 = term<element<0b1110>, monomial<identity, generator<tag<ID>>>>;

using e = multivector<void, term<element<0>, monomial<identity>>>;
using e0 = multivector<void, term<element<0b1>, monomial<identity>>>;
using e1 = multivector<void, term<element<0b10>, monomial<identity>>>;
using e2 = multivector<void, term<element<0b100>, monomial<identity>>>;
using e3 = multivector<void, term<element<0b1000>, monomial<identity>>>;
using e012 = multivector<void, term<element<0b111>, monomial<identity>>>;
using e013 = multivector<void, term<element<0b1011>, monomial<identity>>>;
using e023 = multivector<void, term<element<0b1101>, monomial<identity>>>;
using e123 = multivector<void, term<element<0b1110>, monomial<identity>>>;

template <size_t ID, typename F = float>
struct plane
{
    constexpr static size_t id = ID;
    using type = multivector<void, t0<ID>, t1<ID + 1>, t2<ID + 2>, t3<ID + 3>>;
    F d;
    F x;
    F y;
    F z;
};

template <int D, int X, int Y, int Z>
using iplane = multivector<void,
                           term<element<1>, monomial<multiplier<D>>>,
                           term<element<0b10>, monomial<multiplier<X>>>,
                           term<element<0b100>, monomial<multiplier<Y>>>,
                           term<element<0b1000>, monomial<multiplier<Z>>>>;

template <size_t ID, typename F = float>
struct point
{
    constexpr static size_t id = ID;
    using type = multivector<void, t012<ID>, t013<ID + 1>, t023<ID + 2>, term<element<0b1110>, monomial<identity>>>;
    F x;
    F y;
    F z;
};

template <int X, int Y, int Z>
using ipoint = multivector<void,
                           term<element<0b111>, monomial<multiplier<-Z>>>,
                           term<element<0b1011>, monomial<multiplier<Y>>>,
                           term<element<0b1101>, monomial<multiplier<-X>>>,
                           term<element<0b1110>, monomial<identity>>>;

template <size_t ID, typename F = float>
struct direction
{
    using type = multivector<void, t012<ID>, t013<ID + 1>, t023<ID + 2>>;
    F x;
    F y;
    F z;
};
} // namespace pga
} // namespace gal