#pragma once

#include "finite_algebra.hpp"
#include "ga.hpp"

// The 3D Euclidean Geometric Algebra

namespace gal
{
namespace ega
{
using metric = ::gal::metric<3, 0, 0>;

struct algebra : public ga::algebra<metric>
{
};

// Left contraction
// NOTE: We choose this operator because like the left contraction, >> is right associative
template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator>>(::gal::multivector<void, I...> lhs, ::gal::multivector<void, J...> rhs) noexcept
{
    return ::gal::detail::product<::gal::ga::module<algebra>::contract>(lhs, rhs);
}

template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator^(::gal::multivector<void, I...> lhs, ::gal::multivector<void, J...> rhs) noexcept
{
    return ::gal::detail::product<::gal::ga::module<algebra>::exterior>(lhs, rhs);
}

template <typename... I, typename... J>
[[nodiscard]] constexpr auto operator*(::gal::multivector<void, I...> lhs, ::gal::multivector<void, J...> rhs) noexcept
{
    return ::gal::detail::product<::gal::ga::module<algebra>::geometric>(lhs, rhs);
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

template <size_t ID>
using t = term<element<0>, monomial<identity, generator<tag<ID>>>>;
template <size_t ID>
using t0 = term<element<0b1>, monomial<identity, generator<tag<ID>>>>;
template <size_t ID>
using t1 = term<element<0b10>, monomial<identity, generator<tag<ID>>>>;
template <size_t ID>
using t2 = term<element<0b100>, monomial<identity, generator<tag<ID>>>>;
template <size_t ID>
using t012 = term<element<0b111>, monomial<identity, generator<tag<ID>>>>;

using e = multivector<void, term<element<0>, monomial<identity>>>;
using e0 = multivector<void, term<element<0b1>, monomial<identity>>>;
using e1 = multivector<void, term<element<0b10>, monomial<identity>>>;
using e2 = multivector<void, term<element<0b100>, monomial<identity>>>;
using e01 = multivector<void, term<element<0b11>, monomial<identity>>>;
using e12 = multivector<void, term<element<0b110>, monomial<identity>>>;
using e012 = multivector<void, term<element<0b111>, monomial<identity>>>;

template <size_t ID, typename F = float>
struct vector3
{
    constexpr static size_t id = ID;
    using type = ::gal::multivector<void, t0<ID>, t1<ID + 1>, t2<ID + 2>>;
    F x;
    F y;
    F z;
};

template <size_t ID, typename F = float>
struct point3
{
    constexpr static size_t id = ID;
    using type = ::gal::multivector<void, t0<ID>, t1<ID + 1>, t2<ID + 2>>;
    F x;
    F y;
    F z;
};
} // namespace pga
} // namespace gal