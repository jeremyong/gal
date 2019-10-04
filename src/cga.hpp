#pragma once

#include "ga.hpp"

namespace gal
{
namespace cga
{
// The metric is defined here as the standard Minkowski spacetime. To extract the conformal representations,
// a change of basis is required where o = 1/2 * (e0 + e4) and inf = e4 - e0.
using metric = ::gal::metric<4, 1, 0>;

// The CGA is a graded algebra with 32 basis elements
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

template <typename M>
[[nodiscard]] constexpr auto dual(M input) noexcept
{
    return ::gal::dual<metric>(input);
}

// TODO: provide representations for points, planes, spheres, flats, etc.
} // namespace pga
} // namespace gal