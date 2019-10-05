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

    // Left contraction
    // NOTE: We choose this operator because like the left contraction, >> is right associative
    template <typename... I, typename... J>
    [[nodiscard]] constexpr auto operator>>(::gal::multivector<void, I...> lhs, ::gal::multivector<void, J...> rhs) noexcept
    {
        return detail::product<algebra::contract>(lhs, rhs);
    }

    template <typename... I, typename... J>
    [[nodiscard]] constexpr auto operator^(::gal::multivector<void, I...> lhs, ::gal::multivector<void, J...> rhs) noexcept
    {
        return detail::product<algebra::exterior>(lhs, rhs);
    }

    template <typename... I, typename... J>
    [[nodiscard]] constexpr auto operator*(::gal::multivector<void, I...> lhs, ::gal::multivector<void, J...> rhs) noexcept
    {
        return detail::product<algebra::geometric>(lhs, rhs);
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
} // namespace ega
} // namespace gal