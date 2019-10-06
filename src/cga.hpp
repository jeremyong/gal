#pragma once

#include "ga.hpp"

namespace gal
{
namespace cga
{
    // The metric is defined here as the standard Minkowski spacetime. To extract the conformal representations,
    // a change of basis is required where o = 1/2 * (e3 + e4) and inf = e4 - e3.
    // The element e3 here is the added unit norm basis element and the elements e0, e1, and e2 correspond to
    // the canonical Euclidean R3 basis representation.
    using metric = ::gal::metric<4, 1, 0>;

    // The CGA is a graded algebra with 32 basis elements
    using algebra = ga::algebra<metric>;

    // Left contraction
    // NOTE: We choose this operator because like the left contraction, >> is right associative
    template <typename... I, typename... J>
    [[nodiscard]] constexpr auto operator>>(::gal::multivector<void, I...> lhs, ::gal::multivector<void, J...> rhs) noexcept
    {
        return ::gal::detail::product<algebra::contract>(lhs, rhs);
    }

    template <typename... I, typename... J>
    [[nodiscard]] constexpr auto operator^(::gal::multivector<void, I...> lhs, ::gal::multivector<void, J...> rhs) noexcept
    {
        return ::gal::detail::product<algebra::exterior>(lhs, rhs);
    }

    template <typename... I, typename... J>
    [[nodiscard]] constexpr auto operator*(::gal::multivector<void, I...> lhs, ::gal::multivector<void, J...> rhs) noexcept
    {
        return ::gal::detail::product<algebra::geometric>(lhs, rhs);
    }

    template <typename V, typename T>
    [[nodiscard]] constexpr auto conjugate(V action, T subject) noexcept
    {
        return action * subject * ~action;
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

    template <int X, int Y, int Z>
    using point_t = multivector<void,
        term<element<0b1>, monomial<rational<X>>>,
        term<element<0b10>, monomial<rational<Y>>>,
        term<element<0b100>, monomial<rational<Z>>>,
        term<element<0b1000>, monomial<one_half>, monomial<rational<-(X * X + Y * Y + Z * Z), 2>>>,
        term<element<0b10000>, monomial<one_half>, monomial<rational<X * X + Y * Y + Z * Z, 2>>>>;

    // TODO: provide representations for points, planes, spheres, flats, etc.
    template <typename T = float>
    struct alignas(16) point
    {
        using value_t = T;

        template <typename Q, size_t ID>
        using p2
            = monomial<Q, generator<tag<ID, 0>, degree<2>>, generator<tag<ID, 1>, degree<2>>, generator<tag<ID, 2>, degree<2>>>;

        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b1>, monomial<one, generator<tag<ID, 0>>>>,
                                 term<element<0b10>, monomial<one, generator<tag<ID, 1>>>>,
                                 term<element<0b100>, monomial<one, generator<tag<ID, 2>>>>,
                                 term<element<0b1000>,
                                      monomial<one_half>,
                                      monomial<rational<-1, 2>, generator<tag<ID, 0>, degree<2>>>,
                                      monomial<rational<-1, 2>, generator<tag<ID, 1>, degree<2>>>,
                                      monomial<rational<-1, 2>, generator<tag<ID, 2>, degree<2>>>>,
                                 term<element<0b10000>,
                                      monomial<one_half>,
                                      monomial<one_half, generator<tag<ID, 0>, degree<2>>>,
                                      monomial<one_half, generator<tag<ID, 1>, degree<2>>>,
                                      monomial<one_half, generator<tag<ID, 2>, degree<2>>>>>;

        union
        {
            T x;
            T u;
        };

        union
        {
            T y;
            T v;
        };

        union
        {
            T z;
            T w;
        };

        [[nodiscard]] constexpr const T& operator[](size_t index) const noexcept
        {
            return *(reinterpret_cast<const T*>(this) + index);
        }

        [[nodiscard]] constexpr T& operator[](size_t index) noexcept
        {
            return *(reinterpret_cast<T*>(this) + index);
        }

        template <typename Engine, typename... I>
        [[nodiscard]] constexpr static point<T> convert(const Engine& engine, multivector<void, I...> mv) noexcept
        {
            auto x_e = extract<0b1>(mv);
            auto y_e = extract<0b10>(mv);
            auto z_e = extract<0b100>(mv);

            auto&& [x, y, z] = engine.template evaluate_terms<T>(x_e, y_e, z_e);
            return {x, y, z};
        }
    };
} // namespace cga
} // namespace gal