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
    using algebra = ga::algebra<metric>;

    // Left contraction
    // NOTE: We choose this operator because like the left contraction, >> is right associative
    template <typename... I, typename... J>
    [[nodiscard]] constexpr auto operator>>(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept
    {
        return detail::product<algebra::contract>(lhs, rhs);
    }

    template <typename... I, typename... J>
    [[nodiscard]] constexpr auto operator^(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept
    {
        return detail::product<algebra::exterior>(lhs, rhs);
    }

    template <typename... I, typename... J>
    [[nodiscard]] constexpr auto operator*(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept
    {
        return detail::product<algebra::geometric>(lhs, rhs);
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

    using e    = multivector<void, term<element<0>, monomial<one>>>;
    using e0   = multivector<void, term<element<0b1>, monomial<one>>>;
    using e1   = multivector<void, term<element<0b10>, monomial<one>>>;
    using e2   = multivector<void, term<element<0b100>, monomial<one>>>;
    using e3   = multivector<void, term<element<0b1000>, monomial<one>>>;
    using e012 = multivector<void, term<element<0b111>, monomial<one>>>;
    using e013 = multivector<void, term<element<0b1011>, monomial<one>>>;
    using e023 = multivector<void, term<element<0b1101>, monomial<one>>>;
    using e123 = multivector<void, term<element<0b1110>, monomial<one>>>;
    using e0123 = multivector<void, term<element<0b1111>, monomial<one>>>;

    template <int D, int X, int Y, int Z>
    using plane_t = multivector<void,
                               term<element<1>, monomial<rational<D>>>,
                               term<element<0b10>, monomial<rational<X>>>,
                               term<element<0b100>, monomial<rational<Y>>>,
                               term<element<0b1000>, monomial<rational<Z>>>>;

    template <typename T = float>
    struct alignas(16) plane
    {
        using value_t = T;

        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b1>, monomial<one, generator<tag<ID, 0>>>>,
                                 term<element<0b10>, monomial<one, generator<tag<ID, 1>>>>,
                                 term<element<0b100>, monomial<one, generator<tag<ID, 2>>>>,
                                 term<element<0b1000>, monomial<one, generator<tag<ID, 3>>>>>;
        T d;
        T x;
        T y;
        T z;

        [[nodiscard]] constexpr const T& operator[](size_t index) const noexcept
        {
            return *(reinterpret_cast<const T*>(this) + index);
        }

        [[nodiscard]] constexpr T& operator[](size_t index) noexcept
        {
            return *(reinterpret_cast<T*>(this) + index);
        }

        template <typename Engine, typename... I>
        [[nodiscard]] constexpr static plane<T> convert(Engine& engine, multivector<void, I...> mv) noexcept
        {
            auto d_e = extract<0b1>(mv);
            auto x_e = extract<0b10>(mv);
            auto y_e = extract<0b100>(mv);
            auto z_e = extract<0b1000>(mv);

            auto&& [d, x, y, z] = engine.template evaluate_terms<T>(d_e, x_e, y_e, z_e);

            return {d, x, y, z};
        }
    };

    template <int X, int Y, int Z>
    using point_t = multivector<void,
                               term<element<0b111>, monomial<rational<-Z>>>,
                               term<element<0b1011>, monomial<rational<Y>>>,
                               term<element<0b1101>, monomial<rational<-X>>>,
                               term<element<0b1110>, monomial<one>>>;


    template <typename T = float>
    struct alignas(16) point
    {
        using value_t = T;

        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b111>, monomial<minus_one, generator<tag<ID, 2>>>>,  // -z
                                 term<element<0b1011>, monomial<one, generator<tag<ID, 1>>>>,       // y
                                 term<element<0b1101>, monomial<minus_one, generator<tag<ID, 0>>>>, // -x
                                 term<element<0b1110>, monomial<one>>>;

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
        [[nodiscard]] constexpr static point<T> convert(Engine& engine, multivector<void, I...> mv) noexcept
        {
            auto x_e = -extract<0b1101>(mv);
            auto y_e = extract<0b1011>(mv);
            auto z_e = -extract<0b111>(mv);
            auto c_e = extract<0b1110>(mv);

            auto&& [x, y, z, c] = engine.template evaluate_terms<T>(x_e, y_e, z_e, c_e);

            return {x / c, y / c, z / c};
        }
    };

    // TODO: directions
} // namespace pga
} // namespace gal