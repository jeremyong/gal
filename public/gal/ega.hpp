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

    template <int X, int Y, int Z>
    using point_t = multivector<void,
                                term<element<0b1>, monomial<rational<X>>>,
                                term<element<0b10>, monomial<rational<Y>>>,
                                term<element<0b100>, monomial<rational<Z>>>>;

    template <typename T = float>
    struct alignas(16) point
    {
        constexpr static size_t size = 3;

        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b1>, monomial<one, generator<tag<ID, 0>>>>,
                                 term<element<0b10>, monomial<one, generator<tag<ID, 1>>>>,
                                 term<element<0b100>, monomial<one, generator<tag<ID, 2>>>>>;
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

        GAL_ACCESSORS

        template <typename Engine, typename... I>
        [[nodiscard]] constexpr static point<T> convert(Engine& engine, multivector<void, I...> mv) noexcept
        {
            auto x_e = extract<0b1>(mv);
            auto y_e = extract<0b10>(mv);
            auto z_e = extract<0b100>(mv);

            auto&& [x, y, z] = engine.template evaluate_terms<T>(x_e, y_e, z_e);

            return {x, y, z};
        }
    };

    template <int X, int Y, int Z>
    using rotor_t = multivector<void,
                                term<element<0>, monomial<one, generator<tag<0, 0>>>>,
                                term<element<0b11>, monomial<rational<Z>, generator<tag<0, 1>>>>,   // +z
                                term<element<0b101>, monomial<rational<-Y>, generator<tag<0, 1>>>>, // +y
                                term<element<0b110>, monomial<rational<X>, generator<tag<0, 1>>>>>; // +x

    // A rotor, when conjugated with any geometric entity, rotates it theta radians about its axis
    template <typename T = float, typename AngleT = T>
    struct rotor
    {
        using angle_t = AngleT;
        using value_t = T;
        constexpr static size_t size = 4;

        // theta := ID 0, index 0
        // x := ID 0, index 1
        // y := ID 0, index 2
        // z := ID 0, index 3
        // Cos\theta := ID 0, index 4
        // Sin\theta := ID 0, index 5
        template <size_t ID>
        using type
            = multivector<void,
                          term<element<0>, monomial<one, generator<tag<ID, 4>>>>,
                          term<element<0b11>, monomial<one, generator<tag<ID, 3>>, generator<tag<ID, 5>>>>,        // +z
                          term<element<0b101>, monomial<minus_one, generator<tag<ID, 2>>, generator<tag<ID, 5>>>>, // -y
                          term<element<0b110>, monomial<one, generator<tag<ID, 1>>, generator<tag<ID, 5>>>>>;      // +x

        AngleT theta;
        T x;
        T y;
        T z;

        // As always when doing any normalization operation, NaNs are producible when normalizing
        // vectors of zero length. This is not checked for!
        void normalize() noexcept
        {
            auto l2_inv = T{1} / std::sqrt(x * x + y * y + z * z);
            x = x * l2_inv;
            y = y * l2_inv;
            z = z * l2_inv;
        }

        GAL_ACCESSORS

        template <size_t I>
        [[nodiscard]] constexpr T get_special(std::integral_constant<size_t, I>) const noexcept
        {
            if constexpr (I == 4)
            {
                return std::cos(theta * 0.5);
            }
            else if constexpr (I == 5)
            {
                return std::sin(theta * 0.5);
            }
            else
            {
                static_assert(I == 4 || I == 5, "Unreachable branch");
            }
        }
    };
} // namespace ega
} // namespace gal