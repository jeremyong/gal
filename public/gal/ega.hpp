#pragma once

#include "ga.hpp"

// The 3D Euclidean Geometric Algebra

namespace gal
{
namespace ega
{
    using metric = ::gal::metric<3, 0, 0>;

    using algebra = ga::algebra<metric>;

    GAL_OPERATORS(algebra);

    inline multivector<void, term<element<0>, monomial<one>>> e;
    inline multivector<void, term<element<0b1>, monomial<one>>> e0;
    inline multivector<void, term<element<0b10>, monomial<one>>> e1;
    inline multivector<void, term<element<0b100>, monomial<one>>> e2;
    inline multivector<void, term<element<0b11>, monomial<one>>> e01;
    inline multivector<void, term<element<0b101>, monomial<one>>> e02;
    inline multivector<void, term<element<0b110>, monomial<one>>> e12;
    inline multivector<void, term<element<0b111>, monomial<one>>> e012;

    template <int X, int Y, int Z>
    using point_t = multivector<void,
                                term<element<0b1>, monomial<rational<X>>>,
                                term<element<0b10>, monomial<rational<Y>>>,
                                term<element<0b100>, monomial<rational<Z>>>>;

    template <typename T = float>
    struct alignas(16) point : public entity<T, point<T>, 0b1, 0b10, 0b100>
    {
        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b1>, monomial<one, generator<tag<ID, 0>>>>,
                                 term<element<0b10>, monomial<one, generator<tag<ID, 1>>>>,
                                 term<element<0b100>, monomial<one, generator<tag<ID, 2>>>>>;

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 3;
        }

        point(T x, T y, T z)
            : x{x}
            , y{y}
            , z{z}
        {}

        template <typename T1, size_t... E>
        point(entity<T1, void, E...> const& other)
        {
            x = static_cast<T>(other.template get_by_element<0b1>());
            y = static_cast<T>(other.template get_by_element<0b10>());
            z = static_cast<T>(other.template get_by_element<0b100>());
        }

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
    };

    template <int X, int Y, int Z>
    using rotor_t = multivector<void,
                                term<element<0>, monomial<one, generator<tag<0, 0>>>>,
                                term<element<0b11>, monomial<rational<Z>, generator<tag<0, 1>>>>,   // +z
                                term<element<0b101>, monomial<rational<-Y>, generator<tag<0, 1>>>>, // +y
                                term<element<0b110>, monomial<rational<X>, generator<tag<0, 1>>>>>; // +x

    // A rotor, when conjugated with any geometric entity, rotates it theta radians about its axis
    template <typename T = float, typename AngleT = T>
    struct rotor : public entity<T, rotor<T, AngleT>, 0, 0b11, 0b101, 0b110>
    {
        using angle_t = AngleT;

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

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 4;
        }

        rotor(AngleT theta, T x, T y, T z)
            : theta{theta}
            , x{x}
            , y{y}
            , z{z}
        {}

        AngleT theta;
        T x;
        T y;
        T z;

        // As always when doing any normalization operation, NaNs are producible when normalizing
        // vectors of zero length. This is not checked for!
        void normalize() noexcept
        {
            auto l2_inv = T{1} / std::sqrt(x * x + y * y + z * z);
            x           = x * l2_inv;
            y           = y * l2_inv;
            z           = z * l2_inv;
        }

        template <size_t I>
        [[nodiscard]] T get(std::integral_constant<size_t, I>) const noexcept
        {
            if constexpr (I == 4)
            {
                return std::cos(theta * 0.5);
            }
            else if constexpr (I == 5)
            {
                return std::sin(theta * 0.5);
            }
            // unreachable
        }
    };
} // namespace ega
} // namespace gal