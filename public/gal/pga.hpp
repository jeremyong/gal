#pragma once

#include "ga.hpp"

#include <cmath>

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

    GAL_OPERATORS(algebra);

    inline multivector<void, term<element<0>, monomial<one>>> e;
    inline multivector<void, term<element<0b1>, monomial<one>>> e0;
    inline multivector<void, term<element<0b10>, monomial<one>>> e1;
    inline multivector<void, term<element<0b100>, monomial<one>>> e2;
    inline multivector<void, term<element<0b1000>, monomial<one>>> e3;
    inline multivector<void, term<element<0b11>, monomial<one>>> e01;
    inline multivector<void, term<element<0b101>, monomial<one>>> e02;
    inline multivector<void, term<element<0b1001>, monomial<one>>> e03;
    inline multivector<void, term<element<0b110>, monomial<one>>> e12;
    inline multivector<void, term<element<0b1010>, monomial<one>>> e13;
    inline multivector<void, term<element<0b1100>, monomial<one>>> e23;
    inline multivector<void, term<element<0b111>, monomial<one>>> e012;
    inline multivector<void, term<element<0b1011>, monomial<one>>> e013;
    inline multivector<void, term<element<0b1101>, monomial<one>>> e023;
    inline multivector<void, term<element<0b1110>, monomial<one>>> e123;
    inline multivector<void, term<element<0b1111>, monomial<one>>> e0123;

    template <int D, int X, int Y, int Z>
    using plane_t = multivector<void,
                                term<element<1>, monomial<rational<D>>>,
                                term<element<0b10>, monomial<rational<X>>>,
                                term<element<0b100>, monomial<rational<Y>>>,
                                term<element<0b1000>, monomial<rational<Z>>>>;

    template <typename T = float>
    struct alignas(16) plane : public entity<T, plane<T>, 0b1, 0b10, 0b100, 0b1000>
    {
        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 4;
        }

        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b1>, monomial<one, generator<tag<ID, 0>>>>,
                                 term<element<0b10>, monomial<one, generator<tag<ID, 1>>>>,
                                 term<element<0b100>, monomial<one, generator<tag<ID, 2>>>>,
                                 term<element<0b1000>, monomial<one, generator<tag<ID, 3>>>>>;

        plane(T d, T x, T y, T z)
            : d{d}
            , x{x}
            , y{y}
            , z{z}
        {}

        template <typename T1, size_t... E>
        plane(entity<T1, void, E...> const& other)
        {
            d = static_cast<T>(other.template get_by_element<0b1>());
            x = static_cast<T>(other.template get_by_element<0b10>());
            y = static_cast<T>(other.template get_by_element<0b100>());
            z = static_cast<T>(other.template get_by_element<0b1000>());
        }

        T d;
        T x;
        T y;
        T z;
    };

    template <int X, int Y, int Z>
    using point_t = multivector<void,
                                term<element<0b111>, monomial<rational<-Z>>>,
                                term<element<0b1011>, monomial<rational<Y>>>,
                                term<element<0b1101>, monomial<rational<-X>>>,
                                term<element<0b1110>, monomial<one>>>;

    template <typename T = float>
    struct alignas(16) point : public entity<T, point<T>, 0b111, 0b1011, 0b1101, 0b1110>
    {
        using value_t = T;

        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b111>, monomial<minus_one, generator<tag<ID, 2>>>>,  // -z
                                 term<element<0b1011>, monomial<one, generator<tag<ID, 1>>>>,       // y
                                 term<element<0b1101>, monomial<minus_one, generator<tag<ID, 0>>>>, // -x
                                 term<element<0b1110>, monomial<one>>>;

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 3;
        }

        point(T x, T y, T z)
            : x{x}
            , y{y}
            , z{z}
        {}

        // WARNING: This implicit conversion from an entity does not check if the weight is 0
        template <typename T1, size_t... E>
        point(entity<T1, void, E...> const& other)
        {
            auto c_inv = T{1} / static_cast<T>(other.template get_by_element<0b1110>);
            x          = -static_cast<T>(other.template get_by_element<0b1101>()) * c_inv;
            y          = static_cast<T>(other.template get_by_element<0b1011>()) * c_inv;
            z          = -static_cast<T>(other.template get_by_element<0b111>()) * c_inv;
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
    using direction_t = multivector<void,
                                    term<element<0b111>, monomial<rational<-Z>>>,
                                    term<element<0b1011>, monomial<rational<Y>>>,
                                    term<element<0b1101>, monomial<rational<-X>>>>;

    template <typename T = float>
    struct alignas(16) direction : public entity<T, direction<T>, 0b111, 0b1011, 0b1101>
    {
        using value_t                = T;
        constexpr static size_t size = 3;

        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b111>, monomial<minus_one, generator<tag<ID, 2>>>>,   // -z
                                 term<element<0b1011>, monomial<one, generator<tag<ID, 1>>>>,        // y
                                 term<element<0b1101>, monomial<minus_one, generator<tag<ID, 0>>>>>; // -x

        direction(T x, T y, T z)
            : x{x}
            , y{y}
            , z{z}
        {}

        template <typename T1, size_t... E>
        direction(entity<T1, void, E...> const& other)
        {
            x = -static_cast<T>(other.template get_by_element<0b1101>());
            y = static_cast<T>(other.template get_by_element<0b1011>());
            z = -static_cast<T>(other.template get_by_element<0b111>());
        }

        T x;
        T y;
        T z;
    };

    // Lines in P^3 are defined using Plücker coordinates: https://en.wikipedia.org/wiki/Plücker_coordinates
    // The lines e01, e02, and e03 are the idea lines representing the intersections of e1, e2, and e3 with the ideal
    // plane respectively. The lines e23, e31, and e12 are lines through the origin in the x, y, and z directions
    // respectively. We opt not to provide a 6 coordinate representation for now (join points to construct a line or
    // meet planes)

    template <int X, int Y, int Z>
    using rotor_t = multivector<void,
                                term<element<0>, monomial<one, generator<tag<0, 0>>>>,
                                term<element<0b110>, monomial<rational<Z>, generator<tag<0, 1>>>>,   // +z
                                term<element<0b1010>, monomial<rational<-Y>, generator<tag<0, 1>>>>, // -y
                                term<element<0b1100>, monomial<rational<X>, generator<tag<0, 1>>>>>; // +x

    // A rotor, when conjugated with any geometric entity, rotates it theta radians about its axis
    template <typename T = float, typename AngleT = T>
    struct rotor : public entity<T, rotor<T, AngleT>, 0, 0b110, 0b1010, 0b1100>
    {
        using angle_t = AngleT;
        using value_t = T;
        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 4;
        }

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
                          term<element<0b110>, monomial<one, generator<tag<ID, 3>>, generator<tag<ID, 5>>>>, // +z
                          term<element<0b1010>, monomial<minus_one, generator<tag<ID, 2>>, generator<tag<ID, 5>>>>, // -y
                          term<element<0b1100>, monomial<one, generator<tag<ID, 1>>, generator<tag<ID, 5>>>>>; // +x

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
        [[nodiscard]] constexpr T get(std::integral_constant<size_t, I>) const noexcept
        {
            if constexpr (I == 4)
            {
                return std::cos(theta * 0.5);
            }
            else if constexpr (I == 5)
            {
                return std::sin(theta * 0.5);
            }
        }
    };

    // A translator, when conjugated with any geometric entity, translates it along its axis by the specified distance
    template <typename T = float, typename DistanceT = T>
    struct translator : public entity<T, translator<T, DistanceT>, 0, 0b11, 0b101, 0b1001>
    {
        using distance_t = DistanceT;
        using value_t    = T;
        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 4;
        }

        template <size_t ID>
        using type
            = multivector<void,
                          term<element<0>, monomial<one>>,
                          term<element<0b11>, monomial<minus_one_half, generator<tag<ID, 0>>, generator<tag<ID, 1>>>>, // -x
                          term<element<0b101>, monomial<minus_one_half, generator<tag<ID, 0>>, generator<tag<ID, 2>>>>, // -y
                          term<element<0b1001>, monomial<minus_one_half, generator<tag<ID, 0>>, generator<tag<ID, 3>>>>>; // -z

        translator(DistanceT distance, T x, T y, T z)
            : distance{distance}
            , x{x}
            , y{y}
            , z{z}
        {}

        DistanceT distance;
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
    };
} // namespace pga
} // namespace gal