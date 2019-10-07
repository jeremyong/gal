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
        constexpr static size_t size = 4;

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

        GAL_ACCESSORS

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
        constexpr static size_t size = 3;

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

        GAL_ACCESSORS

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

    template <int X, int Y, int Z>
    using direction_t = multivector<void,
                                    term<element<0b111>, monomial<rational<-Z>>>,
                                    term<element<0b1011>, monomial<rational<Y>>>,
                                    term<element<0b1101>, monomial<rational<-X>>>>;

    template <typename T>
    struct alignas(16) direction
    {
        using value_t = T;
        constexpr static size_t size = 3;

        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b111>, monomial<minus_one, generator<tag<ID, 2>>>>,   // -z
                                 term<element<0b1011>, monomial<one, generator<tag<ID, 1>>>>,        // y
                                 term<element<0b1101>, monomial<minus_one, generator<tag<ID, 0>>>>>; // -x

        T x;
        T y;
        T z;

        GAL_ACCESSORS

        template <typename Engine, typename... I>
        [[nodiscard]] constexpr static direction<T> convert(Engine& engine, multivector<void, I...> mv) noexcept
        {
            auto x_e = -extract<0b1101>(mv);
            auto y_e = extract<0b1011>(mv);
            auto z_e = -extract<0b111>(mv);

            auto&& [x, y, z] = engine.template evaluate_terms<T>(x_e, y_e, z_e);

            return {x, y, z};
        }
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
                          term<element<0b110>, monomial<one, generator<tag<ID, 3>>, generator<tag<ID, 5>>>>, // +z
                          term<element<0b1010>, monomial<minus_one, generator<tag<ID, 2>>, generator<tag<ID, 5>>>>, // -y
                          term<element<0b1100>, monomial<one, generator<tag<ID, 1>>, generator<tag<ID, 5>>>>>; // +x

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

    // A translator, when conjugated with any geometric entity, translates it along its axis by the specified distance
    template <typename T = float, typename DistanceT = T>
    struct translator
    {
        using distance_t = DistanceT;
        using value_t = T;
        constexpr static size_t size = 4;

        template <size_t ID>
        using type
            = multivector<void,
                          term<element<0>, monomial<one>>,
                          term<element<0b11>, monomial<minus_one_half, generator<tag<ID, 0>>, generator<tag<ID, 1>>>>, // -x
                          term<element<0b101>, monomial<minus_one_half, generator<tag<ID, 0>>, generator<tag<ID, 2>>>>, // -y
                          term<element<0b1001>, monomial<minus_one_half, generator<tag<ID, 0>>, generator<tag<ID, 3>>>>>; // -z

        DistanceT distance;
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
    };
} // namespace pga
} // namespace gal