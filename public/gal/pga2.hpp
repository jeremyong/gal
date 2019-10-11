#pragma once

#include "finite_algebra.hpp"
#include "ga.hpp"

#include <cmath>

// Projective Geometric Algebra for representing Euclidean 2-space

namespace gal
{
namespace pga2
{
    // NOTE: the inner product of e0 can be set to +1 or -1 without any change in the algebra's geometric
    // interpretation. Here, we opt to define e0^2 := 1 by convention
    using metric = gal::metric<2, 0, 1>;

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

    // Point at (x, y)
    template <int X, int Y>
    using point_t = multivector<void,
                                term<element<0b11>, monomial<rational<Y>>>,
                                term<element<0b101>, monomial<rational<-X>>>,
                                term<element<0b110>, monomial<one>>>;

    template <typename T = float>
    struct point : public entity<T, point<T>, 0b11, 0b101, 0b110>
    {
        using value_t = T;

        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b11>, monomial<one, generator<tag<ID, 1>>>>,        // y
                                 term<element<0b101>, monomial<minus_one, generator<tag<ID, 0>>>>, // -x
                                 term<element<0b110>, monomial<one>>>;

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 2;
        }

        point(T x, T y)
            : x{x}
            , y{y}
        {}

        // WARNING: This implicit conversion from an entity does not check if the weight is 0
        template <typename T1, size_t... E>
        point(entity<T1, void, E...> const& other)
        {
            auto c_inv = T{1} / static_cast<T>(other.template get_by_element<0b110>);
            x          = -static_cast<T>(other.template get_by_element<0b101>()) * c_inv;
            y          = static_cast<T>(other.template get_by_element<0b11>()) * c_inv;
        }

        union
        {
            T x;
            T u;
            T s;
        };

        union
        {
            T y;
            T v;
            T t;
        };
    };

    // Line with equation ax + by + c = 0
    template <int A, int B, int C>
    using line_t = multivector<void,
                               term<element<0b1>, monomial<rational<C>>>,
                               term<element<0b10>, monomial<rational<A>>>,
                               term<element<0b100>, monomial<rational<B>>>>;

    template <typename T = float>
    struct line : public entity<T, line<T>, 0b1, 0b10, 0b100>
    {
        using value_t = T;

        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b1>, monomial<one, generator<tag<ID, 2>>>>,
                                 term<element<0b10>, monomial<one, generator<tag<ID, 0>>>>,
                                 term<element<0b100>, monomial<one, generator<tag<ID, 1>>>>>;

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 3;
        }

        line(T a, T b, T c)
            : a{a}
            , b{b}
            , c{c}
        {}

        template <typename T1, size_t... E>
        line(entity<T1, void, E...> const& other)
        {
            a = static_cast<T>(other.template get_by_element<0b10>());
            b = static_cast<T>(other.template get_by_element<0b100>());
            c = static_cast<T>(other.template get_by_element<0b1>());
        }

        T a;
        T b;
        T c;
    };

    // Direction pointing toward (x, y) aka an ideal point
    template <int X, int Y>
    using direction_t
        = multivector<void, term<element<0b11>, monomial<rational<Y>>>, term<element<0b101>, monomial<rational<-X>>>>;

    template <typename T>
    struct alignas(16) direction : public entity<T, direction<T>, 0b11, 0b101>
    {
        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b11>, monomial<minus_one, generator<tag<ID, 1>>>>, // y
                                 term<element<0b101>, monomial<one, generator<tag<ID, 0>>>>>;     // -x

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 2;
        }

        direction(T x, T y)
            : x{x}
            , y{y}
        {}

        template <typename T1, size_t... E>
        direction(entity<T1, void, E...> const& other)
        {
            x = -static_cast<T>(other.template get_by_element<0b101>());
            y = static_cast<T>(other.template get_by_element<0b11>());
        }

        T x;
        T y;
    };

    // TODO add rotor and translator

    // Additional helper functions for plane geometry and result verification
    template <typename... T>
    [[nodiscard]] constexpr auto line_slope(multivector<void, T...> line) noexcept
    {
        constexpr auto x = extract<0b10>(line);
        constexpr auto y = extract<0b100>(line);
        // For now this only works on compile time multivectors over a finite field
        static_assert(decltype(x)::size == 1);
        static_assert(decltype(y)::size == 1);
        return -typename decltype(x)::first_t::rational_t{} / typename decltype(y)::first_t::rational_t{};
    }

    template <typename... T>
    [[nodiscard]] constexpr auto cartesian_point(multivector<void, T...> point) noexcept
    {
        constexpr auto x    = -extract<0b101>(point);
        constexpr auto y    = extract<0b11>(point);
        constexpr auto sign = extract<0b110>(point);

        // For now this only works on compile time multivectors over a finite field
        static_assert(decltype(x)::size == 1);
        static_assert(decltype(y)::size == 1);
        static_assert(decltype(sign)::size == 1);

        return std::make_pair(
            typename decltype(x)::first_t::rational_t{} / typename decltype(sign)::first_t::rational_t{},
            typename decltype(y)::first_t::rational_t{} / typename decltype(sign)::first_t::rational_t{});
    }
} // namespace pga2
} // namespace gal