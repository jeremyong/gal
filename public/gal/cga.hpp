#pragma once

#include "ga.hpp"

namespace gal
{
namespace cga
{
    // The metric is defined here as the standard Minkowski spacetime. To extract the conformal representations,
    // a change of basis is required where o = 1/2 * (e3 + e4) and inf = e4 - e3.

    // After the change of basis, we index the elements as follows:
    // e0 := x
    // e1 := y
    // e2 := z
    // e3 := origin
    // e4 := infinity

    struct cga_metric
    {
        using base_metric = metric<4, 1, 0>;
        constexpr static bool is_diagonal = false;
        constexpr static size_t p = 4;
        constexpr static size_t v = 1;
        constexpr static size_t r = 0;
        constexpr static size_t dimension = p + v + r;

        [[nodiscard]] constexpr static bool multi_term_gp(size_t lhs, size_t rhs) noexcept
        {
            bool lhs_off_diag = lhs & 0b11000;
            bool rhs_off_diag = rhs & 0b11000;
            return lhs_off_diag && rhs_off_diag;
        }

        template <typename E, typename... M>
        [[nodiscard]] constexpr static auto diagonalize(term<E, M...>) noexcept
        {
            if constexpr (E::value & 0b11000)
            {
                if constexpr ((E::value & 0b10000) && (E::value & 0b1000))
                {
                    // e_o ^ e_inf = -1 + e_3 ^ e_4
                    return multivector<void, term<element<E::value & 0b111>, monomial<minus_one>>, term<E, M...>>{};
                }
                else if constexpr (E::value & 0b10000)
                {
                    // e_inf = e_4 - e_3
                    constexpr size_t e = E::value ^ 0b11000;
                    return multivector<void, decltype(minus_one{} * term<element<e>, M...>{}), term<E, M...>>{};
                }
                else
                {
                    // e_0 = 1/2 * (e_3 + e_4)
                    constexpr size_t e = E::value ^ 0b11000;
                    return multivector<void,
                                       decltype(one_half{} * term<E, M...>{}),
                                       decltype(one_half{} * term<element<e>, M...>{})>{};
                }
            }
            else
            {
                return multivector<void, term<E, M...>>{};
            }
        }

        template <typename... T>
        [[nodiscard]] constexpr static auto undiagonalize(multivector<void, T...>) noexcept
        {
            if constexpr (sizeof...(T) == 0)
            {
                return multivector<void>{};
            }
            else
            {
                return (undiagonalize(T{}) + ...);
            }
        }

        template <typename E, typename... M>
        [[nodiscard]] constexpr static auto undiagonalize(term<E, M...>) noexcept
        {
            if constexpr (E::value & 0b11000)
            {
                if constexpr ((E::value & 0b10000) && (E::value & 0b1000))
                {
                    // e_3 ^ e_4 = e_o ^ e_inf
                    return multivector<void, term<E, M...>>{};
                }
                else if constexpr (E::value & 0b10000)
                {
                    // e_4 = e_o + 1/2 * e_inf
                    constexpr size_t e = E::value ^ 0b11000;
                    return multivector<void, term<element<e>, M...>, decltype(one_half{} * term<E, M...>{})>{};
                }
                else
                {
                    // e_3 = e_o - 1/2 * e_inf
                    constexpr size_t e = E::value ^ 0b11000;
                    return multivector<void, term<E, M...>, decltype(minus_one_half{} * term<element<e>, M...>{})>{};
                }
            }
            else
            {
                return multivector<void, term<E, M...>>{};
            }
        }

        [[nodiscard]] constexpr static std::pair<int, int> intercept(size_t e, size_t blade) noexcept
        {
            // This is the implementation for a diagonal Cayley table
            if (e == 3)
            {
                if (blade & (1 << 4))
                {
                    return {4, -1};
                }
                else
                {
                    return {-1, 0};
                }
            }
            else if (e == 4)
            {
                if (blade & (1 << 3))
                {
                    return {3, -1};
                }
                else
                {
                    return {-1, 0};
                }
            }
            else if (blade & (1 << e))
            {
                return {e, 1};
            }
            else
            {
                return {-1, 0};
            }
        }
    };

    // The CGA is a graded algebra with 32 basis elements
    using cga_algebra = ga::algebra<cga_metric>;

    inline multivector<void, term<element<0>, monomial<one>>> e;
    inline multivector<void, term<element<0b1>, monomial<one>>> e_x;
    inline multivector<void, term<element<0b10>, monomial<one>>> e_y;
    inline multivector<void, term<element<0b100>, monomial<one>>> e_z;
    inline multivector<void, term<element<0b1000>, monomial<one>>> e_o;
    inline multivector<void, term<element<0b10000>, monomial<one>>> e_inf;

    GAL_OPERATORS(cga_algebra);

    template <int X, int Y, int Z>
    using point_t = multivector<void,
        term<element<0b1>, monomial<rational<X>>>,
        term<element<0b10>, monomial<rational<Y>>>,
        term<element<0b100>, monomial<rational<Z>>>,
        term<element<0b1000>, monomial<one>>,
        term<element<0b10000>, monomial<rational<X * X + Y * Y + Z * Z, 2>>>>;

    template <typename T = float>
    struct alignas(16) point
    {
        using value_t = T;
        constexpr static size_t size = 3;

        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b1>, monomial<one, generator<tag<ID, 0>>>>,
                                 term<element<0b10>, monomial<one, generator<tag<ID, 1>>>>,
                                 term<element<0b100>, monomial<one, generator<tag<ID, 2>>>>,
                                 term<element<0b1000>, monomial<one>>,
                                 term<element<0b10000>,
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

        GAL_ACCESSORS

        template <typename Engine, typename... I>
        [[nodiscard]] constexpr static point<T> convert(const Engine& engine, multivector<void, I...> mv) noexcept
        {
            auto x_e = extract<0b1>(mv);
            auto y_e = extract<0b10>(mv);
            auto z_e = extract<0b100>(mv);
            auto o_e = extract<0b1000>(mv);

            auto&& [o, x, y, z] = engine.template evaluate_terms<T>(o_e, x_e, y_e, z_e);
            return {x / o, y / o, z / o};
        }
    };

    // TODO: provide representations for planes, spheres, flats, etc.
} // namespace cga
} // namespace gal