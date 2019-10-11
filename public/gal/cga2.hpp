#pragma once

#include "ga.hpp"

namespace gal
{
namespace cga2
{
    // The "Compass Ruler Algebra"

    // The metric is defined here as the standard Minkowski spacetime. To extract the conformal representations,
    // a change of basis is required where o = 1/2 * (e2 + e3) and inf = e3 - e2.
    // The element e3 here is the added unit norm basis element and the elements e0 and e2 correspond to
    // the canonical Euclidean R3 basis representation.
    struct cga_metric
    {
        using base_metric                 = metric<3, 1, 0>;
        constexpr static bool is_diagonal = false;
        constexpr static size_t p         = 3;
        constexpr static size_t v         = 1;
        constexpr static size_t r         = 0;
        constexpr static size_t dimension = p + v + r;

        [[nodiscard]] constexpr static bool multi_term_gp(size_t lhs, size_t rhs) noexcept
        {
            bool lhs_off_diag = lhs & 0b1100;
            bool rhs_off_diag = rhs & 0b1100;
            return lhs_off_diag && rhs_off_diag;
        }

        template <typename E, typename... M>
        [[nodiscard]] constexpr static auto diagonalize(term<E, M...>) noexcept
        {
            if constexpr ((E::value & 0b1100) > 0)
            {
                if constexpr (((E::value & 0b1000) > 0) && ((E::value & 0b100) > 0))
                {
                    // e_o ^ e_inf = -1 + e_3 ^ e_4
                    return multivector<void, term<element<E::value & 0b11>, monomial<minus_one>>, term<E, M...>>{};
                }
                else if constexpr ((E::value & 0b1000) > 0)
                {
                    // e_inf = e_4 - e_3
                    constexpr size_t e = E::value ^ 0b1100;
                    return multivector<void, decltype(minus_one{} * term<element<e>, M...>{}), term<E, M...>>{};
                }
                else
                {
                    // e_0 = 1/2 * (e_3 + e_4)
                    constexpr size_t e = E::value ^ 0b1100;
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
            if constexpr ((E::value & 0b1100) > 0)
            {
                if constexpr (((E::value & 0b1000) > 0) && ((E::value & 0b100) > 0))
                {
                    // e_3 ^ e_4 = e_o ^ e_inf
                    return multivector<void, term<E, M...>>{};
                }
                else if constexpr (E::value & 0b1000 > 0)
                {
                    // e_4 = e_o + 1/2 * e_inf
                    constexpr size_t e = E::value ^ 0b1100;
                    return multivector<void, term<element<e>, M...>, decltype(one_half{} * term<E, M...>{})>{};
                }
                else
                {
                    // e_3 = e_o - 1/2 * e_inf
                    constexpr size_t e = E::value ^ 0b1100;
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
            if (e == 2)
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
            else if (e == 3)
            {
                if (blade & (1 << 2))
                {
                    return {2, -1};
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
    inline multivector<void, term<element<0b100>, monomial<one>>> e_o;
    inline multivector<void, term<element<0b1000>, monomial<one>>> e_inf;

    GAL_OPERATORS(algebra);

    template <int X, int Y>
    using point_t = multivector<void,
                                term<element<0b1>, monomial<rational<X>>>,
                                term<element<0b10>, monomial<rational<Y>>>,
                                term<element<0b100>, monomial<one>>,
                                term<element<0b1000>, monomial<rational<X * X + Y * Y, 2>>>>;

    template <typename T = float>
    struct point : public entity<T, point<T>, 0b1, 0b10>
    {
        template <size_t ID>
        using type = multivector<void,
                                 term<element<0b1>, monomial<one, generator<tag<ID, 0>>>>,
                                 term<element<0b10>, monomial<one, generator<tag<ID, 1>>>>,
                                 term<element<0b100>, monomial<one>>,
                                 term<element<0b1000>,
                                      monomial<one_half, generator<tag<ID, 0>, degree<2>>>,
                                      monomial<one_half, generator<tag<ID, 1>, degree<2>>>>>;

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
            auto weight_inv = T{1} / static_cast<T>(other.template get_by_element<0b100>);
            x               = static_cast<T>(other.template get_by_element<0b1>()) * weight_inv;
            y               = static_cast<T>(other.template get_by_element<0b10>()) * weight_inv;
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

    // TODO: provide representations for planes, spheres, flats, etc.
} // namespace cga2
} // namespace gal