#pragma once

#include "entity.hpp"
#include "geometric_algebra.hpp"

namespace gal
{
namespace cga2
{
    // The "Compass Ruler Algebra"

    // The metric is defined here as the standard Minkowski spacetime. To extract the conformal representations,
    // a change of basis is required where o = 1/2 * (e + e-) and inf = e- - e.
    // The element e0 here is the added unit norm basis element and the elements e1 and e2 correspond to
    // the canonical Euclidean R2 basis representation.
    using cra_metric = gal::metric<3, 1, 0>;

    // The CRA is a graded algebra with 16 basis elements
    using cra_algebra = gal::algebra<cra_metric>;

    // 0b1 => e+ extension
    // 0b1000 => e- extension
    namespace detail
    {
        // These tags are needed to provide unique specializations for the expressions for e_o and e_inf
        template <typename T>
        struct e_o_tag
        {};

        template <typename T>
        struct e_inf_tag
        {};
    } // namespace detail
} // namespace cga2


template <typename T>
struct expr<expr_op::identity, mv<cga2::cra_algebra, 0, 2, 2>, cga2::detail::e_o_tag<T>>
{
    using value_t               = T;
    using algebra_t             = cga2::cra_algebra;
    constexpr static expr_op op = expr_op::identity;
    constexpr static mv<cra_algebra, 0, 2, 2> lhs{mv_size{0, 2, 2},
                                                  {},
                                                  {mon{one_half, 0, 0}, mon{one_half, 0, 0}},
                                                  {term{1, 0, 0b1}, term{1, 1, 0b1000}}};
};

template <typename T>
struct expr<expr_op::identity, mv<cga2::cra_algebra, 0, 2, 2>, cga2::detail::e_inf_tag<T>>
{
    using value_t               = T;
    using algebra_t             = cga2::cra_algebra;
    constexpr static expr_op op = expr_op::identity;
    constexpr static mv<cra_algebra, 0, 2, 2> lhs{mv_size{0, 2, 2},
                                                  {},
                                                  {mon{minus_one, 0, 0}, mon{one, 0, 0}},
                                                  {term{1, 0, 0b1}, term{1, 1, 0b1000}}};
};

namespace cga2
{
    template <typename T = float>
    constexpr inline expr<expr_op::identity, mv<cra_algebra, 0, 2, 2>, detail::e_o_tag<T>> e_o;

    template <typename T = float>
    constexpr inline expr<expr_op::identity, mv<cra_algebra, 0, 2, 2>, detail::e_inf_tag<T>> e_inf;

    template <typename T = float>
    union point
    {
        using algebra_t = cga_algebra;
        using value_t   = T;

        std::array<T, 2> data;
        struct
        {
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
        };

        constexpr point(T a, T b) noexcept
            : x{a}
            , y{b}
        {}

        template <uint8_t... E>
        constexpr point(entity<algebra_t, T, E...> in) noexcept
            : data{in.template select<0b10, 0b100>()}
        {}

        [[nodiscard]] constexpr static mv<algebra_t, 4, 6, 4> ie(uint32_t id) noexcept
        {
            return {mv_size{4, 6, 4},
                    {
                        ind{3, 1}, // ind3 = p^2
                        ind{0, 1}, // ind0 = p_x
                        ind{1, 1}, // ind1 = p_y
                        ind{3, 1}  // ind3 = p^2
                    },
                    {
                        mon{one_half, 0, 0, 0},       // 1/2
                        mon{minus_one_half, 1, 0, 1}, // -1/2 p^2
                        mon{one, 1, 1, 1},            // p_x
                        mon{one, 1, 2, 1},            // p_y
                        mon{one_half, 0, 0, 0},       // 1/2
                        mon{one_half, 1, 3, 1}        // 1/2 p^2
                    },
                    {
                        term{2, 0, 0b1},   // (1/2 - 1/2 p^2) e
                        term{1, 2, 0b10},  // p_x
                        term{1, 3, 0b100}, // p_y
                        term{2, 4, 0b1000} // (1/2 + 1/2 p^2) e-
                    }};
        }

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 2;
        }

        [[nodiscard]] constexpr static uint32_t ind_count() noexcept
        {
            return 3;
        }

        [[nodiscard]] constexpr T const& operator[](size_t index) const noexcept
        {
            return data[index];
        }

        [[nodiscard]] constexpr T& operator[](size_t index) noexcept
        {
            return data[index];
        }

        [[nodiscard]] constexpr T get(size_t i) const noexcept
        {
            if (i == 3)
            {
                return x * x + y * y;
            }
            return {};
        }
    };
    // TODO: provide representations for planes, spheres, flats, etc.
} // namespace cga2
} // namespace gal