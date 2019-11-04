#pragma once

#include "engine.hpp"
#include "entity.hpp"
#include "geometric_algebra.hpp"

#include <cmath>

namespace gal
{
namespace cga2
{
    // The "Compass Ruler Algebra"

    // The metric is defined here as the standard Minkowski spacetime. To extract the conformal representations,
    // a change of basis is required where o = 1/2 * (e + e-) and inf = e- - e.
    using cga2_metric = gal::metric<3, 1, 0>;

    // The CRA is a graded algebra with 16 basis elements
    using cga2_algebra = gal::algebra<cga2_metric>;

    // 0b1 => e+ extension
    // 0b1000 => e- extension
    namespace detail
    {
        // These tags are needed to provide unique specializations for the expressions for n_o and n_i
        template <typename T>
        struct n_o_tag
        {};

        template <typename T>
        struct n_i_tag
        {};

        template <typename T>
        struct pseudoscalar_tag
        {};

        template <typename T>
        struct pseudoscalar_inv_tag
        {};
    } // namespace detail
} // namespace cga2

namespace detail
{
    template <>
    constexpr inline bool uses_null_basis<::gal::cga::cga2_algebra> = true;
}

template <typename T>
struct expr<expr_op::identity, mv<cga2::cga2_algebra, 0, 1, 1>, cga2::detail::n_o_tag<T>>
{
    using value_t               = T;
    using algebra_t             = cga2::cga2_algebra;
    constexpr static expr_op op = expr_op::identity;
    constexpr static mv<cga2::cga2_algebra, 0, 1, 1> lhs{mv_size{0, 1, 1}, {}, {mon{one, zero, 0, 0}}, {term{1, 0, 0b1000}}};
};

template <typename T>
struct expr<expr_op::identity, mv<cga2::cga2_algebra, 0, 1, 1>, cga2::detail::n_i_tag<T>>
{
    using value_t               = T;
    using algebra_t             = cga2::cga2_algebra;
    constexpr static expr_op op = expr_op::identity;
    constexpr static mv<cga2::cga2_algebra, 0, 1, 1> lhs{mv_size{0, 1, 1}, {}, {mon{one, zero, 0, 0}}, {term{1, 1, 0b1000}}};
};

template <typename T>
struct expr<expr_op::identity, mv<cga2::cga2_algebra, 0, 1, 1>, cga2::detail::pseudoscalar_tag<T>>
{
    using value_t               = T;
    using algebra_t             = cga2::cga2_algebra;
    constexpr static expr_op op = expr_op::identity;
    constexpr static auto lhs   = cga2::cga2_algebra::pseudoscalar;
};

template <typename T>
struct expr<expr_op::identity, mv<cga::cga2_algebra, 0, 1, 1>, cga2::detail::pseudoscalar_inv_tag<T>>
{
    using value_t               = T;
    using algebra_t             = cga2::cga2_algebra;
    constexpr static expr_op op = expr_op::identity;
    constexpr static auto lhs   = cga2::cga2_algebra::pseudoscalar_inv;
};

namespace cga2
{
    template <typename T = float>
    constexpr inline expr<expr_op::identity, mv<cga2_algebra, 0, 1, 1>, detail::n_o_tag<T>> n_o;

    template <typename T = float>
    constexpr inline expr<expr_op::identity, mv<cga2_algebra, 0, 1, 1>, detail::e_inf_tag<T>> n_i;

    template <typename T = float>
    constexpr inline expr<expr_op::identity, mv<cga2_algebra, 0, 1, 1>, detail::pseudoscalar_tag<T>> ps;

    template <typename T = float>
    constexpr inline expr<expr_op::identity, mv<cga2_algebra, 0, 1, 1>, detail::pseudoscalar_inv_tag<T>> ips;

    template <typename T = float>
    union point
    {
        using algebra_t = cga2_algebra;
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

        template <elem_t... E>
        constexpr point(entity<algebra_t, T, E...> in) noexcept
            : data{in.template select<0b1, 0b10>()}
        {}

        [[nodiscard]] constexpr static mv<algebra_t, 4, 5, 4> ie(uint32_t id) noexcept
        {
            return {mv_size{4, 5, 4},
                    {
                        ind{0, one},    // ind0 = p_x
                        ind{1, one},    // ind1 = p_y
                        ind{0, rat{2}}  // ind3 = p_x^2
                        ind{1, rat{2}}, // ind3 = p_y^2
                    },
                    {
                        mon{one, one, 1, 0},         // p_x
                        mon{one, one, 1, 1},         // p_y
                        mon{one, zero, 0, 0, 0},     // n_o
                        mon{one_half, rat{2}, 1, 2}, // 1/2 p_x^2
                        mon{one_half, rat{2}, 1, 3}  // 1/2 p_y^2
                    },
                    {
                        term{1, 0, 0b1},   // p_x
                        term{1, 1, 0b10},  // p_y
                        term{1, 2, 0b100}, // n_o
                        term{2, 3, 0b1000} // 1/2 p^2 n_i
                    }};
        }

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 2;
        }

        [[nodiscard]] constexpr T const& operator[](size_t index) const noexcept
        {
            return data[index];
        }

        [[nodiscard]] constexpr T& operator[](size_t index) noexcept
        {
            return data[index];
        }
    };
    // TODO: provide representations for planes, spheres, flats, etc.
} // namespace cga2
} // namespace gal