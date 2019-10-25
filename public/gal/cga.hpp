#pragma once

#include "entity.hpp"
#include "geometric_algebra.hpp"

#include <cmath>

namespace gal
{
namespace cga
{
    // The metric is defined here as the standard Minkowski spacetime. To extract the conformal representations,
    // a change of basis is required where o = 1/2 * (e + e-) and inf = e- - e.
    //
    // Internally, the change of basis will occur from the null basis to the natural basis before and after expression
    // evaluation. The elements are ordered such that no and ni (null-basis-origin and null-basis-infinity) come at the
    // end such that the parity of all terms of all blades is unchanged after the change of basis.
    using cga_metric = gal::metric<4, 1, 0>;

    // The CGA is a graded algebra with 32 basis elements
    using cga_algebra = gal::algebra<cga_metric>;

    // 0b1 => e+ extension
    // 0b10000 => e- extension
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
} // namespace cga

namespace detail
{
    template <>
    constexpr inline bool uses_null_basis<::gal::cga::cga_algebra> = true;
}

template <typename T>
struct expr<expr_op::identity, mv<cga::cga_algebra, 0, 1, 1>, cga::detail::n_o_tag<T>>
{
    using value_t               = T;
    using algebra_t             = cga::cga_algebra;
    constexpr static expr_op op = expr_op::identity;
    constexpr static mv<cga::cga_algebra, 0, 1, 1> lhs{mv_size{0, 1, 1}, {}, {mon{one, 0, 0}}, {term{1, 0, 0b1000}}};
};

template <typename T>
struct expr<expr_op::identity, mv<cga::cga_algebra, 0, 1, 1>, cga::detail::n_i_tag<T>>
{
    using value_t               = T;
    using algebra_t             = cga::cga_algebra;
    constexpr static expr_op op = expr_op::identity;
    constexpr static mv<cga::cga_algebra, 0, 1, 1> lhs{mv_size{0, 1, 1}, {}, {mon{one, 0, 0}}, {term{1, 0, 0b10000}}};
};

template <typename T>
struct expr<expr_op::identity, mv<cga::cga_algebra, 0, 1, 1>, cga::detail::pseudoscalar_tag<T>>
{
    using value_t               = T;
    using algebra_t             = cga::cga_algebra;
    constexpr static expr_op op = expr_op::identity;
    constexpr static auto lhs   = cga::cga_algebra::pseudoscalar;
};

template <typename T>
struct expr<expr_op::identity, mv<cga::cga_algebra, 0, 1, 1>, cga::detail::pseudoscalar_inv_tag<T>>
{
    using value_t               = T;
    using algebra_t             = cga::cga_algebra;
    constexpr static expr_op op = expr_op::identity;
    constexpr static auto lhs   = cga::cga_algebra::pseudoscalar_inv;
};

namespace cga
{
    template <typename T = float>
    constexpr inline expr<expr_op::identity, mv<cga_algebra, 0, 1, 1>, detail::n_o_tag<T>> n_o;

    template <typename T = float>
    constexpr inline expr<expr_op::identity, mv<cga_algebra, 0, 1, 1>, detail::n_i_tag<T>> n_i;

    template <typename T = float>
    constexpr inline expr<expr_op::identity, mv<cga_algebra, 0, 1, 1>, detail::pseudoscalar_tag<T>> ps;

    template <typename T = float>
    constexpr inline expr<expr_op::identity, mv<cga_algebra, 0, 1, 1>, detail::pseudoscalar_inv_tag<T>> ips;

    template <typename T = float>
    union point
    {
        using algebra_t = cga_algebra;
        using value_t   = T;

        std::array<T, 3> data;
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

            union
            {
                T z;
                T w;
            };
        };

        constexpr point(T a, T b, T c) noexcept
            : x{a}
            , y{b}
            , z{c}
        {}

        template <uint8_t... E>
        constexpr point(entity<algebra_t, T, E...> in) noexcept
            : data{in.template select<0b1, 0b10, 0b100>()}
        {}

        [[nodiscard]] constexpr static mv<algebra_t, 6, 7, 5> ie(uint32_t id) noexcept
        {
            // A CGA point is represented as no + p + 1/2 p^2 ni
            return {mv_size{6, 7, 5},
                    {
                        ind{id, 1},     // ind0 = p_x
                        ind{id + 1, 1}, // ind1 = p_y
                        ind{id + 2, 1}, // ind2 = p_z
                        ind{id, 2},     // ind3 = p_x^2
                        ind{id + 1, 2}, // ind4 = p_y^2
                        ind{id + 2, 2}, // ind5 = p_z^2
                    },
                    {
                        mon{one, 1, 0, 1},      // p_x
                        mon{one, 1, 1, 1},      // p_y
                        mon{one, 1, 2, 1},      // p_z
                        mon{one, 0, 0, 0},      // no
                        mon{one_half, 1, 3, 2}, // 1/2 p_x^2
                        mon{one_half, 1, 4, 2}, // 1/2 p_y^2
                        mon{one_half, 1, 5, 2}, // 1/2 p_z^2
                    },
                    {
                        term{1, 0, 0b1},    // p_x
                        term{1, 1, 0b10},   // p_y
                        term{1, 2, 0b100},  // p_z
                        term{1, 3, 0b1000}, // no
                        term{3, 4, 0b10000} // 1/2 p^2 ni
                    }};
        }

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 3;
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
            return NAN;
        }
    };
    // TODO: provide representations for planes, spheres, flats, etc.
} // namespace cga
} // namespace gal