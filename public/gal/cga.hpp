#pragma once

#include "engine.hpp"
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
        template <typename T>
        struct n_o_t
        {
            using value_t = T;
            using algebra_t = cga::cga_algebra;
            constexpr static mv<cga::cga_algebra, 0, 1, 1> value{mv_size{0, 1, 1},
                                                                 {},
                                                                 {mon{one, zero, 0, 0}},
                                                                 {term{1, 0, 0b1000}}};
            [[nodiscard]] constexpr static auto ie(uint32_t) noexcept
            {
                return value;
            }
        };

        template <typename T>
        struct n_i_t
        {
            using value_t = T;
            using algebra_t = cga::cga_algebra;
            constexpr static mv<cga::cga_algebra, 0, 1, 1> value{mv_size{0, 1, 1},
                                                                 {},
                                                                 {mon{one, zero, 0, 0}},
                                                                 {term{1, 0, 0b10000}}};
            [[nodiscard]] constexpr static auto ie(uint32_t) noexcept
            {
                return value;
            }
        };

        template <typename T>
        struct ps_t
        {
            using value_t = T;
            using algebra_t = cga::cga_algebra;
            [[nodiscard]] constexpr static auto ie(uint32_t) noexcept
            {
                return ::gal::detail::to_null_basis(cga::cga_algebra::pseudoscalar);
            }
        };

        template <typename T>
        struct ips_t
        {
            using value_t = T;
            using algebra_t = cga::cga_algebra;
            [[nodiscard]] constexpr static auto ie(uint32_t) noexcept
            {
                return ::gal::detail::to_null_basis(cga::cga_algebra::pseudoscalar_inv);
            }
        };
    } // namespace detail
} // namespace cga

namespace detail
{
    template <>
    constexpr inline bool uses_null_basis<::gal::cga::cga_algebra> = true;
}

namespace cga
{
    template <typename T = float>
    constexpr inline ::gal::detail::expr_id<detail::n_o_t<T>, 0> n_o;

    template <typename T = float>
    constexpr inline ::gal::detail::expr_id<detail::n_i_t<T>, 0> n_i;

    template <typename T = float>
    constexpr inline ::gal::detail::expr_id<detail::ps_t<T>, 0> ps;

    template <typename T = float>
    constexpr inline ::gal::detail::expr_id<detail::ips_t<T>, 0> ips;

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
                        ind{id, rat{1}},     // ind0 = p_x
                        ind{id + 1, rat{1}}, // ind1 = p_y
                        ind{id + 2, rat{1}}, // ind2 = p_z
                        ind{id, rat{2}},     // ind3 = p_x^2
                        ind{id + 1, rat{2}}, // ind4 = p_y^2
                        ind{id + 2, rat{2}}, // ind5 = p_z^2
                    },
                    {
                        mon{one, one, 1, 0},         // p_x
                        mon{one, one, 1, 1},         // p_y
                        mon{one, one, 1, 2},         // p_z
                        mon{one, zero, 0, 0},        // no
                        mon{one_half, rat{2}, 1, 3}, // 1/2 p_x^2
                        mon{one_half, rat{2}, 1, 4}, // 1/2 p_y^2
                        mon{one_half, rat{2}, 1, 5}, // 1/2 p_z^2
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
} // namespace cga
} // namespace gal