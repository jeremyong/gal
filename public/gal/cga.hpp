#pragma once

#include "engine.hpp"
#include "entity.hpp"
#include "geometric_algebra.hpp"

#include <cmath>

namespace gal
{
namespace cga
{
    // The metric is defined here as the standard Minkowski spacetime. To extract the conformal
    // representations, a change of basis is required where o = 1/2 * (e + e-) and inf = e- - e.
    //
    // Internally, the change of basis will occur from the null basis to the natural basis before
    // and after expression evaluation. The elements are ordered such that no and ni
    // (null-basis-origin and null-basis-infinity) come at the end such that the parity of all terms
    // of all blades is unchanged after the change of basis.
    using cga_metric = gal::metric<4, 1, 0>;

    // The CGA is a graded algebra with 32 basis elements
    using cga_algebra = gal::algebra<cga_metric>;

    // 0b1 => e+ extension
    // 0b10000 => e- extension

    constexpr detail::rpne<cga_algebra, 1> operator"" _e1(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b1;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<cga_algebra, 1> operator"" _e2(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b10;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<cga_algebra, 1> operator"" _e3(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b100;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    // Before change-of-basis:
    // 0b1000 => no
    // 0b10000 => ni

    constexpr detail::rpne<cga_algebra, 1> operator"" _no(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b1000;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<cga_algebra, 1> operator"" _ni(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b10000;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<cga_algebra, 1> operator"" _ps(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b11111;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<cga_algebra, 1> operator"" _ips(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b11111;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(-n), 1}};
    }
} // namespace cga

namespace detail
{
    template <>
    constexpr inline bool uses_null_basis<::gal::cga::cga_algebra> = true;
}

namespace cga
{
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

        template <elem_t... E>
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
