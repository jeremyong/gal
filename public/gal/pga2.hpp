#pragma once

#include "engine.hpp"
#include "entity.hpp"
#include "geometric_algebra.hpp"

#include <cmath>

// Projective Geometric Algebra for representing Euclidean 2-space

namespace gal
{
namespace pga2
{
    // NOTE: the inner product of e0 can be set to +1 or -1 without any change in the algebra's
    // geometric interpretation. Here, we opt to define e0^2 := 1 by convention
    using pga2_metric = ::gal::metric<2, 0, 1>;

    // PGA2 is a graded algebra with 8 basis elements
    using pga2_algebra = gal::algebra<pga2_metric>;

    constexpr detail::rpne<pga2_algebra, 1> operator"" _e0(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b1;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<pga2_algebra, 1> operator"" _e1(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b10;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<pga2_algebra, 1> operator"" _e2(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b100;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<pga2_algebra, 1> operator"" _e01(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b11;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<pga2_algebra, 1> operator"" _e02(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b101;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<pga2_algebra, 1> operator"" _e12(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b110;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<pga2_algebra, 1> operator"" _e012(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b111;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<pga2_algebra, 1> operator"" _ps(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b111;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    template <typename T = float>
    union line
    {
        using algebra_t               = pga2_algebra;
        using value_t                 = T;
        constexpr static bool is_dual = true;

        std::array<T, 3> data;
        // Line coordinates in PGA2 satisfy the equation 0 = ax + by + c = 0
        struct
        {
            T d;
            T x;
            T y;
        };

        [[nodiscard]] constexpr static auto ie(uint32_t id) noexcept
        {
            return gal::detail::construct_ie<algebra_t>(
                id,
                std::make_integer_sequence<width_t, 3>{},
                std::integer_sequence<elem_t, 0b1, 0b10, 0b100>{});
        }

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 3;
        }

        constexpr line(T d, T x, T y) noexcept
            : d{d}
            , x{x}
            , y{y}
        {}

        template <elem_t... E>
        constexpr line(entity<pga2_algebra, T, E...> in) noexcept
            : data{in.template select<0b1, 0b10, 0b100>()}
        {}

        [[nodiscard]] constexpr T const& operator[](size_t index) const noexcept
        {
            return data[index];
        }

        [[nodiscard]] constexpr T& operator[](size_t index) noexcept
        {
            return data[index];
        }
    };

    template <typename T = float>
    union point
    {
        using algebra_t               = pga2_algebra;
        using value_t                 = T;
        constexpr static bool is_dual = true;

        std::array<T, 2> data;
        struct
        {
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

        // Like planes, points are represented dually as the intersection of three planes
        [[nodiscard]] constexpr static mv<algebra_t, 2, 3, 3> ie(uint32_t id) noexcept
        {
            return {mv_size{2, 3, 3},
                    {
                        ind{id + 1, 1}, // -y
                        ind{id, 1}      // x
                    },
                    {mon{minus_one, one, 0, 1}, // -y
                     mon{one, one, 1, 1},       // x
                     mon{one, zero, 0, 0}},     // point at origin
                    {
                        term{1, 0, 0b11},  // -y * e01
                        term{1, 1, 0b101}, // x * e02
                        term{1, 2, 0b110}  // e12
                    }};
        }

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 2;
        }

        constexpr point(T x, T y) noexcept
            : x{x}
            , y{y}
        {}

        template <elem_t... E>
        constexpr point(entity<pga2_algebra, T, E...> in) noexcept
            : x{in.template select<0b101>()}
            , y{in.template select<0b11>()}
        {
            auto w_inv = T{1} / in.template select<0b110>();
            x          = w_inv;
            y          = -w_inv;
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

    template <typename T = float>
    union vector
    {
        using algebra_t               = pga2_algebra;
        using value_t                 = T;
        constexpr static bool is_dual = true;

        std::array<T, 3> data;
        struct
        {
            T x;
            T y;
            T z;
        };

        // Like planes, points are represented dually as the intersection of three planes
        [[nodiscard]] constexpr static mv<algebra_t, 3, 3, 3> ie(uint32_t id) noexcept
        {
            return {mv_size{3, 4, 4},
                    {
                        ind{id + 2, 1}, // -z
                        ind{id + 1, 1}, // y
                        ind{id, 1}      // -x
                    },
                    {
                        mon{minus_one, one, 0, 1}, // -z
                        mon{one, one, 1, 1},       // y
                        mon{minus_one, one, 2, 1}  // -x
                    },
                    {
                        term{1, 0, 0b111},  // -z * e012
                        term{1, 1, 0b1011}, // y * e013
                        term{1, 2, 0b1101}  // x * e023
                    }};
        }

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 3;
        }

        constexpr vector(T x, T y, T z) noexcept
            : x{x}
            , y{y}
            , z{z}
        {}

        template <elem_t... E>
        constexpr vector(entity<pga2_algebra, T, E...> in) noexcept
            : data{in.template select<0b1101, 0b1011, 0b111>()}
        {
            z = -z;
            x = -x;
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
} // namespace pga2
} // namespace gal
