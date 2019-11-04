#pragma once

#include "engine.hpp"
#include "entity.hpp"
#include "geometric_algebra.hpp"

#include <cmath>

// The 3D Euclidean Geometric Algebra

namespace gal
{
namespace ega
{
    using ega_metric = gal::metric<3, 0, 0>;

    using ega_algebra = gal::algebra<ega_metric>;

    constexpr detail::rpne<ega_algebra, 1> operator"" _e0(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b1;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<ega_algebra, 1> operator"" _e1(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b10;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<ega_algebra, 1> operator"" _e2(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b100;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<ega_algebra, 1> operator"" _e01(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b11;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<ega_algebra, 1> operator"" _e02(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b101;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<ega_algebra, 1> operator"" _e12(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b110;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<ega_algebra, 1> operator"" _e012(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b111;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<ega_algebra, 1> operator"" _ps(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b111;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<ega_algebra, 1> operator"" _ips(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b111;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(-n), 1}};
    }

    template <typename T = float>
    union vector
    {
        using algebra_t = ega_algebra;
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

        [[nodiscard]] constexpr static auto ie(uint32_t id) noexcept
        {
            return ::gal::detail::construct_ie<algebra_t>(
                id,
                std::make_integer_sequence<width_t, 3>{},
                std::integer_sequence<elem_t, 0b1, 0b10, 0b100>{});
        }

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 3;
        }

        constexpr vector(T a, T b, T c) noexcept
            : x{a}
            , y{b}
            , z{c}
        {}

        template <elem_t... E>
        constexpr vector(entity<ega_algebra, T, E...> in) noexcept
            : data{in.template select<0b1, 0b10, 0b100>()}
        {}

        void normalize() noexcept
        {
            auto l2_inv = T{1} / std::sqrt(x * x + y * y + z * z);
            x           = x * l2_inv;
            y           = y * l2_inv;
            z           = z * l2_inv;
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

    template <typename T>
    union rotor
    {
        using algebra_t = ega_algebra;
        using value_t   = T;

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 5;
        }

        std::array<T, 5> data;
        struct
        {
            T cos_theta;
            T sin_theta;
            T x;
            T y;
            T z;
        };

        constexpr rotor(T theta, T x, T y, T z) noexcept
            : cos_theta{std::cos(T{0.5} * theta)}
            , sin_theta{std::sin(T{0.5} * theta)}
            , x{x}
            , y{y}
            , z{z}
        {}

        // As always when doing any normalization operation, NaNs are producible when normalizing
        // vectors of zero length. This is not checked for!
        void normalize() noexcept
        {
            auto l2_inv = T{1} / std::sqrt(x * x + y * y + z * z);
            x           = x * l2_inv;
            y           = y * l2_inv;
            z           = z * l2_inv;
        }

        // Cos\theta := ID 0
        // Sin\theta := ID 1
        // x := ID 2
        // y := ID 3
        // z := ID 4
        [[nodiscard]] constexpr static mv<ega_algebra, 8, 4, 4> ie(uint32_t id) noexcept
        {
            return {mv_size{7, 4, 4},
                    {ind{id, one},     // cos(t/2)
                     ind{id + 1, one}, // z * sin(t/2)
                     ind{id + 4, one},
                     ind{id + 1, one}, // -y * sin(t/2)
                     ind{id + 3, one},
                     ind{id + 1, one}, // x * sin(t/2)
                     ind{id + 2, one}},
                    {mon{one, one, 1, 0},
                     mon{minus_one, rat{2}, 2, 1},
                     mon{one, rat{2}, 2, 3},
                     mon{minus_one, rat{2}, 2, 5}},
                    {term{1, 0, 0}, term{1, 1, 0b11}, term{1, 2, 0b101}, term{1, 3, 0b110}}};
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
} // namespace ega
} // namespace gal
