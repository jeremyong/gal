#pragma once

#include "engine.hpp"
#include "entity.hpp"
#include "geometric_algebra.hpp"
#include "pga.hpp"

#include <cmath>

// 3D Euclidean Vector Space Geometric Algebra

namespace gal
{
namespace vga
{
    using vga_metric = gal::metric<3, 0, 0>;

    using vga_algebra = gal::algebra<vga_metric>;

    constexpr detail::rpne<vga_algebra, 1> operator"" _e0(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b1;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<vga_algebra, 1> operator"" _e1(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b10;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<vga_algebra, 1> operator"" _e2(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b100;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<vga_algebra, 1> operator"" _e01(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b11;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<vga_algebra, 1> operator"" _e02(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b101;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<vga_algebra, 1> operator"" _e12(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b110;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<vga_algebra, 1> operator"" _e012(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b111;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<vga_algebra, 1> operator"" _ps(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b111;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(n), 1}};
    }

    constexpr detail::rpne<vga_algebra, 1> operator"" _ips(unsigned long long n)
    {
        uint32_t op = detail::c_scalar + 0b111;
        return {{detail::node{op, op}}, 1, rat{static_cast<num_t>(-n), 1}};
    }

    template <typename T = float>
    union vector
    {
        using algebra_t = vga_algebra;
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

        GAL_NODISCARD constexpr static auto ie(uint32_t id) noexcept
        {
            return ::gal::detail::construct_ie<algebra_t>(
                id,
                std::make_integer_sequence<width_t, 3>{},
                std::integer_sequence<elem_t, 0b1, 0b10, 0b100>{});
        }

        GAL_NODISCARD constexpr static size_t size() noexcept
        {
            return 3;
        }

        constexpr vector(T a, T b, T c) noexcept
            : x{a}
            , y{b}
            , z{c}
        {}

        template <elem_t... E>
        constexpr vector(entity<vga_algebra, T, E...> in) noexcept
            : data{in.template select<0b1, 0b10, 0b100>()}
        {}

        void normalize() noexcept
        {
            auto l2_inv = T{1} / std::sqrt(x * x + y * y + z * z);
            x           = x * l2_inv;
            y           = y * l2_inv;
            z           = z * l2_inv;
        }

        GAL_NODISCARD constexpr T const& operator[](size_t index) const noexcept
        {
            return data[index];
        }

        GAL_NODISCARD constexpr T& operator[](size_t index) noexcept
        {
            return data[index];
        }
    };

    // A "point" is not typically a recognized quantity in vector spaces as points live on the
    // affine plane. This class represents the point as a vector and provides conversions from other
    // algebras. For example, given a dual projectivized point entity, the constructor will
    // automatically undualize and homogenize the entity. Given a CGA point, the constructor will
    // read off the point coordinates.
    template <typename T>
    union point
    {
        using algebra_t = vga_algebra;
        using value_t   = T;

        constexpr static mv<gal::pga::pga_algebra, 3, 4, 4>
        ie(gal::pga::pga_algebra, uint32_t id) noexcept
        {
            return {mv_size{3, 4, 4},
                    {
                        ind{id + 2, one}, // -z
                        ind{id + 1, one}, // y
                        ind{id, one}      // -x
                    },
                    {mon{minus_one, one, 1, 0}, // -z
                     mon{one, one, 1, 1},       // y
                     mon{minus_one, one, 1, 2}, // -x
                     mon{one, zero, 0, 0}},
                    {
                        term{1, 0, 0b111},  // -z * e012
                        term{1, 1, 0b1011}, // y * e013
                        term{1, 2, 0b1101}, // -x * e023
                        term{1, 3, 0b1110}  // e123
                    }};
        }

        GAL_NODISCARD constexpr static size_t size() noexcept
        {
            return 3;
        }

        T data[3];
        struct
        {
            T x;
            T y;
            T z;
        };

        constexpr point(T x, T y, T z) noexcept
            : x{x}
            , y{y}
            , z{z}
        {}

        // Convert projectivized point
        template <elem_t... E>
        constexpr point(entity<gal::pga::pga_algebra, T, E...> const& other) noexcept
            : data{other.template select<0b1101>(),
                   other.template select<0b1011>(),
                   other.template select<0b111>()}
        {
            auto w_inv = T{1} / other.template select<0b1110>();
            x *= -w_inv;
            y *= w_inv;
            z *= -w_inv;
        }

        GAL_NODISCARD constexpr operator pga::point<T>() const noexcept
        {
            // Dualize and concatenate projective origin
            return {-z, y, -x, T{1}};
        }

        GAL_NODISCARD constexpr T const& operator[](size_t index) const noexcept
        {
            return data[index];
        }

        GAL_NODISCARD constexpr T& operator[](size_t index) noexcept
        {
            return data[index];
        }
    };

    template <typename T>
    union rotor
    {
        using algebra_t = vga_algebra;
        using value_t   = T;

        GAL_NODISCARD constexpr static size_t size() noexcept
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
        GAL_NODISCARD constexpr static mv<vga_algebra, 8, 4, 4> ie(uint32_t id) noexcept
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

        GAL_NODISCARD constexpr T const& operator[](size_t index) const noexcept
        {
            return data[index];
        }

        GAL_NODISCARD constexpr T& operator[](size_t index) noexcept
        {
            return data[index];
        }
    };

    template <typename L, typename... Data>
    auto compute(L lambda, Data const&... input)
    {
        return ::gal::detail::compute<::gal::vga::vga_algebra>(lambda, input...);
    }

    template <typename... Data>
    using evaluate = ::gal::detail::evaluate<gal::vga::vga_algebra, Data...>;
} // namespace vga
} // namespace gal
