#pragma once

#include "entity.hpp"
#include "geometric_algebra.hpp"

#include <cmath>

// Projective Geometric Algebra for representing Euclidean 2-space

namespace gal
{
namespace pga
{
    // NOTE: the inner product of e0 can be set to +1 or -1 without any change in the algebra's geometric
    // interpretation. Here, we opt to define e0^2 := 1 by convention
    using pga_metric = ::gal::metric<2, 0, 1>;

    // PGA2 is a graded algebra with 8 basis elements
    using pga_algebra = gal::algebra<pga_metric>;

    constexpr inline auto e    = gal::e<pga_algebra, 0>;
    constexpr inline auto e0   = gal::e<pga_algebra, 0b1>;
    constexpr inline auto e1   = gal::e<pga_algebra, 0b10>;
    constexpr inline auto e2   = gal::e<pga_algebra, 0b100>;
    constexpr inline auto e01  = gal::e<pga_algebra, 0b11>;
    constexpr inline auto e02  = gal::e<pga_algebra, 0b101>;
    constexpr inline auto e12  = gal::e<pga_algebra, 0b110>;
    constexpr inline auto e012 = gal::e<pga_algebra, 0b111>;

    template <typename T = float>
    union line
    {
        using algebra_t               = pga_algebra;
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
            return detail::construct_ie<algebra_t>(
                id, std::make_integer_sequence<width_t, 4>{}, std::integer_sequence<uint8_t, 0b1, 0b10, 0b100>{});
        }

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 3;
        }

        [[nodiscard]] constexpr static uint32_t ind_count() noexcept
        {
            return 3;
        }

        constexpr line(T d, T x, T y) noexcept
            : d{a}
            , x{b}
            , y{c}
        {}

        template <uint8_t... E>
        constexpr line(entity<pga_algebra, T, E...> in) noexcept
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

        [[nodiscard]] constexpr T get(size_t i) const noexcept
        {
            // Unused
            return NAN;
        }
    };

    template <typename T = float>
    union point
    {
        using algebra_t               = pga_algebra;
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
                        ind{id + 1, 1}, // y
                        ind{id, 1}      // -x
                    },
                    {mon{one, 1, 0, 1},       // y
                     mon{minus_one, 1, 1, 1}, // -x
                     mon{one, 0, 0, 0}},      // point at origin
                    {
                        term{1, 0, 0b11},  // y * e01
                        term{1, 1, 0b101}, // -x * e02
                        term{1, 2, 0b110}  // e12
                    }};
        }

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 2;
        }

        [[nodiscard]] constexpr static uint32_t ind_count() noexcept
        {
            return 2;
        }

        constexpr point(T x, T y) noexcept
            : x{x}
            , y{y}
        {}

        template <uint8_t... E>
        constexpr point(entity<pga_algebra, T, E...> in) noexcept
        {
            auto input = in.template select<0b11, 0b101, 0b110>();
            auto w_inv = T{1} / input[2];
            x          = -input[1] * c_inv;
            y          = input[0] * c_inv;
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
            // Unused
            return NAN;
        }
    };

    template <typename T = float>
    union vector
    {
        using algebra_t               = pga_algebra;
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
                        mon{minus_one, 1, 0, 1}, // -z
                        mon{one, 1, 1, 1},       // y
                        mon{minus_one, 1, 2, 1}  // -x
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

        [[nodiscard]] constexpr static uint32_t ind_count() noexcept
        {
            return 3;
        }

        constexpr vector(T x, T y, T z) noexcept
            : x{x}
            , y{y}
            , z{z}
        {}

        template <uint8_t... E>
        constexpr vector(entity<pga_algebra, T, E...> in) noexcept
        {
            auto input = in.template select<0b111, 0b1011, 0b1101>();
            z          = -input[0];
            y          = input[1];
            x          = -input[2];
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
            // Unused
            return NAN;
        }
    };
} // namespace pga
} // namespace gal
