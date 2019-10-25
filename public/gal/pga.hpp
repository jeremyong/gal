#pragma once

#include "entity.hpp"
#include "geometric_algebra.hpp"

#include <cmath>

// The Projective Geometric Algebra for representing Euclidean 3-space
// NOTE: In comments, we often write "the PGA" to mean literally "the projective geometric algebra."

namespace gal
{
namespace pga
{
    // NOTE: the inner product of e0 can be set to +1 or -1 without any change in the algebra's geometric
    // interpretation. Here, we opt to define e0^2 := 1 by convention
    using pga_metric = ::gal::metric<3, 0, 1>;

    // The PGA is a graded algebra with 16 basis elements
    using pga_algebra = gal::algebra<pga_metric>;

    constexpr inline auto e = gal::e<pga_algebra, 0>;
    constexpr inline auto e0 = gal::e<pga_algebra, 0b1>;
    constexpr inline auto e1 = gal::e<pga_algebra, 0b10>;
    constexpr inline auto e2 = gal::e<pga_algebra, 0b100>;
    constexpr inline auto e3 = gal::e<pga_algebra, 0b1000>;
    constexpr inline auto e01 = gal::e<pga_algebra, 0b11>;
    constexpr inline auto e02 = gal::e<pga_algebra, 0b101>;
    constexpr inline auto e03 = gal::e<pga_algebra, 0b1001>;
    constexpr inline auto e12 = gal::e<pga_algebra, 0b110>;
    constexpr inline auto e13 = gal::e<pga_algebra, 0b1010>;
    constexpr inline auto e23 = gal::e<pga_algebra, 0b1100>;
    constexpr inline auto e012 = gal::e<pga_algebra, 0b111>;
    constexpr inline auto e013 = gal::e<pga_algebra, 0b1011>;
    constexpr inline auto e023 = gal::e<pga_algebra, 0b1101>;
    constexpr inline auto e123 = gal::e<pga_algebra, 0b1110>;
    constexpr inline auto e1234 = gal::e<pga_algebra, 0b1111>;

    template <typename T = float>
    union plane
    {
        using algebra_t = pga_algebra;
        using value_t   = T;
        constexpr static bool is_dual = true;

        std::array<T, 4> data;
        struct
        {
            T d;
            T x;
            T y;
            T z;
        };

        [[nodiscard]] constexpr static auto ie(uint32_t id) noexcept
        {
            return detail::construct_ie<algebra_t>(
                id, std::make_integer_sequence<width_t, 4>{}, std::integer_sequence<uint8_t, 0b1, 0b10, 0b100, 0b1000>{});
        }

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 4;
        }

        [[nodiscard]] constexpr static uint32_t ind_count() noexcept
        {
            return 4;
        }

        constexpr plane(T d, T x, T y, T z) noexcept
            : d{d}
            , x{x}
            , y{y}
            , z{z}
        {}

        template <uint8_t... E>
        constexpr plane(entity<pga_algebra, T, E...> in) noexcept
            : data{in.template select<0b1, 0b10, 0b100, 0b1000>()}
        {
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
    union point
    {
        using algebra_t = pga_algebra;
        using value_t   = T;
        constexpr static bool is_dual = true;

        std::array<T, 3> data;
        struct
        {
            T x;
            T y;
            T z;
        };

        // Like planes, points are represented dually as the intersection of three planes
        [[nodiscard]] constexpr static mv<algebra_t, 3, 4, 4> ie(uint32_t id) noexcept
        {
            return {mv_size{3, 4, 4},
                    {
                        ind{id + 2, 1}, // -z
                        ind{id + 1, 1}, // y
                        ind{id, 1}      // -x
                    },
                    {mon{minus_one, 1, 0, 1}, // -z
                     mon{one, 1, 1, 1},       // y
                     mon{minus_one, 1, 2, 1}, // -x
                     mon{one, 0, 0, 0}},
                    {
                        term{1, 0, 0b111},  // -z * e012
                        term{1, 1, 0b1011}, // y * e013
                        term{1, 2, 0b1101}, // x * e023
                        term{1, 3, 0b1110}  // e123
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

        constexpr point(T x, T y, T z) noexcept
            : x{x}
            , y{y}
            , z{z}
        {}

        template <uint8_t... E>
        constexpr point(entity<pga_algebra, T, E...> in) noexcept
            : data{}
        {
            auto input = in.template select<0b111, 0b1011, 0b1101, 0b1110>();
            auto w_inv = T{1} / input[3];
            z = -input[0] * w_inv;
            y = input[1] * w_inv;
            x = -input[2] * w_inv;
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
        using algebra_t = pga_algebra;
        using value_t   = T;
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
            : data{}
        {
            auto input = in.template select<0b111, 0b1011, 0b1101>();
            z = -input[0];
            y = input[1];
            x = -input[2];
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

    // Lines in P^3 are defined using Plücker coordinates: https://en.wikipedia.org/wiki/Plücker_coordinates
    // The lines e01, e02, and e03 are the ideal lines representing the intersections of e1, e2, and e3 with the ideal
    // plane respectively. The lines e23, e31, and e12 are lines through the origin in the x, y, and z directions
    // respectively. We opt not to provide a 6 coordinate representation for now (join points to construct a line or
    // meet planes)

    // TODO: It isn't great that we cache the cos and sin of the rotor since storing 5 elements prevents a more natural
    // tightly-packed alignment
    template <typename T = float>
    union rotor
    {
        using algebra_t = pga_algebra;
        using value_t = T;

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 5;
        }

        [[nodiscard]] constexpr static uint32_t ind_count() noexcept
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
        // A rotation t around a line is given by the expression cos(t/2) + sin(t/2)(l_x + l_y + l_z)
        [[nodiscard]] constexpr static mv<pga_algebra, 8, 4, 4> ie(uint32_t id) noexcept
        {
            return {mv_size{7, 4, 4},
                    {ind{id, 1},     // cos(t/2)
                     ind{id + 1, 1}, // z * sin(t/2)
                     ind{id + 4, 1},
                     ind{id + 1, 1}, // -y * sin(t/2)
                     ind{id + 3, 1},
                     ind{id + 1, 1}, // x * sin(t/2)
                     ind{id + 2, 1}},
                    {mon{one, 1, 0, 1}, mon{one, 1, 1, 1}, mon{minus_one, 2, 2, 2}, mon{one, 2, 3, 2}},
                    {term{1, 0, 0b0}, term{1, 1, 0b110}, term{1, 2, 0b1010}, term{1, 3, 0b1100}}};
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

    template <typename T = float>
    union translator
    {
        using algebra_t = pga_algebra;
        using value_t = T;

        [[nodiscard]] constexpr static size_t size() noexcept
        {
            return 4;
        }

        [[nodiscard]] constexpr static uint32_t ind_count() noexcept
        {
            return 4;
        }

        std::array<T, 4> data;
        struct
        {
            T d;
            T x;
            T y;
            T z;
        };

        constexpr translator(T d, T x, T y, T z) noexcept
            : d{d}
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

        // A translation along a line with distance d is given by the expression 1 + d/2(P_inf)
        [[nodiscard]] constexpr static mv<pga_algebra, 6, 4, 4> ie(uint32_t id) noexcept
        {
            return {mv_size{6, 4, 4},
                    {ind{id, 1}, // d * l_x
                     ind{id + 1, 1},
                     ind{id, 1}, // d * l_y
                     ind{id + 1, 1},
                     ind{id, 1}, // d * l_z
                     ind{id + 2, 1}},
                    {
                        mon{one, 0, 0, 0},            // 1
                        mon{minus_one_half, 2, 1, 2}, // -1/2 * l_x
                        mon{minus_one_half, 2, 2, 2}, // -1/2 * l_y
                        mon{minus_one_half, 2, 3, 2}  // -1/2 * l_z
                    },
                    {term{1, 0, 0b0}, term{1, 1, 0b11}, term{1, 2, 0b101}, term{1, 3, 0b1001}}};
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
} // namespace pga
} // namespace gal