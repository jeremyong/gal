#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>

#include "opt.hpp"

#ifdef _MSC_VER
#    include <intrin.h>
#endif

namespace gal
{
// right-to-left binary exponentiation
template <typename T, int N, int D>
[[nodiscard]] GAL_FORCE_INLINE constexpr T
pow(T s, std::integral_constant<int, N>, std::integral_constant<int, D>) noexcept
{
    if constexpr (D > 1)
    {
        return std::pow(s, T{N} / T{D});
    }
    else if constexpr (N == 1)
    {
        return s;
    }
    else if constexpr (N == -1)
    {
        return T{1} / s;
    }
    else if constexpr (N < 0)
    {
        return T{1}
               / ::gal::pow(s, std::integral_constant<int, -N>{}, std::integral_constant<int, 1>{});
    }
    else if constexpr (N == 2)
    {
        return s * s;
    }
    else if constexpr (N == 3)
    {
        return s * s * s;
    }
    else if constexpr (N == 4)
    {
        auto s2 = s * s;
        return s2 * s2;
    }
    else if constexpr (N == 5)
    {
        auto s2 = s * s;
        return s2 * s2 * s;
    }
    else if constexpr (N == 6)
    {
        auto s2 = s * s;
        auto s3 = s2 * s;
        return s3 * s3;
    }
    else if constexpr (N == 7)
    {
        auto s2 = s * s;
        auto s4 = s2 * s2;
        return s4 * s2 * s;
    }
    else if constexpr (N == 8)
    {
        auto s2 = s * s;
        auto s4 = s2 * s2;
        return s4 * s4;
    }
    else
    {
        T temp1{1};
        T temp2{s};
        int e = N;
        while (e > 1)
        {
            if ((e & 1) == 1)
            {
                temp1 *= temp2;
            }
            temp2 *= temp2;
            e = e >> 1;
        }
        return temp1 * temp2;
    }
}

template <typename T>
[[nodiscard]] constexpr T next_pow_2(T s) noexcept
{
    // This implementation will change when the size of the input changes
    if constexpr (sizeof(T) == 4)
    {
        if (s == 0)
        {
            return 0;
        }

        --s;
        s |= s >> 1;
        s |= s >> 2;
        s |= s >> 4;
        s |= s >> 8;
        s |= s >> 16;
        return ++s;
    }
    else if constexpr (sizeof(T) == 8)
    {
        if (s == 0)
        {
            return 0;
        }

        --s;
        s |= s >> 1;
        s |= s >> 2;
        s |= s >> 4;
        s |= s >> 8;
        s |= s >> 16;
        s |= s >> 32;
        return ++s;
    }
}

// Rotate an input a single byte to the right.
[[nodiscard]] constexpr uint32_t rot32(uint32_t input) noexcept
{
    return (input >> 8) | ((input & 0xff) << 24);
}

// Precondition: input is expected to be greater than 0 or results are undefined
[[nodiscard]] constexpr uint32_t leading_set_index(uint32_t input) noexcept
{
#if defined(__clang) || defined(__GNUG__)
    return 31 - __builtin_clz(input);
#elif defined(_MSC_VER)
    uint32_t leading_zero = 0;
    if (_BitScanReverse(&leading_zero, input))
    {
        return 31 - leading_zero;
    }
    else
    {
        return 32;
    }
#endif
}

[[nodiscard]] constexpr uint32_t pop_count(uint32_t input) noexcept
{
#if defined(__clang) || defined(__GNUG__)
    return __builtin_popcount(input);
#elif defined(_MSC_VER)
    return __popcnt(input);
#endif
}

using num_t = int32_t;
// Even though the denominator will always be kept greater than zero as an invariant, it is
// convenient for several operations to use the signed quantity as an intermediate quantity.
using den_t = int32_t;

// The module we work with is attached to the field of rational numbers.
// The numerator and denominator are left as ints (even though D > 0 is an invariant) so the
// compiler can help detect overflows
struct rat
{
    num_t num = 0;
    den_t den = 1;

    [[nodiscard]] constexpr bool is_zero() const noexcept
    {
        return num == 0;
    }

    [[nodiscard]] constexpr bool is_unit() const noexcept
    {
        return num == den;
    }

    [[nodiscard]] constexpr rat reciprocal() const noexcept
    {
        return {den, num};
    }

    [[nodiscard]] constexpr rat negation() const noexcept
    {
        return {-num, den};
    }

    constexpr rat& operator/=(rat other) noexcept
    {
        this->operator*=(other.reciprocal());
        return *this;
    }

    constexpr rat& operator*=(rat other) noexcept
    {
        num *= other.num;
        den *= other.den;
        if (num == 0)
        {
            den = 1;
        }
        else if (den > 1)
        {
            auto gcd = std::gcd(num, den);
            if (gcd > 1)
            {
                num /= gcd;
                den /= gcd;
            }
        }
        return *this;
    }

    constexpr rat& operator+=(rat other) noexcept
    {
        num = num * other.den + other.num * den;
        if (num == 0)
        {
            den = 1;
        }
        else
        {
            den = den * other.den;
            if (den > 1)
            {
                auto gcd = std::gcd(num, den);
                if (gcd > 1)
                {
                    num /= gcd;
                    den /= gcd;
                }
            }
        }
        return *this;
    }

    template <typename T>
    [[nodiscard]] constexpr operator T() const noexcept
    {
        static_assert(std::is_floating_point_v<T>,
                      "Attempting to cast a rational number to non-floating-point type.");
        return static_cast<T>(num) / static_cast<T>(den);
    }
};

// Common rational types used for brevity
constexpr inline rat one{1, 1};
constexpr inline rat minus_one{-1, 1};
constexpr inline rat one_half{1, 2};
constexpr inline rat minus_one_half{-1, 2};
constexpr inline rat zero{0, 1};

namespace detail
{
    [[nodiscard]] constexpr num_t abs(num_t in) noexcept
    {
        return in < 0 ? -in : in;
    }

    [[nodiscard]] constexpr double abs(double in) noexcept
    {
        return in < 0.0 ? -in : in;
    }

    [[nodiscard]] constexpr rat overflow_gate(rat in) noexcept
    {
        // As expressions expand, there may be cases where terms becoming vanishingly small or N and
        // D become great enough to risk overflow. This function is a pure function which
        // deterministically nudges the result so that compilation stays fast without sacrificing
        // too much precision.

        // The goal is to at least match the precision that would have been afforded with floating
        // point precision for the range of numbers dealt with.

        // If we call the error introduced e, let's choose e to be 1e-7. This is only slightly worse
        // than the precision afforded by floating point numbers in the single digit ranges. Any
        // fraction of the form p/q is the mediant of two fractions which bound it from above and
        // below. A possible choice is floor(p/2)/ceil(q/2) and ceil(p/2)/floor(q/2). Suppose p is
        // even. Then floor(p/2) = ceil(p/2) = p/2. Assuming p/q is irreducible, this implies that
        // floor(q/2) = (q-1)/2 and ceil(q/2) = (q + 1)/2. The error of the lower bound then is p/q
        // - p/(q
        // + 1) = p/(q^2 + 1). Conversely, ceil(p/2)/floor(q/2) for even p reduces to p/(q - 1) and
        // the error is given by p/(q^2 - 1). For large q, the error then is approximately p/q^2.

        // We need to account for erroneously introducing systematic bias into the expression. This
        // occurs if we consistently perturb the denominator up or down. To balance things out, the
        // direction of the perturbation will be based on the parity of the second least significant
        // bit of the denominator. This is by no means perfect given that the distribution of
        // rationals encountered is itself biased, but cheap to compute and difficult to beat
        // without introducing unreasonable amounts of complexity and compilation time.

        if (in.den < (1 << 10))
        {
            auto gcd = std::gcd(in.num, in.den);
            if (gcd > 1)
            {
                return {in.num / gcd, in.den / gcd};
            }
            else
            {
                return in;
            }
        }
        else
        {
            auto gcd = std::gcd(in.num, in.den);
            if (gcd > 1)
            {
                return {in.num / gcd, in.den / gcd};
            }
            else
            {
                double n_d     = static_cast<double>(in.num);
                double d_d     = static_cast<double>(in.den);
                double epsilon = abs(n_d / d_d / d_d);
                double frac    = n_d / d_d;

                if (abs(frac) < 1e-7)
                {
                    return zero;
                }
                else if (epsilon < 1e-7)
                {
                    // The rational is the mediant of two fractions with smaller denominators. Pick
                    // one of them based on the parity.
                    if (in.num % 2 == 1)
                    {
                        // Perturb to the left-side of the mediant
                        auto n2  = (in.num - 1) / 2;
                        auto d2  = in.den / 2;
                        auto gcd = std::gcd(n2, d2);
                        if (gcd > 1)
                        {
                            return {n2 / gcd, d2 / gcd};
                        }
                        else
                        {
                            return {n2, d2};
                        }
                    }
                    else
                    {
                        // Perturb to the right-side of the mediant
                        auto n2  = in.num / 2;
                        auto d2  = (in.den - 1) / 2;
                        auto gcd = std::gcd(n2, d2);
                        if (gcd > 1)
                        {
                            return {n2 / gcd, d2 / gcd};
                        }
                        else
                        {
                            return {n2, d2};
                        }
                    }
                }
                else
                {
                    return in;
                }
            }
        }
    }
} // namespace detail

[[nodiscard]] constexpr bool operator==(rat lhs, rat rhs) noexcept
{
    auto gcd1 = std::gcd(lhs.num, lhs.den);
    auto gcd2 = std::gcd(rhs.den, rhs.den);
    return lhs.num / gcd1 == rhs.num / gcd2 && lhs.den / gcd1 == rhs.den / gcd2;
}

[[nodiscard]] constexpr bool operator!=(rat lhs, int rhs) noexcept
{
    return lhs.den != 1 || lhs.num != rhs;
}

[[nodiscard]] constexpr bool operator!=(rat lhs, rat rhs) noexcept
{
    return lhs.num != rhs.num || lhs.den != rhs.den;
}

[[nodiscard]] constexpr bool operator<(rat lhs, rat rhs) noexcept
{
    return (lhs.num * rhs.den) < (rhs.num * lhs.den);
}

[[nodiscard]] constexpr bool operator>(rat lhs, rat rhs) noexcept
{
    return (lhs.num * rhs.den) > (rhs.num * lhs.den);
}

[[nodiscard]] constexpr rat operator-(rat in) noexcept
{
    return in.negation();
}

[[nodiscard]] constexpr rat operator*(int lhs, rat rhs) noexcept
{
    if (lhs == 0 || rhs.num == 0)
    {
        return zero;
    }
    else
    {
        return {rhs.num * lhs, rhs.den};
    }
}

[[nodiscard]] constexpr rat operator*(rat lhs, rat rhs) noexcept
{
    rat product{lhs.num * rhs.num, lhs.den * rhs.den};
    if (product.num == 0)
    {
        return zero;
    }
    else if (product.den > 1)
    {
        // Periodically reduce the fraction if possible to prevent overflow
        return detail::overflow_gate(product);
    }
    else
    {
        return product;
    }
}

[[nodiscard]] constexpr rat operator/(rat lhs, int rhs) noexcept
{
    if (rhs < 0)
    {
        lhs.num = -lhs.num;
        lhs.den *= -rhs;
    }
    else
    {
        lhs.den *= rhs;
    }
    return detail::overflow_gate(lhs);
}

[[nodiscard]] constexpr rat operator/(rat lhs, rat rhs) noexcept
{
    rat quotient{lhs.num * rhs.den, rhs.num * lhs.den};
    if (quotient.num == 0)
    {
        return zero;
    }
    else if (detail::abs(quotient.den) > 1)
    {
        if (quotient.den < 0)
        {
            return detail::overflow_gate(rat{-quotient.num, -quotient.den});
        }
        else
        {
            return detail::overflow_gate(quotient);
        }
    }
    else
    {
        if (quotient.den < 0)
        {
            return {-quotient.num, -quotient.den};
        }
        else
        {
            return quotient;
        }
    }
}

[[nodiscard]] constexpr rat operator+(rat lhs, rat rhs) noexcept
{
    int n = lhs.num * rhs.den + rhs.num * lhs.den;
    if (n == 0)
    {
        return zero;
    }
    else
    {
        int d = lhs.den * rhs.den;
        if (d > 1)
        {
            return detail::overflow_gate(rat{n, d});
        }
        else
        {
            return {n, d};
        }
    }
}

[[nodiscard]] constexpr auto operator-(rat lhs, rat rhs) noexcept
{
    return lhs + -rhs;
}

// This is a relatively specialized operation that returns the next integer that has an even number
// of bits set. The input is assumed to also have an even number of bits set.
[[nodiscard]] constexpr uint32_t next_even(uint32_t in) noexcept
{
    ++in;

    auto pop = pop_count(in);
    while (pop % 2 == 1)
    {
        ++in;
        pop = pop_count(in);
    }
    return in;
}
} // namespace gal
