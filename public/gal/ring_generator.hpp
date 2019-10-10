#pragma once

#include <cstddef>
#include <limits>
#include <numeric>
#include <tuple>

namespace gal
{
template <typename T1, typename T2>
struct equals
{
    constexpr static bool value = false;
};

template <typename T1, typename T2>
struct less
{
    constexpr static bool value = false;
};

// A tag identifies a geometric object and the index of the parameter that defines it.
// For example, a point in E3 might by associated with a unique tuple of (x, y, z), where
// each coordinate is associated with the tags <ID, 0>, <ID, 1>, and <ID, 2> where ID is
// the identifier of the point in question.
// A tag with id ~0ull compares != to all other tags
template <size_t S = ~0ull, size_t I = 0>
struct tag
{
    constexpr static size_t id  = S;
    constexpr static size_t index = I;
    constexpr static bool untagged = S == ~0ull;
};

template <size_t S1, size_t I1, size_t S2, size_t I2>
struct equals<tag<S1, I1>, tag<S2, I2>>
{
    constexpr static bool value = S1 != ~0ull && S2 != ~0ull && S1 == S2 && I1 == I2;
};

// Careful. This trait is implemented with the expectation that equality is checked first
template <size_t S1, size_t I1, size_t S2, size_t I2>
struct less<tag<S1, I1>, tag<S2, I2>>
{
    constexpr static bool value = S1 < S2 || (S1 == S2 && I1 < I2);
};

template <int D>
struct degree
{
    constexpr static int value = D;
};

// Although the multivector space will ultimately be defined over a field, we decompose the field into the
// product of scalars (essentially factoring out a free module). The free module over the ring of
// integers has the nice property that we can condense computation by performing arithmetic exactly
// at compiole time.
// The free module we factor out is a bimodule (i.e. there is no preference for left or right
// multiplication by the scalar).
// A generator encodes its degree in the monomial, as well as its identifier (if available).
// NOTE: "Degree" here is meant in the sense of a polynimal (e.g. x^2 has degree 2).
// We permit negative degrees to allow expressing linear combinations of nth-roots as well.
// The order of a generator is the positive integer n such that g^k = 0 for all k >= n
// The dual unit in particular has order 2. An order of "0" here, by convention, refers to
// an infinite order.
template <typename Tag, typename Degree = degree<1>, size_t Order = 0>
struct generator
{
    using tag_t                   = Tag;
    constexpr static int degree   = Degree::value;
    constexpr static bool is_zero = false;
};

// The derived generator is an expression that is dependent on other generators and indicates to the engine that the
// result contained should be memoized
template <typename T>
struct derived_generator
{
    T value;
};

// This struct exists purely as a tag to indicate that a generator has been annihilated
struct zero_generator
{
    constexpr static bool is_zero = true;
};

// NOTE: the tag of ~0 - 1 is RESERVED by the library for the dual unit. This ensures that the
// dual unit (which has a high chance of extinguishing the monomial it's a factor of) comes
// first in the factor ordering.
using dual_generator = generator<tag<~0ull - 1>, degree<1>, 2>;

// The generators that make up a monomial are weakly ordered based on the source identifiers.
// If all generators are identified (ID != ~0ull), the ordering becomes a total order.
template <typename T1, typename D1, size_t O1, typename T2, typename D2, size_t O2>
struct equals<generator<T1, D1, O1>, generator<T2, D2, O2>>
{
    constexpr static bool value = equals<T1, T2>::value && D1::value == D2::value && O1 == O2;
};

template <typename T1, typename D1, size_t O1, typename T2, typename D2, size_t O2>
struct less<generator<T1, D1, O1>, generator<T2, D2, O2>>
{
    constexpr static bool value = less<T1, T2>::value || (equals<T1, T2>::value && D1::value < D2::value);
};

// The module we work with is attached to the field of rational numbers.
// The numerator and denominator are left as ints (even though D > 0 is an invariant) so the compiler can help detect
// overflows
template <int N, int D = 1>
struct rational
{
    static_assert(D != 0, "Indeterminant rational encountered");
    constexpr static int num      = N;
    constexpr static int den      = D;
    constexpr static bool is_zero = N == 0;

    using neg_t = rational<-N, D>;

    template <typename T>
    [[nodiscard]] constexpr static T convert() noexcept
    {
        return static_cast<T>(N) / static_cast<T>(D);
    }
};

// Common rational types used for brevity
using one            = rational<1>;
using minus_one      = rational<-1>;
using one_half       = rational<1, 2>;
using minus_one_half = rational<-1, 2>;
using zero           = rational<0>;

namespace detail
{
    template <int N, int D>
    [[nodiscard]] constexpr auto overflow_gate(rational<N, D> in)
    {
        // As expressions expand, there may be cases where terms becoming vanishingly small or N and D become great
        // enough to risk overflow. This function is a pure function which deterministically nudges the result so that
        // compilation stays fast without sacrificing too much precision.

        // The goal is to at least match the precision that would have been afforded with floating point precision for
        // the range of numbers dealt with.

        // If we call the error introduced e, let's choose e to be 1e-7. This is only slightly worse than the precision
        // afforded by floating point numbers in the single digit ranges. Any fraction of the form p/q is the mediant of
        // two fractions which bound it from above and below. A possible choice is floor(p/2)/ceil(q/2) and
        // ceil(p/2)/floor(q/2). Suppose p is even. Then floor(p/2) = ceil(p/2) = p/2. Assuming p/q is irreducible, this
        // implies that floor(q/2) = (q-1)/2 and ceil(q/2) = (q + 1)/2. The error of the lower bound then is p/q - p/(q
        // + 1) = p/(q^2 + 1). Conversely, ceil(p/2)/floor(q/2) for even p reduces to p/(q - 1) and the error is given
        // by p/(q^2 - 1). For large q, the error then is approximately p/q^2.

        // We need to account for erroneously introducing systematic bias into the expression. This occurs if we
        // consistently perturb the denominator up or down. To balance things out, the direction of the perturbation
        // will be based on the parity of the second least significant bit of the denominator. This is by no means
        // perfect given that the distribution of rationals encountered is itself biased, but cheap to compute and
        // difficult to beat without introducing unreasonable amounts of complexity and compilation time.

        if constexpr (D < (1 << 10))
        {
            constexpr auto gcd = std::gcd(N, D);
            if constexpr (gcd > 1)
            {
                return rational<N / gcd, D / gcd>{};
            }
            else
            {
                return in;
            }
        }
        else
        {
            constexpr auto gcd = std::gcd(N, D);
            if constexpr (gcd > 1)
            {
                return rational<N / gcd, D / gcd>{};
            }
            else
            {
                constexpr double n_d     = static_cast<double>(N);
                constexpr double d_d     = static_cast<double>(D);
                constexpr double epsilon = std::abs(n_d / d_d / d_d);
                constexpr double frac    = n_d / d_d;

                if constexpr (std::abs(frac) < 1e-7)
                {
                    return zero{};
                }
                else if constexpr (epsilon < 1e-7)
                {
                    // The rational is the mediant of two fractions with smaller denominators. Pick one of them based on
                    // the parity.
                    if constexpr (N % 2 == 1)
                    {
                        // Perturb to the left-side of the mediant
                        constexpr auto N2 = (N - 1) / 2;
                        constexpr auto D2 = D / 2;
                        constexpr auto gcd = std::gcd(N2, D2);
                        if constexpr (gcd > 1)
                        {
                            return rational<N2 / gcd, D2 / gcd>{};
                        }
                        else
                        {
                            return rational<N2, D2>{};
                        }
                    }
                    else
                    {
                        // Perturb to the right-side of the mediant
                        constexpr auto N2 = N / 2;
                        constexpr auto D2 = (D - 1) / 2;
                        constexpr auto gcd = std::gcd(N2, D2);
                        if constexpr (gcd > 1)
                        {
                            return rational<N2 / gcd, D2 / gcd>{};
                        }
                        else
                        {
                            return rational<N2, D2>{};
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

template <int N1, int D1, int N2, int D2>
[[nodiscard]] constexpr auto operator==(rational<N1, D1>, rational<N2, D2>) noexcept
{
    constexpr auto gcd1 = std::gcd(N1, D1);
    constexpr auto gcd2 = std::gcd(N2, D2);
    return N1 / gcd1 == N2 / gcd2 && D1 / gcd1 == D2 / gcd2;
}

template <int N, int D>
[[nodiscard]] constexpr auto operator-(rational<N, D>) noexcept
{
    return rational<-N, D>{};
}

template <int N1, int D1, int N2, int D2>
[[nodiscard]] constexpr auto operator*(rational<N1, D1>, rational<N2, D2>)noexcept
{
    constexpr auto N = N1 * N2;
    constexpr auto D = D1 * D2;
    if constexpr (N == 0)
    {
        return zero{};
    }
    else if constexpr (D > 1)
    {
        // Periodically reduce the fraction if possible to prevent overflow
        return detail::overflow_gate(rational<N, D>{});
    }
    else
    {
        return rational<N, D>{};
    }
}

template <int N1, int D1, int N2, int D2>
[[nodiscard]] constexpr auto operator/(rational<N1, D1>, rational<N2, D2>) noexcept
{
    constexpr auto N = N1 * D2;
    constexpr auto D = N2 * D1;
    if constexpr (N == 0)
    {
        return zero{};
    }
    else if constexpr (D > 1)
    {
        if constexpr (D < 0)
        {
            return detail::overflow_gate(rational<-N, -D>{});
        }
        else
        {
            return detail::overflow_gate(rational<N, D>{});
        }
    }
    else
    {
        if constexpr (D < 0)
        {
            return rational<-N, -D>{};
        }
        else
        {
            return rational<N, D>{};
        }
    }
}

template <int N1, int D1, int N2, int D2>
[[nodiscard]] constexpr auto operator+(rational<N1, D1>, rational<N2, D2>) noexcept
{
    constexpr auto N = N1 * D2 + N2 * D1;
    if constexpr (N == 0)
    {
        return zero{};
    }
    else
    {
        constexpr auto D = D1 * D2;
        if constexpr (D > 1)
        {
            return detail::overflow_gate(rational<N, D>{});
        }
        else
        {
            return rational<N, D>{};
        }
    }
}

template <int N1, int D1, int N2, int D2>
[[nodiscard]] constexpr auto operator-(rational<N1, D1> lhs, rational<N2, D2> rhs) noexcept
{
    return lhs + -rhs;
}

} // namespace gal
