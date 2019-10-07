#pragma once

#include <cstddef>
#include <limits>
#include <numeric>

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
    constexpr static bool value = S1 < S2 || I1 < I2;
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
    constexpr static bool value = less<T1, T2>::value || D1::value < D2::value;
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
using one       = rational<1>;
using minus_one = rational<-1>;
using one_half  = rational<1, 2>;
using minus_one_half = rational<-1, 2>;
using zero      = rational<0>;

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
    else if constexpr (D != 1 && N >= (std::numeric_limits<int>::max() >> 2))
    {
        // Periodically reduce the fraction if possible to prevent overflow
        constexpr auto gcd = std::gcd(N, D);
        return rational<N / gcd, D / gcd>{};
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
    else if constexpr (D != 1 && N >= (std::numeric_limits<int>::max() >> 2))
    {
        constexpr auto gcd = std::gcd(N, D);
        if constexpr (D < 0)
        {
            return rational<-N / gcd, -D / gcd>{};
        }
        else
        {
            return rational<N / gcd, D / gcd>{};
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
        if constexpr (D != 1 && N >= (std::numeric_limits<int>::max() >> 2))
        {
            constexpr auto gcd = std::gcd(N, D);
            return rational<N / gcd, D / gcd>{};
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
