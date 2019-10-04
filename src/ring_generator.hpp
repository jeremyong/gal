#pragma once

#include <cstddef>

namespace gal
{
template <typename T1, typename T2> struct equals
{
    constexpr static bool value = false;
};

template <typename T1, typename T2> struct less
{
    constexpr static bool value = false;
};

// A tag with id ~0ull compares != to all other tags
template <size_t S = ~0ull> struct tag
{
    constexpr static size_t value = S;
    constexpr static bool untagged = S == ~0ull;
};

template <size_t S1, size_t S2> struct equals<tag<S1>, tag<S2>>
{
    constexpr static bool value = S1 != ~0ull && S2 != ~0ull && S1 == S2;
};

template <size_t S1, size_t S2> struct less<tag<S1>, tag<S2>>
{
    constexpr static bool value = S1 < S2;
};

template <int D> struct degree
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
template <typename Tag, typename Degree = degree<1>, size_t Order = 0> struct generator
{
    using tag_t = Tag;
    constexpr static int degree = Degree::value;
    constexpr static bool is_zero = false;
};

// This struct exists purely as a tag to indicate that a generator has been annihilated
struct zero_generator
{
    constexpr static bool is_zero = true;
};

// NOTE: the tag of 0 is RESERVED by the library for the dual unit. This ensures that the
// dual unit (which has a high chance of extinguishing the monomial it's a factor of) comes
// first in the factor ordering.
using dual_generator = generator<tag<0>, degree<1>, 2>;

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
} // namespace gal
