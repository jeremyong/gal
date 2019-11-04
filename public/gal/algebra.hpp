#pragma once

#include "algorithm.hpp"
#include "numeric.hpp"

#include <array>

namespace gal
{
using width_t = std::uint_fast32_t;
using elem_t  = std::uint_fast8_t;

namespace detail
{
    enum ind_constant_id : width_t
    {
        ind_constant_start = ~0u - 128,
        pi_ind             = ~0u - 128,
        e_ind,
    };

    template <typename T>
    constexpr inline T ind_constants[] = {M_PI, 2.71828182845904523536};
} // namespace detail

// Although the multivector space will ultimately be defined over a field, we decompose the
// field into the product of scalars (essentially factoring out a free module). The free module
// over the ring of integers has the nice property that we can condense computation by
// performing arithmetic exactly at compile time. The free module we factor out is a bimodule
// (i.e. there is no preference for left or right multiplication by the scalar). An
// indeterminate encodes its degree in the monomial, as well as its identifier (if available).
// NOTE: "Degree" here is meant in the sense of a polynomial/monomial degree (e.g. x^2 has
// degree 2 and x*y^2*z has degree 3). We permit negative degrees to allow expressing linear
// combinations of nth-roots as well. The order of an indeterminate is the positive integer n
// such that g^k = 0 for all k >= n The dual unit in particular has order 2. An order of "0"
// here, by convention, refers to an infinite order.
struct ind
{
    width_t id{~0u};
    rat degree;
};

// The indeterminates that make up a monomial are weakly ordered based on the source identifiers.
// If all indeterminates are identified (ID != ~0ull), the ordering becomes a total order.
[[nodiscard]] constexpr bool operator==(ind lhs, ind rhs) noexcept
{
    return lhs.id == rhs.id && lhs.degree == rhs.degree;
}

[[nodiscard]] constexpr bool operator!=(ind lhs, ind rhs) noexcept
{
    return lhs.id != rhs.id || lhs.degree != rhs.degree;
}

[[nodiscard]] constexpr bool operator<(ind lhs, ind rhs) noexcept
{
    return lhs.id < rhs.id || (lhs.id == rhs.id && lhs.degree < rhs.degree);
}

[[nodiscard]] constexpr ind operator^(ind lhs, int rhs) noexcept
{
    lhs.degree *= rat{rhs};
    return lhs;
}

// TRICK
// The shorter type names are intentional for producing less verbose type signatures in
// compile-errors and diagnostics.
struct mon
{
    rat q;
    rat degree; // Sum of all exponents of indeterminates
    width_t count      = 0;
    width_t ind_offset = 0;
};

// Temporal structure useful for sorting monomials in place
struct mon_view
{
    mon m;
    ind* ind_begin = nullptr;
};

[[nodiscard]] constexpr bool operator==(mon_view const& lhs, mon_view const& rhs) noexcept
{
    if (lhs.m.degree != rhs.m.degree || lhs.m.count != rhs.m.count)
    {
        return false;
    }
    else
    {
        for (width_t i = 0; i != lhs.m.count; ++i)
        {
            auto const& lhs_ind = *(lhs.ind_begin + i);
            auto const& rhs_ind = *(rhs.ind_begin + i);
            if (lhs_ind != rhs_ind)
            {
                return false;
            }
        }

        return true;
    }
}

[[nodiscard]] constexpr bool operator<(mon_view const& lhs, mon_view const& rhs) noexcept
{
    if (lhs.m.degree < rhs.m.degree)
    {
        return true;
    }
    else if (lhs.m.degree > rhs.m.degree)
    {
        return false;
    }
    else
    {
        auto min_ind = std::min(lhs.m.count, rhs.m.count);
        for (width_t i = 0; i != min_ind; ++i)
        {
            auto const& lhs_ind = *(lhs.ind_begin + i);
            auto const& rhs_ind = *(rhs.ind_begin + i);
            if (lhs_ind < rhs_ind)
            {
                return true;
            }
            else if (rhs_ind < lhs_ind)
            {
                return false;
            }
        }

        if (lhs.m.count == rhs.m.count)
        {
            // The monomials compare exactly equal
            return false;
        }
        else if (min_ind == lhs.m.count)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

struct term
{
    width_t count      = 0;
    width_t mon_offset = 0;
    uint32_t element   = 0;
};

[[nodiscard]] constexpr bool operator<(term const& lhs, term const& rhs) noexcept
{
    return lhs.element < rhs.element;
}

[[nodiscard]] constexpr bool operator==(term const& lhs, term const& rhs) noexcept
{
    return lhs.element == rhs.element;
}

// B := indicates if the size has been initialized
// I := number of indeterminates
// M := number of monomials
// T := number of terms
template <bool B = false, width_t I = 0, width_t M = 0, width_t T = 0>
struct mv_size_t
{
    constexpr static bool initialized  = B;
    constexpr static width_t ind_size  = I;
    constexpr static width_t mon_size  = M;
    constexpr static width_t term_size = T;
};

constexpr inline mv_size_t<> null_size_v;

struct mv_size
{
    width_t ind  = 0;
    width_t mon  = 0;
    width_t term = 0;
};

struct mon_it
{
    ind* ind_begin_;
    mon* mon_;

    [[nodiscard]] constexpr mon& operator*() noexcept
    {
        return *mon_;
    }

    [[nodiscard]] constexpr mon* operator->() noexcept
    {
        return mon_;
    }

    // NOTE: it is not checked that these two iterators refer to the same multivector
    [[nodiscard]] constexpr size_t operator-(mon_it const& other) const noexcept
    {
        return mon_ - other.mon_;
    }

    [[nodiscard]] constexpr size_t operator-(mon const* other) const noexcept
    {
        return mon_ - other;
    }

    [[nodiscard]] constexpr bool operator==(mon_it const& other) const noexcept
    {
        return ind_begin_ == other.ind_begin_ && mon_ == other.mon_;
    }

    [[nodiscard]] constexpr bool operator!=(mon_it const& other) const noexcept
    {
        return ind_begin_ != other.ind_begin_ || mon_ != other.mon_;
    }

    constexpr mon_it& operator++() noexcept
    {
        ++mon_;
        return *this;
    }

    constexpr mon_it operator++(int) noexcept
    {
        return {ind_begin_, mon_++};
    }

    constexpr mon_it& operator+=(size_t jump) noexcept
    {
        mon_ += jump;
        return *this;
    }

    [[nodiscard]] constexpr ind* begin() noexcept
    {
        return ind_begin_ + mon_->ind_offset;
    }

    [[nodiscard]] constexpr ind* end() noexcept
    {
        return ind_begin_ + mon_->ind_offset + mon_->count;
    }
};

struct const_mon_it
{
    ind const* ind_begin_;
    mon const* mon_;

    [[nodiscard]] constexpr mon const& operator*() const noexcept
    {
        return *mon_;
    }

    [[nodiscard]] constexpr mon const* operator->() const noexcept
    {
        return mon_;
    }

    // NOTE: it is not checked that these two iterators refer to the same multivector
    [[nodiscard]] constexpr size_t operator-(const_mon_it const& other) const noexcept
    {
        return mon_ - other.mon_;
    }

    [[nodiscard]] constexpr size_t operator-(mon const* other) const noexcept
    {
        return mon_ - other;
    }

    [[nodiscard]] constexpr bool operator==(const_mon_it const& other) const noexcept
    {
        return ind_begin_ == other.ind_begin_ && mon_ == other.mon_;
    }

    [[nodiscard]] constexpr bool operator!=(const_mon_it const& other) const noexcept
    {
        return ind_begin_ != other.ind_begin_ || mon_ != other.mon_;
    }

    constexpr const_mon_it& operator++() noexcept
    {
        ++mon_;
        return *this;
    }

    constexpr const_mon_it operator++(int) noexcept
    {
        return {ind_begin_, mon_++};
    }

    constexpr const_mon_it& operator+=(size_t jump) noexcept
    {
        mon_ += jump;
        return *this;
    }

    [[nodiscard]] constexpr ind const* cbegin() const noexcept
    {
        return ind_begin_ + mon_->ind_offset;
    }

    [[nodiscard]] constexpr ind const* cend() const noexcept
    {
        return ind_begin_ + mon_->ind_offset + mon_->count;
    }
};

struct term_it
{
    ind* ind_begin_;
    mon* mon_begin_;
    term* term_;

    [[nodiscard]] constexpr term& operator*() const noexcept
    {
        return *term_;
    }

    [[nodiscard]] constexpr term* operator->() const noexcept
    {
        return term_;
    }

    // NOTE: it is not checked that these two iterators refer to the same multivector
    [[nodiscard]] constexpr size_t operator-(term_it const& other) const noexcept
    {
        return term_ - other.term_;
    }

    [[nodiscard]] constexpr bool operator==(term_it const& other) const noexcept
    {
        return term_ == other.term_;
    }

    [[nodiscard]] constexpr bool operator!=(term_it const& other) const noexcept
    {
        return term_ != other.term_;
    }

    constexpr term_it& operator++() noexcept
    {
        ++term_;
        return *this;
    }

    constexpr term_it operator++(int) noexcept
    {
        return {ind_begin_, mon_begin_, term_++};
    }

    constexpr term_it& operator+=(size_t jump) noexcept
    {
        term_ += jump;
        return *this;
    }

    [[nodiscard]] constexpr mon_it begin() noexcept
    {
        return {ind_begin_, mon_begin_ + term_->mon_offset};
    }

    [[nodiscard]] constexpr mon_it end() noexcept
    {
        return {ind_begin_, mon_begin_ + term_->mon_offset + term_->count};
    }
};

struct const_term_it
{
    ind const* ind_begin_;
    mon const* mon_begin_;
    term const* term_;

    [[nodiscard]] constexpr term const& operator*() const noexcept
    {
        return *term_;
    }

    [[nodiscard]] constexpr term const* operator->() const noexcept
    {
        return term_;
    }

    // NOTE: it is not checked that these two iterators refer to the same multivector
    [[nodiscard]] constexpr size_t operator-(const_term_it const& other) const noexcept
    {
        return term_ - other.term_;
    }

    [[nodiscard]] constexpr bool operator==(const_term_it const& other) const noexcept
    {
        return term_ == other.term_;
    }

    [[nodiscard]] constexpr bool operator!=(const_term_it const& other) const noexcept
    {
        return term_ != other.term_;
    }

    constexpr const_term_it& operator++() noexcept
    {
        ++term_;
        return *this;
    }

    constexpr const_term_it operator++(int) noexcept
    {
        return {ind_begin_, mon_begin_, term_++};
    }

    constexpr const_term_it& operator+=(size_t jump) noexcept
    {
        term_ += jump;
        return *this;
    }

    [[nodiscard]] constexpr const_mon_it cbegin() const noexcept
    {
        return {ind_begin_, mon_begin_ + term_->mon_offset};
    }

    [[nodiscard]] constexpr const_mon_it cend() const noexcept
    {
        return {ind_begin_, mon_begin_ + term_->mon_offset + term_->count};
    }
};

template <width_t IndMax, width_t MonMax, width_t TermMax>
struct dense_mv
{
    std::array<std::pair<term, std::array<std::pair<mon, std::array<ind, IndMax>>, MonMax>>, TermMax> data;
};

enum class mv_op : uint32_t
{
    id,
    sin,
    cos,
    tan,
    sqrt,
};

// Multivector representation, intended to be a compile-time representation
// A := Algebra
// I := Indeterminate capacity
// M := Monomial capacity
// T := Term capacity
// The mv struct supports nested iteration. For example:
//
//     for (auto&& t : mv)
//     {
//         // t is a term in mv
//         for (auto&& m : t)
//         {
//             // m is a monomial in t
//             for (auto&& i : m)
//             {
//                 // i is an indeterminate in m
//             }
//         }
//     }
// In the snippet above, dereferencing any of t, m, or g will result in the reference to the object
// in question.
template <typename A, width_t I, width_t M, width_t T>
struct mv
{
    using algebra_t = A;

    [[nodiscard]] constexpr static width_t ind_capacity() noexcept
    {
        return I;
    }

    [[nodiscard]] constexpr static width_t mon_capacity() noexcept
    {
        return M;
    }

    [[nodiscard]] constexpr static width_t term_capacity() noexcept
    {
        return T;
    }

    mv_size size;
    std::array<ind, I> inds;
    std::array<mon, M> mons;
    std::array<term, T> terms;
    mv_op o{mv_op::id};

    // For evaluation, it is often convenient to fully evaluate a multivector and refer to it later.
    // This creates a contracted form referencing the appropriate elements of this multivector
    // (which will presumably be evaluated first).
    constexpr mv<A, T, T, T> create_ref(uint32_t id) const noexcept
    {
        mv<A, T, T, T> out{};
        out.size = mv_size{T, T, T};
        for (width_t i = 0; i != T; ++i)
        {
            out.inds[i]  = ind{id + i, one};
            out.mons[i]  = mon{one, one, 1, i};
            out.terms[i] = term{1, i, terms[i].element};
        }
        return out;
    }

    constexpr void scale(rat q) noexcept
    {
        for (size_t i = 0; i != size.mon; ++i)
        {
            mons[i].q = q * mons[i].q;
        }
    }

    // The transcendental operations here are applied to the final reified value.
    constexpr void sin(rat q) noexcept
    {
        o = mv_op::sin;
        scale(q);
    }

    constexpr void cos(rat q) noexcept
    {
        o = mv_op::cos;
        scale(q);
    }

    constexpr void tan(rat q) noexcept
    {
        o = mv_op::tan;
        scale(q);
    }

    constexpr void sqrt(rat q) noexcept
    {
        o = mv_op::sqrt;
        scale(q);
        // TODO:
        // Check if we can take the sqrt of the scalar multiplier

        // WARNING: this function is only defined for a multivector with a single indeterminate
        // value for a single term and monomial.
        // inds[0].degree *= one_half;
        // mons[0].degree *= one_half;
    }

    // Select a single component and emit it as a scalar
    constexpr mv<A, I, M, 1> operator[](elem_t e) const noexcept
    {
        mv<A, I, M, 1> out{};
        for (auto it = cbegin(); it != cend(); ++it)
        {
            if (it->element == e)
            {
                out.push(it, one, 0);
                out.size.term = 1;
                return out;
            }
        }
        return out;
    }

    // Push a term onto this multivector with a scaling factor and change of element
    constexpr void push(const_term_it it, rat scale, elem_t e) noexcept
    {
        auto out_mons_it = mons.begin() + size.mon;

        for (auto mon_it = it.cbegin(); mon_it != it.cend(); ++mon_it)
        {
            auto out_inds_it = inds.begin() + size.ind;
            for (auto ind_it = mon_it.cbegin(); ind_it != mon_it.cend(); ++ind_it)
            {
                *out_inds_it++ = *ind_it;
            }
            *out_mons_it            = *mon_it;
            out_mons_it->q          = out_mons_it->q * scale;
            out_mons_it->ind_offset = size.ind;
            ++out_mons_it;
            size.ind = static_cast<width_t>(out_inds_it - inds.begin());
        }

        terms[size.term]            = *it;
        terms[size.term].element    = e;
        terms[size.term].mon_offset = size.mon;
        size.mon                    = static_cast<width_t>(out_mons_it - mons.begin());
        ++size.term;
    }

    template <width_t I2, width_t M2, width_t T2>
    [[nodiscard]] constexpr auto resize() const noexcept
    {
        if constexpr (I == I2 && M == M2 && T == T2)
        {
            return *this;
        }
        else
        {
            mv<A, I2, M2, T2> out{};
            out.size = size;
            for (size_t i = 0; i != out.size.ind; ++i)
            {
                out.inds[i] = inds[i];
            }
            for (size_t i = 0; i != out.size.mon; ++i)
            {
                out.mons[i] = mons[i];
            }
            for (size_t i = 0; i != out.size.term; ++i)
            {
                out.terms[i] = terms[i];
            }
            return out;
        }
    }

    // Reshape the multivector to be fully dense in the arrays to simplify compilation of the final
    // reduction
    template <width_t IndMax, width_t MonMax, width_t TermMax>
    [[nodiscard]] constexpr auto regularize() const noexcept
    {
        dense_mv<IndMax, MonMax, TermMax> out;

        for (width_t i = 0; i != size.term; ++i)
        {
            auto const& term = terms[i];
            auto& out_term   = out.data[i];
            out_term.first   = term;

            for (width_t j = term.mon_offset; j != term.mon_offset + term.count; ++j)
            {
                auto const& mon = mons[j];
                auto& out_mon   = out_term.second[j - term.mon_offset];
                out_mon.first   = mon;

                for (width_t k = mon.ind_offset; k != mon.ind_offset + mon.count; ++k)
                {
                    out_mon.second[k - mon.ind_offset] = inds[k];
                }
            }
        }

        return out;
    }

    [[nodiscard]] constexpr mv_size extent() const noexcept
    {
        // Determine the maximum number of indeterminates across all monomials and the maximum
        // number of monomials across all terms

        width_t mon_count = 0;
        width_t ind_count = 0;
        for (auto term_it = cbegin(); term_it != cend(); ++term_it)
        {
            if (term_it->count > mon_count)
            {
                mon_count = term_it->count;
            }

            for (auto mon_it = term_it.cbegin(); mon_it != term_it.cend(); ++mon_it)
            {
                if (mon_it->count > ind_count)
                {
                    ind_count = mon_it->count;
                }
            }
        }
        return {ind_count, mon_count, size.term};
    }

    [[nodiscard]] constexpr const_term_it cbegin() const noexcept
    {
        return {inds.data(), mons.data(), terms.data()};
    }

    [[nodiscard]] constexpr term_it begin() noexcept
    {
        return {inds.data(), mons.data(), terms.data()};
    }

    [[nodiscard]] constexpr const_term_it cend() const noexcept
    {
        return {inds.data(), mons.data(), terms.data() + size.term};
    }

    [[nodiscard]] constexpr term_it end() noexcept
    {
        return {inds.data(), mons.data(), terms.data() + size.term};
    }
};

namespace detail
{
    template <typename A, width_t I1, width_t M1, width_t T1, width_t I2, width_t M2, width_t T2>
    constexpr auto sum(mv<A, I1, M1, T1> const& lhs, mv<A, I2, M2, T2> const& rhs) noexcept
    {
        constexpr size_t ind_size  = I1 + I2;
        constexpr size_t mon_size  = M1 + M2;
        constexpr size_t term_size = T1 + T2;
        mv<A, ind_size, mon_size, term_size> out{};

        auto lhs_it      = lhs.cbegin();
        auto rhs_it      = rhs.cbegin();
        auto lhs_end     = lhs.cend();
        auto rhs_end     = rhs.cend();
        auto out_ind_it  = out.inds.begin();
        auto out_mon_it  = out.mons.begin();
        auto out_term_it = out.terms.begin();

        while (true)
        {
            if (lhs_it == lhs_end && rhs_it == rhs_end)
            {
                // Base case
                out.size = mv_size{static_cast<width_t>(out_ind_it - out.inds.begin()),
                                   static_cast<width_t>(out_mon_it - out.mons.begin()),
                                   static_cast<width_t>(out_term_it - out.terms.begin())};
                return out;
            }
            else if (lhs_it == lhs_end || rhs_it == rhs_end)
            {
                auto& it  = lhs_it == lhs_end ? rhs_it : lhs_it;
                auto& end = lhs_it == lhs_end ? rhs_end : lhs_end;

                for (; it != end; ++it)
                {
                    *out_term_it            = *it;
                    out_term_it->mon_offset = out_mon_it - out.mons.begin();
                    ++out_term_it;

                    for (auto m_it = it.cbegin(); m_it != it.cend(); ++m_it)
                    {
                        *out_mon_it            = *m_it;
                        out_mon_it->ind_offset = out_ind_it - out.inds.begin();
                        ++out_mon_it;

                        for (auto i_it = m_it.cbegin(); i_it != m_it.cend(); ++i_it)
                        {
                            *out_ind_it++ = *i_it;
                        }
                    }
                }
                out.size = mv_size{static_cast<width_t>(out_ind_it - out.inds.begin()),
                                   static_cast<width_t>(out_mon_it - out.mons.begin()),
                                   static_cast<width_t>(out_term_it - out.terms.begin())};
                return out;
            }
            else if (rhs_it->element == lhs_it->element)
            {
                auto lhs_mon_it  = lhs_it.cbegin();
                auto lhs_mon_end = lhs_it.cend();
                auto rhs_mon_it  = rhs_it.cbegin();
                auto rhs_mon_end = rhs_it.cend();
                auto mon_cursor  = out_mon_it;

                while (true) // Graded lexicographic comparison
                {
                    if (lhs_mon_it == lhs_mon_end && rhs_mon_it == rhs_mon_end)
                    {
                        break;
                    }
                    else if (lhs_mon_it == lhs_mon_end || rhs_mon_it == rhs_mon_end)
                    {
                        auto& it  = lhs_mon_it == lhs_mon_end ? rhs_mon_it : lhs_mon_it;
                        auto& end = lhs_mon_it == lhs_mon_end ? rhs_mon_end : lhs_mon_end;
                        for (; it != end; ++it)
                        {
                            *out_mon_it            = *it;
                            out_mon_it->ind_offset = out_ind_it - out.inds.begin();
                            ++out_mon_it;
                            for (auto i_it = it.cbegin(); i_it != it.cend(); ++i_it)
                            {
                                *out_ind_it++ = *i_it;
                            }
                        }
                        break;
                    }
                    else
                    {
                        // The monomials need to be evaluated for equality
                        if (lhs_mon_it->degree != rhs_mon_it->degree)
                        {
                            auto& it
                                = lhs_mon_it->degree < rhs_mon_it->degree ? lhs_mon_it : rhs_mon_it;
                            *out_mon_it            = *it;
                            out_mon_it->ind_offset = out_ind_it - out.inds.begin();
                            ++out_mon_it;
                            for (auto i_it = it.cbegin(); i_it != it.cend(); ++i_it)
                            {
                                *out_ind_it++ = *i_it;
                            }
                            ++it;
                        }
                        else
                        {
                            // Monomial degrees match so a lexicographic comparison is needed
                            auto lhs_ind_it  = lhs_mon_it.cbegin();
                            auto lhs_ind_end = lhs_mon_it.cend();
                            auto rhs_ind_it  = rhs_mon_it.cbegin();
                            auto rhs_ind_end = rhs_mon_it.cend();
                            while (true)
                            {
                                if (lhs_ind_it == lhs_ind_end && rhs_ind_it == rhs_ind_end)
                                {
                                    // Important! When adding, we adopt the lhs rational scaling
                                    // factor so added terms from the rhs need to be divided by the
                                    // rhs rational scaling factor.
                                    auto q = lhs_mon_it->q + rhs_mon_it->q;

                                    if (!q.is_zero())
                                    {
                                        // The monomials compare equal and do not cancel
                                        *out_mon_it            = *lhs_mon_it;
                                        out_mon_it->ind_offset = out_ind_it - out.inds.begin();
                                        out_mon_it->q          = q;
                                        ++out_mon_it;

                                        for (auto i_it = lhs_mon_it.cbegin();
                                             i_it != lhs_mon_it.cend();
                                             ++i_it)
                                        {
                                            *out_ind_it++ = *i_it;
                                        }
                                    }
                                    ++lhs_mon_it;
                                    ++rhs_mon_it;
                                    break;
                                }
                                else if (lhs_ind_it == lhs_ind_end || rhs_ind_it == rhs_ind_end)
                                {
                                    auto& it = lhs_ind_it == lhs_ind_end ? lhs_mon_it : rhs_mon_it;
                                    *out_mon_it            = *it;
                                    out_mon_it->ind_offset = out_ind_it - out.inds.begin();
                                    ++out_mon_it;

                                    for (auto i_it = it.cbegin(); i_it != it.cend(); ++i_it)
                                    {
                                        *out_ind_it++ = *i_it;
                                    }
                                    ++it;
                                    break;
                                }
                                else
                                {
                                    auto const& lhs_ind = *lhs_ind_it;
                                    auto const& rhs_ind = *rhs_ind_it;
                                    if (lhs_ind == rhs_ind)
                                    {
                                        ++lhs_ind_it;
                                        ++rhs_ind_it;
                                    }
                                    else
                                    {
                                        auto& it    = lhs_ind < rhs_ind ? lhs_mon_it : rhs_mon_it;
                                        *out_mon_it = *it;
                                        out_mon_it->ind_offset = out_ind_it - out.inds.begin();
                                        ++out_mon_it;
                                        for (auto i_it = it.cbegin(); i_it != it.cend(); ++i_it)
                                        {
                                            *out_ind_it++ = *i_it;
                                        }
                                        ++it;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                } // Graded lexicographic comparison

                width_t mon_count = static_cast<width_t>(out_mon_it - mon_cursor);
                if (mon_count != 0)
                {
                    *out_term_it            = *lhs_it;
                    out_term_it->mon_offset = mon_cursor - out.mons.begin();
                    out_term_it->count      = mon_count;
                    ++out_term_it;
                }
                ++lhs_it;
                ++rhs_it;
            }
            else
            {
                auto& it        = lhs_it->element < rhs_it->element ? lhs_it : rhs_it;
                auto& end       = lhs_it->element < rhs_it->element ? lhs_end : rhs_end;
                auto mon_cursor = out_mon_it;

                for (auto mon_it = it.cbegin(); mon_it != it.cend(); ++mon_it)
                {
                    *out_mon_it            = *mon_it;
                    out_mon_it->ind_offset = out_ind_it - out.inds.begin();
                    ++out_mon_it;
                    for (auto i_it = mon_it.cbegin(); i_it != mon_it.cend(); ++i_it)
                    {
                        *out_ind_it++ = *i_it;
                    }
                }

                *out_term_it            = *it++;
                out_term_it->mon_offset = static_cast<width_t>(mon_cursor - out.mons.begin());
                ++out_term_it;
            }
        }
    }

    template <typename A, width_t I, width_t M, width_t T>
    [[nodiscard]] constexpr auto scale(rat lhs, mv<A, I, M, T> rhs) noexcept
    {
        rhs.scale(lhs);
        return rhs;
    }

    template <typename A, width_t I, width_t M, width_t T>
    [[nodiscard]] constexpr auto shift(rat lhs, mv<A, I, M, T> const& rhs) noexcept
    {
        mv<A, 0, 1, 1> addend{mv_size{0, 1, 1}, {}, {mon{lhs, zero, 0, 0}}, {term{1, 0, 0}}};
        return sum(rhs, addend);
    }

    template <typename T>
    [[nodiscard]] constexpr auto negate(T in) noexcept
    {
        for (auto it = in.begin(); it != in.end(); ++it)
        {
            for (auto& mon : it)
            {
                mon.q = -mon.q;
            }
        }
        return in;
    }

    // Specialized op for dividing a multivector by a scalar (in multivector form).
    template <typename A, width_t I, width_t M, width_t T, width_t I2, width_t M2, width_t T2>
    constexpr auto divide(mv<A, I, M, T> const& lhs, mv<A, I2, M2, T2> const& rhs, rat q) noexcept
    {
        static_assert(M2 != 0 && T2 != 0, "Divide by zero detected!");
        mon m = rhs.mons[0];
        if constexpr (I2 == 0)
        {
            // Divide each monomial of the LHS by the RHS monomial scaling constant
            auto out = lhs;
            for (auto it = out.begin(); it != out.end(); ++it)
            {
                for (auto mon_it = it.begin(); mon_it != it.end(); ++mon_it)
                {
                    mon_it->q = q * mon_it->q / m.q;
                }
            }
            return out;
        }
        else if constexpr (I2 == 1)
        {
            ind d    = rhs.inds[0];
            d.degree = d.degree.negation();
            mv<A, I + M, M, T> out;
            auto out_terms_it = out.terms.begin();
            auto out_mons_it  = out.mons.begin();
            auto out_inds_it  = out.inds.begin();

            // Multiply each monomial of the lhs by the reciprocal of the rhs.
            for (auto it = lhs.cbegin(); it != lhs.cend(); ++it)
            {
                for (auto mon_it = it.cbegin(); mon_it != it.cend(); ++mon_it)
                {
                    width_t ind_offset = out_inds_it - out.inds.begin();
                    bool canceled      = false;
                    bool written       = false;
                    for (auto ind_it = mon_it.cbegin(); ind_it != mon_it.cend(); ++ind_it)
                    {
                        if (!written)
                        {
                            if (d.id == ind_it->id)
                            {
                                rat deg = d.degree + ind_it->degree;
                                if (!deg.is_zero())
                                {
                                    canceled = true;
                                }
                                else
                                {
                                    *out_inds_it++ = ind{d.id, deg};
                                }
                                written = true;
                            }
                            else if (d.id < ind_it->id)
                            {
                                *out_inds_it++ = d;
                                *out_inds_it++ = *ind_it;
                                written        = true;
                            }
                            else
                            {
                                *out_inds_it++ = *ind_it;
                            }
                        }
                        else
                        {
                            *out_inds_it++ = *ind_it;
                        }
                    }

                    if (!written)
                    {
                        *out_inds_it++ = d;
                    }

                    rat deg        = mon_it->degree - d.degree;
                    width_t count  = canceled ? mon_it->count - 1 : mon_it->count + 1;
                    *out_mons_it++ = mon{mon_it->q / m.q, deg, count, ind_offset};
                }

                *out_terms_it++ = *it;
            }

            out.size.term = out_terms_it - out.terms.begin();
            out.size.mon  = out_mons_it - out.mons.begin();
            out.size.ind  = out_inds_it - out.inds.begin();

            out.scale(q);
            return out;
        }
    }

    template <typename T>
    [[nodiscard]] constexpr auto reverse(T in) noexcept
    {
        for (auto it = in.begin(); it != in.end(); ++it)
        {
            if (it->element == 0)
            {
                continue;
            }

            auto grade  = pop_count(it->element);
            auto parity = grade * (grade - 1) / 2;
            if (parity % 2 == 1)
            {
                for (auto& mon : it)
                {
                    mon.q = -mon.q;
                }
            }
        }
        return in;
    }

    [[nodiscard]] constexpr std::pair<uint32_t, int>
    poincare_complement(elem_t element, elem_t dim) noexcept
    {
        uint32_t complement = ((1 << dim) - 1) ^ element;

        uint32_t swaps = 0;
        uint32_t grade = pop_count(element);
        while (element > 0)
        {
            if ((element & 1) == 0)
            {
                swaps += grade;
            }
            else
            {
                --grade;
            }
            element >>= 1;
        }

        return {complement, swaps % 2 == 0 ? 1 : -1};
    }

    [[nodiscard]] constexpr auto collate(term* const start,
                                         term* const end,
                                         mon_view* const mon_start,
                                         ind* const ind_start,
                                         term* out_terms_it,
                                         mon* out_mons_it,
                                         ind* out_inds_it,
                                         mv_size& out_size) noexcept
    {
        auto out_terms_begin = out_terms_it;
        auto out_mons_begin  = out_mons_it;
        auto out_inds_begin  = out_inds_it;
        // We need to collate terms that refer to the same element together into one
        for (auto term_it = start; term_it != end;)
        {
            auto mon_cursor = out_mons_it;

            auto next = term_it;
            for (; next != end; ++next)
            {
                if (next->element == term_it->element)
                {
                    for (auto mon_it = mon_start + next->mon_offset;
                         mon_it != mon_start + next->mon_offset + next->count;
                         ++mon_it)
                    {
                        *out_mons_it++ = mon_it->m;
                    }
                }
                else
                {
                    // Non-collatable term
                    break;
                }
            }

            // Now, the monomials need to be sorted and collated
            sort(mon_cursor, out_mons_it, [ind_begin = ind_start](auto&& lhs, auto&& rhs) {
                return mon_view{lhs, ind_begin + lhs.ind_offset}
                       < mon_view{rhs, ind_begin + rhs.ind_offset};
            });

            // We *overwrite* monomials in the final output as we collate them (conserves memory)
            auto mon_writer = mon_cursor;

            for (auto mon_it = mon_cursor; mon_it != out_mons_it;)
            {
                rat q         = mon_it->q;
                auto next_mon = mon_it + 1;
                mon_view current_mon{*mon_it, ind_start + mon_it->ind_offset};
                for (; next_mon != out_mons_it; ++next_mon)
                {
                    if (current_mon == mon_view{*next_mon, ind_start + next_mon->ind_offset})
                    {
                        // Coincident monomials need to be added together
                        q = q + next_mon->q;
                    }
                    else
                    {
                        // Non-collatable monomial
                        break;
                    }
                }

                if (!q.is_zero())
                {
                    // CAREFUL! The mon_writer above may overwrite what mon_it points to
                    *mon_writer++ = mon{q,
                                        mon_it->degree,
                                        mon_it->count,
                                        static_cast<width_t>(out_inds_it - out_inds_begin)};

                    // Copy over indeterminates referenced by the monomial
                    for (auto ind_it = ind_start + current_mon.m.ind_offset;
                         ind_it != ind_start + current_mon.m.ind_offset + mon_it->count;
                         ++ind_it)
                    {
                        *out_inds_it++ = *ind_it;
                    }
                }

                mon_it = next_mon;
            }

            width_t mon_count = static_cast<width_t>(mon_writer - mon_cursor);
            if (mon_count > 0)
            {
                *out_terms_it++ = term{
                    mon_count, static_cast<width_t>(mon_cursor - out_mons_begin), term_it->element};
            }
            out_mons_it = mon_writer;
            // Set the iterator to the term after a repeated sequence of terms of the same element.
            term_it = next;
        }

        out_size = mv_size{static_cast<width_t>(out_inds_it - out_inds_begin),
                           static_cast<width_t>(out_mons_it - out_mons_begin),
                           static_cast<width_t>(out_terms_it - out_terms_begin)};
    }

    template <typename T>
    [[nodiscard]] constexpr auto poincare_dual(T const& in) noexcept
    {
        // Because the dual is not order preserving, we write the contents to a new multivector
        T out{};
        out.size = in.size;

        // Under the poincare dual map, the order of the terms will reverse
        auto out_term_it = out.terms.begin();
        auto out_mon_it  = out.mons.begin();
        auto out_ind_it  = out.inds.begin();

        for (auto it = in.cbegin(); it != in.cend(); ++it)
        {
            auto [g, parity] = poincare_complement(it->element, T::algebra_t::metric_t::dimension);
            *out_term_it++ = term{it->count, static_cast<width_t>(out_mon_it - out.mons.begin()), g};
            for (auto mon_it = it.cbegin(); mon_it != it.cend(); ++mon_it)
            {
                *out_mon_it++ = mon{parity * mon_it->q,
                                    mon_it->degree,
                                    mon_it->count,
                                    static_cast<width_t>(out_ind_it - out.inds.begin())};
                for (auto ind_it = mon_it.cbegin(); ind_it != mon_it.cend(); ++ind_it)
                {
                    *out_ind_it++ = *ind_it;
                }
            }
        }

        std::array<mon_view, T::mon_capacity()> mon_views;
        for (width_t i = 0; i != out.size.mon; ++i)
        {
            mon_views[i] = mon_view{out.mons[i], out.inds.begin()};
        }

        sort(out.terms.begin(), out.terms.begin() + out.size.term);

        T collated{};
        collate(out.terms.begin(),
                out.terms.begin() + out.size.term,
                mon_views.begin(),
                out.inds.begin(),
                collated.terms.begin(),
                collated.mons.begin(),
                collated.inds.begin(),
                collated.size);
        return collated;
    }

    // Given a specified product operation, compute the product between the lhs and the rhs.
    // If the size is not yet initialized, compute the size that would result from the
    // multiplication. Multiplication is always done left-to-right. P := product operation between
    // basis elements (returns a pair of a multiplier and target element)
    template <typename P, typename A, width_t I1, width_t M1, width_t T1, width_t I2, width_t M2, width_t T2>
    [[nodiscard]] constexpr auto
    product(P, mv<A, I1, M1, T1> const& lhs, mv<A, I2, M2, T2> const& rhs) noexcept
    {
        // The total number of terms conservatively is O(n*m) where n is the number of terms in the
        // lhs and m is the number of terms in the rhs. Note that this applies to both the number of
        // indeterminates and the number of monomials.
        constexpr width_t term_size = T1 * T2;
        constexpr width_t mon_size  = M1 * M2;
        constexpr width_t ind_size  = M1 * I2 + M2 * I1;

        // The monomials and indeterminates start out unsorted so we place them in temporary storage
        // first before the final sort-on-copy.
        std::array<ind, ind_size> temp_inds;
        std::array<mon_view, mon_size> temp_mons;
        std::array<term, term_size> temp_terms;
        auto temp_inds_it  = temp_inds.begin();
        auto temp_mons_it  = temp_mons.begin();
        auto temp_terms_it = temp_terms.begin();

        // It is NOT the case that the product will move distinct term pairs to distinct term
        // elements. For example, the inner or contraction product can easily collapse multiple
        // values to a scalar quantity.

        for (auto lhs_it = lhs.cbegin(); lhs_it != lhs.cend(); ++lhs_it)
        {
            for (auto rhs_it = rhs.cbegin(); rhs_it != rhs.cend(); ++rhs_it)
            {
                auto&& [element, multiplier] = P::product(lhs_it->element, rhs_it->element);
                if (multiplier != 0)
                {
                    auto mon_cursor = temp_mons_it;

                    // Multiply the polynomials of the lhs term and rhs term, scale it by the
                    // multiplier, and accumulate it into out. We do not bother to sort OR reduce as
                    // this will happen in the final pass.
                    for (auto lhs_mon = lhs_it.cbegin(); lhs_mon != lhs_it.cend(); ++lhs_mon)
                    {
                        for (auto rhs_mon = rhs_it.cbegin(); rhs_mon != rhs_it.cend(); ++rhs_mon)
                        {
                            // Merge the indeterminates of the lhs and rhs monomials.
                            auto lhs_ind_it  = lhs_mon.cbegin();
                            auto lhs_ind_end = lhs_mon.cend();
                            auto rhs_ind_it  = rhs_mon.cbegin();
                            auto rhs_ind_end = rhs_mon.cend();
                            auto ind_cursor  = temp_inds_it;
                            rat degree;

                            while (true)
                            {
                                if (lhs_ind_it == lhs_ind_end && rhs_ind_it == rhs_ind_end)
                                {
                                    break;
                                }
                                else if (lhs_ind_it == lhs_ind_end || rhs_ind_it == rhs_ind_end)
                                {
                                    auto& it = lhs_ind_it == lhs_ind_end ? rhs_ind_it : lhs_ind_it;
                                    auto& it_end
                                        = lhs_ind_it == lhs_ind_end ? rhs_ind_end : lhs_ind_end;
                                    for (; it != it_end; ++it)
                                    {
                                        degree += it->degree;
                                        *temp_inds_it++ = *it;
                                    }
                                    break;
                                }
                                else
                                {
                                    auto const& lhs_ind = *lhs_ind_it;
                                    auto const& rhs_ind = *rhs_ind_it;
                                    if (lhs_ind.id == rhs_ind.id)
                                    {
                                        rat next_degree = lhs_ind.degree + rhs_ind.degree;
                                        degree += next_degree;
                                        if (next_degree != 0)
                                        {
                                            // TODO: handle dual numbers
                                            *temp_inds_it++ = ind{lhs_ind.id, next_degree};
                                        }
                                        ++lhs_ind_it;
                                        ++rhs_ind_it;
                                    }
                                    else if (lhs_ind.id < rhs_ind.id)
                                    {
                                        *temp_inds_it++ = lhs_ind;
                                        degree += lhs_ind.degree;
                                        ++lhs_ind_it;
                                    }
                                    else
                                    {
                                        *temp_inds_it++ = rhs_ind;
                                        degree += rhs_ind.degree;
                                        ++rhs_ind_it;
                                    }
                                }
                            }

                            *temp_mons_it++
                                = mon_view{mon{rat{multiplier * lhs_mon->q * rhs_mon->q},
                                               degree,
                                               static_cast<width_t>(temp_inds_it - ind_cursor),
                                               static_cast<width_t>(ind_cursor - temp_inds.begin())},
                                           ind_cursor};
                        }
                    }

                    auto mon_count   = static_cast<width_t>(temp_mons_it - mon_cursor);
                    *temp_terms_it++ = term{
                        mon_count, static_cast<width_t>(mon_cursor - temp_mons.begin()), element};
                }
            }
        }

        // The product between terms is not necessarily order-preserving so we need to both sort
        // terms and monomials
        sort(temp_terms.begin(), temp_terms_it);

        mv<A, ind_size, mon_size, term_size> out{};
        collate(temp_terms.begin(),
                temp_terms_it,
                temp_mons.begin(),
                temp_inds.begin(),
                out.terms.begin(),
                out.mons.begin(),
                out.inds.begin(),
                out.size);
        return out;
    }

    template <typename A, width_t I, width_t M, width_t T>
    constexpr mv<A, I, M, 1> component(mv<A, I, M, T> const& in, elem_t e) noexcept
    {
        return in[e];
    }

    template <typename A, width_t I, width_t M, width_t T, size_t N>
    [[nodiscard]] constexpr auto
    extract(mv<A, I, M, T> const& in, std::array<elem_t, N> const& elements) noexcept
    {
        mv<A, I, M, N> out{};
        auto element = elements.begin();
        auto term    = in.cbegin();

        while (true)
        {
            if (element == elements.end() || term == in.cend())
            {
                return out;
            }

            if (term->element > *element)
            {
                ++element;
            }
            else if (term->element < *element)
            {
                ++term;
            }
            else
            {
                out.push(term++, one, *element);
            }
        }
    }
} // namespace detail

// Convenience template variable for making basis elements
template <typename A, uint32_t G, int N = 1, int D = 1>
constexpr inline mv<A, 0, 1, 1> e{mv_size{0, 1, 1}, {}, {mon{rat{N, D}, zero, 0, 0}}, {term{1, 0, G}}};

} // namespace gal
