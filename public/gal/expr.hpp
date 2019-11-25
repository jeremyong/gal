#pragma once

// expr.hpp
// Defines operators for constructing expression template trees as an intermediate structure for
// computation.
//
// Node checksums:
// Each expr node exposes a static constant corresponding to a checksum of itself along with its
// children. Checksums are computed recursively by summing child node checksums to produce the
// parent checksum. Afterwards, the operation represented by the node is written into the result.

#include "algebra.hpp"
#include "crc.hpp"
#include "null_algebra.hpp"
#include "tuple.hpp"

#include <type_traits>

#ifdef GAL_DEBUG
#    include <sstream>
#    include <string>
#endif

namespace gal
{
namespace detail
{
    enum op : uint32_t
    {
        // Utility ops
        op_id,  // (00) Specifies a leaf of the expression tree referring to an input
        op_cse, // (01) This is used during DFA (data flow analysis) for referring to computed
                // intermediate values
        op_se, // (02) Indicate the start of a subexpression in the subexpression list. This is used
               // for polyadic operators such as sum and the exterior product
        op_noop, // (03) Noop used to indicate that current arg on the stack should be evaluated as
                 // a temporary

        // Algebraic ops
        op_rev,  // (04) Reversion
        op_pd,   // (05) Poincare dual
        op_exp,  // (06) Exponential map
        op_log,  // (07) Logarithmic map
        op_sum,  // (08) Sum
        op_gp,   // (09) Geometric product
        op_ep,   // (10) Exterior product
        op_lc,   // (11) Left contraction
        op_sip,  // (12) Symmetric inner product
        op_comp, // (13) Extract an individual component as a scalar
        op_div,  // (14) Division by scalar indeterminant

        // Scalar ops (these ops are applied term by term)
        op_sqrt, // (15) Square root
        op_sin,  // (16) Sine
        op_cos,  // (17) Cosine
        op_tan,  // (18) Tangent

        // Constants
        c_zero = 1 << 16, // multivector that is exactly zero
        c_const_start,
        c_pi = c_const_start, // The constant pi
        c_e,                  // The Euler constant
        c_const_end,
        c_scalar, // scalar component of a multivector
                  // c_scalar + E corresponds to the unit constant element E
                  // For example, c_scalar + (1 << (dim - 1)) represents the pseudoscalar
    };

    // RPN node
    struct node
    {
        uint32_t o = 0;

        // For an expression reference, the checksum is used as a placeholder to skip to a
        // subsequent node in the expression.
        crc_t checksum = 0;

        // Ex is an overloaded member variable (unions are not modifiable in a constexpr context)
        //
        // For op summands, the ex corresponds to the number of nodes contained in the summand (not
        // including this one).
        //
        // For exterior product factors, ex holds the number of operands (useful for computing the
        // merge parity).
        //
        // For the sum and exterior operators, ex refers to the number of arguments.
        //
        // For the id op (used for an input with an identifier), ex is the index of the input in the
        // input list.
        //
        // For extractions, ex refers to the element we wish to extract.
        //
        // For subexpression references (during DFA), ex is an index to the subexpression.
        width_t ex = 0;

        // This rational scaling constant can be used for sum subexpressions. If this is zero for a
        // subexpression, the subexpression is part of an exterior product.
        //
        // For sum and exterior product nodes, the denominator here contains the total length of the
        // expression that constitutes the subexpression up to and including the op_sum/op_ep.
        rat q{0, 0};

        constexpr bool operator==(node const& other) const
        {
            return o == other.o && checksum == other.checksum && ex == other.ex;
        }

        constexpr bool operator!=(node const& other) const
        {
            return !(*this == other);
        }
    };

    // Reverse-Polish notation expression
    // A := Algebra which defines dimensionality, metric tensor, product evaluation functions
    // S := Size
    template <typename A, width_t S>
    struct rpne
    {
        using algebra_t                   = A;
        constexpr static width_t capacity = S;

        // To simplify processing, we arrange this structure so that each node maps to exactly one
        // factor
        node nodes[S];
        // The count here can easily diverge from the capacity due to CSE and term elimination
        width_t count = 0;
        rat q{1, 1};

        constexpr auto operator[](elem_t e) noexcept
        {
            rpne<A, S + 1> out;
            out.append(*this);
            // Using operator[] successfully is not expected behavior, so we don't bother computing
            // the checksum in a fancy manner.
            out.nodes[S] = node{op_comp, back().checksum + e, e};
            out.count    = S + 1;
            out.q        = q;
            return out;
        }

        constexpr node* begin() noexcept
        {
            return nodes;
        }

        constexpr node const* begin() const noexcept
        {
            return nodes;
        }

        constexpr node* end() noexcept
        {
            return nodes + count;
        }

        constexpr node const* end() const noexcept
        {
            return nodes + count;
        }

        constexpr node pop() noexcept
        {
            return nodes[count--];
        }

        constexpr void pop(width_t n) noexcept
        {
            count -= n;
        }

        template <width_t S1>
        constexpr void append(rpne<A, S1> const& src) noexcept
        {
            for (width_t i = 0; i != src.count; ++i)
            {
                nodes[count++] = src.nodes[i];
            }
        }

        constexpr void append(node const& in) noexcept
        {
            nodes[count++] = in;
        }

        template <width_t S1>
        constexpr void append_omit_ends(rpne<A, S1> const& src) noexcept
        {
            for (width_t i = 0; i != src.count - 1; ++i)
            {
                nodes[count++] = src.nodes[i];
            }
        }

        constexpr void append_cse(width_t i) noexcept
        {
            // The index to the ref passed here is 1-indexed (disambiguates it from a not-found
            // result).
            nodes[count++] = node{op_cse, 0, i - 1, zero};
        }

        constexpr void append_noop() noexcept
        {
            nodes[count++] = node{op_noop, 0, 0, zero};
        }

        template <width_t S1>
        constexpr void append_as_se(rpne<A, S1> const& src) noexcept
        {
            nodes[count++] = node{op_se, 0, src.count, src.q};
            for (width_t i = 0; i != src.count; ++i)
            {
                nodes[count++] = src.nodes[i];
            }
        }

        constexpr void append(op o) noexcept
        {
            nodes[count++] = {o, crc32((o << 8) + back().checksum)};
        }

        constexpr void append(op o, crc_t c) noexcept
        {
            nodes[count++] = {o, c};
        }

        constexpr node& back() noexcept
        {
            return nodes[count - 1];
        }

        constexpr node const& back() const noexcept
        {
            return nodes[count - 1];
        }

        template <width_t S1>
        constexpr bool operator==(rpne<A, S1> const& other) const noexcept
        {
            if (count != other.count)
            {
                return false;
            }

            for (width_t i = 0; i != count; ++i)
            {
                if (nodes[i] != other.nodes[i])
                {
                    return false;
                }
            }

            return true;
        }
    };

    template <typename A, size_t S, typename D, typename... Ds>
    constexpr auto rpne_entities(std::array<rpne<A, 1>, S>& out, uint32_t current_id, size_t i) noexcept
    {
        out[i].nodes[0] = node{op_id, current_id, i};
        out[i].count    = 1;

        if constexpr (sizeof...(Ds) > 0)
        {
            if (std::is_floating_point_v<D>)
            {
                return rpne_entities<A, S, Ds...>(out, current_id + 1, i + 1);
            }
            else
            {
                return rpne_entities<A, S, Ds...>(out, current_id + D::size(), i + 1);
            }
        }
        else
        {
            if constexpr (std::is_floating_point_v<D>)
            {
                return ::gal::make_pair(current_id + 1, i + 1);
            }
            else
            {
                return ::gal::make_pair(current_id + D::size(), i + 1);
            }
        }
    }

    template <typename A, typename... D>
    constexpr auto rpne_entities() noexcept
    {
        std::array<rpne<A, 1>, sizeof...(D)> out;
        auto next = rpne_entities<A, sizeof...(D), D...>(out, 0, 0);
        return ::gal::make_pair(out, next);
    }

    template <typename A>
    constexpr auto rpne_from_constant(num_t n, den_t d) noexcept
    {
        rpne<A, 1> out;
        out.nodes[0] = node{c_scalar, c_scalar};
        out.q        = rat{n, d};
        out.count    = 1;
        return out;
    }

    template <typename T, size_t... I>
    constexpr auto rpne_total_capacity(T const&, std::index_sequence<I...>) noexcept
    {
        return (std::decay_t<decltype(static_cast<T*>(nullptr)->template get<I>())>::capacity + ...);
    }

    // Concatenate a variadic number of expressions in the input tuple as subexpressions
    template <auto const& input>
    constexpr auto rpne_concat() noexcept
    {
        using type = std::decay_t<decltype(input)>;

        if constexpr (is_tuple_v<type>)
        {
            using A                    = typename decltype(input.template get<0>())::algebra_t;
            constexpr width_t capacity = input.apply([](auto const&... args) {
                return (std::decay_t<decltype(args)>::capacity + ...);
            }) + static_cast<width_t>(std::decay_t<decltype(input)>::size());

            return input.apply([capacity](auto const&... args) {
                rpne<A, capacity> out;
                (out.append_as_se(args), ...);
                return out;
            });
        }
        else
        {
            return input;
        }
    }

    struct rpne_constant
    {
        width_t id;
        rat q = one;

        template <typename A>
        constexpr rpne<A, 1> convert() const noexcept
        {
            return {{node{c_pi, c_pi}}, 1, q};
        }
    };
} // namespace detail

constexpr inline detail::rpne_constant PI{detail::pi_ind};
constexpr inline detail::rpne_constant E{detail::e_ind};

template <typename A, width_t S1, width_t S2>
constexpr detail::rpne<A, S1 + S2 + 3> operator+(detail::rpne<A, S1> lhs, detail::rpne<A, S2> rhs)
{
    using detail::op_sum;
    using detail::rpne;

    // Conservative size accounts for 2 summand ops, 1 sum op, and the inputs.
    rpne<A, S1 + S2 + 3> out;

    // Early out if either the lhs or the rhs is exactly zero
    if (lhs.count == 0 || rhs.count == 0)
    {
        if (lhs.count == 0)
        {
            out.append(rhs);
        }
        else
        {
            out.append(lhs);
        }
        return out;
    }

    // If the trailing operation of the lhs is a summation, we can move the summation op to the
    // right of the rhs. By the same token, if the trailing operation on the rhs is a summation,
    // that can be coalesced as well.

    auto const& lhs_final = lhs.back();
    auto const& rhs_final = rhs.back();

    bool lhs_sum = lhs_final.o == op_sum;
    bool rhs_sum = rhs_final.o == op_sum;

    width_t lhs_index = 0;
    width_t rhs_index = 0;

    uint32_t summand_count = 0;

    if (lhs_sum && rhs_sum)
    {
        // Merge the lhs and rhs summands
        while (true)
        {
            // Note that we subtract 1 from the upper limit to avoid the final sum op
            if (lhs_index == lhs.count - 1)
            {
                // Copy remaining rhs nodes. (exclude final sum op)
                while (rhs_index != rhs.count - 1)
                {
                    auto rhs_summand = rhs.nodes[rhs_index++];
                    ++summand_count;
                    rhs_summand.q *= rhs.q;
                    out.nodes[out.count++] = rhs_summand;
                    for (width_t i = 0; i != rhs_summand.ex; ++i)
                    {
                        out.nodes[out.count++] = rhs.nodes[rhs_index++];
                    }
                }
                break;
            }
            else if (rhs_index == rhs.count - 1)
            {
                // Copy remaining lhs nodes.
                while (lhs_index != lhs.count - 1)
                {
                    auto lhs_summand = lhs.nodes[lhs_index++];
                    ++summand_count;
                    lhs_summand.q *= rhs.q;
                    out.nodes[out.count++] = lhs_summand;
                    for (width_t i = 0; i != lhs_summand.ex; ++i)
                    {
                        out.nodes[out.count++] = lhs.nodes[lhs_index++];
                    }
                }
                break;
            }

            auto const& lhs_summand = lhs.nodes[lhs_index];
            auto const& rhs_summand = rhs.nodes[rhs_index];

            detail::crc_t lhs_checksum = lhs_summand.checksum;
            detail::crc_t rhs_checksum = rhs_summand.checksum;

            if (lhs_checksum == rhs_checksum && lhs_summand.ex == rhs_summand.ex)
            {
                bool equal = true;
                // Check if the two summands are in fact equivalent
                for (width_t i = 0; i != lhs_summand.ex; ++i)
                {
                    if (lhs.nodes[lhs_index + 1 + i] != rhs.nodes[rhs_index + 1 + i])
                    {
                        equal = false;
                        break;
                    }
                }

                if (equal)
                {
                    rat q = lhs_summand.q * lhs.q + rhs_summand.q * rhs.q;
                    if (!q.is_zero())
                    {
                        auto& s    = out.nodes[out.count++];
                        s.q        = q;
                        s.o        = detail::op_se;
                        s.checksum = lhs_checksum;
                        s.ex       = lhs_summand.ex;
                        ++summand_count;

                        for (width_t i = 0; i != lhs_summand.ex; ++i)
                        {
                            out.nodes[out.count++] = lhs.nodes[lhs_index + 1 + i];
                        }
                    }

                    lhs_index += 1 + lhs_summand.ex;
                    rhs_index += 1 + rhs_summand.ex;

                    continue;
                }
            }

            // We let the if statement above *fall through* to the cases below in case of checksum
            // equality but term inequality

            if (lhs_checksum < rhs_checksum)
            {
                out.nodes[out.count] = lhs.nodes[lhs_index++];
                out.nodes[out.count++].q *= lhs.q;
                for (width_t i = 0; i != lhs_summand.ex; ++i)
                {
                    out.nodes[out.count++] = lhs.nodes[lhs_index++];
                }
            }
            else
            {
                out.nodes[out.count] = rhs.nodes[rhs_index++];
                out.nodes[out.count++].q *= rhs.q;
                for (width_t i = 0; i != rhs_summand.ex; ++i)
                {
                    out.nodes[out.count++] = rhs.nodes[rhs_index++];
                }
            }
            ++summand_count;
        }
        // Merge of sums complete
    }
    else if (lhs_sum || rhs_sum)
    {
        // Indicator for whether the "non-sum" term has been merged
        bool written              = false;
        rat sum_q                 = lhs_sum ? lhs.q : rhs.q;
        rat non_sum_q             = lhs_sum ? rhs.q : lhs.q;
        auto const* sum_nodes     = lhs_sum ? lhs.nodes : rhs.nodes;
        auto sum_count            = lhs_sum ? lhs.count : rhs.count;
        auto const* non_sum_nodes = lhs_sum ? rhs.nodes : lhs.nodes;
        auto non_sum_count        = lhs_sum ? rhs.count : lhs.count;
        auto& sum_index           = lhs_sum ? lhs_index : rhs_index;
        auto& non_sum_index       = lhs_sum ? rhs_index : lhs_index;
        auto const& non_sum_final = lhs_sum ? rhs_final : lhs_final;

        while (true)
        {
            if (sum_index == sum_count - 1)
            {
                if (!written)
                {
                    out.nodes[out.count++] = detail::node{
                        detail::op_se, non_sum_final.checksum, non_sum_count, non_sum_q};
                    ++summand_count;
                    for (width_t i = 0; i != non_sum_count; ++i)
                    {
                        out.nodes[out.count + i] = non_sum_nodes[i];
                    }
                    out.count += non_sum_count;
                }
                break;
            }

            auto& summand = sum_nodes[sum_index];

            if (!written && summand.checksum == non_sum_final.checksum && summand.ex == non_sum_count)
            {
                // Compare summand equality to non-sum
                bool equal = true;
                for (width_t i = 0; i != non_sum_count; ++i)
                {
                    if (sum_nodes[sum_index + 1 + i] != non_sum_nodes[i])
                    {
                        equal = false;
                        break;
                    }
                }

                if (equal)
                {
                    rat q   = summand.q * sum_q + non_sum_q;
                    written = true;
                    if (q.is_zero())
                    {
                        // The terms cancel
                        sum_index += 1 + summand.ex;
                    }
                    else
                    {
                        out.nodes[out.count++]
                            = detail::node{detail::op_se, non_sum_final.checksum, non_sum_count, q};
                        ++summand_count;
                        for (width_t i = 0; i != non_sum_count; ++i)
                        {
                            out.nodes[out.count + i] = non_sum_nodes[i];
                        }
                        out.count += non_sum_count;
                        sum_index += 1 + summand.ex;
                    }
                    continue;
                }
            }

            if (!written && summand.checksum > non_sum_final.checksum)
            {
                // Insert the non-sum into the summand list
                written = true;
                out.nodes[out.count++]
                    = detail::node{detail::op_se, non_sum_final.checksum, non_sum_count, non_sum_q};
                for (width_t i = 0; i != non_sum_count; ++i)
                {
                    out.nodes[out.count + i] = non_sum_nodes[i];
                }
                out.count += non_sum_count;
            }
            else
            {
                out.nodes[out.count] = summand;
                out.nodes[out.count++].q *= sum_q;
                for (width_t i = 0; i != summand.ex; ++i)
                {
                    out.nodes[out.count + i] = sum_nodes[sum_index + 1 + i];
                }
                sum_index += summand.ex + 1;
                out.count += summand.ex;
            }
            ++summand_count;
        }
    }
    else
    {
        // Neither the lhs nor the rhs are sums. Add both as new summands.

        auto& final_lhs = lhs.back();
        auto& final_rhs = rhs.back();

        if (final_lhs.checksum == final_rhs.checksum && lhs.count == rhs.count)
        {
            if (lhs == rhs)
            {
                // No summation is necessary.
                rat q = lhs.q + rhs.q;
                if (q.is_zero())
                {
                    out.count = 0;
                    return out;
                }
                else
                {
                    for (width_t i = 0; i != lhs.count; ++i)
                    {
                        out.nodes[i] = lhs.nodes[i];
                    }
                    out.q     = q;
                    out.count = lhs.count;
                    return out;
                }
            }
        }

        if (final_lhs.checksum < final_rhs.checksum)
        {
            out.nodes[out.count++]
                = detail::node{detail::op_se, final_lhs.checksum, lhs.count, lhs.q};
            out.append(lhs);
            out.nodes[out.count++]
                = detail::node{detail::op_se, final_rhs.checksum, rhs.count, rhs.q};
            out.append(rhs);
        }
        else
        {
            out.nodes[out.count++]
                = detail::node{detail::op_se, final_rhs.checksum, rhs.count, rhs.q};
            out.append(rhs);
            out.nodes[out.count++]
                = detail::node{detail::op_se, final_lhs.checksum, lhs.count, lhs.q};
            out.append(lhs);
        }
        summand_count = 2;
    }

    // Do a final pass to see if the sum has been collapsed to a single term
    if (summand_count == 1)
    {
        // Copy-erase the initial subexpression op
        out.q = out.nodes[0].q;
        for (width_t i = 0; i != out.count - 1; ++i)
        {
            out.nodes[i] = out.nodes[i + 1];
        }
        --out.count;
        return out;
    }

    // Combine all summand CRCs into a pseudo bloom filter
    width_t i       = 0;
    detail::crc_t c = 0;
    while (i < out.count)
    {
        auto const& summand = out.nodes[i];
        // The lowest 5 bits of the summand map directly to a bit in the crc
        c |= 1 << (summand.checksum & 31);
        i += 1 + summand.ex;
    }

    auto& sum_node = out.nodes[out.count++];
    sum_node.o     = op_sum;
    // NOTE: we don't bother applying a crc function here as the checksum acts as a simple bloom
    // filter
    sum_node.checksum = c;
    sum_node.ex       = summand_count;
    sum_node.q.den    = out.count;
    return out;
}

template <typename A, width_t S>
constexpr auto operator*(int n, detail::rpne<A, S> rhs)
{
    if (n == 0 || rhs.count == 0)
    {
        rhs.count = 0;
        return rhs;
    }
    else
    {
        rhs.q *= rat{n, 1};
        return rhs;
    }
}

constexpr auto operator*(int n, detail::rpne_constant c) noexcept
{
    c.q *= rat{n, 1};
    return c;
}

constexpr auto operator*(detail::rpne_constant c, int n) noexcept
{
    c.q *= rat{n, 1};
    return c;
}

constexpr auto operator/(detail::rpne_constant c, int d) noexcept
{
    c.q *= rat{1, d};
    return c;
}

template <typename A, width_t S>
constexpr auto operator-(detail::rpne<A, S> in)
{
    if (in.count > 0)
    {
        in.q *= minus_one;
    }
    return in;
}

template <typename A, width_t S1, width_t S2>
constexpr auto operator-(detail::rpne<A, S1> const& lhs, detail::rpne<A, S2> const& rhs)
{
    return lhs + -rhs;
}

template <typename A, width_t S>
constexpr auto operator/(detail::rpne<A, S> lhs, int d)
{
    if (lhs.count > 0)
    {
        lhs.q *= rat{1, d};
    }
    return lhs;
}

template <typename A, width_t S1, width_t S2>
constexpr auto operator/(detail::rpne<A, S1> const& lhs, detail::rpne<A, S2> const& rhs)
{
    detail::rpne<A, S1 + S2 + 1> out;
    out.append(lhs);
    out.append(rhs);
    out.append(detail::op_div);
    out.back().q = lhs.q / rhs.q;
    out.q        = one;
    return out;
}

template <typename A, width_t S>
constexpr auto operator/(detail::rpne_constant const& lhs, detail::rpne<A, S> const& rhs)
{
    return lhs.template convert<A>() / rhs;
}

template <typename A, width_t S>
constexpr auto operator/(detail::rpne<A, S> const& lhs, detail::rpne_constant const& rhs)
{
    return lhs / rhs.template convert<A>();
}

template <typename A, width_t S>
constexpr auto operator+(int n, detail::rpne<A, S> const& rhs)
{
    return detail::rpne_from_constant<A>(n, 1) + rhs;
}

template <typename A, width_t S>
constexpr auto operator-(int n, detail::rpne<A, S> const& rhs)
{
    return detail::rpne_from_constant<A>(n, 1) + -rhs;
}

template <typename A, width_t S>
constexpr auto operator+(detail::rpne<A, S> const& lhs, int n)
{
    return detail::rpne_from_constant<A>(n, 1) + lhs;
}

template <typename A, width_t S>
constexpr auto operator-(detail::rpne<A, S> const& lhs, int n)
{
    return detail::rpne_from_constant<A>(-n, 1) + lhs;
}

template <typename A, width_t S>
constexpr auto operator~(detail::rpne<A, S> const& in)
{
    detail::rpne<A, S + 1> out;
    if (in.count == 0)
    {
        return out;
    }

    out.append(in);

    if (in.back().o == detail::op_rev)
    {
        --out.count;
    }
    else
    {
        out.append(detail::op_rev);
    }
    return out;
}

template <typename A, width_t S>
constexpr auto operator!(detail::rpne<A, S> const& in)
{
    detail::rpne<A, S + 1> out;
    if (in.count == 0)
    {
        return out;
    }

    out.append(in);
    out.append(detail::op_pd);
    return out;
}

template <typename A, width_t S>
constexpr auto sqrt(detail::rpne<A, S> const& in)
{
    detail::rpne<A, S + 1> out;
    out.append(in);
    out.append(detail::op_sqrt);
    out.back().q = in.q;
    out.q        = one;
    return out;
}

template <typename A, width_t S>
constexpr auto sin(detail::rpne<A, S> const& in)
{
    detail::rpne<A, S + 1> out;
    out.append(in);
    out.append(detail::op_sin);
    out.back().q = in.q;
    out.q        = one;
    return out;
}

template <typename A, width_t S>
constexpr auto cos(detail::rpne<A, S> const& in)
{
    detail::rpne<A, S + 1> out;
    out.append(in);
    out.append(detail::op_cos);
    out.back().q = in.q;
    out.q        = one;
    return out;
}

template <typename A, width_t S>
constexpr auto tan(detail::rpne<A, S> const& in)
{
    detail::rpne<A, S + 1> out;
    out.append(in);
    out.append(detail::op_tan);
    out.back().q = in.q;
    out.q        = one;
    return out;
}

// Closed-form exponential function
// NOTE: results are *undefined* when the arg is not a bivector
template <typename A, width_t S>
constexpr auto exp(detail::rpne<A, S> const& in)
{
    detail::rpne<A, S + 1> out;
    out.append(in);
    out.append(detail::op_exp);
    return out;
}

// Closed-form logarithmic function
// NOTE: results are *undefined* when the arg is not a member of the even subalgebra
template <typename A, width_t S>
constexpr auto log(detail::rpne<A, S> const& in)
{
    detail::rpne<A, S + 1> out;
    out.append(in);
    out.append(detail::op_log);
    return out;
}

template <typename A, width_t S1, width_t S2>
constexpr auto operator*(detail::rpne<A, S1> const& lhs, detail::rpne<A, S2> const& rhs)
{
    detail::rpne<A, S1 + S2 + 1> out;
    if (lhs.count == 0)
    {
        out.append(rhs);
        return out;
    }
    else if (rhs.count == 0)
    {
        out.append(lhs);
        return out;
    }

    auto q = lhs.q * rhs.q;
    if (q.is_zero())
    {
        return out;
    }

    out.q = q;
    out.append(lhs);
    out.append(rhs);
    out.append(detail::op_gp);
    out.nodes[out.count - 1].checksum
        = detail::crc32(lhs.back().checksum + ~rhs.back().checksum + (detail::op_gp << 8));
    return out;
}

template <typename A, width_t S>
constexpr auto operator*(detail::rpne<A, S> const& lhs, detail::rpne_constant const& rhs)
{
    return lhs * rhs.template convert<A>();
}

template <typename A, width_t S>
constexpr auto operator*(detail::rpne_constant const& lhs, detail::rpne<A, S> const& rhs)
{
    return rhs * lhs.template convert<A>();
}

template <typename A, width_t S1, width_t S2>
constexpr auto operator^(detail::rpne<A, S1> const& lhs, detail::rpne<A, S2> const& rhs)
{
    detail::rpne<A, S1 + S2 + 3> out;

    if (lhs.count == 0 || rhs.count == 0 || rhs.q.is_zero() || lhs.q.is_zero())
    {
        return out;
    }

    bool lhs_ep = lhs.back().o == detail::op_ep;
    bool rhs_ep = rhs.back().o == detail::op_ep;

    width_t lhs_index = 0;
    width_t rhs_index = 0;

    uint32_t factor_count = 0;

    if (lhs_ep && rhs_ep)
    {
        factor_count = lhs.back().ex + rhs.back().ex;
        out.append_omit_ends(lhs);
        out.append_omit_ends(rhs);
    }
    else if (lhs_ep || rhs_ep)
    {
        if (lhs_ep)
        {
            factor_count = lhs.back().ex + 1;
            // Omit final op
            out.append_omit_ends(lhs);
            out.nodes[out.count++]
                = detail::node{detail::op_se, rhs.back().checksum, rhs.count, zero};
            out.append(rhs);
        }
        else
        {
            factor_count = rhs.back().ex + 1;
            out.nodes[out.count++]
                = detail::node{detail::op_se, lhs.back().checksum, lhs.count, zero};
            out.append(lhs);
            out.append_omit_ends(rhs);
        }
        // Merged non-ep factor into ep
    }
    else
    {
        // Neither the lhs nor the rhs is an exterior product.
        factor_count = 2;

        // Check first for equivalency
        auto lhs_checksum = lhs.back().checksum;
        auto rhs_checksum = rhs.back().checksum;

        out.nodes[out.count++] = detail::node{detail::op_se, lhs_checksum, lhs.count, zero};
        out.append(lhs);
        out.nodes[out.count++] = detail::node{detail::op_se, rhs_checksum, rhs.count, zero};
        out.append(rhs);
        out.q = lhs.q * rhs.q;
    }

    // Loop through all ep factors to compute the bloom filter
    width_t i       = 1;
    detail::crc_t c = 0;
    while (i < out.count)
    {
        auto const& ep_f = out.nodes[i];
        c |= 1 << (ep_f.checksum & 31);
        i += 1 + ep_f.ex;
    }

    auto& ep_node    = out.nodes[out.count++];
    ep_node.o        = detail::op_ep;
    ep_node.checksum = c;
    ep_node.ex       = factor_count;
    ep_node.q.den    = out.count;
    return out;
}

template <typename A, width_t S1, width_t S2>
constexpr auto operator&(detail::rpne<A, S1> const& lhs, detail::rpne<A, S2> const& rhs)
{
    return !(!lhs ^ !rhs);
}

template <typename A, width_t S1, width_t S2>
constexpr auto operator>>(detail::rpne<A, S1> const& lhs, detail::rpne<A, S2> const& rhs)
{
    detail::rpne<A, S1 + S2 + 1> out;
    if (lhs.count == 0)
    {
        out.append(rhs);
        return out;
    }
    else if (rhs.count == 0)
    {
        out.append(lhs);
        return out;
    }

    auto q = lhs.q * rhs.q;
    if (q.is_zero())
    {
        return out;
    }

    out.q = q;
    out.append(lhs);
    out.append(rhs);
    out.append(detail::op_lc);
    out.nodes[out.count - 1].checksum
        = detail::crc32(lhs.back().checksum + ~rhs.back().checksum + (detail::op_lc << 8));
    return out;
}

template <typename A, width_t S1, width_t S2>
constexpr auto operator|(detail::rpne<A, S1> const& lhs, detail::rpne<A, S2> const& rhs)
{
    detail::rpne<A, S1 + S2 + 1> out;
    if (lhs.count == 0)
    {
        out.append(rhs);
        return out;
    }
    else if (rhs.count == 0)
    {
        out.append(lhs);
        return out;
    }

    auto q = lhs.q * rhs.q;
    if (q.is_zero())
    {
        return out;
    }

    out.q = q;
    if (lhs.back().checksum < rhs.back().checksum)
    {
        out.append(lhs);
        out.append(rhs);
    }
    else
    {
        out.append(rhs);
        out.append(lhs);
    }
    out.append(detail::op_sip);
    out.nodes[out.count - 1].checksum
        = detail::crc32(lhs.back().checksum + rhs.back().checksum + (detail::op_sip << 8));
    return out;
}

template <typename A, width_t S1, width_t S2>
constexpr auto operator%(detail::rpne<A, S1> const& lhs, detail::rpne<A, S2> const& rhs)
{
    return rhs * lhs * (~rhs);
}

template <typename A, width_t S1, width_t S2>
constexpr auto scalar_product(detail::rpne<A, S1> const& lhs, detail::rpne<A, S2> const& rhs)
{
    detail::rpne<A, S1 + S2 + 2> out;
    if (lhs.count == 0)
    {
        out.append(rhs);
        return out;
    }
    else if (rhs.count == 0)
    {
        out.append(lhs);
        return out;
    }

    auto q = lhs.q * rhs.q;
    if (q.is_zero())
    {
        return out;
    }

    out.q = q;
    if (lhs.back().checksum < rhs.back().checksum)
    {
        out.append(lhs);
        out.append(rhs);
    }
    else
    {
        out.append(rhs);
        out.append(lhs);
    }

    // Implement the scalar product as the scalar component of the symmetric inner product

    out.append(detail::op_sip);
    out.nodes[out.count - 1].checksum
        = detail::crc32(lhs.back().checksum + rhs.back().checksum + (detail::op_sip << 8));

    out.append(detail::op_comp);
    out.nodes[out.count - 1].ex       = 0;
    out.nodes[out.count - 1].checksum = out.nodes[out.count - 2].checksum;
    return out;
}

#ifdef GAL_DEBUG
template <typename A, width_t S>
std::string to_string(detail::rpne<A, S> const& in, bool show_checksums = false, bool index = false)
{
    using namespace gal::detail;

    std::stringstream str;

    str << in.count << " node(s): " << in.q.num << '/' << in.q.den << '*';
    if (index)
    {
        str << '\n';
    }
    for (width_t i = 0; i != in.count; ++i)
    {
        if (index)
        {
            str << '[' << i << "] - ";
        }

        auto const& n = in.nodes[i];
        switch (n.o)
        {
        case op_id:
            str << "id(" << n.checksum << ") ";
            break;
        case op_cse:
            str << "cse(" << n.ex << ") ";
            break;
        case op_se:
            str << "se(" << n.q.num << '/' << n.q.den << " count: " << n.ex << ") ";
            break;
        case op_noop:
            str << "noop ";
            break;
        case op_rev:
            if (show_checksums)
            {
                str << "~(" << n.checksum << ") ";
            }
            else
            {
                str << "~ ";
            }
            break;
        case op_pd:
            if (show_checksums)
            {
                str << "!(" << n.checksum << ") ";
            }
            else
            {
                str << "! ";
            }
            break;
        case op_sum:
            if (show_checksums)
            {
                str << "+(" << n.ex << ", " << n.checksum << ") ";
            }
            else
            {
                str << "+(" << n.ex << ") ";
            }
            break;
        case op_gp:
            if (show_checksums)
            {
                str << "*(" << n.checksum << ") ";
            }
            else
            {
                str << "* ";
            }
            break;
        case op_ep:
            if (show_checksums)
            {
                str << "^(" << n.ex << ", " << n.checksum << ") ";
            }
            else
            {
                str << "^ ";
            }
            break;
        case op_lc:
            if (show_checksums)
            {
                str << "<<(" << n.checksum << ") ";
            }
            else
            {
                str << "<< ";
            }
            break;
        case op_sip:
            if (show_checksums)
            {
                str << "|(" << n.checksum << ") ";
            }
            else
            {
                str << "| ";
            }
            break;
        case op_exp:
            if (show_checksums)
            {
                str << "exp(" << n.checksum << ") ";
            }
            else
            {
                str << "exp ";
            }
            break;
        case op_log:
            if (show_checksums)
            {
                str << "log(" << n.checksum << ") ";
            }
            else
            {
                str << "log ";
            }
            break;
        case op_comp:
            str << "[" << n.ex << "] ";
            break;
        case op_div:
            if (show_checksums)
            {
                str << "/(" << n.checksum << ") ";
            }
            else
            {
                str << "/"
                    << "(" << n.q.num << '/' << n.q.den << ") ";
            }
            break;
        case op_sqrt:
            str << "sqrt"
                << "(" << n.q.num << '/' << n.q.den << ") ";
            break;
        case op_sin:
            str << "sin"
                << "(" << n.q.num << '/' << n.q.den << ") ";
            break;
        case op_cos:
            str << "cos"
                << "(" << n.q.num << '/' << n.q.den << ") ";
            break;
        case op_tan:
            str << "tan"
                << "(" << n.q.num << '/' << n.q.den << ") ";
            break;
        case c_zero:
            str << "0 ";
            break;
        case c_pi:
            str << "pi ";
            break;
        case c_e:
            str << "e ";
            break;
        default:
            uint32_t c = n.o - c_scalar;
            str << "e_" << c << ' ';
            break;
        }
        if (index)
        {
            str << '\n';
        }
    }
    str << std::endl;

    return str.str();
}
#endif
} // namespace gal
