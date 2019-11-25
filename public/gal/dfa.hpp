#pragma once

// dfa.hpp
// Implementation file for performing data flow analysis on expression trees.

#include "algebra.hpp"
#include "algorithm.hpp"
#include "expr.hpp"
#include "stack.hpp"
#include "tuple.hpp"

#include <array>

namespace gal
{
template <typename A, width_t I1, width_t M1, width_t T1, width_t I2, width_t M2, width_t T2>
constexpr auto operator+(mv<A, I1, M1, T1> const& lhs, mv<A, I2, M2, T2> const& rhs)
{
    return detail::sum(lhs, rhs);
}

template <typename A, width_t I1, width_t M1, width_t T1, width_t I2, width_t M2, width_t T2>
constexpr auto operator^(mv<A, I1, M1, T1> const& lhs, mv<A, I2, M2, T2> const& rhs) noexcept
{
    return detail::product(typename A::exterior{}, lhs, rhs);
}

namespace detail
{
    // DFA refers to "data flow analysis," one of the main ideas from compiler engineering loosely
    // applied here. The primary purpose of the routines in this file is to convert an expression
    // tree into a reduced form with common subexpressions factored out from the original.

    // Count the number of op_id nodes in an expression
    template <typename A, width_t S>
    constexpr width_t rpn_id_count(rpne<A, S> const& exp) noexcept
    {
        width_t count = 0;
        for (auto const& node : exp)
        {
            if (node.o == op_id)
            {
                ++count;
            }
        }
        return count;
    }

    // Given an expression and a known count of id nodes, return an array of the ids of such nodes
    template <typename A, width_t S, width_t C>
    constexpr auto rpn_ids(rpne<A, S> const& exp, std::integral_constant<width_t, C>) noexcept
    {
        std::array<uint32_t, C> ids{};
        std::array<uint32_t, C> indices{};
        size_t i = 0;
        for (auto const& node : exp)
        {
            if (node.o == op_id)
            {
                ids[i]       = node.checksum;
                indices[i++] = node.ex;
                if (i == C)
                {
                    break;
                }
            }
        }
        return ::gal::pair{ids, indices};
    }

    // Functor to retrieve a tuple of multivectors given an array of ids
    template <typename A, auto const& Ids, auto const& Indices, typename... D>
    struct rpn_inputs
    {
        template <size_t... I>
        constexpr auto operator()(std::index_sequence<I...>) noexcept
        {
            return make_tuple(ie<nth_element<Indices[I], D...>>(Ids[I])...);
        }

        template <typename I>
        constexpr auto ie(uint32_t id) noexcept
        {
            if constexpr (std::is_floating_point_v<I>)
            {
                return mv<A, 1, 1, 1>{
                    mv_size{1, 1, 1}, {ind{id, one}}, {mon{one, one, 1, 0}}, {term{1, 0, 0}}};
            }
            else
            {
                return I::ie(id);
            }
        }
    };

    template <typename T>
    struct mv_scaler
    {
        rat q = zero;

        constexpr auto operator()(T const& m) noexcept
        {
            return pair{m.first, scale(q, m.second)};
        }
    };

    template <typename T>
    struct rpn_temp
    {
        T ie;
        op o;
        crc_t checksum;
        uint32_t id;
    };

    template <typename T>
    rpn_temp(T&&, op, crc_t, uint32_t)->rpn_temp<std::decay_t<T>>;

    template <typename I, typename T, typename R>
    struct rpn_result
    {
        // Tuple of input multivectors in indeterminate form
        I inputs;
        // Tuple of rpn_temp values
        T temps;
        // Tuple of output multivectors in indeterminate form
        R results;
        // Number of temporary indeterminates needed
        uint32_t id_count;
    };

    struct cse
    {
        crc_t crc           = 0;
        width_t offset      = 0;
        width_t count       = 0;
        width_t refs        = 0;
        width_t final_index = 0;
        bool required       = false;
    };

    template <width_t C>
    struct cses
    {
        cse ses[C];
        width_t count = 0;
    };

    template <typename A, width_t C>
    constexpr bool
    compare_se(cses<C> const& known, rpne<A, C> const& exp, cse const& lhs, cse const& rhs) noexcept
    {
        if (lhs.crc != rhs.crc || lhs.count != rhs.count)
        {
            return false;
        }

        for (width_t i = 0; i != lhs.count;)
        {
            node const* lhs_n = &exp.nodes[lhs.offset + i];
            node const* rhs_n = &exp.nodes[rhs.offset + i];

            if (lhs_n->o == op_cse && rhs_n->o == op_cse)
            {
                if (lhs_n->ex == rhs_n->ex)
                {
                    i = lhs_n->checksum - lhs.offset;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                if (lhs_n->o == op_cse)
                {
                    while (lhs_n->o == op_cse)
                    {
                        width_t offset = known.ses[lhs_n->ex].offset;
                        lhs_n          = &exp.nodes[offset];
                    }
                }

                if (rhs_n->o == op_cse)
                {
                    while (rhs_n->o == op_cse)
                    {
                        width_t offset = known.ses[rhs_n->ex].offset;
                        rhs_n          = &exp.nodes[offset];
                    }
                }

                if (*lhs_n != *rhs_n)
                {
                    return false;
                }
                ++i;
            }
        }

        return true;
    }

    // Given a subexpression, increment its counter if it was registered previously. Otherwise,
    // append it.
    template <typename A, width_t C>
    constexpr width_t register_se(rpne<A, C> const& exp,
                                  cses<C>& known,
                                  crc_t crc,
                                  width_t offset,
                                  width_t count,
                                  bool required) noexcept
    {
        cse next{crc, offset, count, 1, 0, required};
        for (width_t i = 0; i != known.count; ++i)
        {
            auto& known_se = known.ses[i];
            if (compare_se(known, exp, known_se, next))
            {
                ++known_se.refs;
                known_se.required = required || known_se.required;
                return i + 1;
            }
        }

        known.ses[known.count++] = next;
        return 0;
    }

    template <typename A, width_t C>
    constexpr static auto rpn_reshape(rpne<A, C> expr) noexcept
    {
        // Unroll the rpn expression. As args are pushed onto the stack, check against an array
        // of known subexpressions. If there is a match, increment the occurence counter.
        // Otherwise, append it to the list. Finally, reshape the input expression to refer to
        // intermediate results.

        // Track known subexpressions here.
        cses<C> known{};

        // Re-express the expression with CSEs factored out
        // If no cses are found, the output will match the input. A cse is guaranteed to have
        // length 2 or more, so finding one will decrease the length of the output.
        // However, the output length can *increase* when no CSEs are found but tracendentals are
        // present. Each transcendental produces an additional op to lift it out of the primary
        // expression so in the worst case scenario, 2 * C additional ops are necessary.
        rpne<A, 3 * C> out{};

        stack<width_t, C> offsets{};

        width_t se_refs[C] = {};
        width_t se_count   = 0;

        for (width_t i = 0; i != expr.count; ++i)
        {
            auto& n = expr.nodes[i];

            switch (n.o)
            {
            case op_comp:
                // Component selection does not constitute a CSE candidate
                break;
            case op_rev:
                // fallthrough
            case op_pd: {
                width_t offset = offsets.peek();
                width_t se_i = register_se(expr, known, n.checksum, offset, i - offset + 1, false);
                if (se_i > 0)
                {
                    node const& prev = expr.nodes[offset];
                    if (prev.o == op_cse)
                    {
                        --known.ses[prev.ex].refs;
                    }

                    node& ref    = expr.nodes[offset];
                    ref.o        = op_cse;
                    ref.checksum = i + 1;
                    ref.ex       = se_i - 1;
                    if (se_count > 0)
                    {
                        // Track the subexpression index if we have an active addend or factor of a
                        // polyadic operator accumulating.
                        se_refs[se_count - 1] = se_i;
                    }
                }
                break;
            }
            case op_div: {
                // Division requires its second operand to be extracted.

                // Mark the last known se as required IF
                // 1. The prior se is not a constant AND
                // 2. The prior se is not a cse

                if (i > 2)
                {
                    width_t last_offset = offsets.peek();
                    node const& prev    = expr.nodes[last_offset];
                    if (last_offset == i - 2 && expr.nodes[i - 1].o != op_comp)
                    {
                        known.ses[known.count - 1].required = true;
                    }
                    else if (last_offset < i - 2 && prev.o != op_cse && prev.checksum != i)
                    {
                        known.ses[known.count - 1].required = true;
                    }
                }

                offsets.pop();
                width_t offset = offsets.peek();
                width_t se_i   = register_se(expr, known, n.checksum, offset, i - offset + 1, true);
                if (se_i > 0)
                {
                    // Unlike the unary and polyadic ops, we don't decrement the ref count of a
                    // preceding common subexpression
                    node& ref    = expr.nodes[offset];
                    ref.o        = op_cse;
                    ref.checksum = i + 1;
                    ref.ex       = se_i - 1;
                    if (se_count > 0)
                    {
                        se_refs[se_count - 1] = se_i;
                    }
                }
                break;
            }
            case op_sqrt:
            case op_sin:
            case op_cos:
            case op_tan:
            case op_exp:
                // fallthrough
            case op_log: {
                width_t offset = offsets.peek();

                // Check if the argument has been extracted. If not, mark it as required.
                if (i > 2)
                {
                    node const& prev = expr.nodes[offset];
                    if (offset == i - 2 && expr.nodes[i - 1].o != op_comp)
                    {
                        known.ses[known.count - 1].required = true;
                    }
                    else if (offset < i - 2 && !(prev.o == op_cse && prev.checksum == i))
                    {
                        known.ses[known.count - 1].required = true;
                    }
                }

                width_t se_i = register_se(expr, known, n.checksum, offset, i - offset + 1, true);
                if (se_i > 0)
                {
                    // Unlike the unary and polyadic ops, we don't decrement the ref count of a
                    // preceding common subexpression
                    node& ref    = expr.nodes[offset];
                    ref.o        = op_cse;
                    ref.checksum = i + 1;
                    ref.ex       = se_i - 1;
                    if (se_count > 0)
                    {
                        se_refs[se_count - 1] = se_i;
                    }
                }
                break;
            }
            case op_gp:
                // fallthrough
            case op_lc:
                // fallthrough
            case op_sip: {
                width_t last_offset = offsets.pop();
                width_t offset      = offsets.peek();
                width_t se_i = register_se(expr, known, n.checksum, offset, i - offset + 1, false);
                if (se_i > 0)
                {
                    // If this binary op produces a common subexpression, decrement the ref count of
                    // its argument subexpressions
                    node& second = expr.nodes[last_offset];
                    if (second.o == op_cse)
                    {
                        --known.ses[second.ex].refs;
                    }
                    node& first = expr.nodes[offset];
                    if (first.o == op_cse)
                    {
                        --known.ses[first.ex].refs;
                    }

                    node& ref    = expr.nodes[offset];
                    ref.o        = op_cse;
                    ref.checksum = i + 1;
                    ref.ex       = se_i - 1;
                    if (se_count > 0)
                    {
                        se_refs[se_count - 1] = se_i;
                    }
                }
                break;
            }
            case op_sum:
                // fallthrough
            case op_ep: {
                // TODO: register variadic ops with a different function that compares
                // subsets
                // This offset refers to the subexpression op just before the first operand.
                offsets.pop(2 * n.ex - 1);
                width_t offset = offsets.peek();
                width_t se_i = register_se(expr, known, n.checksum, offset, i - offset + 1, false);
                if (se_i > 0)
                {
                    // Decrement the ref counts of all constituent addends/factors
                    for (width_t j = se_count - n.ex; j != se_count; ++j)
                    {
                        if (se_refs[j] > 0)
                        {
                            --known.ses[se_refs[j] - 1].refs;
                        }
                    }

                    // NOTE: It is *critical* that these references are reset to zero, since zero
                    // has a special meaning (CSE substitution for this summand did not occur).
                    for (width_t j = se_count - n.ex; j != se_count; ++j)
                    {
                        se_refs[j] = 0;
                    }

                    se_count -= n.ex;

                    node& ref    = expr.nodes[offset];
                    ref.o        = op_cse;
                    ref.checksum = i + 1;
                    ref.ex       = se_i - 1;
                }
                else
                {
                    for (width_t j = se_count - n.ex; j != se_count; ++j)
                    {
                        se_refs[j] = 0;
                    }

                    se_count -= n.ex;
                }

                // The acrobatics here are needed to support nested expressions containing polyadic
                // ops.
                if (se_i > 0 && se_count > 0)
                {
                    se_refs[se_count - 1] = se_i;
                }
                break;
            }
            case op_se:
                ++se_count;
                offsets.push(i);
                break;
            case op_id:
                // fallthrough
            default:
                offsets.push(i);
                break;
            }
        }

        // Identify all candidates that are required or have a refcount greater than one and append
        // them to the output. Also, we use this counter to track the CSE candidates that are
        // ultimately not written to get the final indexing right.
        width_t cse_count = 0;

        // During subexpression substitution, the length of factors and addends can decrease.
        stack<pair<node*, width_t>, C> se_stack{};

        for (width_t i = 0; i != known.count; ++i)
        {
            cse& c = known.ses[i];
            if (c.refs > 1 || c.required)
            {
                c.final_index = cse_count++;
                se_stack.clear();

                for (width_t j = c.offset; j != c.offset + c.count;)
                {
                    node& n = expr.nodes[j];
                    out.append(n);

                    if (n.o == op_se)
                    {
                        se_stack.push(make_pair(&out.back(), out.count));
                    }
                    else if (n.o == op_sum || n.o == op_ep)
                    {
                        se_stack.push(make_pair(&out.back(), out.count));

                        // All subexpressions on the stack can now have their lengths adjusted
                        for (width_t k = se_stack.count - n.ex - 1; k != se_stack.count - 1; ++k)
                        {
                            auto [se, offset]   = se_stack[k];
                            width_t next_offset = se_stack[k + 1].second;
                            se->ex              = next_offset - offset - 1;
                        }

                        auto first_se    = se_stack.pop(n.ex + 1);
                        out.back().q.den = out.count - first_se.second - 1;
                    }

                    if (n.o == op_cse)
                    {
                        // If we found a common subexpression within a common subexpression, skip
                        // nodes that should not be written.
                        j = n.checksum;
                    }
                    else
                    {
                        ++j;
                    }
                }
                // Use a noop as a separator
                out.append_noop();

                node& n = expr.nodes[c.offset];

                n.o        = op_cse;
                n.checksum = c.offset + c.count;
                n.ex       = i;
            }
        }

        // Perform a pass to fix up collapsed subexpression indices
        for (width_t i = 0; i != out.count; ++i)
        {
            node& n = out.nodes[i];
            if (n.o == op_cse)
            {
                n.ex = known.ses[n.ex].final_index;
            }
        }

        // The final result concatenates the common subexpressions with the reshaped expression.
        se_stack.clear();

        for (width_t i = 0; i != expr.count;)
        {
            auto const& n = expr.nodes[i];
            out.append(n);

            // The branch below is needed to fix subexpression lengths (which shorten on CSE
            // contraction).
            if (n.o == op_se)
            {
                se_stack.push(make_pair(&out.back(), out.count));
            }
            else if (n.o == op_sum || n.o == op_ep)
            {
                se_stack.push(make_pair(&out.back(), out.count));

                // All subexpressions on the stack can now have their lengths adjusted
                for (width_t k = se_stack.count - n.ex - 1; k != se_stack.count - 1; ++k)
                {
                    auto [se, offset]   = se_stack[k];
                    width_t next_offset = se_stack[k + 1].second;
                    se->ex              = next_offset - offset - 1;
                }

                auto first_se    = se_stack.pop(n.ex + 1);
                out.back().q.den = out.count - first_se.second - 1;
            }

            if (n.o == op_cse)
            {
                // Fix up the CSE index (unused expressions are skipped).
                out.back().ex = known.ses[n.ex].final_index;
                i             = n.checksum;
            }
            else
            {
                ++i;
            }
        }

        out.q = expr.q;
        return out;
    }

    template <typename I, typename T, typename R>
    rpn_result(I, T, R, uint32_t)->rpn_result<I, T, R>;

    template <typename I, typename T, typename M>
    struct rpn_state
    {
        I inputs;
        T temps;
        M args;
        uint32_t id_count;
    };

    template <typename I, typename T, typename M>
    rpn_state(I, T, M, uint32_t)->rpn_state<I, T, M>;

    // We thread the expression and evaluation index through as type parameters here to permit
    // heterogeneous return types.
    template <auto const& exp, width_t i, width_t l, auto const& State>
    struct rpn_ctx
    {
        using algebra_t = typename std::decay_t<decltype(exp)>::algebra_t;

        constexpr static auto const& n = exp.nodes[i];

        // Unused except when evaluating subexpressions
        constexpr static rpn_state temp_state{State.inputs, State.temps, tuple<>{}, State.id_count};

        constexpr static auto compute_state() noexcept
        {
            if constexpr (n.o == op_id)
            {
                auto pop = State.inputs.pop();

                if constexpr (uses_null_basis<algebra_t>)
                {
                    return rpn_state{
                        pop.second,
                        State.temps,
                        State.args.push(make_pair(n.checksum, to_natural_basis(pop.first))),
                        State.id_count};
                }
                else
                {
                    return rpn_state{pop.second,
                                     State.temps,
                                     State.args.push(make_pair(n.checksum, pop.first)),
                                     State.id_count};
                }
            }
            else if constexpr (n.o == op_cse)
            {
                auto const& temp = State.temps.template get<n.ex>();
                return rpn_state{State.inputs,
                                 State.temps,
                                 State.args.push(make_pair(n.checksum, temp.ie.create_ref(temp.id))),
                                 State.id_count};
            }
            else if constexpr (n.o == op_se)
            {
                // Evaluate up to the subexpression limit and push that as an arg
                constexpr auto& sub_state = rpn_ctx<exp, i + 1, i + 1 + n.ex, temp_state>::state;
                auto next_arg             = sub_state.args.template get<0>().second;
                if constexpr (n.q.num != 0)
                {
                    // This is a summand, so apply the scaling factor
                    next_arg.scale(n.q);
                    return rpn_state{sub_state.inputs,
                                     sub_state.temps,
                                     State.args.push(make_pair(n.checksum, next_arg)),
                                     State.id_count};
                }
                else
                {
                    return rpn_state{sub_state.inputs,
                                     sub_state.temps,
                                     State.args.push(make_pair(n.checksum, next_arg)),
                                     State.id_count};
                }
            }
            else if constexpr (n.o == op_noop)
            {
                // Force the current arg on the stack to be evaluated as a temporary.
                constexpr auto pop = State.args.pop();
                return rpn_state{State.inputs,
                                 State.temps.append(rpn_temp{
                                     pop.first.second, op_noop, n.checksum, State.id_count}),
                                 pop.second,
                                 State.id_count + pop.first.second.size.term};
            }
            else if constexpr (n.o == op_rev)
            {
                auto args = State.args;
                // Reverse the first argument (top of the stack). The type of the
                // multivector does not change upon reversion.
                args.template set<0>(make_pair(n.checksum, reverse(args.template get<0>().second)));
                return rpn_state{State.inputs, State.temps, args, State.id_count};
            }
            else if constexpr (n.o == op_pd)
            {
                auto args = State.args;
                args.template set<0>(
                    make_pair(n.checksum, poincare_dual(args.template get<0>().second)));
                return rpn_state{State.inputs, State.temps, args, State.id_count};
            }
            else if constexpr (n.o == op_sum)
            {
                constexpr auto split = State.args.template split<n.ex>();
                constexpr auto sum
                    = split.first.apply([](auto... args) { return (args.second + ...); });
                return rpn_state{
                    State.inputs,
                    State.temps,
                    split.second.push(make_pair(
                        n.checksum, sum.template resize<sum.size.ind, sum.size.mon, sum.size.term>())),
                    State.id_count};
            }
            else if constexpr (n.o == op_gp)
            {
                constexpr auto split = State.args.template split<2>();
                constexpr auto gp    = product(typename algebra_t::geometric{},
                                            split.first.template get<1>().second,
                                            split.first.template get<0>().second);
                return rpn_state{
                    State.inputs,
                    State.temps,
                    split.second.push(make_pair(
                        n.checksum, gp.template resize<gp.size.ind, gp.size.mon, gp.size.term>())),
                    State.id_count};
            }
            else if constexpr (n.o == op_ep)
            {
                constexpr auto split = State.args.template split<n.ex>();
                constexpr auto ep
                    = split.first.apply_reverse([](auto... args) { return (args.second ^ ...); });
                return rpn_state{
                    State.inputs,
                    State.temps,
                    split.second.push(make_pair(
                        n.checksum, ep.template resize<ep.size.ind, ep.size.mon, ep.size.term>())),
                    State.id_count};
            }
            else if constexpr (n.o == op_lc)
            {
                constexpr auto split = State.args.template split<2>();
                constexpr auto lc    = product(typename algebra_t::contract{},
                                            split.first.template get<1>().second,
                                            split.first.template get<0>().second);
                return rpn_state{
                    State.inputs,
                    State.temps,
                    split.second.push(make_pair(
                        n.checksum, lc.template resize<lc.size.ind, lc.size.mon, lc.size.term>())),
                    State.id_count};
            }
            else if constexpr (n.o == op_sip)
            {
                constexpr auto split = State.args.template split<2>();
                constexpr auto sip   = product(typename algebra_t::symmetric_inner{},
                                             split.first.template get<1>().second,
                                             split.first.template get<0>().second);
                return rpn_state{
                    State.inputs,
                    State.temps,
                    split.second.push(make_pair(
                        n.checksum, sip.template resize<sip.size.ind, sip.size.mon, sip.size.term>())),
                    State.id_count};
            }
            else if constexpr (n.o == op_div)
            {
                constexpr auto split = State.args.template split<2>();
                constexpr auto div   = divide(
                    split.first.template get<1>().second, split.first.template get<0>().second, n.q);
                return rpn_state{State.inputs,
                                 State.temps,
                                 split.second.push(make_pair(n.checksum, div)),
                                 State.id_count};
            }
            else if constexpr (n.o == op_sqrt)
            {
                auto args                    = State.args;
                args.template get<0>().first = n.checksum;
                args.template get<0>().second.sqrt(n.q);
                return rpn_state{State.inputs, State.temps, args, State.id_count};
            }
            else if constexpr (n.o == op_sin)
            {
                auto args                    = State.args;
                args.template get<0>().first = n.checksum;
                args.template get<0>().second.sin(n.q);
                return rpn_state{State.inputs, State.temps, args, State.id_count};
            }
            else if constexpr (n.o == op_cos)
            {
                auto args                    = State.args;
                args.template get<0>().first = n.checksum;
                args.template get<0>().second.cos(n.q);
                return rpn_state{State.inputs, State.temps, args, State.id_count};
            }
            else if constexpr (n.o == op_tan)
            {
                auto args                    = State.args;
                args.template get<0>().first = n.checksum;
                args.template get<0>().second.tan(n.q);
                return rpn_state{State.inputs, State.temps, args, State.id_count};
            }
            else if constexpr (n.o == op_exp)
            {
                // Transcendental functions require operations performed on reified results
                auto pop = State.args.pop();

                // The result of the exp function is an element of the even subalgebra
                auto exp_result = algebra_t::even_mv(State.id_count);

                return rpn_state{State.inputs,
                                 State.temps.append(rpn_temp{
                                     pop.first.second, op_exp, n.checksum, State.id_count}),
                                 pop.second.push(make_pair(n.checksum, exp_result)),
                                 State.id_count + decltype(exp_result)::ind_capacity()};
            }
            else if constexpr (n.o == op_log)
            {
                // Transcendental functions require operations performed on reified results
                auto pop = State.args.pop();

                // The result of the log function is a bivector
                auto log_result = algebra_t::bivector_mv(State.id_count);

                return rpn_state{State.inputs,
                                 State.temps.append(rpn_temp{
                                     pop.first.second, op_log, n.checksum, State.id_count}),
                                 pop.second.push(make_pair(n.checksum, log_result)),
                                 State.id_count + decltype(log_result)::ind_capacity()};
            }
            else if constexpr (n.o == op_comp)
            {
                constexpr auto pop  = State.args.pop();
                constexpr auto comp = pop.first.second[n.ex];
                return rpn_state{
                    State.inputs,
                    State.temps,
                    pop.second.push(make_pair(
                        pop.first.first,
                        comp.template resize<comp.size.ind, comp.size.mon, comp.size.term>())),
                    State.id_count};
            }
            else if constexpr (n.o == c_zero)
            {
                mv<algebra_t, 0, 0, 0> z{mv_size{0, 0, 0}, {}, {}, {}};
                return rpn_state{State.inputs,
                                 State.temps,
                                 State.args.push(make_pair(n.checksum, z)),
                                 State.id_count};
            }
            else if constexpr (n.o < c_const_end)
            {
                mv<algebra_t, 1, 1, 1> c{mv_size{1, 1, 1},
                                         {ind{n.o - c_const_start + ind_constant_start, one}},
                                         {mon{one, one, 1, 0}},
                                         {term{1, 0, 0}}};
                return rpn_state{State.inputs,
                                 State.temps,
                                 State.args.push(make_pair(n.checksum, c)),
                                 State.id_count};
            }
            else
            {
                auto g = n.o - c_scalar;
                mv<algebra_t, 0, 1, 1> c{
                    mv_size{0, 1, 1}, {}, {mon{one, zero, 0, 0}}, {term{1, 0, g}}};
                return rpn_state{State.inputs,
                                 State.temps,
                                 State.args.push(make_pair(n.checksum, c)),
                                 State.id_count};
            }
        }

        constexpr static width_t next_i() noexcept
        {
            if constexpr (n.o == op_se)
            {
                return i + 1 + n.ex;
            }
            else
            {
                return i + 1;
            }
        }

        constexpr static auto inner_state = compute_state();

        constexpr static auto& state = rpn_ctx<exp, next_i(), l, inner_state>::state;
    };

    // Termination condition
    template <auto const& expr, width_t l, auto const& State>
    struct rpn_ctx<expr, l, l, State>
    {
        constexpr static auto finalize_state() noexcept
        {
            auto out = State;
            if constexpr (l == expr.count)
            {
                // When finalizing an expression, only at this point do we account for the
                // outermost scaling factor. Subexpression scaling factors are stored on the
                // SE node.
                // out.args.template mutate<mv_scaler>(expr.q);
            }
            return out;
        }

        // constexpr static auto state = finalize_state();
        constexpr static auto const& state = State;
    };

    // TODO: We need a way to subdivide an expression to prevent unbounded term growth
} // namespace detail
} // namespace gal
