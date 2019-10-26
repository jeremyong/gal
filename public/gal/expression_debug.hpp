#pragma once

#include "expression.hpp"

namespace gal
{
// This function is like `gal::reify` but it doesn't shrink multivectors using constexpr and is thus amenable for
// runtime debugging
template <typename exp_t>
[[nodiscard]] auto debug_reify() noexcept
{
    // Recursive function which compiles an expression into a reification table suitable for further processing needed
    // to evaluate all polynomial terms.

    // Base case
    if constexpr (exp_t::op == expr_op::identity)
    {
        if constexpr (detail::uses_null_basis<typename exp_t::algebra_t>)
        {
            return detail::to_natural_basis(exp_t::lhs);
        }
        else
        {
            return exp_t::lhs;
        }
    }
    else if constexpr (exp_t::op == expr_op::negate)
    {
        return detail::negate(debug_reify<typename exp_t::lhs_t>());
    }
    else if constexpr (exp_t::op == expr_op::reverse)
    {
        return detail::reverse(debug_reify<typename exp_t::lhs_t>());
    }
    else if constexpr (exp_t::op == expr_op::poincare_dual)
    {
        return detail::poincare_dual(debug_reify<typename exp_t::lhs_t>());
    }
    else if constexpr (exp_t::op == expr_op::clifford_conjugate)
    {
        // TODO: implement me
        // return ::gal::clifford_conjugate(debug_reify(T1{}));
    }
    else if constexpr (exp_t::op == expr_op::shift)
    {
        return detail::shift(exp_t::rhs_t::q(), debug_reify<typename exp_t::lhs_t>());
    }
    else if constexpr (exp_t::op == expr_op::scale)
    {
        return detail::scale(exp_t::rhs_t::q(), debug_reify<typename exp_t::lhs_t>());
    }
    else if constexpr (exp_t::op == expr_op::extract)
    {
        return detail::extract(reify<typename exp_t::lhs_t>(), exp_t::elements);
    }
    else if constexpr (exp_t::op == expr_op::select)
    {
        // TODO: implement me
    }
    else
    {
        auto lhs = debug_reify<typename exp_t::lhs_t>();
        auto rhs = debug_reify<typename exp_t::rhs_t>();

        if constexpr (exp_t::op == expr_op::sum)
        {
            return detail::sum(lhs, rhs);
        }
        else if constexpr (exp_t::op == expr_op::difference)
        {
            auto n_rhs = detail::negate(rhs);
            return detail::sum(lhs, n_rhs);
        }
        else if constexpr (exp_t::op == expr_op::geometric)
        {
            return detail::product(typename decltype(lhs)::algebra_t::geometric{}, lhs, rhs);
        }
        else if constexpr (exp_t::op == expr_op::sandwich)
        {
            // rhs * lhs * ~rhs
            auto rhs_reverse = detail::reverse(rhs);
            auto temp        = detail::product(typename decltype(lhs)::algebra_t::geometric{}, lhs, rhs_reverse);
            return detail::product(typename decltype(lhs)::algebra_t::geometric{}, rhs, temp);
        }
        else if constexpr (exp_t::op == expr_op::exterior)
        {
            return detail::product(typename decltype(lhs)::algebra_t::exterior{}, lhs, rhs);
        }
        else if constexpr (exp_t::op == expr_op::regressive)
        {
            // TODO: check if both lhs and rhs are dual
            auto lhs_dual = detail::poincare_dual(lhs);
            auto rhs_dual = detail::poincare_dual(rhs);
            auto out = detail::product(typename decltype(lhs)::algebra_t::exterior{}, lhs_dual, rhs_dual);
            return detail::poincare_dual(out);
        }
        else if constexpr (exp_t::op == expr_op::contract)
        {
            return detail::product(typename decltype(lhs)::algebra_t::contract{}, lhs, rhs);
        }
        else if constexpr (exp_t::op == expr_op::symmetric_inner)
        {
            return detail::product(typename decltype(lhs)::algebra_t::symmetric_inner{}, lhs, rhs);
        }
    }
}
} // namespace gal
