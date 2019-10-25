#pragma once

#include "algebra.hpp"
#include "null_algebra.hpp"

#include <type_traits>

namespace gal
{
enum class expr_op
{
    //////////////////////
    // Unary operations //
    //////////////////////

    identity, // Special op used to inject the first operand into an expression template (expression-leaves)
    negate,
    reverse,
    poincare_dual,
    clifford_conjugate,

    ///////////////////////
    // Binary operations //
    ///////////////////////

    // Scalar operations
    scalar_sum,
    scalar_product,

    // Arithmetric
    sum,
    difference,

    // General products
    geometric,
    sandwich,

    // Coordinate-free products
    exterior,
    regressive,

    // Metric products
    contract,
    symmetric_inner,

    // Projection operators
    extract, // Extract a single multivector term corresponding to a variadic list of basis elements
    select,  // Select all terms corresponding to a specified grade

    ///////////////////////////////////
    // Special non-linear operations //
    ///////////////////////////////////

    fast_exp, // Closed form exponential
};

template <int num, int den = 1>
struct frac_t
{
    [[nodiscard]] constexpr static rat q() noexcept
    {
        return {num, den};
    }
};

template <int num, int den = 1>
constexpr inline frac_t<num, den> frac;

// If O is a unary operation, T2 should be void
template <expr_op O, typename T1, typename T2 = void>
struct expr
{
    using value_t               = typename T1::value_t;
    using algebra_t             = typename T1::algebra_t;
    constexpr static expr_op op = O;
    using lhs_t                 = T1;
    using rhs_t                 = T2;
};

template <expr_op O, typename T>
struct expr<O, T, void>
{
    using value_t               = typename T::value_t;
    using algebra_t             = typename T::algebra_t;
    constexpr static expr_op op = O;
    using lhs_t                 = T;
};

template <typename T, uint32_t ID>
struct expr<expr_op::identity, T, std::integral_constant<uint32_t, ID>>
{
    constexpr static expr_op op = expr_op::identity;
    using value_t               = typename T::value_t;
    using algebra_t             = typename T::algebra_t;
    constexpr static auto lhs   = T::ie(ID);
};

template <typename T, uint8_t... N>
struct expr<expr_op::extract, T, std::integer_sequence<uint8_t, N...>>
{
    using value_t                                               = typename T::value_t;
    using algebra_t                                             = typename T::algebra_t;
    constexpr static expr_op op                                 = expr_op::extract;
    using lhs_t                                                 = T;
    constexpr static std::array<uint8_t, sizeof...(N)> elements = {N...};
};

struct fast_exp_tag
{};

template <typename T>
struct expr<expr_op::fast_exp, T, fast_exp_tag>
{
    using value_t               = typename T::value_t;
    using algebra_t             = typename T::algebra_t;
    constexpr static expr_op op = expr_op::fast_exp;
    using lhs_t                 = T;
};

template <typename T, uint8_t G>
struct expr<expr_op::select, T, std::integral_constant<uint8_t, G>>
{
    using value_t                  = typename T::value_t;
    using algebra_t                = typename T::algebra_t;
    constexpr static expr_op op    = expr_op::select;
    using lhs_t                    = T;
    constexpr static uint8_t grade = G;
};

// TODO: Support scaling operation by a constant

template <expr_op O, typename T1, typename T2>
[[nodiscard]] constexpr auto operator~(expr<O, T1, T2>) noexcept
{
    return expr<expr_op::reverse, expr<O, T1, T2>, void>{};
}

template <expr_op O, typename T1, typename T2>
[[nodiscard]] constexpr auto operator-(expr<O, T1, T2>) noexcept
{
    return expr<expr_op::negate, expr<O, T1, T2>, void>{};
}

template <expr_op O, typename T1, typename T2>
[[nodiscard]] constexpr auto operator!(expr<O, T1, T2>) noexcept
{
    return expr<expr_op::poincare_dual, expr<O, T1, T2>, void>{};
}

template <expr_op O1, typename T1, typename T2, expr_op O2, typename S1, typename S2>
[[nodiscard]] constexpr auto operator+(expr<O1, T1, T2>, expr<O2, S1, S2>) noexcept
{
    return expr<expr_op::sum, expr<O1, T1, T2>, expr<O2, S1, S2>>{};
}

template <expr_op O, typename T1, typename T2, int F1, int F2>
[[nodiscard]] constexpr auto operator+(frac_t<F1, F2>, expr<O, T1, T2>)
{
    return expr<expr_op::scalar_sum, expr<O, T1, T2>, frac_t<F1, F2>>{};
}

template <expr_op O, typename T1, typename T2, int F1, int F2>
[[nodiscard]] constexpr auto operator+(expr<O, T1, T2>, frac_t<F1, F2>)
{
    return expr<expr_op::scalar_sum, expr<O, T1, T2>, frac_t<F1, F2>>{};
}

template <expr_op O, typename T1, typename T2, int F1, int F2>
[[nodiscard]] constexpr auto operator-(frac_t<F1, F2>, expr<O, T1, T2>)
{
    return expr<expr_op::scalar_sum, expr<expr_op::negate, expr<O, T1, T2>>, frac_t<F1, F2>>{};
}

template <expr_op O, typename T1, typename T2, int F1, int F2>
[[nodiscard]] constexpr auto operator-(expr<O, T1, T2>, frac_t<F1, F2>)
{
    return expr<expr_op::scalar_sum, expr<O, T1, T2>, frac_t<-F1, F2>>{};
}

template <expr_op O, typename T1, typename T2, int F1, int F2>
[[nodiscard]] constexpr auto operator*(expr<O, T1, T2>, frac_t<F1, F2>)
{
    return expr<expr_op::scalar_product, expr<O, T1, T2>, frac_t<F1, F2>>{};
}

template <expr_op O, typename T1, typename T2, int F1, int F2>
[[nodiscard]] constexpr auto operator*(frac_t<F1, F2>, expr<O, T1, T2>)
{
    return expr<expr_op::scalar_product, expr<O, T1, T2>, frac_t<F1, F2>>{};
}

template <expr_op O, typename T1, typename T2, int F1, int F2>
[[nodiscard]] constexpr auto operator/(expr<O, T1, T2>, frac_t<F1, F2>)
{
    return expr<expr_op::scalar_product, expr<O, T1, T2>, frac_t<F2, F1>>{};
}

template <expr_op O1, typename T1, typename T2, expr_op O2, typename S1, typename S2>
[[nodiscard]] constexpr auto operator-(expr<O1, T1, T2>, expr<O2, S1, S2>) noexcept
{
    return expr<expr_op::difference, expr<O1, T1, T2>, expr<O2, S1, S2>>{};
}

template <expr_op O1, typename T1, typename T2, expr_op O2, typename S1, typename S2>
[[nodiscard]] constexpr auto operator*(expr<O1, T1, T2>, expr<O2, S1, S2>)noexcept
{
    return expr<expr_op::geometric, expr<O1, T1, T2>, expr<O2, S1, S2>>{};
}

template <expr_op O1, typename T1, typename T2, expr_op O2, typename S1, typename S2>
[[nodiscard]] constexpr auto operator%(expr<O1, T1, T2>, expr<O2, S1, S2>) noexcept
{
    return expr<expr_op::sandwich, expr<O1, T1, T2>, expr<O2, S1, S2>>{};
}

template <expr_op O1, typename T1, typename T2, expr_op O2, typename S1, typename S2>
[[nodiscard]] constexpr auto operator^(expr<O1, T1, T2>, expr<O2, S1, S2>) noexcept
{
    return expr<expr_op::exterior, expr<O1, T1, T2>, expr<O2, S1, S2>>{};
}

template <expr_op O1, typename T1, typename T2, expr_op O2, typename S1, typename S2>
[[nodiscard]] constexpr auto operator&(expr<O1, T1, T2>, expr<O2, S1, S2>)noexcept
{
    return expr<expr_op::regressive, expr<O1, T1, T2>, expr<O2, S1, S2>>{};
}

template <expr_op O1, typename T1, typename T2, expr_op O2, typename S1, typename S2>
[[nodiscard]] constexpr auto operator>>(expr<O1, T1, T2>, expr<O2, S1, S2>) noexcept
{
    return expr<expr_op::contract, expr<O1, T1, T2>, expr<O2, S1, S2>>{};
}

template <expr_op O1, typename T1, typename T2, expr_op O2, typename S1, typename S2>
[[nodiscard]] constexpr auto operator|(expr<O1, T1, T2>, expr<O2, S1, S2>) noexcept
{
    return expr<expr_op::symmetric_inner, expr<O1, T1, T2>, expr<O2, S1, S2>>{};
}

// This operation is ONLY well defined for bivectors and implements a closed form solution by decomposing the bivector
// into a euclidean and ideal line in order to normalize properly
template <expr_op O, typename T1, typename T2>
[[nodiscard]] constexpr auto fast_exp(expr<O, T1, T2>) noexcept
{
    return expr<expr_op::fast_exp, expr<O, T1, T2>, fast_exp_tag>{};
}

template <uint8_t... E>
struct extract
{
    template <expr_op O, typename T1, typename T2>
    [[nodiscard]] constexpr auto operator()(expr<O, T1, T2>) noexcept
    {
        return expr<expr_op::extract, expr<O, T1, T2>, std::integer_sequence<uint8_t, E...>>{};
    }
};

template <expr_op O, typename T1, typename T2>
[[nodiscard]] constexpr auto select(expr<O, T1, T2>, uint32_t grade) noexcept
{
    return expr<expr_op::select, expr<O, T1, T2>, void>{grade};
}

template <typename exp_t>
[[nodiscard]] constexpr auto reify() noexcept
{
    // Recursive function which compiles an expression into a reification table suitable for further processing needed
    // to evaluate all polynomial terms.

    // Base case
    if constexpr (exp_t::op == expr_op::identity)
    {
        if constexpr (detail::uses_null_basis<typename exp_t::algebra_t>)
        {
            constexpr auto out = detail::to_natural_basis(exp_t::lhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }
        else
        {
            return exp_t::lhs;
        }
    }
    else if constexpr (exp_t::op == expr_op::negate)
    {
        return detail::negate(reify<typename exp_t::lhs_t>());
    }
    else if constexpr (exp_t::op == expr_op::reverse)
    {
        return detail::reverse(reify<typename exp_t::lhs_t>());
    }
    else if constexpr (exp_t::op == expr_op::poincare_dual)
    {
        return detail::poincare_dual(reify<typename exp_t::lhs_t>());
    }
    else if constexpr (exp_t::op == expr_op::clifford_conjugate)
    {
        // TODO: implement me
        // return ::gal::clifford_conjugate(reify(T1{}));
    }
    else if constexpr (exp_t::op == expr_op::scalar_sum)
    {
        return detail::scalar_sum(exp_t::rhs_t::q(), reify<typename exp_t::lhs_t>());
    }
    else if constexpr (exp_t::op == expr_op::scalar_product)
    {
        return detail::scalar_product(exp_t::rhs_t::q(), reify<typename exp_t::lhs_t>());
    }
    else if constexpr (exp_t::op == expr_op::extract)
    {
        constexpr auto out = detail::extract(reify<typename exp_t::lhs_t>(), exp_t::elements);
        return out.template resize<out.size.ind, out.size.mon, out.size.term>();
    }
    else if constexpr (exp_t::op == expr_op::select)
    {
        // TODO: implement me
    }
    else if constexpr (exp_t::op == expr_op::fast_exp)
    {
        constexpr auto out = detail::fast_exp(reify<typename exp_t::lhs_t>());
        return out.template resize<out.size.ind, out.size.mon, out.size.term>();
    }
    else // Binary operation
    {
        constexpr auto lhs = reify<typename exp_t::lhs_t>();
        constexpr auto rhs = reify<typename exp_t::rhs_t>();

        if constexpr (exp_t::op == expr_op::sum)
        {
            constexpr auto out = detail::sum(lhs, rhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }
        else if constexpr (exp_t::op == expr_op::difference)
        {
            constexpr auto n_rhs = detail::negate(rhs);
            constexpr auto out   = detail::sum(lhs, n_rhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }
        else if constexpr (exp_t::op == expr_op::geometric)
        {
            constexpr auto out = detail::product(typename decltype(lhs)::algebra_t::geometric{}, lhs, rhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }
        else if constexpr (exp_t::op == expr_op::sandwich)
        {
            // rhs * lhs * ~rhs
            constexpr auto rhs_reverse = detail::reverse(rhs);
            constexpr auto temp  = detail::product(typename decltype(lhs)::algebra_t::geometric{}, lhs, rhs_reverse);
            constexpr auto temp2 = detail::product(typename decltype(lhs)::algebra_t::geometric{},
                                                   rhs,
                                                   temp.template resize<temp.size.ind, temp.size.mon, temp.size.term>());
            return temp2.template resize<temp2.size.ind, temp2.size.mon, temp2.size.term>();
        }
        else if constexpr (exp_t::op == expr_op::exterior)
        {
            constexpr auto out = detail::product(typename decltype(lhs)::algebra_t::exterior{}, lhs, rhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }
        else if constexpr (exp_t::op == expr_op::regressive)
        {
            // TODO: check if both lhs and rhs are dual
            constexpr auto lhs_dual = detail::poincare_dual(lhs);
            constexpr auto rhs_dual = detail::poincare_dual(rhs);
            constexpr auto out = detail::product(typename decltype(lhs)::algebra_t::exterior{}, lhs_dual, rhs_dual);
            return detail::poincare_dual(out.template resize<out.size.ind, out.size.mon, out.size.term>());
        }
        else if constexpr (exp_t::op == expr_op::contract)
        {
            constexpr auto out = detail::product(typename decltype(lhs)::algebra_t::contract{}, lhs, rhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }
        else if constexpr (exp_t::op == expr_op::symmetric_inner)
        {
            constexpr auto out = detail::product(typename decltype(lhs)::algebra_t::symmetric_inner{}, lhs, rhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }
    }
}
} // namespace gal
