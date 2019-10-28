#pragma once

#include "algebra.hpp"
#include "null_algebra.hpp"

#include <type_traits>

namespace gal
{
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

namespace detail
{
    template <typename T>
    struct expr
    {};

    struct id_t;

    // Identity
    template <typename T, uint32_t ID>
    struct expr_id : public expr<expr_id<T, ID>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            if constexpr (detail::uses_null_basis<typename T::algebra_t>)
            {
                constexpr auto out = detail::to_natural_basis(T::ie(ID));
                return out.template resize<out.size.ind, out.size.mon, out.size.term>();
            }
            else
            {
                return T::ie(ID);
            }
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            if constexpr (detail::uses_null_basis<typename T::algebra_t>)
            {
                return detail::to_natural_basis(T::ie(ID));
            }
            else
            {
                return T::ie(ID);
            }
        }
#endif
    };

    // Negate
    template <typename T>
    struct expr_neg : public expr<expr_neg<T>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            return detail::negate(T::reify());
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            return detail::negate(T::reify_debug());
        }
#endif
    };

    // Reverse
    template <typename T>
    struct expr_rev : public expr<expr_rev<T>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            return detail::reverse(T::reify());
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            return detail::reverse(T::reify_debug());
        }
#endif
    };

    // Poincare dual
    template <typename T>
    struct expr_pd : public expr<expr_pd<T>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            return detail::poincare_dual(T::reify());
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            return detail::poincare_dual(T::reify_debug());
        }
#endif
    };

    // Shift by a constant value
    template <typename T, int Num, int Den>
    struct expr_sh : public expr<expr_sh<T, Num, Den>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            return detail::shift(rat{Num, Den}, T::reify());
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            return detail::shift(rat{Num, Den}, T::reify_debug());
        }
#endif
    };

    // Scale by a constant value
    template <typename T, int Num, int Den>
    struct expr_sc : public expr<expr_sc<T, Num, Den>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            return detail::scale(rat{Num, Den}, T::reify());
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            return detail::scale(rat{Num, Den}, T::reify_debug());
        }
#endif
    };

    // Multivector sum
    template <typename L, typename R>
    struct expr_sum : public expr<expr_sum<L, R>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            constexpr auto lhs = L::reify();
            constexpr auto rhs = R::reify();
            constexpr auto out = detail::sum(lhs, rhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            auto lhs = L::reify_debug();
            auto rhs = R::reify_debug();
            return detail::sum(lhs, rhs);
        }
#endif
    };

    // Multivector difference
    template <typename L, typename R>
    struct expr_diff : public expr<expr_diff<L, R>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            constexpr auto lhs     = L::reify();
            constexpr auto rhs     = R::reify();
            constexpr auto neg_rhs = detail::negate(rhs);
            constexpr auto out     = detail::sum(lhs, neg_rhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            auto lhs = L::reify_debug();
            auto rhs = R::reify_debug();
            auto neg_rhs = detail::negate(rhs);
            return detail::sum(lhs, neg_rhs);
        }
#endif
    };

    // Geometric product
    template <typename L, typename R>
    struct expr_gp : public expr<expr_gp<L, R>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            constexpr auto lhs = L::reify();
            constexpr auto rhs = R::reify();
            constexpr typename decltype(lhs)::algebra_t::geometric gp{};
            constexpr auto out = detail::product(gp, lhs, rhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            auto lhs = L::reify_debug();
            auto rhs = R::reify_debug();
            constexpr typename decltype(lhs)::algebra_t::geometric gp{};
            return detail::product(gp, lhs, rhs);
        }
#endif
    };

    // Sandwich operator
    template <typename L, typename R>
    struct expr_sw : public expr<expr_sw<L, R>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            constexpr auto lhs     = L::reify();
            constexpr auto rhs     = R::reify();
            constexpr auto rhs_rev = detail::reverse(rhs);
            constexpr typename decltype(lhs)::algebra_t::geometric gp{};
            constexpr auto tmp  = detail::product(gp, lhs, rhs_rev);
            constexpr auto tmp2 = tmp.template resize<tmp.size.ind, tmp.size.mon, tmp.size.term>();
            constexpr auto out  = detail::product(gp, rhs, tmp2);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            auto lhs = L::reify_debug();
            auto rhs = R::reify_debug();
            auto rhs_rev = detail::reverse(rhs);
            constexpr typename decltype(lhs)::algebra_t::geometric gp{};
            auto tmp = detail::product(gp, lhs, rhs_rev);
            return detail::product(gp, rhs, tmp);
        }
#endif
    };

    // Exterior product
    template <typename L, typename R>
    struct expr_ep : public expr<expr_ep<L, R>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            constexpr auto lhs = L::reify();
            constexpr auto rhs = R::reify();
            constexpr typename decltype(lhs)::algebra_t::exterior ep{};
            constexpr auto out = detail::product(ep, lhs, rhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            auto lhs = L::reify_debug();
            auto rhs = R::reify_debug();
            constexpr typename decltype(lhs)::algebra_t::exterior ep{};
            return detail::product(ep, lhs, rhs);
        }
#endif
    };

    // Regressive product
    template <typename L, typename R>
    struct expr_rp : public expr<expr_rp<L, R>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            constexpr auto lhs = detail::poincare_dual(L::reify());
            constexpr auto rhs = detail::poincare_dual(R::reify());
            constexpr typename decltype(lhs)::algebra_t::exterior ep{};
            constexpr auto out = detail::product(ep, lhs, rhs);
            return detail::poincare_dual(out.template resize<out.size.ind, out.size.mon, out.size.term>());
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            auto lhs = detail::poincare_dual(L::reify_debug());
            auto rhs = detail::poincare_dual(R::reify_debug());
            constexpr typename decltype(lhs)::algebra_t::exterior ep{};
            return detail::poincare_dual(detail::product(ep, lhs, rhs));
        }
#endif
    };

    // Left contraction
    template <typename L, typename R>
    struct expr_lc : public expr<expr_lc<L, R>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            constexpr auto lhs = L::reify();
            constexpr auto rhs = R::reify();
            constexpr typename decltype(lhs)::algebra_t::contract lc{};
            constexpr auto out = detail::product(lc, lhs, rhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            auto lhs = L::reify_debug();
            auto rhs = R::reify_debug();
            constexpr typename decltype(lhs)::algebra_t::contract lc{};
            return detail::product(lc, lhs, rhs);
        }
#endif
    };

    // Symmetric inner product
    template <typename L, typename R>
    struct expr_sip : public expr<expr_sip<L, R>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            constexpr auto lhs = L::reify();
            constexpr auto rhs = R::reify();
            constexpr typename decltype(lhs)::algebra_t::symmetric_inner si{};
            constexpr auto out = detail::product(si, lhs, rhs);
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            auto lhs = L::reify_debug();
            auto rhs = R::reify_debug();
            constexpr typename decltype(lhs)::algebra_t::symmetric_inner si{};
            return detail::product(si, lhs, rhs);
        }
#endif
    };

    // Scalar product
    template <typename L, typename R>
    struct expr_sp : public expr<expr_sp<L, R>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            constexpr auto lhs = L::reify();
            constexpr auto rhs = R::reify();
            using algebra_t    = typename decltype(lhs)::algebra_t;
            constexpr typename decltype(lhs)::algebra_t::symmetric_inner si{};
            constexpr auto out = detail::extract(detail::product(si, lhs, rhs), {0});
            return out.template resize<out.size.ind, out.size.mon, 1>();
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            auto lhs = L::reify_debug();
            auto rhs = R::reify_debug();
            constexpr typename decltype(lhs)::algebra_t::symmetric_inner si{};
            return detail::extract(detail::product(si, lhs, rhs), {0});
        }
#endif
    };

    // Extract
    template <typename T, uint8_t... N>
    struct expr_ex : public expr<expr_ex<T, N...>>
    {
        [[nodiscard]] constexpr static auto reify() noexcept
        {
            constexpr auto out = detail::extract(T::reify(), std::array<uint8_t, sizeof...(N)>{N...});
            return out.template resize<out.size.ind, out.size.mon, out.size.term>();
        }

#ifdef GAL_DEBUG
        [[nodiscard]] static auto reify_debug() noexcept
        {
            return detail::extract(T::reify_debug(), std::array<uint8_t, sizeof...(N)>{N...});
        }
#endif
    };
} // namespace detail

template <uint8_t... E>
struct extract
{
    template <typename T>
    [[nodiscard]] constexpr auto operator()(detail::expr<T>) noexcept
    {
        return detail::expr_ex<T, E...>{};
    }
};

template <typename T>
[[nodiscard]] constexpr detail::expr_rev<T> operator~(detail::expr<T> const&) noexcept
{
    return {};
}

template <typename T>
[[nodiscard]] constexpr detail::expr_neg<T> operator-(detail::expr<T> const&) noexcept
{
    return {};
}

template <typename T>
[[nodiscard]] constexpr detail::expr_pd<T> operator!(detail::expr<T> const&) noexcept
{
    return {};
}

template <typename T1, typename T2>
[[nodiscard]] constexpr detail::expr_sum<T1, T2> operator+(detail::expr<T1> const&, detail::expr<T2> const&) noexcept
{
    return {};
}

template <typename T, int N, int D>
[[nodiscard]] constexpr auto operator+(frac_t<N, D>, detail::expr<T> const&) noexcept
{
    if constexpr (D < 0)
    {
        return detail::expr_sh<T, -N, -D>{};
    }
    else
    {
        return detail::expr_sh<T, N, D>{};
    }
}

template <typename T, int N, int D>
[[nodiscard]] constexpr auto operator+(detail::expr<T> const&, frac_t<N, D>) noexcept
{
    if constexpr (D < 0)
    {
        return detail::expr_sh<T, -N, -D>{};
    }
    else
    {
        return detail::expr_sh<T, N, D>{};
    }
}

template <typename T1, typename T2>
[[nodiscard]] constexpr detail::expr_diff<T1, T2> operator-(detail::expr<T1> const&, detail::expr<T2> const&) noexcept
{
    return {};
}

template <typename T, int N, int D>
[[nodiscard]] constexpr auto operator-(frac_t<N, D> const& lhs, detail::expr<T> const& rhs) noexcept
{
    return lhs + -rhs;
}

template <typename T, int N, int D>
[[nodiscard]] constexpr auto operator-(detail::expr<T> const&, frac_t<N, D>) noexcept
{
    if constexpr (D < 0)
    {
        return detail::expr_sh<T, N, -D>{};
    }
    else
    {
        return detail::expr_sh<T, -N, D>{};
    }
}

template <typename T, int N, int D>
[[nodiscard]] constexpr auto operator*(frac_t<N, D>, detail::expr<T> const&)noexcept
{
    if constexpr (D < 0)
    {
        return detail::expr_sc<T, -N, -D>{};
    }
    else
    {
        return detail::expr_sc<T, N, D>{};
    }
}

template <typename T, int N, int D>
[[nodiscard]] constexpr auto operator*(detail::expr<T> const&, frac_t<N, D>)noexcept
{
    if constexpr (D < 0)
    {
        return detail::expr_sc<T, -N, -D>{};
    }
    else
    {
        return detail::expr_sc<T, N, D>{};
    }
}

template <typename T, int N, int D>
[[nodiscard]] constexpr auto operator/(detail::expr<T> const&, frac_t<N, D>) noexcept
{
    if constexpr (N < 0)
    {
        return detail::expr_sc<T, -D, -N>{};
    }
    else
    {
        return detail::expr_sc<T, D, N>{};
    }
}

template <typename T1, typename T2>
[[nodiscard]] constexpr detail::expr_gp<T1, T2> operator*(detail::expr<T1> const&, detail::expr<T2> const&)noexcept
{
    return {};
}

template <typename T1, typename T2>
[[nodiscard]] constexpr detail::expr_sw<T1, T2> operator%(detail::expr<T1> const&, detail::expr<T2> const&) noexcept
{
    return {};
}

template <typename T1, typename T2>
[[nodiscard]] constexpr detail::expr_ep<T1, T2> operator^(detail::expr<T1> const&, detail::expr<T2> const&) noexcept
{
    return {};
}

template <typename T1, typename T2>
[[nodiscard]] constexpr detail::expr_rp<T1, T2> operator&(detail::expr<T1> const&, detail::expr<T2> const&)noexcept
{
    return {};
}

template <typename T1, typename T2>
[[nodiscard]] constexpr detail::expr_lc<T1, T2> operator>>(detail::expr<T1> const&, detail::expr<T2> const&) noexcept
{
    return {};
}

template <typename T1, typename T2>
[[nodiscard]] constexpr detail::expr_sip<T1, T2> operator|(detail::expr<T1> const&, detail::expr<T2> const&) noexcept
{
    return {};
}

template <typename T1, typename T2>
[[nodiscard]] constexpr detail::expr_sp<T1, T2> scalar_product(detail::expr<T1> const&, detail::expr<T2> const&) noexcept
{
    return {};
}
} // namespace gal
