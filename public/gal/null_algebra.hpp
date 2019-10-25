#pragma once

#include "algebra.hpp"

// To compute using a null-basis, the geometric product behaves differently in the sense that the geometric product of
// two basis elements will not necessarily result in a single term. In particular, no * ni = -1 + no^ni.

// This file contains routines for converting an ie (indeterminate form) multivector to and from the null basis from the
// natural basis.

namespace gal
{
namespace detail
{
    // Null basis specialization (e.g. true in the case of conformal geometric algebra)
    template <typename A>
    constexpr inline bool uses_null_basis = false;

    // By convention, the last two generators are always the null elements, with the point at infinity coming last.
    // The size estimate of the returned vector is conservative so the result is stored at compile time, it is recommend
    // that `shrink` be invoked.
    template <typename A, width_t I, width_t M, width_t T>
    [[nodiscard]] constexpr auto to_null_basis(mv<A, I, M, T> const& in) noexcept
    {
        mv<A, I, M, T> lhs{};

        mv<A, 2 * I, 2 * M, 2 * T> temp{};

        // The second to last element we denote ep, and the last element we denote en
        auto ep = 1 << (A::metric_t::dimension - 2);
        auto en = 1 << (A::metric_t::dimension - 1);
        auto enp = ep ^ en;

        for (auto it = in.cbegin(); it != in.cend(); ++it)
        {
            if (it->element & enp)
            {
                if (it->element & en && it->element & ep)
                {
                    // ep ^ en = no ^ ni
                    // No change needed
                    lhs.push(it, one, it->element);
                }
                else if (it->element & ep)
                {
                    // ep = no - 1/2 ni
                    temp.push(it, one, it->element);
                    temp.push(it, minus_one_half, it->element ^ enp);
                }
                else // if (it->element & en)
                {
                    // en = no + 1/2 ni
                    temp.push(it, one, it->element ^ enp);
                    temp.push(it, one_half, it->element);
                }
            }
            else
            {
                lhs.push(it, one, it->element);
            }
        }

        // The temp mv now needs to be sorted so that coincident terms can be collated
        sort(temp.terms.begin(), temp.terms.begin() + temp.size.term);

        // Transform the array of monomials to monomial views so they can be independently sorted as well
        std::array<mon_view, 2 * M> mon_views;
        for (width_t i = 0; i != temp.size.mon; ++i)
        {
            mon_views[i] = mon_view{temp.mons[i], temp.inds.begin()};
        }

        mv<A, 2 * I, 2 * M, 2 * T> rhs{};
        collate(temp.terms.begin(),
                temp.terms.begin() + temp.size.term,
                mon_views.begin(),
                temp.inds.begin(),
                rhs.terms.begin(),
                rhs.mons.begin(),
                rhs.inds.begin(),
                rhs.size);
        return sum(lhs, rhs);
    }

    template <typename A, width_t I, width_t M, width_t T>
    [[nodiscard]] constexpr auto to_natural_basis(mv<A, I, M, T> const& in) noexcept
    {
        mv<A, I, M, T> lhs{};

        mv<A, 2 * I, 2 * M, 2 * T> temp{};

        // The second to last element we denote no, and the last element we denote ni
        auto no = 1 << (A::metric_t::dimension - 2);
        auto ni = 1 << (A::metric_t::dimension - 1);
        auto noi = no ^ ni;

        for (auto it = in.cbegin(); it != in.cend(); ++it)
        {
            if (it->element & noi)
            {
                if (it->element & no && it->element & ni)
                {
                    // ep ^ en = no ^ ni
                    // No change needed
                    lhs.push(it, one, it->element);
                }
                else if (it->element & no)
                {
                    // no = 1/2 ep + 1/2 en
                    temp.push(it, one_half, it->element);
                    temp.push(it, one_half, it->element ^ noi);
                }
                else // if (it->element & ni)
                {
                    // ni = en - ep
                    temp.push(it, minus_one, it->element ^ noi);
                    temp.push(it, one, it->element);
                }
            }
            else
            {
                lhs.push(it, one, it->element);
            }
        }

        // The temp mv now needs to be sorted so that coincident terms can be collated
        sort(temp.terms.begin(), temp.terms.begin() + temp.size.term);

        // Transform the array of monomials to monomial views so they can be independently sorted as well
        std::array<mon_view, 2 * M> mon_views;
        for (width_t i = 0; i != temp.size.mon; ++i)
        {
            mon_views[i] = mon_view{temp.mons[i], temp.inds.begin()};
        }

        mv<A, 2 * I, 2 * M, 2 * T> rhs{};
        collate(temp.terms.begin(),
                temp.terms.begin() + temp.size.term,
                mon_views.begin(),
                temp.inds.begin(),
                rhs.terms.begin(),
                rhs.mons.begin(),
                rhs.inds.begin(),
                rhs.size);
        return sum(lhs, rhs);

    }
}
}