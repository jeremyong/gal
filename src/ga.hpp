#pragma once

#include "finite_algebra.hpp"
#include "utility.hpp"

namespace gal
{
template <typename E>
[[nodiscard]] constexpr static auto grade() noexcept
{
    return count_bits(E::value);
}

namespace ga
{
    template <typename Metric>
    struct algebra
    {
        template <typename... T>
        using term_t = term<T...>;

        struct contract
        {
            constexpr static bool has_order_preserving_product = false;

            template <typename E1, typename E2, typename... M1, typename... M2>
            [[nodiscard]] constexpr static auto product(term_t<E1, M1...> lhs, term_t<E2, M2...> rhs) noexcept
            {
                constexpr auto contraction = contract_product<E1::value, E2::value>();
                if constexpr (contraction.second == 0)
                {
                    return multivector<void>{};
                }
                else
                {
                    using element_t = element<contraction.first>;
                    using term_t    = decltype(typename scale<contraction.second, term_t<element_t, M1...>>::type{}
                                            * term_t<element_t, M2...>{});
                    return multivector<void, term_t>{};
                }
            }

            template <size_t E1, size_t E2>
            [[nodiscard]] constexpr static std::pair<size_t, int> contract_product() noexcept
            {
                if constexpr (E1 == 0)
                {
                    return {E2, 1};
                }
                else if constexpr (E1 > E2)
                {
                    return {0, 0};
                }
                else
                {
                    size_t lhs = E1;
                    size_t rhs = E2;
                    int parity = 1;
                    while (lhs > 0)
                    {
                        auto lhs_element = leading_index(lhs);
                        if (rhs & (1 << lhs_element))
                        {
                            auto en = Metric::element_norm(lhs_element);
                            if (en == 0)
                            {
                                return {0, 0};
                            }

                            if (count_bits(rhs & ((1 << lhs_element) - 1)))
                            {
                                parity *= -1;
                            }
                            lhs &= ~(1 << lhs_element);
                            rhs &= ~(1 << lhs_element);

                            // Account for negative elements in the metric signature
                            parity *= en;
                        }
                        else
                        {
                            return {0, 0};
                        }
                    }

                    return {rhs, parity};
                }
            }
        };

        struct exterior
        {
            constexpr static bool has_order_preserving_product = false;

            template <typename E1, typename E2, typename... M1, typename... M2>
            [[nodiscard]] constexpr static auto product(term_t<E1, M1...> lhs, term_t<E2, M2...> rhs) noexcept
            {
                constexpr auto wedge = exterior_product<E1::value, E2::value>();
                if constexpr (wedge.second == 0)
                {
                    return multivector<void>{};
                }
                else
                {
                    using element_t = element<wedge.first>;
                    using term_t    = decltype(typename scale<wedge.second, term_t<element_t, M1...>>::type{}
                                            * term_t<element_t, M2...>{});
                    return multivector<void, term_t>{};
                }
            }

            template <size_t E1, size_t E2>
            [[nodiscard]] constexpr static std::pair<size_t, int> exterior_product() noexcept
            {
                if constexpr (E1 == 0)
                {
                    return {E2, 1};
                }
                else if constexpr (E2 == 0)
                {
                    return {E1, 1};
                }
                else
                {
                    constexpr auto intersection = E1 & E2;
                    if constexpr (intersection != 0)
                    {
                        return {0, 0};
                    }
                    else
                    {
                        constexpr size_t element = E1 | E2;
                        size_t swaps             = 0;
                        size_t lhs               = E1;
                        size_t rhs               = E2;

                        while (lhs > 0)
                        {
                            auto lhs_element = leading_index(lhs);
                            swaps += count_bits(rhs & ((1 << lhs_element) - 1));
                            lhs &= ~(1 << lhs_element);
                            rhs |= 1 << lhs_element;
                        }
                        return {element, swaps % 2 == 0 ? 1 : -1};
                    }
                }
            }
        };

        struct geometric
        {
            constexpr static bool has_order_preserving_product = false;

            template <typename E1, typename E2, typename... M1, typename... M2>
            [[nodiscard]] constexpr static auto product(term_t<E1, M1...> lhs, term_t<E2, M2...> rhs) noexcept
            {
                constexpr auto gp = geometric_product<E1::value, E2::value>();
                if constexpr (gp.second == 0)
                {
                    return multivector<void>{};
                }
                else
                {
                    using element_t = element<gp.first>;
                    using term_t    = decltype(typename scale<gp.second, term_t<element_t, M1...>>::type{}
                                            * term_t<element_t, M2...>{});
                    return multivector<void, term_t>{};
                }
            }

            template <size_t E1, size_t E2>
            [[nodiscard]] constexpr static std::pair<size_t, int> geometric_product() noexcept
            {
                if constexpr (E1 == 0)
                {
                    return {E2, 1};
                }
                else if constexpr (E2 == 0)
                {
                    return {E1, 1};
                }
                else
                {
                    constexpr size_t element = E1 ^ E2;
                    size_t swaps             = 0;
                    size_t lhs               = E1;
                    size_t rhs               = E2;

                    while (lhs > 0)
                    {
                        auto lhs_element = leading_index(lhs);
                        swaps += count_bits(rhs & ((1 << lhs_element) - 1));
                        lhs &= ~(1 << lhs_element);
                        if (rhs & (1 << lhs_element))
                        {
                            auto en = Metric::element_norm(lhs_element);
                            if (en == 0)
                            {
                                return {0, 0};
                            }
                            else if (en == -1)
                            {
                                ++swaps;
                            }
                            rhs &= ~(1 << lhs_element);
                        }
                        else
                        {
                            rhs |= 1 << lhs_element;
                        }
                    }
                    return {element, swaps % 2 == 0 ? 1 : -1};
                }
            }
        };
    };

    template <size_t Dim, size_t E>
    [[nodiscard]] constexpr int poincare_complement_parity()
    {
        int swaps = 0;
        size_t e  = E;
        int i     = 0;
        int grade = count_bits(E);
        while (e > 0)
        {
            if ((e & 1) == 0)
            {
                swaps += grade;
            }
            else
            {
                --grade;
            }
            ++i;
            e = e >> i;
        }
        if ((swaps & 1) == 1)
        {
            return -1;
        }
        else
        {
            return 1;
        }
    }

    template <size_t Dim, typename E, typename... As>
    [[nodiscard]] constexpr auto poincare_complement(term<E, As...>)
    {
        constexpr size_t complement = ((1 << Dim) - 1) ^ E::value;
        // We require that the concatenation of E and its complement form an even permutation
        // of the basis element sequence.
        constexpr int parity = poincare_complement_parity<Dim, E::value>();
        if constexpr (parity == 1)
        {
            return term<element<complement>, As...>{};
        }
        else
        {
            return typename scale<-1, term<element<complement>, As...>>::type{};
        }
    }
} // namespace ga

// The geometric algebra also leverages additional operations
// - reversion
// - grade involution
// - grade selection
// - dual (contraction onto the unit pseudoscalar)

// The following operations are general regardless of metric so we keep them in the upper namespace
// Term reversion
template <typename E, typename... M>
[[nodiscard]] constexpr auto operator~(term<E, M...> t) noexcept
{
    constexpr auto g    = grade<E>();
    constexpr bool flip = (g * (g - 1) / 2) % 2 == 1;
    if constexpr (flip)
    {
        return -t;
    }
    else
    {
        return t;
    }
}

// General multivector reversion
template <typename... T>
[[nodiscard]] constexpr auto operator~(multivector<void, T...>) noexcept
{
    return multivector<void, decltype(~T{})...>{};
}

template <typename Metric>
struct pseudoscalar
{
    constexpr static multivector<void, term<element<(1 << Metric::dimension) - 1>, monomial<one>>> value{};
    constexpr static auto inverse = ~value;
};

template <typename Metric, typename M>
[[nodiscard]] constexpr auto polarity_dual(M input) noexcept
{
    return input >> pseudoscalar<Metric>::inverse;
}

// Poincaré dual (does not require a non-degenerate metric)
template <typename Metric, typename... T>
[[nodiscard]] constexpr auto dual(multivector<void, T...>) noexcept
{
    if constexpr (sizeof...(T) == 0)
    {
        return multivector<void>{};
    }
    else
    {
        return ((multivector<void, decltype(ga::poincare_complement<Metric::dimension>(T{}))>{}) + ...);
    }
}

namespace detail
{
    template <size_t E, typename T, typename... Ts>
    [[nodiscard]] constexpr auto extract(multivector<void, T, Ts...>) noexcept
    {
        if constexpr (E == T::element_t::value)
        {
            return T{};
        }
        else if constexpr (sizeof...(Ts) == 0)
        {
            return term<element<0>, monomial<zero>>{};
        }
        else
        {
            return extract<E>(multivector<void, Ts...>{});
        }
    }
} // namespace detail

// Extract a term containing a given basis element from a multivector
template <size_t E, typename... T>
[[nodiscard]] constexpr auto extract(multivector<void, T...> mv) noexcept
{
    if constexpr (sizeof...(T) == 0)
    {
        return term<element<0>, monomial<zero>>{};
    }
    else
    {
        return detail::extract<E>(mv);
    }
}
} // namespace gal
