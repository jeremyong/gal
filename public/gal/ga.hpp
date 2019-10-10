#pragma once

#include "finite_algebra.hpp"
#include "utility.hpp"

#include <array>

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
        using metric_t = Metric;

        struct inner
        {
            constexpr static bool has_order_preserving_product = true;

            template <typename E1, typename E2, typename... M1, typename... M2>
            [[nodiscard]] constexpr static auto product(term<E1, M1...> lhs, term<E2, M2...> rhs) noexcept
            {
                constexpr auto inner = inner_product<E1::value, E2::value>();
                if constexpr (inner.second == 0)
                {
                    return multivector<void>{};
                }
                else
                {
                    using element_t = element<inner.first>;
                    if constexpr (inner.second == 1)
                    {
                        using term_t = decltype(term<element_t, M1...>{} * term<element_t, M2...>{});
                        return multivector<void, term_t>{};
                    }
                    else
                    {
                        using term_t = decltype(rational<inner.second>{} * term<element_t, M1...>{}
                                                * term<element_t, M2...>{});
                        return multivector<void, term_t>{};
                    }
                }
            }

            template <size_t E1, size_t E2>
            [[nodiscard]] constexpr static std::pair<size_t, int> inner_product() noexcept
            {
                if constexpr (E1 == 0 || E2 == 0)
                {
                    return {0, 0};
                }
                else
                {
                    constexpr auto gp = geometric::template geometric_product<E1, E2>();
                    constexpr auto desired_grade = static_cast<int>(E2) - static_cast<int>(E1);
                    if constexpr ((gp.first == desired_grade || gp.first == -desired_grade) && gp.second != 0)
                    {
                        return gp;
                    }
                    else
                    {
                        return {0, 0};
                    }
                }
            }
        };

        struct contract
        {
            constexpr static bool has_order_preserving_product = false;

            template <typename E1, typename E2, typename... M1, typename... M2>
            [[nodiscard]] constexpr static auto product(term<E1, M1...> lhs, term<E2, M2...> rhs) noexcept
            {
                constexpr auto contraction = contract_product<E1::value, E2::value>();
                if constexpr (contraction.second == 0)
                {
                    return multivector<void>{};
                }
                else
                {
                    using element_t = element<contraction.first>;
                    if constexpr (contraction.second == 1)
                    {
                        using term_t = decltype(term<element_t, M1...>{} * term<element_t, M2...>{});
                        return multivector<void, term_t>{};
                    }
                    else
                    {
                        using term_t = decltype(rational<contraction.second>{} * term<element_t, M1...>{}
                                                * term<element_t, M2...>{});
                        return multivector<void, term_t>{};
                    }
                }
            }

            template <size_t E1, size_t E2>
            [[nodiscard]] constexpr static std::pair<size_t, int> contract_product() noexcept
            {
                if constexpr (E1 == 0)
                {
                    return {E2, 1};
                }
                else if constexpr (count_bits(E1) > count_bits(E2))
                {
                    return {0, 0};
                }
                else
                {
                    size_t lhs = E1;
                    size_t rhs = E2;
                    int swaps = 0;
                    while (lhs > 0)
                    {
                        auto lhs_element = leading_index(lhs);
                        auto [index, dot] = Metric::intercept(lhs_element, rhs);
                        if (index == -1 || dot == 0)
                        {
                            return {0, 0};
                        }
                        else
                        {
                            swaps += count_bits(rhs & ((1 << index) - 1));
                            lhs &= ~(1 << lhs_element);
                            rhs &= ~(1 << index);
                            if (dot == -1)
                            {
                                ++swaps;
                            }
                        }
                    }

                    return {rhs, swaps % 2 == 0 ? 1 : -1};
                }
            }
        };

        struct exterior
        {
            constexpr static bool has_order_preserving_product = false;

            template <typename E1, typename E2, typename... M1, typename... M2>
            [[nodiscard]] constexpr static auto product(term<E1, M1...> lhs, term<E2, M2...> rhs) noexcept
            {
                constexpr auto wedge = exterior_product<E1::value, E2::value>();
                if constexpr (wedge.second == 0)
                {
                    return multivector<void>{};
                }
                else
                {
                    using element_t = element<wedge.first>;
                    if constexpr (wedge.second == 1)
                    {
                        using term_t = decltype(term<element_t, M1...>{} * term<element_t, M2...>{});
                        return multivector<void, term_t>{};
                    }
                    else
                    {
                        using term_t = decltype(rational<wedge.second>{} * term<element_t, M1...>{}
                                                * term<element_t, M2...>{});
                        return multivector<void, term_t>{};
                    }
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
            [[nodiscard]] constexpr static auto product(term<E1, M1...> lhs, term<E2, M2...> rhs) noexcept
            {
                if constexpr (!Metric::is_diagonal)
                {
                    if constexpr (Metric::multi_term_gp(E1::value, E2::value))
                    {
                        auto lhs2 = Metric::diagonalize(lhs);
                        auto rhs2 = Metric::diagonalize(rhs);
                        auto p = detail::product<typename algebra<typename Metric::base_metric>::geometric>(lhs2, rhs2);
                        return Metric::undiagonalize(p);
                    }
                    else
                    {
                        return product_diagonal(lhs, rhs);
                    }
                }
                else
                {
                    return product_diagonal(lhs, rhs);
                }
            }

            template <typename E1, typename E2, typename... M1, typename... M2>
            [[nodiscard]] constexpr static auto product_diagonal(term<E1, M1...> lhs, term<E2, M2...> rhs) noexcept
            {
                constexpr auto gp = geometric_product<E1::value, E2::value>();
                if constexpr (gp.second == 0)
                {
                    return multivector<void>{};
                }
                else
                {
                    using element_t = element<gp.first>;
                    if constexpr (gp.second == 1)
                    {
                        using term_t = decltype(term<element_t, M1...>{} * term<element_t, M2...>{});
                        return multivector<void, term_t>{};
                    }
                    else
                    {
                        using term_t
                            = decltype(rational<gp.second>{} * term<element_t, M1...>{} * term<element_t, M2...>{});
                        return multivector<void, term_t>{};
                    }
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
                        auto [index, dot] = Metric::intercept(lhs_element, rhs);
                        if (index == -1)
                        {
                            // Exterior
                            swaps += count_bits(rhs & ((1 << lhs_element) - 1));
                            rhs |= 1 << lhs_element;
                        }
                        else
                        {
                            // Contract
                            swaps += count_bits(rhs & ((1 << index) - 1));
                            if (dot == 0)
                            {
                                return {0, 0};
                            }
                            else if (dot == -1)
                            {
                                ++swaps;
                            }
                            rhs &= ~(1 << index);
                        }
                        
                        lhs &= ~(1 << lhs_element);
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
            e = e >> 1;
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
            return minus_one{} * term<element<complement>, As...>{};
        }
    }

} // namespace ga

// The geometric algebra also leverages additional operations
// - reversion
// - grade involution
// - grade selection
// - dual (contraction onto the unit pseudoscalar)

// Term reversion
template <typename E, typename... M>
[[nodiscard]] constexpr auto reverse(term<E, M...> t) noexcept
{
    constexpr auto g    = grade<E>();
    constexpr int flips = (g * (g - 1) / 2);
    if constexpr (flips % 2 == 1)
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
[[nodiscard]] constexpr auto reverse(multivector<void, T...>) noexcept
{
    return multivector<void, decltype(reverse(T{}))...>{};
}

template <typename Metric>
struct pseudoscalar
{
    constexpr static multivector<void, term<element<(1 << Metric::dimension) - 1>, monomial<one>>> value{};
    constexpr static multivector<
        void,
        term<element<(1 << Metric::dimension) - 1>,
             monomial<std::conditional_t<((Metric::dimension * (Metric::dimension - 1) / 2) + Metric::v) % 2 == 0, one, minus_one>>>>
        inverse{};
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

// Filter the term that matches the element E
template <size_t E, typename... T>
[[nodiscard]] constexpr auto filter(multivector<void, T...> mv) noexcept
{
    if constexpr (sizeof...(T) == 0)
    {
        return mv;
    }
    else
    {
        return (std::conditional_t<T::element_t::value == E, multivector<void>, multivector<void, T>>{} + ...);
    }
}

// Project a multivector onto a subset of the graded-basis provided in the template parameter sequence
template <size_t... E, typename... T>
[[nodiscard]] constexpr auto component_select(multivector<void, T...> mv) noexcept
{
    static_assert(sizeof...(E) > 0, "Can't select 0 terms from a multivector");
    return multivector<void, decltype(extract<E>(mv))...>{};
}

// Remove all terms in a multivector that are elements contained in the template parameter sequence
template <size_t E, size_t... Es, typename... T>
[[nodiscard]] constexpr auto component_filter(multivector<void, T...> mv) noexcept
{
    if constexpr (sizeof...(Es) == 0)
    {
        return filter<E>(mv);
    }
    else
    {
        return component_filter<Es...>(filter<E>(mv));
    }
}

#define GAL_OPERATORS(Algebra) \
    template <typename... I, typename... J> \
    [[nodiscard]] constexpr auto operator|(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept \
    {\
        return ::gal::detail::product<Algebra::inner>(lhs, rhs);\
    }\
    template <typename... I, typename... J> \
    [[nodiscard]] constexpr auto operator>>(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept \
    {\
        return ::gal::detail::product<Algebra::contract>(lhs, rhs);\
    }\
    template <typename... I, typename... J>\
    [[nodiscard]] constexpr auto operator^(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept\
    {\
        return ::gal::detail::product<Algebra::exterior>(lhs, rhs);\
    }\
    template <typename... I, typename... J>\
    [[nodiscard]] constexpr auto operator*(multivector<void, I...> lhs, multivector<void, J...> rhs) noexcept\
    {\
        return ::gal::detail::product<Algebra::geometric>(lhs, rhs);\
    }\
    template <typename... I>\
    [[nodiscard]] constexpr auto operator~(multivector<void, I...> lhs) noexcept\
    {\
        return ::gal::reverse(lhs);\
    }\
    template <typename V, typename T>\
    [[nodiscard]] constexpr auto conjugate(V action, T subject) noexcept\
    {\
        return action * subject * ::gal::reverse(action);\
    }\
    template <typename... I>\
    [[nodiscard]] constexpr auto operator!(multivector<void, I...> input) noexcept\
    {\
        return dual<Algebra::metric_t>(input);\
    }\
    template <typename M1, typename M2>\
    [[nodiscard]] constexpr auto operator&(M1 lhs, M2 rhs) noexcept\
    {\
        return !(!lhs ^ !rhs);\
    }\
    template <typename T = float> using scalar = ::gal::scalar<T>;\
    using ::gal::simplify;\
    using pseudoscalar = ::gal::pseudoscalar<Algebra::metric_t>

#define GAL_ACCESSORS \
        [[nodiscard]] constexpr const T& operator[](size_t index) const noexcept\
        {\
            return *(reinterpret_cast<const T*>(this) + index);\
        }\
        [[nodiscard]] constexpr T& operator[](size_t index) noexcept\
        {\
            return *(reinterpret_cast<T*>(this) + index);\
        }\
        template <size_t I>\
        [[nodiscard]] constexpr auto get() const noexcept\
        {\
            if constexpr (I < size)\
            {\
                return *(reinterpret_cast<const T*>(this) + I);\
            }\
            else\
            {\
                return get_special(std::integral_constant<size_t, I>{});\
            }\
        }

template <typename T = float>
struct scalar
{
    constexpr static size_t size = 1;

    template <size_t ID>
    using type = multivector<void, term<element<0>, monomial<one, generator<tag<ID, 0>>>>>;

    T data;

    [[nodiscard]] constexpr operator T() const noexcept
    {
        return data;
    }

    GAL_ACCESSORS

    template <typename Engine, typename... I>
    [[nodiscard]] constexpr static scalar<T> convert(const Engine& engine, multivector<void, I...> mv) noexcept
    {
        auto s_e = extract<0>(mv);
        return {engine.template evaluate<T>(s_e)};
    }
};

namespace detail
{
    template <size_t ID, size_t... I, size_t... E>
    [[nodiscard]] constexpr static auto
    compute_type(std::integral_constant<size_t, ID>, std::index_sequence<I...>, std::index_sequence<E...>) noexcept
    {
        return multivector<void, term<element<E>, monomial<one, generator<tag<ID, I>>>>...>{};
    }
}

template <typename T, size_t... E>
struct entity
{
    constexpr static size_t size = sizeof...(E);

    template <size_t ID>
    using type = decltype(detail::compute_type(std::integral_constant<size_t, ID>{},
                                               std::make_index_sequence<sizeof...(E)>{},
                                               std::index_sequence<E...>{}));

    std::array<T, size> data;

    GAL_ACCESSORS

    template <typename Engine, typename... I>
    [[nodiscard]] constexpr static entity convert(const Engine& engine, multivector<void, I...> mv) noexcept
    {
        return {engine.template evaluate<T>(extract<E>(mv))...};
    }
};

namespace detail
{
    template <typename T, typename... E>
    [[nodiscard]] constexpr static auto compute_entity(multivector<void, E...>) noexcept
    {
        return entity<T, E::element_t::value...>{};
    }
}

// Convenience type for situations where all values need to be extracted
// USE SPARINGLY (or for debugging only)
template <typename Metric, typename T = float>
struct dense_vector
{
    constexpr static size_t size = (1 << Metric::dimension);

    std::array<T, size> data;

    GAL_ACCESSORS

    template <typename Engine, typename... I>
    [[nodiscard]] constexpr static dense_vector convert(const Engine& engine, multivector<void, I...> mv) noexcept
    {
        return convert(engine, mv, std::make_index_sequence<size>{});
    }

private:
    template <typename Engine, typename M, size_t... I>
    [[nodiscard]] constexpr static dense_vector convert(const Engine& engine, M mv, std::index_sequence<I...>) noexcept
    {
        return {engine.template evaluate<T>(extract<I>(mv))...};
    }
};
} // namespace gal
