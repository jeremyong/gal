#pragma once

#include "finite_algebra.hpp"

#include <cmath>

namespace gal
{
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

    // Returns the position of a term with the corresponding element in the multivector.
    // Returns -1 if the element is not found.
    template <size_t E, typename... Ts, size_t... I>
    [[nodiscard]] constexpr int element_index(multivector<void, Ts...>, std::index_sequence<I...>) noexcept
    {
        if constexpr (sizeof...(I) == 0)
        {
            return -1;
        }
        else
        {
            return ((Ts::element_t::value == E ? I + 1 : 0) + ...) - 1;
        }
    }

    template <size_t I, size_t... Is>
    [[nodiscard]] constexpr size_t first_nonzero() noexcept
    {
        if constexpr (I == 0)
        {
            if constexpr (sizeof...(Is) == 0)
            {
                return 0;
            }
            else
            {
                return first_nonzero<Is...>();
            }
        }
        else
        {
            return I;
        }
    }

    template <size_t Index, typename... Ts, size_t... I>
    [[nodiscard]] constexpr int generator_index(multivector<void, Ts...>, std::index_sequence<I...>) noexcept
    {
        constexpr size_t index = first_nonzero((Ts::first_t::tag_t::index == Index ? I + 1 : 0)...);
        if constexpr (index == 0)
        {
            return -1;
        }
        else
        {
            return index - 1;
        }
    }

    template <size_t ID, size_t... I, size_t... E>
    [[nodiscard]] constexpr auto
        compute_type(std::integral_constant<size_t, ID>, std::index_sequence<I...>, std::index_sequence<E...>) noexcept
    {
        return multivector<void, term<element<E>, monomial<one, generator<tag<ID, I>>>>...>{};
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

// This is required so that all template instantiations of entities can be coalesced in a single array of references
struct base_entity {};

template <typename T, typename S, size_t... E>
struct entity : public base_entity
{
    using value_t = T;

    [[nodiscard]] constexpr static size_t size() noexcept
    {
        return S::size();
    }

    entity()
    {}

    [[nodiscard]] constexpr const T& operator[](size_t index) const noexcept
    {
        return *(reinterpret_cast<const T*>(this) + index);
    }

    [[nodiscard]] constexpr T& operator[](size_t index) noexcept
    {
        return *(reinterpret_cast<T*>(this) + index);
    }

    template <size_t I>
    [[nodiscard]] constexpr auto get() const noexcept
    {
        if constexpr (I < size())
        {
            return *(reinterpret_cast<const T*>(this) + I);
        }
        else
        {
            return S::template get(std::integral_constant<size_t, I>{});
        }
    }
};

template <typename T, size_t... E>
struct entity<T, void, E...> : public base_entity
{
    using value_t  = T;
    using elements = std::index_sequence<E...>;

    template <size_t ID>
    using type = decltype(detail::compute_type(std::integral_constant<size_t, ID>{},
                                               std::make_index_sequence<sizeof...(E)>{},
                                               std::index_sequence<E...>{}));

    [[nodiscard]] constexpr static size_t size() noexcept
    {
        return sizeof...(E);
    }

    entity()
        : data{}
    {}

    std::array<T, size()> data;

    [[nodiscard]] constexpr const T& operator[](size_t index) const noexcept
    {
        return *(reinterpret_cast<const T*>(this) + index);
    }

    [[nodiscard]] constexpr T& operator[](size_t index) noexcept
    {
        return *(reinterpret_cast<T*>(this) + index);
    }

    template <size_t I>
    [[nodiscard]] constexpr auto get() const noexcept
    {
        static_assert(I < size(), "Overflow detected while access compile time value in entity");
        return *(reinterpret_cast<const T*>(this) + I);
    }

    template <size_t Element>
    [[nodiscard]] constexpr auto get_by_element() const noexcept
    {
        constexpr auto index = detail::element_index<Element>(type<0>{}, std::make_index_sequence<size()>{});
        if constexpr (index == -1)
        {
            return 0;
        }
        else
        {
            return get<index>();
        }
    }
};

// Helper for scalar entities that can be implicitly casted to and from the field type
template <typename T = float>
struct scalar : public entity<T, void, 0>
{
    template <size_t ID>
    using type = multivector<void, term<element<0>, monomial<one, generator<tag<ID, 0>>>>>;

    [[nodiscard]] constexpr static size_t size() noexcept
    {
        return 1;
    }

    scalar()
        : value{}
    {}

    scalar(T value)
        : value{value}
    {}

    template <typename T2, size_t... E>
    scalar(const entity<T2, void, E...>& other)
    {
        value = other.template get_by_element<0>();
    }

    T value;

    [[nodiscard]] constexpr operator T() const noexcept
    {
        return value;
    }
};

namespace detail
{
    // Helper operator for flattening the monomials across multiple terms. The element of the term is irrelevant
    template <typename E, typename... M1, typename... M2>
    [[nodiscard]] constexpr auto operator&(std::tuple<M1...>, term<E, M2...>)noexcept
    {
        return std::tuple<M1..., M2...>{};
    }

    // Helper function for flattening generators across multiple monomials
    template <typename Q, typename... G1, typename... G2>
    [[nodiscard]] constexpr auto operator&(std::tuple<G1...>, monomial<Q, G2...>)noexcept
    {
        return std::tuple<G1..., G2...>{};
    }

    // Flattens rational coefficients
    template <typename M, typename... Q1>
    [[nodiscard]] constexpr auto operator|(std::tuple<Q1...>, M) noexcept
    {
        return std::tuple<Q1..., typename M::rational_t>{};
    }

    // T... := Parameter pack of terms extracted from a multivector
    template <typename... T>
    [[nodiscard]] constexpr auto flattened_monomials(T...)
    {
        return (std::tuple<>{} & ... & T{});
    }

    // M... := Flattened tuple of monomials
    template <typename... M>
    [[nodiscard]] constexpr auto flattened_generators(std::tuple<M...>)
    {
        return (std::tuple<>{} & ... & M{});
    }

    // M... := Flattened tuple of monomials
    template <typename... M>
    [[nodiscard]] constexpr auto flattened_rationals(std::tuple<M...>)
    {
        return (std::tuple<>{} | ... | M{});
    }

    template <typename F, typename... T>
    [[nodiscard]] constexpr auto compute_entity(multivector<void, T...>) noexcept
    {
        return entity<F, void, T::element_t::value...>{};
    }

    template <typename E, typename T>
    struct reifier
    {};

    template <typename F, size_t... E, typename... T>
    struct reifier<entity<F, void, E...>, multivector<void, T...>>
    {
        using entity_t               = entity<F, void, E...>;
        using mv_t                   = multivector<void, T...>;
        using flattened_monomials_t  = decltype(detail::flattened_monomials(extract<E>(mv_t{})...));
        using flattened_generators_t = decltype(detail::flattened_generators(flattened_monomials_t{}));
        using flattened_rationals_t  = decltype(detail::flattened_rationals(flattened_monomials_t{}));

        constexpr static size_t monomial_count  = std::tuple_size_v<flattened_monomials_t>;
        constexpr static size_t generator_count = std::tuple_size_v<flattened_generators_t>;

        constexpr static std::array<int, entity_t::size()> monomial_counts
            = {decltype(extract<E>(mv_t{}))::size...};
        constexpr static auto generator_counts = std::apply(
            [](auto... args) -> std::array<int, monomial_count> { return {decltype(args)::gen_count()...}; },
            flattened_monomials_t{});
        constexpr static auto rationals = std::apply(
            [](auto... args) -> std::array<std::pair<int, int>, monomial_count> {
                return {std::make_pair(decltype(args)::rational_t::num, decltype(args)::rational_t::den)...};
            },
            flattened_monomials_t{});
        constexpr static auto generators = std::apply(
            [](auto... args) -> std::array<std::tuple<int, int, int>, generator_count> {
                return {std::make_tuple(
                    decltype(args)::tag_t::id, decltype(args)::tag_t::index, decltype(args)::degree)...};
            },
            flattened_generators_t{});
    };

    template <typename, typename>
    struct entity_type
    {};

    template <typename F, typename... T>
    struct entity_type<F, multivector<void, T...>>
    {
        using type = entity<F, void, T::element_t::value...>;
    };
} // namespace detail
} // namespace gal