#pragma once

#include "entity.hpp"

#ifdef GAL_DEBUG
#include "expression_debug.hpp"
#endif

#include <cmath>
#include <tuple>

namespace gal
{
namespace detail
{
    template <typename D, typename... Ds, typename... Out, uint ID>
    [[nodiscard]] constexpr static auto ies(std::tuple<Out...> out, std::integral_constant<uint, ID> id) noexcept
    {
        if constexpr (sizeof...(Ds) == 0)
        {
            return std::tuple_cat(out, std::make_tuple(expr<expr_op::identity, D, decltype(id)>{}));
        }
        else
        {
            return ies<Ds...>(std::tuple_cat(out, std::make_tuple(expr<expr_op::identity, D, decltype(id)>{})),
                              std::integral_constant<uint, ID + D::ind_count()>{});
        }
    }

    template <typename... Out>
    [[nodiscard]] constexpr static auto ies(std::tuple<Out...> out) noexcept
    {
        return out;
    }

    template <typename T>
    struct is_tuple
    {
        constexpr static bool value = false;
    };

    template <typename... T>
    struct is_tuple<std::tuple<T...>>
    {
        constexpr static bool value = true;
    };

    template <typename T>
    inline constexpr bool is_tuple_v = is_tuple<T>::value;

    // The indeterminate value is either a pointer to an entity's value or an evaluated expression
    template <typename T>
    struct ind_value
    {
        union
        {
            T const* pointer;
            T value;
        };

        bool is_value;

        [[nodiscard]] constexpr T operator*() const noexcept
        {
            return is_value ? value : *pointer;
        }
    };

    template <typename T, typename D, typename... Ds>
    constexpr static void fill(T* out, D const& datum, Ds const&... data) noexcept
    {
        for (size_t i = 0; i != D::size(); ++i)
        {
            auto& iv    = *(out + i);
            iv.pointer  = &datum[i];
            iv.is_value = false;
        }

        for (size_t i = D::size(); i != D::ind_count(); ++i)
        {
            auto& iv    = *(out + i);
            iv.value    = datum.get(i);
            iv.is_value = true;
        }

        if constexpr (sizeof...(Ds) > 0)
        {
            fill(out + D::ind_count(), data...);
        }
    }

    template <typename, auto const&, width_t, typename>
    struct cmon
    {};

    template <typename F, auto const& ie, width_t Index, size_t... I>
    struct cmon<F, ie, Index, std::index_sequence<I...>>
    {
        template <size_t N>
        constexpr static F value(std::array<ind_value<F>, N> const& data) noexcept
        {
            constexpr auto m = ie.mons[Index];
            if constexpr (m.q.is_zero())
            {
                return {0};
            }
            else if constexpr (sizeof...(I) == 0)
            {
                return static_cast<F>(m.q);
            }
            else
            {
                return static_cast<F>(m.q)
                       * (::gal::pow(*data[ie.inds[m.ind_offset + I].id],
                                     ie.inds[m.ind_offset + I].degree.num,
                                     ie.inds[m.ind_offset + I].degree.den)
                          * ...);
            }
        }
    };

    template <typename, auto const&, size_t, typename>
    struct cterm
    {};

    template <typename F, auto const& ie, size_t Offset, size_t... I>
    struct cterm<F, ie, Offset, std::index_sequence<I...>>
    {
        template <size_t N>
        constexpr static F value(std::array<ind_value<F>, N> const& data) noexcept
        {
            if constexpr (sizeof...(I) == 0)
            {
                return {0};
            }
            else
            {
                return (cmon<F, ie, Offset + I, std::make_index_sequence<ie.mons[Offset + I].count>>::value(data)
                        + ...);
            }
        }
    };

    template <auto const& ie, typename F, typename A, size_t N, size_t... I>
    [[nodiscard]] constexpr static auto
    compute_entity(std::array<ind_value<F>, N> const& data, std::index_sequence<I...>) noexcept
    {
        using entity_t = entity<A, F, ie.terms[I].element...>;
        return entity_t{cterm<F, ie, ie.terms[I].mon_offset, std::make_index_sequence<ie.terms[I].count>>::value(data)...};
    }

    template <typename A, typename V, typename T, typename D>
    [[nodiscard]] static auto finalize_entity(D const& data)
    {
        constexpr static auto reified = reify<T>();
        if constexpr (detail::uses_null_basis<A>)
        {
            constexpr static auto null_conversion = detail::to_null_basis(reified);
            return compute_entity<null_conversion, V, A>(data, std::make_index_sequence<null_conversion.size.term>());
        }
        else
        {
            return compute_entity<reified, V, A>(data, std::make_index_sequence<reified.size.term>());
        }
    }
} // namespace detail

template <typename... Data>
struct evaluate
{
    template <typename L>
    [[nodiscard]] constexpr auto operator()(L&& lambda) noexcept
    {
        constexpr auto ies = detail::ies<Data...>(std::tuple<>{}, std::integral_constant<uint, 0>{});
        using ie_result_t  = decltype(std::apply(lambda, ies));
        return reify<ie_result_t>();
    }

#ifdef GAL_DEBUG
    // Non-constexpr variant of the main evaluation operator for runtime debugging
    template <typename L>
    [[nodiscard]] auto debug(L&& lambda) noexcept
    {
        static auto ies        = detail::ies<Data...>(std::tuple<>{}, std::integral_constant<uint, 0>{});
        static auto expression = std::apply(lambda, ies);
        using ie_result_t      = decltype(std::apply(lambda, ies));
        return debug_reify<ie_result_t>();
    }
#endif
};

template <typename L, typename... Data>
[[nodiscard]] static auto compute(L&& lambda, Data const&... input) noexcept
{
    constexpr auto ies = detail::ies<Data...>(std::tuple<>{}, std::integral_constant<uint, 0>{});
    using ie_result_t  = decltype(std::apply(lambda, ies));
    // Produce a lookup table keyed to the indeterminate id mapping to a union containing either an entity property or
    // an evaluated property
    if constexpr (detail::is_tuple_v<ie_result_t>)
    {
        if constexpr (std::tuple_size_v<ie_result_t> != 0)
        {
            using value_t   = typename std::tuple_element_t<0, ie_result_t>::value_t;
            using algebra_t = typename std::tuple_element_t<0, ie_result_t>::algebra_t;

            std::array<detail::ind_value<value_t>, (Data::ind_count() + ...)> data{};
            detail::fill(data.data(), input...);

            return std::apply(
                [&](auto&&... args) {
                    return std::make_tuple(
                        detail::finalize_entity<algebra_t, value_t, std::decay_t<decltype(args)>>(data)...);
                },
                ie_result_t{});
        }
    }
    else
    {
        using value_t   = typename ie_result_t::value_t;
        using algebra_t = typename ie_result_t::algebra_t;

        std::array<detail::ind_value<value_t>, (Data::ind_count() + ...)> data{};
        detail::fill(data.data(), input...);
        return detail::finalize_entity<algebra_t, value_t, ie_result_t>(data);
    }
}
} // namespace gal
