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

    template <auto const& ie, int num, int den, width_t Offset, typename Seq>
    struct ie_tree_mon
    {};

    template <auto const& ie, int num, int den, width_t Offset, size_t... I>
    struct ie_tree_mon<ie, num, den, Offset, std::index_sequence<I...>>
    {
        constexpr static rat q{num, den};
        constexpr static auto inds
            = std::make_tuple(std::make_pair(std::integral_constant<width_t, ie.inds[Offset + I].id>{},
                                             std::integral_constant<int, ie.inds[Offset + I].degree>{})...);
    };

    template <auto const& ie, size_t Offset, typename Seq>
    struct ie_tree_term
    {};

    template <auto const& ie, size_t Offset, size_t... I>
    struct ie_tree_term<ie, Offset, std::index_sequence<I...>>
    {
        constexpr static std::tuple<ie_tree_mon<ie,
                                                ie.mons[Offset + I].q.num,
                                                ie.mons[Offset + I].q.den,
                                                ie.mons[Offset + I].ind_offset,
                                                std::make_index_sequence<ie.mons[Offset + I].count>>...>
            children{};
    };

    // Data := flattened array of inputs
    template <auto const& ie, typename F, typename A, typename Data, size_t... I>
    [[nodiscard]] constexpr static auto compute_entity_typed(Data const& data, std::index_sequence<I...>) noexcept
    {
        using entity_t = entity<A, F, ie.terms[I].element...>;
        // NOTE: this should be flattened to arrays of constant size ideally
        constexpr std::tuple<ie_tree_term<ie, ie.terms[I].mon_offset, std::make_index_sequence<ie.terms[I].count>>...> tree{};
        return std::apply(
            [&data](auto&&... t) {
                return entity_t{std::apply(
                    [&](auto&&... m) {
                        if constexpr (sizeof...(m) == 0)
                        {
                            return F{0};
                        }
                        else
                        {
                            return ((static_cast<F>(m.q)
                                     * std::apply(
                                         [&](auto&&... i) {
                                             if constexpr (sizeof...(i) == 0)
                                             {
                                                 return F{1};
                                             }
                                             else
                                             {
                                                 return (::gal::pow(*data[i.first], i.second) * ...);
                                             }
                                         },
                                         m.inds))
                                    + ...);
                        }
                    },
                    t.children)...};
            },
            tree);
    }

    template <auto const& dmv, typename F, typename A, typename Data, size_t... I>
    [[nodiscard]] constexpr static auto compute_entity(Data const& data, std::index_sequence<I...>) noexcept
    {
        using entity_t = entity<A, F, dmv.data[I].first.element...>;
        // NOTE: this should be flattened to arrays of constant size ideally
        return std::apply(
            [&data](auto&&... t) {
                return entity_t{std::apply(
                    [&](auto&&... m) {
                        if constexpr (sizeof...(m) == 0)
                        {
                            return F{0};
                        }
                        else
                        {
                            return ((m.first.q.num == 0
                                         ? 0
                                         : (static_cast<F>(m.first.q)
                                            * std::apply(
                                                [&](auto&&... i) {
                                                    if constexpr (sizeof...(i) == 0)
                                                    {
                                                        return F{1};
                                                    }
                                                    else
                                                    {
                                                        return ((i.id == ~0u ? 1 : ::gal::pow(*data[i.id], i.degree))
                                                                * ...);
                                                    }
                                                },
                                                m.second)))
                                    + ...);
                        }
                    },
                    t.second)...};
            },
            dmv.data);
    }

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

    template <typename A, typename V, typename T, typename D>
    [[nodiscard]] static auto finalize_entity_typed(D const& data)
    {
        constexpr static auto reified = reify<T>();
        if constexpr (detail::uses_null_basis<A>)
        {
            constexpr static auto null_conversion = detail::to_null_basis(reified);
            return compute_entity_typed<null_conversion, V, A>(data, std::make_index_sequence<null_conversion.size.term>());
        }
        else
        {
            return compute_entity_typed<reified, V, A>(data, std::make_index_sequence<reified.size.term>());
        }
    }

    template <typename A, typename V, typename T, typename D>
    [[nodiscard]] static auto finalize_entity(D const& data)
    {
        constexpr auto reified = reify<T>();
        if constexpr (detail::uses_null_basis<A>)
        {
            constexpr auto null_conversion = detail::to_null_basis(reified);
            constexpr auto extent          = null_conversion.extent();
            constexpr static auto regularized
                = null_conversion.template regularize<extent.ind, extent.mon, extent.term>();
            return compute_entity<regularized, V, A>(data, std::make_index_sequence<null_conversion.size.term>());
        }
        else
        {
            constexpr auto extent             = reified.extent();
            constexpr static auto regularized = reified.template regularize<extent.ind, extent.mon, extent.term>();
            return compute_entity<regularized, V, A>(data, std::make_index_sequence<reified.size.term>());
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
        if constexpr (std::tuple_size_v<ie_result_t>> 0)
        {
            using value_t   = typename std::tuple_element_t<0, ie_result_t>::value_t;
            using algebra_t = typename std::tuple_element_t<0, ie_result_t>::algebra_t;

            std::array<detail::ind_value<value_t>, (Data::ind_count() + ...)> data{};
            detail::fill(data.data(), input...);

            return std::apply(
                [&](auto&&... args) {
                    return std::make_tuple(
                        detail::finalize_entity_typed<algebra_t, value_t, std::decay_t<decltype(args)>>(data)...);
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
        return detail::finalize_entity_typed<algebra_t, value_t, ie_result_t>(data);
    }
}
} // namespace gal
