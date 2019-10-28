#pragma once

#include "entity.hpp"

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
            return std::tuple_cat(out, std::make_tuple(expr_id<D, ID>{}));
        }
        else
        {
            return ies<Ds...>(std::tuple_cat(out, std::make_tuple(expr_id<D, ID>{})),
                              std::integral_constant<uint, ID + D::size()>{});
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
        T const* pointer;

        [[nodiscard]] constexpr T operator*() const noexcept
        {
            return *pointer;
        }
    };

    template <typename T, typename D, typename... Ds>
    GAL_FORCE_INLINE constexpr static void fill(T* out, D const& datum, Ds const&... data) noexcept
    {
        if constexpr (D::size() > 0)
        {
            for (size_t i = 0; i != D::size(); ++i)
            {
                auto& iv   = *(out + i);
                iv.pointer = &datum[i];
            }
        }

        if constexpr (sizeof...(Ds) > 0)
        {
            fill(out + D::size(), data...);
        }
    }

    template <typename, auto const&, width_t, typename>
    struct cmon
    {};

    template <typename F, auto const& ie, width_t Index, size_t... I>
    struct cmon<F, ie, Index, std::index_sequence<I...>>
    {
        template <size_t N>
        GAL_FORCE_INLINE constexpr static F value(std::array<ind_value<F>, N> const& data) noexcept
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
                                     std::integral_constant<int, ie.inds[m.ind_offset + I].degree.num>{},
                                     std::integral_constant<int, ie.inds[m.ind_offset + I].degree.den>{})
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
        GAL_FORCE_INLINE constexpr static F value(std::array<ind_value<F>, N> const& data) noexcept
        {
            return (cmon<F, ie, Offset + I, std::make_index_sequence<ie.mons[Offset + I].count>>::value(data) + ...);
        }
    };

    template <auto const& ie, typename F, typename A, size_t N, size_t... I>
    [[nodiscard]] GAL_FORCE_INLINE constexpr static inline auto
    compute_entity(std::array<ind_value<F>, N> const& data, std::index_sequence<I...>) noexcept
    {
        using entity_t = entity<A, F, ie.terms[I].element...>;
        if constexpr (sizeof...(I) == 0)
        {
            return entity_t{};
        }
        else
        {
            return entity_t{
                cterm<F, ie, ie.terms[I].mon_offset, std::make_index_sequence<ie.terms[I].count>>::value(data)...};
        }
    }

    template <typename A, typename V, typename T, typename D>
    [[nodiscard]] GAL_FORCE_INLINE static inline auto finalize_entity(D const& data)
    {
        constexpr static auto reified = T::reify();
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

    template <typename T, typename... Ts>
    struct reflect_first
    {
        using value_t = typename T::value_t;
        using algebra_t = typename T::algebra_t;
    };
} // namespace detail

template <typename... Data>
struct evaluate
{
    template <typename L>
    [[nodiscard]] constexpr auto operator()(L&& lambda) noexcept
    {
        constexpr auto ies = detail::ies<Data...>(std::tuple<>{}, std::integral_constant<uint, 0>{});
        using ie_result_t  = decltype(std::apply(lambda, ies));
        return ie_result_t::reify();
    }

#ifdef GAL_DEBUG
    // Non-constexpr variant of the main evaluation operator for runtime debugging
    template <typename L>
    [[nodiscard]] auto debug(L&& lambda) noexcept
    {
        static auto ies        = detail::ies<Data...>(std::tuple<>{}, std::integral_constant<uint, 0>{});
        static auto expression = std::apply(lambda, ies);
        using ie_result_t      = decltype(std::apply(lambda, ies));
        return ie_result_t::reify_debug();
    }
#endif
};

namespace detail
{
    template <typename D, typename... Ds>
    struct category
    {
        using value_t   = typename D::value_t;
        using algebra_t = typename D::algebra_t;
    };
}

template <typename L, typename... Data>
[[nodiscard]] static auto compute(L&& lambda, Data const&... input) noexcept
{
    constexpr auto ies = detail::ies<Data...>(std::tuple<>{}, std::integral_constant<uint, 0>{});
    using cat          = detail::category<Data...>;
    using value_t      = typename cat::value_t;
    using algebra_t    = typename cat::algebra_t;
    using ie_result_t  = decltype(std::apply(lambda, ies));

    // Produce a lookup table keyed to the indeterminate id mapping to a union containing either an entity property or
    // an evaluated property
    if constexpr (detail::is_tuple_v<ie_result_t>)
    {
        if constexpr (std::tuple_size_v<ie_result_t> != 0)
        {
            std::array<detail::ind_value<value_t>, (Data::size() + ...)> data{};
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
        std::array<detail::ind_value<value_t>, (Data::size() + ...)> data{};
        detail::fill(data.data(), input...);
        return detail::finalize_entity<algebra_t, value_t, ie_result_t>(data);
    }
}
} // namespace gal
