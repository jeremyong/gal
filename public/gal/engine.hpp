#pragma once

#include "dfa.hpp"
#include "entity.hpp"

#include <cmath>
#include <tuple>
#include <type_traits>

namespace gal
{
namespace detail
{
    // The indeterminate value is either a pointer to an entity's value or an evaluated expression
    template <typename T>
    struct ind_value
    {
        union
        {
            T value;
            T const* pointer;
        };
        bool is_pointer;

        [[nodiscard]] GAL_FORCE_INLINE constexpr T operator*() const noexcept
        {
            return is_pointer ? *pointer : value;
        }

        GAL_FORCE_INLINE constexpr ind_value<T>& operator=(T v) noexcept
        {
            is_pointer = false;
            value      = v;
            return *this;
        }
    };

    template <typename T, typename D, typename... Ds>
    GAL_FORCE_INLINE constexpr static void fill(T* out, D const& datum, Ds const&... data) noexcept
    {
        if constexpr (std::is_floating_point_v<D>)
        {
            auto& iv      = *out;
            iv.is_pointer = true;
            iv.pointer    = &datum;
        }
        else if constexpr (D::size() > 0)
        {
            for (size_t i = 0; i != D::size(); ++i)
            {
                auto& iv      = *(out + i);
                iv.is_pointer = true;
                iv.pointer    = &datum[i];
            }
        }

        if constexpr (sizeof...(Ds) > 0)
        {
            if constexpr (std::is_floating_point_v<D>)
            {
                fill(out + 1, data...);
            }
            else
            {
                fill(out + D::size(), data...);
            }
        }
    }

    template <typename, auto const&, width_t, typename>
    struct cmon
    {};

    template <typename F, mv_op Op>
    GAL_FORCE_INLINE constexpr F apply_mv_op(F in)
    {
        if constexpr (Op == mv_op::id)
        {
            return in;
        }
        else if constexpr (Op == mv_op::sin)
        {
            return std::sin(in);
        }
        else if constexpr (Op == mv_op::cos)
        {
            return std::cos(in);
        }
        else if constexpr (Op == mv_op::tan)
        {
            return std::tan(in);
        }
        else if constexpr (Op == mv_op::sqrt)
        {
            return std::sqrt(in);
        }
    }

    template <typename F, auto const& ie, width_t Index, size_t... I>
    struct cmon<F, ie, Index, std::index_sequence<I...>>
    {
        GAL_FORCE_INLINE constexpr static F data_value(ind_value<F> const* data, width_t id) noexcept
        {
            if (id >= ind_constant_start)
            {
                return ind_constants<F>[id - ind_constant_start];
            }
            else
            {
                return *data[id];
            }
        }

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
                return apply_mv_op<F, ie.o>(static_cast<F>(m.q));
            }
            else
            {
                return apply_mv_op<F, ie.o>(
                    static_cast<F>(m.q)
                    * (::gal::pow(data_value(data.data(), ie.inds[m.ind_offset + I].id),
                                  std::integral_constant<int, ie.inds[m.ind_offset + I].degree.num>{},
                                  std::integral_constant<int, ie.inds[m.ind_offset + I].degree.den>{})
                       * ...));
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
            return (
                cmon<F, ie, Offset + I, std::make_index_sequence<ie.mons[Offset + I].count>>::value(data)
                + ...);
        }
    };

    template <auto const& ie, typename F, typename A, size_t N, size_t... I>
    GAL_FORCE_INLINE constexpr static inline auto
    compute_entity(std::array<ind_value<F>, N> const& data, F q, std::index_sequence<I...>) noexcept
    {
        if constexpr (sizeof...(I) == 0)
        {
            return entity<A, F>{};
        }
        else
        {
            using entity_t = entity<A, F, ie.terms[I].element...>;
            return entity_t{(
                q
                * cterm<F, ie, ie.terms[I].mon_offset, std::make_index_sequence<ie.terms[I].count>>::value(
                    data))...};
        }
    }

    template <auto const& ie, auto o, typename F, typename A, size_t N, size_t... I>
    GAL_FORCE_INLINE constexpr static inline void
    compute_temp(std::array<ind_value<F>, N>& data, std::index_sequence<I...>, size_t offset) noexcept
    {
        if constexpr (sizeof...(I) == 0)
        {
            return;
        }
        else
        {
            if constexpr (o == op_exp)
            {
                using entity_t = entity<A, F, ie.terms[I].element...>;
                typename special_entities<A, F>::line_t line{entity_t{
                    cterm<F, ie, ie.terms[I].mon_offset, std::make_index_sequence<ie.terms[I].count>>::value(
                        data)...}};
                auto motor = line.exp();
                for (size_t i = 0; i != decltype(motor)::size(); ++i)
                {
                    data[offset + i] = motor[i];
                }
            }
            else if constexpr (o == op_log)
            {
                using entity_t = entity<A, F, ie.terms[I].element...>;
                typename special_entities<A, F>::motor_t motor{entity_t{
                    cterm<F, ie, ie.terms[I].mon_offset, std::make_index_sequence<ie.terms[I].count>>::value(
                        data)...}};
                auto line = static_cast<decltype(motor)>(motor.log()).bivector();
                for (size_t i = 0; i != decltype(line)::size(); ++i)
                {
                    data[offset + i] = line[i];
                }
            }
            else
            {
                ((data[offset + I]
                  = cterm<F, ie, ie.terms[I].mon_offset, std::make_index_sequence<ie.terms[I].count>>::value(
                      data)),
                 ...);
            }
        }
    }

    template <typename A, typename V, auto const& temps, typename D, size_t I>
    GAL_FORCE_INLINE static inline void finalize_temps(D& data, std::integral_constant<size_t, I>)
    {
        if constexpr (I == std::decay_t<decltype(temps)>::size())
        {
            return;
        }
        else
        {
            constexpr static auto ie = temps.template get<I>().ie;
            constexpr static auto o  = temps.template get<I>().o;
            constexpr static auto id = temps.template get<I>().id;
            if constexpr (detail::uses_null_basis<A>)
            {
                constexpr static auto null_conversion = detail::to_null_basis(ie);
                compute_temp<null_conversion, o, V, A>(
                    data, std::make_index_sequence<null_conversion.size.term>{}, id);
            }
            else
            {
                compute_temp<ie, o, V, A>(data, std::make_index_sequence<ie.size.term>(), id);
            }

            if constexpr (I + 1 != std::decay_t<decltype(temps)>::size())
            {
                finalize_temps<A, V, temps>(data, std::integral_constant<size_t, I + 1>{});
            }
        }
    }

    template <typename A, typename V, auto const& result, typename D>
    GAL_FORCE_INLINE static inline auto finalize_entity(D const& data, V q)
    {
        if constexpr (detail::uses_null_basis<A>)
        {
            constexpr static auto null_conversion = detail::to_null_basis(result);
            return compute_entity<null_conversion, V, A>(
                data, q, std::make_index_sequence<null_conversion.size.term>());
        }
        else
        {
            return compute_entity<result, V, A>(
                data, q, std::make_index_sequence<result.size.term>());
        }
    }

    template <typename A, typename V, auto const& results, size_t I>
    struct result_finalizer
    {
        constexpr static auto ie = results.template get<I>().second;

        template <typename D>
        GAL_FORCE_INLINE constexpr static auto compute(D const& data, V q) noexcept
        {
            return finalize_entity<A, V, ie>(data, q);
        }
    };

    template <typename A, typename V, auto const& results>
    struct finalize_entities
    {
        constexpr static size_t size = std::decay_t<decltype(results)>::size();
        template <typename D>
        GAL_FORCE_INLINE constexpr static auto execute(D const& data, V q) noexcept
        {
            return execute_impl(data, q, std::decay_t<decltype(results)>::indices);
        }

        template <typename D, size_t... I>
        GAL_FORCE_INLINE constexpr static auto
        execute_impl(D const& data, V q, std::index_sequence<I...>) noexcept
        {
            return std::make_tuple(result_finalizer<A, V, results, size - I - 1>::compute(data, q)...);
        }
    };

    template <typename D, typename... Ds>
    struct category
    {
        using value_t   = typename D::value_t;
        using algebra_t = typename D::algebra_t;
    };

    template <typename D>
    constexpr size_t data_size() noexcept
    {
        if constexpr (std::is_floating_point_v<D>)
        {
            return 1;
        }
        else
        {
            return D::size();
        }
    }
} // namespace detail

template <typename... Data>
struct evaluate
{
#ifdef GAL_DEBUG
    // Produce the RPN expression
    template <typename L>
    [[nodiscard]] static auto rpnf(L lambda) noexcept
    {
        using A                          = typename detail::category<Data...>::algebra_t;
        constexpr static auto entities   = detail::rpne_entities<A, Data...>();
        constexpr static auto expression = std::apply(lambda, entities.first);
        return detail::rpne_concat<expression>();
    }

    // Produce an RPN expression with CSEs factored out
    template <typename L>
    [[nodiscard]] static auto rpnf_reshaped(L lambda) noexcept
    {
        using A                          = typename detail::category<Data...>::algebra_t;
        constexpr static auto entities   = detail::rpne_entities<A, Data...>();
        constexpr static auto expression = std::apply(lambda, entities.first);
        return detail::rpn_reshape(detail::rpne_concat<expression>());
    }

    template <typename L>
    static auto input_map(L lambda) noexcept
    {
        using A                          = typename detail::category<Data...>::algebra_t;
        constexpr static auto entities   = detail::rpne_entities<A, Data...>();
        constexpr static auto expression = std::apply(lambda, entities.first);
        constexpr static auto rpn        = detail::rpne_concat<expression>();
        constexpr static auto id_count   = detail::rpn_id_count(rpn);
        return detail::rpn_ids(rpn, std::integral_constant<width_t, id_count>{});
    }

    // Produce a multivector
    template <typename L>
    [[nodiscard]] static auto ie(L lambda) noexcept
    {
        using A                          = typename detail::category<Data...>::algebra_t;
        constexpr static auto entities   = detail::rpne_entities<A, Data...>();
        constexpr static auto expression = std::apply(lambda, entities.first);
        constexpr static auto rpn        = detail::rpne_concat<expression>();
        constexpr static auto id_count   = detail::rpn_id_count(rpn);
        constexpr static auto reshaped   = detail::rpn_reshape(rpn);
        constexpr static auto flattened
            = detail::rpn_ids(rpn, std::integral_constant<width_t, id_count>{});
        constexpr static auto ids     = flattened.first;
        constexpr static auto indices = flattened.second;
        constexpr static auto inputs
            = detail::rpn_inputs<A, ids, indices, Data...>{}(std::make_index_sequence<ids.size()>{});
        // The inputs are now multivector ies.
        constexpr static detail::rpn_state input_state{
            inputs, tuple<>{}, tuple<>{}, entities.second.first};
        auto result = detail::rpn_ctx<reshaped, 0, reshaped.count, input_state>::state;
        return result.args.template get<0>().second;
    }

    template <typename L>
    [[nodiscard]] static auto ie_reshaped(L lambda) noexcept
    {
        using A                          = typename detail::category<Data...>::algebra_t;
        constexpr static auto entities   = detail::rpne_entities<A, Data...>();
        constexpr static auto expression = std::apply(lambda, entities.first);
        constexpr static auto rpn        = detail::rpne_concat<expression>();
        constexpr static auto id_count   = detail::rpn_id_count(rpn);
        constexpr static auto reshaped   = detail::rpn_reshape(rpn);
        constexpr static auto flattened
            = detail::rpn_ids(reshaped, std::integral_constant<width_t, id_count>{});
        constexpr static auto ids     = flattened.first;
        constexpr static auto indices = flattened.second;
        constexpr static auto inputs
            = detail::rpn_inputs<A, ids, indices, Data...>{}(std::make_index_sequence<ids.size()>{});
        // The inputs are now multivector ies.
        constexpr static detail::rpn_state input_state{
            inputs, tuple<>{}, tuple<>{}, entities.second.first};
        auto result = detail::rpn_ctx<reshaped, 0, reshaped.count, input_state>::state;
        return result.args.template get<0>().second;
    }
#endif
};

template <typename L, typename... Data>
static inline auto compute(L lambda, Data const&... input) noexcept
{
    // At least one input is needed to infer the metric and value type
    static_assert(sizeof...(Data) > 0, "Compute contexts without any inputs are not permitted");

    using A = typename detail::category<Data...>::algebra_t;
    using V = typename detail::category<Data...>::value_t;

    // See expr.hpp
    // Governs compiling an expression from a user-supplied lambda in RPN form.
    constexpr static auto entities   = detail::rpne_entities<A, Data...>();
    constexpr static auto expression = std::apply(lambda, entities.first);
    constexpr static auto rpn        = detail::rpne_concat<expression>();

    // Extract common subexpressions and dependent args.
    constexpr static auto reshaped  = detail::rpn_reshape(rpn);
    constexpr static V scale_factor = static_cast<V>(reshaped.q);

    constexpr static auto id_count = detail::rpn_id_count(reshaped);
    constexpr static auto flattened
        = detail::rpn_ids(reshaped, std::integral_constant<width_t, id_count>{});
    constexpr static auto ids     = flattened.first;
    constexpr static auto indices = flattened.second;
    constexpr static auto inputs
        = detail::rpn_inputs<A, ids, indices, Data...>{}(std::make_index_sequence<ids.size()>{});

    // The inputs are now multivector ies.
    constexpr static detail::rpn_state input_state{
        inputs, tuple<>{}, tuple<>{}, entities.second.first};
    constexpr static auto const& processed
        = detail::rpn_ctx<reshaped, 0, reshaped.count, input_state>::state;
    constexpr static auto temps = processed.temps;

    // All temporaries need to be evaluated in order
    std::array<detail::ind_value<V>,
               ((detail::data_size<Data>() + ...) + processed.id_count - entities.second.first)>
        data{};
    detail::fill(data.data(), input...);
    // Evaluate temporaries which will be appended to the data array as value types
    detail::finalize_temps<A, V, temps>(data, std::integral_constant<size_t, 0>{});

    if constexpr (decltype(processed.args)::size() == 0)
    {
        return entity<A, V>{};
    }
    else if constexpr (decltype(processed.args)::size() == 1)
    {
        constexpr static auto result_ie = processed.args.template get<0>().second;

        return detail::finalize_entity<A, V, result_ie>(data, scale_factor);
    }
    else
    {
        // Pack each returned entity into a tuple
        constexpr static auto args = processed.args;
        return detail::finalize_entities<A, V, args>::execute(data, scale_factor);
    }
}
} // namespace gal
