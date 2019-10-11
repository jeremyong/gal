#pragma once

#include "ga.hpp"

#include <cmath>

namespace gal
{
// An engine is used to execute the computation expressed in a multivector expression template.
// The default engine evaluates the expression immediately on the CPU using a specified precision.
// TODO: SPIR-V engine
// TODO: SIMD engine
// TODO: constraint evaluation

// Memoization trie used to store and retrieve results of subexpressions that repeat during a computation
// The number of combinations possible of multiplications between N factors scales as the factorial of N, so it is
// unreasonable to use a dense representation for memoized products. However, a total ordering exists between the
// generators so intermediate products are amenable to storage in a trie-like structure.

// TODO: AT THE MOMENT (10/06/2019) this optimization does not yet seem necessary because the compiler correctly
// memoizes results, even without -ffast-math. However, depending on benchmarks, this may need to be revisted in the
// feature.
template <typename... T>
struct mem_trie
{};

namespace detail
{
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

} // namespace detail

enum class scalar_type
{
    single_precision,
    double_precision,
};

// F := field over which computation is done (i.e. float)
template <typename F = float>
class engine
{
public:
    template <typename L, typename... Data>
    [[nodiscard]] inline constexpr static auto compute(L&& lambda, Data const&... source_data) noexcept
    {
        std::array<base_entity const*, sizeof...(Data)> const data{&source_data...};
        std::array<scalar_type, sizeof...(Data)> const type_sizes{
            (std::is_same_v<typename std::decay_t<decltype(source_data)>::value_t, float>
                 ? scalar_type::single_precision
                 : scalar_type::double_precision)...};

        constexpr auto inputs = types<Data...>(std::make_index_sequence<sizeof...(Data)>{});

        // It is important that this lambda NOT be explicitly invoked as it is a constexpr lambda
        using result_t = std::decay_t<decltype(std::apply(std::forward<L>(lambda), inputs))>;

        if constexpr (detail::is_tuple_v<result_t>)
        {
            return reify_entities(data, type_sizes, std::make_index_sequence<std::tuple_size_v<result_t>>{}, result_t{});
        }
        else
        {
            return reify_entity(data, type_sizes, result_t{});
        }
    }

private:
    template <typename Data, typename E, typename Sizes, size_t... I>
    [[nodiscard]] inline constexpr static auto
    reify_entities(Data const& data, Sizes const& sizes, std::index_sequence<I...>, E) noexcept
    {
        return std::make_tuple(reify_entity(data, sizes, std::tuple_element_t<I, E>{})...);
    }

    template <typename Data, typename Sizes, typename E>
    [[nodiscard]] inline constexpr static auto reify_entity(Data const& data, Sizes const& sizes, E) noexcept
    {
        using entity_t = typename detail::entity_type<F, E>::type;
        entity_t out;

        using reifier          = detail::reifier<entity_t, E>;
        size_t monomial_index  = 0;
        size_t generator_index = 0;

        for (size_t i = 0; i != reifier::monomial_counts.size(); ++i)
        {
            auto const monomial_count = reifier::monomial_counts[i];
            F accum{0};

            for (size_t j = monomial_index; j != monomial_index + monomial_count; ++j)
            {
                auto const generator_count = reifier::generator_counts[j];
                auto const rational        = reifier::rationals[j];
                F product                  = static_cast<F>(rational.first) / static_cast<F>(rational.second);

                for (size_t k = generator_index; k != generator_index + generator_count; ++k)
                {
                    auto const generator = reifier::generators[k];
                    auto const id        = std::get<0>(generator);
                    auto const index     = std::get<1>(generator);
                    auto const degree    = std::get<2>(generator);
                    auto const& datum    = *data[id];
                    // TODO: come up with a more robust way of negotiating the scalar types here
                    auto const& gen = sizes[id] == scalar_type::single_precision
                                          ? *(reinterpret_cast<float const*>(&datum) + index)
                                          : *(reinterpret_cast<double const*>(&datum) + index);

                    product *= std::pow(gen, degree);
                }
                accum += product;

                generator_index += generator_count;
            }

            out[i] = accum;
            monomial_index += monomial_count;
        }

        return out;
    }

    template <typename... Data, size_t... IDs>
    [[nodiscard]] inline constexpr static auto types(std::index_sequence<IDs...>) noexcept
    {
        return std::tuple<typename Data::template type<IDs>...>();
    }
};
} // namespace gal
