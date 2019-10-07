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
{
};

template <typename... I>
class engine
{
public:
    constexpr engine(const I&... inputs) noexcept
        : data{inputs...}
    {
        static_assert(sizeof...(I) > 0, "An engine without any inputs cannot do any useful computation");
    }

    template <typename Result, typename T>
    [[nodiscard]] constexpr auto compute(T&& lambda) const noexcept
    {
        const auto inputs = types(std::make_index_sequence<sizeof...(I)>{});
        const auto result = std::apply(std::forward<T>(lambda), inputs);
        return Result::convert(*this, result);
    }

    template <typename T, typename... Ts>
    [[nodiscard]] constexpr auto evaluate_terms(Ts... terms) const noexcept
    {
        return std::tuple(evaluate<T>(terms)...);
    }

    template <typename T, typename E, typename... M>
    [[nodiscard]] constexpr T evaluate(term<E, M...>) const noexcept
    {
        if constexpr (sizeof...(M) == 0)
        {
            return {};
        }
        else
        {
            return (evaluate<T>(M{}) + ...);
        }
    }

private:
    template <typename T, typename Q, typename... G>
    [[nodiscard]] constexpr T evaluate(monomial<Q, G...>) const noexcept
    {
        constexpr auto q = Q::template convert<T>();
        if constexpr (sizeof...(G) == 0)
        {
            return q;
        }
        else
        {
            return q * ((evaluate<T>(G{})) * ...);
        }
    }

    template <typename T, typename Tag, typename Degree, size_t Order>
    [[nodiscard]] constexpr T evaluate(generator<Tag, Degree, Order>) const noexcept
    {
        static_assert(!Tag::untagged, "A generator in the field field of an expression is unlabeled");
        // TODO: leverage sub-expression memoization table to accelerate compile times and unoptimized codegen
        if constexpr (Degree::value == 0)
        {
            return T{1};
        }
        else
        {
            const auto value = get<T, Tag::id, Tag::index>();
            // constexpr bool is_derived
                // = std::is_same<typename std::decay<decltype(value)>::type, derived_generator<T>>::value;
            // TODO: memoize derived results
            if constexpr (Degree::value == 1)
            {
                return value;
            }
            else
            {
                // TODO: permit substitution of a different power function for more exotic types
                return std::pow(value, Degree::value);
            }
        }
    }

    template <typename T, size_t ID, size_t Index>
    [[nodiscard]] constexpr auto get() const noexcept
    {
        return std::get<ID>(data).template get<Index>();
    }

    template <size_t... N>
    [[nodiscard]] constexpr auto types(std::index_sequence<N...>) const noexcept
    {
        return std::tuple<typename I::template type<N>...>();
    }

    std::tuple<const I&...> data;
};
} // namespace gal