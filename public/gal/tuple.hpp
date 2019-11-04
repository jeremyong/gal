#pragma once

#include "algorithm.hpp"

#include <cstddef>
#include <utility>

namespace gal
{
// Our own tuple and pair types (supports additional operations and compiles faster)
// The author apologizes for needing to introduce these types but they are lightweight and will not
// add any executable bloat. The *reason* this tuple compiles so much faster is because for our
// use-case, we only need to handle value-semantics (the elements of the tuple are presumed not to
// be references).

template <typename T1, typename T2>
struct pair
{
    T1 first{};
    T2 second{};
};

template <typename T1, typename T2>
pair(T1&&, T2 &&)->pair<std::decay_t<T1>, std::decay_t<T2>>;

template <typename T1, typename T2>
constexpr pair<T1, T2> make_pair(T1 v1, T2 v2)
{
    return {v1, v2};
}

template <size_t I, typename T>
struct tuple_node
{
    T value;
};

template <typename S, typename... T>
struct tuple_impl
{};

template <size_t... I, typename... T>
struct tuple_impl<std::index_sequence<I...>, T...> : public tuple_node<I, T>...
{};

template <typename... T>
struct tuple : public tuple_impl<std::make_index_sequence<sizeof...(T)>, T...>
{
    static constexpr size_t size()
    {
        return sizeof...(T);
    }

    static constexpr std::make_index_sequence<sizeof...(T)> indices{};

    // Returns an index to the first occurence of an element of type E in the tuple.
    template <typename E>
    static constexpr auto find_t()
    {
        return find_t_impl<E, 0>();
    }

    template <size_t I>
    constexpr auto get() const
    {
        return static_cast<tuple_node<I, detail::nth_element<I, T...>> const*>(this)->value;
    }

    template <size_t I>
    constexpr auto& get()
    {
        return static_cast<tuple_node<I, detail::nth_element<I, T...>>*>(this)->value;
    }

    constexpr auto back() const
    {
        return get<sizeof...(T) - 1>();
    }

    constexpr auto& back()
    {
        return get<sizeof...(T) - 1>();
    }

    template <size_t I, typename V>
    constexpr void set(V&& v)
    {
        static_cast<tuple_node<I, detail::nth_element<I, T...>>*>(this)->value = std::forward<V>(v);
    }

    // The push, append, and pop operations do not mutate because the type signature is necessarily
    // different
    template <typename V>
    constexpr auto push(V&& v) const
    {
        using type = std::decay_t<V>;
        tuple<type, T...> out{};
        out.template set<0>(std::forward<V>(v));
        push_impl(out, std::make_index_sequence<sizeof...(T)>{});
        return out;
    }

    template <typename V>
    constexpr auto append(V&& v) const
    {
        using type = std::decay_t<V>;
        tuple<T..., type> out{};
        out.template set<sizeof...(T)>(std::forward<V>(v));
        append_impl(out, std::make_index_sequence<sizeof...(T)>{});
        return out;
    }

    // Return a pair containing the popped element and a tuple of the remaining elements
    constexpr auto pop() const
    {
        return pop_impl(std::make_index_sequence<sizeof...(T) - 1>{});
    }

    // Returns a pair of the tuple partitioned in two at the specified partition point (the first
    // tuple returned will have I elements)
    template <size_t I>
    constexpr auto split() const
    {
        return split_n(std::make_index_sequence<I>{}, std::make_index_sequence<sizeof...(T) - I>{});
    }

    // Similar to split but the element at the index supplied is omitted from the returned tuples.
    template <size_t I>
    constexpr auto split_drop() const
    {
        return split_n(
            std::make_index_sequence<I>{}, std::make_index_sequence<sizeof...(T) - I - 1>{});
    }

    // This is a somewhat specialized operation that splits the tuple on the first occurence of the
    // supplied type. The element that first matched the type is consumed in the operation.
    template <typename E>
    constexpr auto split_at_type() const
    {
        constexpr auto index = find_t<E>();
        return split_drop<index>();
    }

    template <typename L>
    constexpr auto apply(L l) const
    {
        return apply_impl(std::make_index_sequence<sizeof...(T)>{}, l);
    }

    template <typename L>
    constexpr auto apply_reverse(L l) const
    {
        return apply_reverse_impl(std::make_index_sequence<sizeof...(T)>{}, l);
    }

    template <template <typename> typename F, typename... A>
    constexpr void mutate(A&&... args)
    {
        mutate_impl<F>(std::make_index_sequence<sizeof...(T)>{}, std::forward<A>(args)...);
    }

private:
    template <typename E, size_t I>
    static constexpr auto find_t_impl()
    {
        if constexpr (std::is_same_v<detail::nth_element<I, T...>, E>)
        {
            return I;
        }
        else
        {
            static_assert(
                I + 1 < sizeof...(T), "Attempted to locate a type in a tuple that does not exist.");
            return find_t_impl<E, I + 1>();
        }
    }

    template <typename S, size_t... I>
    constexpr auto push_impl(S& out, std::index_sequence<I...>) const
    {
        (out.template set<I + 1>(get<I>()), ...);
    }

    template <typename S, size_t... I>
    constexpr auto append_impl(S& out, std::index_sequence<I...>) const
    {
        (out.template set<I>(get<I>()), ...);
    }

    template <size_t... I>
    constexpr auto pop_impl(std::index_sequence<I...>) const
    {
        if constexpr (sizeof...(I) == 0)
        {
            return make_pair(get<0>(), tuple<>{});
        }
        else
        {
            tuple<decltype(get<1 + I>())...> second{get<1 + I>()...};
            return make_pair(get<0>(), second);
        }
    }

    template <size_t... I1, size_t... I2>
    constexpr auto split_n(std::index_sequence<I1...>, std::index_sequence<I2...>) const
    {
        if constexpr (sizeof...(I1) == 0)
        {
            if constexpr (sizeof...(I1) == 0)
            {
                return make_pair(tuple<>{}, tuple<>{});
            }
            else
            {
                tuple<decltype(get<I2>())...> second{get<I2>()...};
                return make_pair(tuple<>{}, second);
            }
        }
        else if constexpr (sizeof...(I2) == 0)
        {
            tuple<decltype(get<I1>())...> first{get<I1>()...};
            return make_pair(first, tuple<>{});
        }
        else
        {
            tuple<decltype(get<I1>())...> first{get<I1>()...};
            tuple<decltype(get<sizeof...(I1) + I2>())...> second{get<sizeof...(I1) + I2>()...};
            return make_pair(first, second);
        }
    }

    template <typename L, size_t... I>
    constexpr auto apply_impl(std::index_sequence<I...>, L l) const
    {
        return l(get<I>()...);
    }

    template <typename L, size_t... I>
    constexpr auto apply_reverse_impl(std::index_sequence<I...>, L l) const
    {
        return l(get<sizeof...(I) - I - 1>()...);
    }

    template <template <typename> typename F, size_t... I, typename... A>
    constexpr auto mutate_impl(std::index_sequence<I...>, A&&... args)
    {
        (set<I>(F<decltype(get<I>())>{std::forward<A>(args)...}(get<I>())), ...);
    }
};

template <typename... T>
tuple(T&&...)->tuple<std::decay_t<T>...>;

template <typename... T>
constexpr auto make_tuple(T&&... args)
{
    return tuple<std::decay_t<T>...>{std::forward<T>(args)...};
}

template <typename T>
struct is_tuple
{
    constexpr static bool value = false;
};

template <typename... T>
struct is_tuple<tuple<T...>>
{
    constexpr static bool value = true;
};

template <typename T>
constexpr inline bool is_tuple_v = is_tuple<T>::value;

template <typename T>
struct is_pair
{
    constexpr static bool value = false;
};

template <typename... T>
struct is_pair<pair<T...>>
{
    constexpr static bool value = true;
};

template <typename T>
constexpr inline bool is_pair_v = is_pair<T>::value;
} // namespace gal
