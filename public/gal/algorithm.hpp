#pragma once

#include <cstddef>
#include <utility>

namespace gal
{
namespace detail
{
    // Retrieve the nth type from a type list
    template <size_t I, typename T>
    struct indexed
    {
        using type = T;
    };

    template <typename Is, typename... Ts>
    struct indexer;

    template <size_t... Is, typename... Ts>
    struct indexer<std::index_sequence<Is...>, Ts...> : indexed<Is, Ts>...
    {};

    template <size_t I, typename T>
    constexpr static indexed<I, T> select(indexed<I, T>);

    // Uses ADL to select the correct type
    template <size_t I, typename... Ts>
    using nth_element =
        typename decltype(select<I>(indexer<std::index_sequence_for<Ts...>, Ts...>{}))::type;

    template <typename T>
    constexpr void swap(T& lhs, T& rhs) noexcept
    {
        if (&lhs == &rhs)
        {
            return;
        }

        T tmp = lhs;
        lhs   = rhs;
        rhs   = tmp;
    }

    // Needed for the time being because std::sort is not yet declared constexpr
    template <typename T>
    constexpr void sort(T first, T last) noexcept
    {
        if (last - first < 2)
        {
            // Don't bother sorting an empty or single element range
            return;
        }
        else if (last - first < 3)
        {
            if (*(first + 1) < *first)
            {
                swap(*first, *(first + 1));
            }
        }
        else if (last - first < 4)
        {
            // For a short 3 element run, do a few quick comparisons and swaps to terminate sooner
            if (*(first + 1) < *first)
            {
                if (*(first + 2) < *first)
                {
                    if (*(first + 2) < *(first + 1))
                    {
                        // 3 2 1
                        swap(*first, *(first + 2));
                    }
                    else
                    {
                        // 3 1 2
                        swap(*first, *(first + 1));
                        swap(*(first + 1), *(first + 2));
                    }
                }
                else
                {
                    // 2 1 3
                    swap(*first, *(first + 1));
                }
            }
            else
            {
                if (*(first + 2) < *first)
                {
                    // 2 3 1
                    swap(*first, *(first + 1));
                    swap(*first, *(first + 2));
                }
                else
                {
                    if (*(first + 2) < *(first + 1))
                    {
                        // 1 3 2
                        swap(*(first + 1), *(first + 2));
                    }
                }
            }
        }
        else
        {
            // NOTE: This is NOT the most efficient qsort implementation but it is implemented this
            // way for simplicity. If it is a bottleneck, a faster sort will be implemented in the
            // future. Use the last as the pivot
            auto pivot  = *(last - 1);
            auto cursor = first;
            for (auto it = first; it != last - 1; ++it)
            {
                if (*it < pivot)
                {
                    swap(*cursor++, *it);
                }
            }
            swap(*(last - 1), *cursor);
            sort(first, cursor);
            sort(cursor + 1, last);
        }
    }

    // The third argument is a comparator
    template <typename T, typename L>
    constexpr void sort(T first, T last, L&& less) noexcept
    {
        if (last - first < 2)
        {
            // Don't bother sorting an empty or single element range
            return;
        }
        else if (last - first < 3)
        {
            if (less(*(first + 1), *first))
            {
                swap(*first, *(first + 1));
            }
        }
        else if (last - first < 4)
        {
            // For a short 3 element run, do a few quick comparisons and swaps to terminate sooner
            if (less(*(first + 1), *first))
            {
                if (less(*(first + 2), *first))
                {
                    if (less(*(first + 2), *(first + 1)))
                    {
                        // 3 2 1
                        swap(*first, *(first + 2));
                    }
                    else
                    {
                        // 3 1 2
                        swap(*first, *(first + 1));
                        swap(*(first + 1), *(first + 2));
                    }
                }
                else
                {
                    // 2 1 3
                    swap(*first, *(first + 1));
                }
            }
            else
            {
                if (less(*(first + 2), *first))
                {
                    // 2 3 1
                    swap(*first, *(first + 1));
                    swap(*first, *(first + 2));
                }
                else
                {
                    if (less(*(first + 2), *(first + 1)))
                    {
                        // 1 3 2
                        swap(*(first + 1), *(first + 2));
                    }
                }
            }
        }
        else
        {
            // NOTE: This is NOT the most efficient qsort implementation but it is implemented this
            // way for simplicity. If it is a bottleneck, a faster sort will be implemented in the
            // future. Use the last as the pivot
            auto pivot  = *(last - 1);
            auto cursor = first;
            for (auto it = first; it != last - 1; ++it)
            {
                if (less(*it, pivot))
                {
                    swap(*cursor++, *it);
                }
            }
            swap(*(last - 1), *cursor);
            sort(first, cursor, less);
            sort(cursor + 1, last, less);
        }
    }
} // namespace detail
} // namespace gal
