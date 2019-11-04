#pragma once

#include <cstddef>

namespace gal
{
namespace detail
{
    // Dead-simple constexpr stack
    template <typename T, size_t C>
    struct stack
    {
        // Initialization is required for constexpr-contexts
        T data[C]    = {};
        size_t count = 0;

        constexpr bool empty() const noexcept
        {
            return count == 0;
        }

        constexpr void clear() noexcept
        {
            count = 0;
        }

        constexpr T operator[](size_t i) const noexcept
        {
            return data[i];
        }

        constexpr T& operator[](size_t i) noexcept
        {
            return data[i];
        }

        constexpr void push(T in) noexcept
        {
            data[count++] = in;
        }

        constexpr T pop() noexcept
        {
            return data[count--];
        }

        constexpr T pop(size_t n) noexcept
        {
            count -= n;
            return data[count];
        }

        constexpr T& peek() noexcept
        {
            return data[count - 1];
        }

        constexpr T peek() const noexcept
        {
            return data[count - 1];
        }

        constexpr T peek(size_t i) const noexcept
        {
            return data[count - i];
        }
    };
} // namespace detail
} // namespace gal
