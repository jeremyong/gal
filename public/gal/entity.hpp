#pragma once

#include "expression.hpp"

namespace gal
{
namespace detail
{
    template <typename A, width_t... N, uint8_t... E>
    [[nodiscard]] constexpr static auto
    construct_ie(uint32_t id, std::integer_sequence<width_t, N...>, std::integer_sequence<uint8_t, E...>) noexcept
    {
        constexpr size_t count         = sizeof...(E);

        return mv<A, count, count, count>{
            mv_size{count, count, count}, {ind{id + N, 1, 0}...}, {mon{one, 1, N, 1}...}, {term{1, N, E}...}};
    }

    template <uint8_t T, uint8_t... E>
    [[nodiscard]] constexpr uint8_t index_of()
    {
        // This is ONE-indexed to distinguish between identifying the zero-th element versus not finding the element.
        constexpr uint8_t result = 1 + ((T < E ? 1 : 0) + ...);
        static_assert(result > 0, "Attempted to select a component of an entity that does not exist!");
        return result - 1;
    }

} // namespace detail

// All entities are expected to provide a static function to retrieve the entity's indeterminate expression given an
// entity id. See the implementation for a scalar below for an example.
// This generic entity is computed by a GAL engine and all entities are expected to be convertible from an entity.
template <typename A, typename T, uint8_t... E>
struct entity
{
    using algebra_t = A;
    using value_t   = T;
    constexpr static std::array<uint8_t, sizeof...(E)> elements{E...};

    std::array<T, sizeof...(E)> data_;

    // NOTE: in GAL code, `ie` refers always to "indeterminate expression"
    [[nodiscard]] constexpr static auto ie(uint32_t id) noexcept
    {
        return detail::construct_ie<A>(
            id, std::make_integer_sequence<width_t, sizeof...(E)>{}, std::integer_sequence<uint8_t, E...>{});
    }

    [[nodiscard]] constexpr static size_t size() noexcept
    {
        return sizeof...(E);
    }

    [[nodiscard]] constexpr static uint32_t ind_count() noexcept
    {
        return sizeof...(E);
    }

    template <uint8_t... S>
    [[nodiscard]] constexpr std::array<T, sizeof...(S)> select() noexcept
    {
        return {data_[detail::index_of<S, E...>()]...};
    }

    [[nodiscard]] constexpr T const* data() const noexcept
    {
        return data_.data();
    }

    [[nodiscard]] constexpr auto begin() noexcept
    {
        return data_.begin();
    }

    [[nodiscard]] constexpr auto end() noexcept
    {
        return data_.end();
    }

    [[nodiscard]] constexpr const T& operator[](size_t index) const noexcept
    {
        return data_[index];
    }

    [[nodiscard]] constexpr T& operator[](size_t index) noexcept
    {
        return data_[index];
    }

    [[nodiscard]] constexpr T get(size_t) const noexcept
    {
        // Unreachable
        return {};
    }
};

template <typename A, typename T>
struct scalar
{
    using algebra_t = A;
    using value_t   = T;

    [[nodiscard]] constexpr static size_t size() noexcept
    {
        return 1;
    }

    [[nodiscard]] constexpr static uint32_t ind_count() noexcept
    {
        return 1;
    }

    // NOTE: in GAL code, `ie` refers always to "indeterminate expression"
    [[nodiscard]] constexpr static mv<A, 1, 1, 1> ie(uint32_t id) noexcept
    {
        return {mv_size{1, 1, 1}, {ind{id, 1}}, {mon{one, 1, 0, 1}}, {term{1, 0, 0}}};
    }

    [[nodiscard]] constexpr T const* data() const noexcept
    {
        return &value;
    }

    [[nodiscard]] constexpr operator T() const noexcept
    {
        return value;
    }

    [[nodiscard]] constexpr T const& operator[](size_t) const noexcept
    {
        return value;
    }

    [[nodiscard]] constexpr T& operator[](size_t) noexcept
    {
        return value;
    }

    [[nodiscard]] constexpr T get(size_t) const noexcept
    {
        // Unreachable
        return {};
    }

    T value;
};
} // namespace gal
