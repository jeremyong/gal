#pragma once

#include "expr.hpp"

namespace gal
{
namespace detail
{
    // This struct should be specialized for each algebra to supply describe which structs/unions
    // represent various special entities (e.g. lines, motors, etc.)
    template <typename A, typename V>
    struct special_entities
    {};

    template <typename A, width_t... N, elem_t... E>
    GAL_NODISCARD constexpr auto construct_ie(uint32_t id,
                                              std::integer_sequence<width_t, N...>,
                                              std::integer_sequence<elem_t, E...>) noexcept
    {
        constexpr size_t count = sizeof...(E);

        return mv<A, count, count, count>{mv_size{count, count, count},
                                          {ind{id + N, one}...},
                                          {mon{one, one, 1, N}...},
                                          {term{1, N, E}...}};
    }

    template <elem_t E, int16_t I, elem_t E1, elem_t... Es>
    GAL_NODISCARD constexpr int16_t indexof() noexcept
    {
        if constexpr (E == E1)
        {
            return I;
        }
        else if constexpr (sizeof...(Es) == 0)
        {
            return -1;
        }
        else
        {
            return indexof<E, I + 1, Es...>();
        }
    }

    template <elem_t... E, size_t... N>
    GAL_NODISCARD constexpr std::array<int16_t, sizeof...(N)>
        construct_element_lut(std::index_sequence<N...>) noexcept
    {
        return {indexof<N, 0, E...>()...};
    }
} // namespace detail

// All entities are expected to provide a static function to retrieve the entity's indeterminate
// expression given an entity id. See the implementation for a scalar below for an example. This
// generic entity is computed by a GAL engine and all entities are expected to be convertible from
// an entity.
template <typename A, typename T, elem_t... E>
struct entity
{
    using algebra_t = A;
    using value_t   = T;
    constexpr static std::array<elem_t, sizeof...(E)> elements{E...};
    constexpr static std::array<int16_t, (1 << algebra_t::metric_t::dimension)> element_lut
        = detail::construct_element_lut<E...>(
            std::make_index_sequence<(1 << algebra_t::metric_t::dimension)>{});

    std::array<T, sizeof...(E)> data_;

    // NOTE: in GAL code, `ie` refers always to "indeterminate expression"
    GAL_NODISCARD constexpr static auto ie(uint32_t id) noexcept
    {
        return detail::construct_ie<A>(id,
                                       std::make_integer_sequence<width_t, sizeof...(E)>{},
                                       std::integer_sequence<elem_t, E...>{});
    }

    GAL_NODISCARD constexpr static size_t size() noexcept
    {
        return sizeof...(E);
    }

    template <elem_t... S>
    GAL_NODISCARD constexpr auto select() const noexcept
    {
        if constexpr (sizeof...(S) == 1)
        {
            return select(S...);
        }
        else
        {
            return std::array<T, sizeof...(S)>{(element_lut[S] == -1 ? 0 : data_[element_lut[S]])...};
        }
    }

    GAL_NODISCARD constexpr T select(elem_t e) const noexcept
    {
        auto index = element_lut[e];
        return (index == -1 ? T{} : data_[index]);
    }

    GAL_NODISCARD constexpr T* select(elem_t e) noexcept
    {
        auto index = element_lut[e];
        return (index == -1 ? nullptr : &data_[index]);
    }

    GAL_NODISCARD constexpr T const* data() const noexcept
    {
        return data_.data();
    }

    GAL_NODISCARD constexpr auto begin() noexcept
    {
        return data_.begin();
    }

    GAL_NODISCARD constexpr auto end() noexcept
    {
        return data_.end();
    }

    GAL_NODISCARD constexpr const T& operator[](size_t index) const noexcept
    {
        return data_[index];
    }

    GAL_NODISCARD constexpr T& operator[](size_t index) noexcept
    {
        return data_[index];
    }
};

template <typename A, typename T>
struct scalar
{
    using algebra_t = A;
    using value_t   = T;

    GAL_NODISCARD constexpr static size_t size() noexcept
    {
        return 1;
    }

    // NOTE: in GAL code, `ie` refers always to "indeterminate expression"
    GAL_NODISCARD constexpr static mv<A, 1, 1, 1> ie(uint32_t id) noexcept
    {
        return {mv_size{1, 1, 1}, {ind{id, one}}, {mon{one, one, 1, 0}}, {term{1, 0, 0}}};
    }

    GAL_NODISCARD constexpr T const* data() const noexcept
    {
        return &value;
    }

    GAL_NODISCARD constexpr operator T() const noexcept
    {
        return value;
    }

    GAL_NODISCARD constexpr T const& operator[](size_t) const noexcept
    {
        return value;
    }

    GAL_NODISCARD constexpr T& operator[](size_t) noexcept
    {
        return value;
    }

    T value;
};

namespace detail
{
    template <typename T>
    struct is_scalar
    {
        constexpr static bool value = false;
    };

    template <typename A, typename T>
    struct is_scalar<scalar<A, T>>
    {
        constexpr static bool value = true;
    };

    template <typename T>
    constexpr inline bool is_scalar_v = is_scalar<T>::value;
} // namespace detail
} // namespace gal
