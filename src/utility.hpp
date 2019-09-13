#pragma once

#include <cstdint>

#define GAL_RESTRICT __restrict

namespace gal
{
#if defined(__clang__) || defined(__GNUG__)

[[nodiscard]] constexpr int count_bits64(uint64_t input) noexcept
{
    return __builtin_popcountll(input);
}

[[nodiscard]] constexpr int count_bits(uint32_t input) noexcept
{
    return __builtin_popcount(input);
}

[[nodiscard]] constexpr int count_trailing_zeros64(uint64_t input) noexcept
{
    return __builtin_ctzll(input);
}

[[nodiscard]] constexpr int count_trailing_zeros(uint32_t input) noexcept
{
    return __builtin_ctz(input);
}

[[nodiscard]] constexpr int count_leading_zeros(uint32_t input) noexcept
{
    return __builtin_clz(input);
}

// Returns the index of the most significant set bit (0-indexed)
[[nodiscard]] constexpr int leading_index(uint64_t input) noexcept
{
    return 63 - __builtin_clzll(input);
}

// Returns the index of the most significant set bit (0-indexed)
[[nodiscard]] constexpr int leading_index(uint32_t input) noexcept
{
    return 31 - __builtin_clz(input);
}


#elif defined(_MSC_VER)

// TODO: test MSVC

#include <intrin.h>
#pragma intrinsic(_BitScanReverse)
#pragma intrinsic(_BitScanReverse64)
#pragma intrinsic(_BitScanForward)
#pragma intrinsic(_BitScanForward64)

constexpr int count_bits64(uint64_t input) noexcept
{
    // MSVC popcnt returns an unsigned int
    return static_cast<int>(__popcnt64(input));
}

constexpr int count_bits(uint32_t input) noexcept
{
    // MSVC popcnt returns an unsigned int
    return static_cast<int>(__popcnt(input));
}

constexpr int count_trailing_zeros64(uint64_t input) noexcept
{
    unsigned long index;
    if (_BitScanForward64(&index, input))
    {
        return index;
    }
    else
    {
        return 0;
    }
}

constexpr int count_trailing_zeros(uint32_t input) noexcept
{
    unsigned long index;
    if (_BitScanForward(&index, input))
    {
        return index;
    }
    else
    {
        return 0;
    }
}

constexpr int count_leading_zeros_complement(uint32_t input) noexcept
{
    unsigned long index;
    if (_BitScanReverse(&index, input))
    {
        return index;
    }
    else
    {
        return 0;
    }
}

#endif
}