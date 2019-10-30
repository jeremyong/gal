#pragma once

#include "log.hpp"

#include <SDL.h>
#include <cassert>
#include <vulkan/vulkan.h>

[[nodiscard]] constexpr bool check(VkResult result) noexcept
{
    if (result == VK_SUCCESS)
    {
        return true;
    }
    return false;
}

[[nodiscard]] constexpr bool check(VkBool32 result) noexcept
{
    if (result == VK_TRUE)
    {
        return true;
    }
    return false;
}

[[nodiscard]] constexpr bool check(SDL_bool result) noexcept
{
    if (result == SDL_TRUE)
    {
        return true;
    }
    return false;
}

[[nodiscard]] constexpr bool check(bool result) noexcept
{
    return result;
}

#define ENSURE(condition, ...) \
    if (!check(condition))     \
    {                          \
        ERROR(__VA_ARGS__);    \
        assert(false);         \
    }
