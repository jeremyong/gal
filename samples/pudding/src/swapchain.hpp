#pragma once

#include "context.hpp"

#include <image_view.hpp>
#include <render_target.hpp>
#include <vulkan/vulkan.h>
#include <cstdint>
#include <optional>

struct SDL_Window;

class swapchain
{
public:
    swapchain(SDL_Window& window);

    void create();

    std::optional<std::pair<pd::render_target const&, VkSemaphore>> acquire();
    void present(VkSemaphore signal);

private:
    // Resets the swapchain in the event of an invalidation (e.g. window resize)
    void reset();

    SDL_Window* window_;
    VkSwapchainKHR swapchain_ = VK_NULL_HANDLE;
    VkSurfaceKHR surface_;
    VkSurfaceFormatKHR format_;
    VkPresentModeKHR mode_;
    std::array<VkSemaphore, 3> acquire_sems_;
    std::array<VkImage, 3> images_;
    std::array<pd::image_view, 3> views_;
    std::array<pd::render_target, 3> rts_;

    VkSurfaceCapabilitiesKHR surface_caps_;

    int width_;
    int height_;
    uint32_t active_index_;
    uint32_t image_count_;
};