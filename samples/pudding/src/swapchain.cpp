#include "swapchain.hpp"

#include <SDL.h>
#include <SDL_vulkan.h>
#include <vector>

#include "context.hpp"
#include "ensure.hpp"
#include "log.hpp"

swapchain::swapchain(SDL_Window& window)
{
    window_ = &window;
    ENSURE(SDL_Vulkan_CreateSurface(window_, context.instance, &surface_), "Unable to create vulkan surface");

    VkBool32 surface_supported;
    vkGetPhysicalDeviceSurfaceSupportKHR(
        context.physical_device, context.present_queue_family, surface_, &surface_supported);
    ENSURE(surface_supported, "The surface associated with a created swapchain is unsupported");

    // https://vulkan.lunarg.com/doc/view/1.0.26.0/linux/vkspec.chunked/ch29s05.html
    // Query available presentation modes and swapchain surface formats
    uint32_t format_count;
    vkGetPhysicalDeviceSurfaceFormatsKHR(context.physical_device, surface_, &format_count, nullptr);
    ENSURE(format_count > 0, "Unable to find any suitable surface formats");
    std::vector<VkSurfaceFormatKHR> formats;
    formats.resize(format_count);
    vkGetPhysicalDeviceSurfaceFormatsKHR(context.physical_device, surface_, &format_count, formats.data());

    if (format_count == 1 && formats[0].format == VK_FORMAT_UNDEFINED)
    {
        // This indicates that all formats are supported
        format_.format = VK_FORMAT_R8G8B8A8_UNORM;
        format_.colorSpace = VK_COLOR_SPACE_SRGB_NONLINEAR_KHR;
    }
    else
    {
        for (auto& format : formats)
        {
            if (format.format == VK_FORMAT_B8G8R8A8_UNORM && format.colorSpace == VK_COLOR_SPACE_SRGB_NONLINEAR_KHR)
            {
                // Not all formats are supported but we found the one we want
                format_.format = VK_FORMAT_B8G8R8A8_UNORM;
                format_.colorSpace = VK_COLOR_SPACE_SRGB_NONLINEAR_KHR;
                break;
            }
        }
    }

    if (format_.format == VK_FORMAT_UNDEFINED)
    {
        // We weren't able to find the format that we wanted so use the first one available instead
        format_ = formats[0];
        WARN("Preferred surface format not found; falling back to {}:{} instead", format_.format, format_.colorSpace);
    }

    uint32_t mode_count;
    vkGetPhysicalDeviceSurfacePresentModesKHR(context.physical_device, surface_, &mode_count, nullptr);
    std::vector<VkPresentModeKHR> modes;
    modes.resize(mode_count);
    vkGetPhysicalDeviceSurfacePresentModesKHR(context.physical_device, surface_, &mode_count, modes.data());

    // Check if the desired presentation mode is available
    VkPresentModeKHR desired = context.settings.vsync ? VK_PRESENT_MODE_FIFO_KHR : VK_PRESENT_MODE_MAILBOX_KHR;
    bool fifo_available = false;
    bool mailbox_available = false;
    bool relaxed_fifo_available = false;
    for (auto const& mode : modes)
    {
        if (mode == VK_PRESENT_MODE_FIFO_KHR)
        {
            fifo_available = true;
        }

        if (mode == VK_PRESENT_MODE_MAILBOX_KHR)
        {
            mailbox_available = true;
        }

        if (mode == VK_PRESENT_MODE_FIFO_RELAXED_KHR)
        {
            relaxed_fifo_available = true;
        }

        if (mode == desired)
        {
            // Found the mode we want
            mode_ = mode;
        }
    }

    if (mode_ != desired)
    {
        if (desired == VK_PRESENT_MODE_MAILBOX_KHR)
        {
            // Users turning off vsync want a render rate that exceeds the refresh rate.
            if (relaxed_fifo_available)
            {
                mode_ = VK_PRESENT_MODE_FIFO_RELAXED_KHR;
            }
            else
            {
                mode_ = VK_PRESENT_MODE_IMMEDIATE_KHR;
            }
        }
        else
        {
            if (mailbox_available)
            {
                // Prioritize present modes that do not exhibit tearing
                mode_ = VK_PRESENT_MODE_MAILBOX_KHR;
            }
            else if (relaxed_fifo_available)
            {
                mode_ = VK_PRESENT_MODE_FIFO_RELAXED_KHR;
            }
            else
            {
                mode_ = VK_PRESENT_MODE_IMMEDIATE_KHR;
            }
        }
    }

    switch (mode_)
    {
    case VK_PRESENT_MODE_FIFO_KHR:
        DEBUG("Creating swapchain using FIFO presentation");
        break;
    case VK_PRESENT_MODE_FIFO_RELAXED_KHR:
        DEBUG("Creating swapchain using FIFO relaxed presentation");
        break;
    case VK_PRESENT_MODE_MAILBOX_KHR:
        DEBUG("Creating swapchain using mailbox presentation");
        break;
    case VK_PRESENT_MODE_IMMEDIATE_KHR:
        DEBUG("Creating swapchain using immediate presentation");
        break;
    default:
        break;
    }

    create();
    DEBUG("Swapchain of size {} created", image_count_);

    for (auto& sem : acquire_sems_)
    {
        VkSemaphoreCreateInfo sem_info{};
        sem_info.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;
        ENSURE(vkCreateSemaphore(context.device, &sem_info, nullptr, &sem), "Unable to create vk semaphore");
    }
}

void swapchain::create()
{
    // Query surface capabilities
    vkGetPhysicalDeviceSurfaceCapabilitiesKHR(context.physical_device, surface_, &surface_caps_);

    uint8_t desired_image_count = context.settings.vsync ? 2 : 3;

    if (surface_caps_.maxImageCount < desired_image_count)
    {
        DEBUG("Unable to accommodate requested swapchain image size of {}, dropping to {}",
              desired_image_count, surface_caps_.maxImageCount);
        desired_image_count = surface_caps_.maxImageCount;
    }

    // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkSwapchainCreateInfoKHR.html
    VkSwapchainCreateInfoKHR swapchain_info{};
    swapchain_info.sType = VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR;
    swapchain_info.surface = surface_;
    swapchain_info.minImageCount = desired_image_count;
    swapchain_info.imageFormat = format_.format;
    swapchain_info.imageColorSpace = format_.colorSpace;

    SDL_GetWindowSize(window_, &width_, &height_);
    if (surface_caps_.currentExtent.width == std::numeric_limits<uint32_t>::max())
    {
        // Set swapchain extent to the size of the window
        swapchain_info.imageExtent.width = width_;
        swapchain_info.imageExtent.height = height_;
    }
    else
    {
        swapchain_info.imageExtent = surface_caps_.currentExtent;
    }
    width_ = swapchain_info.imageExtent.width;
    height_ = swapchain_info.imageExtent.height;

    // Not supporting VR at this time
    swapchain_info.imageArrayLayers = 1;
    swapchain_info.imageUsage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
    swapchain_info.preTransform = surface_caps_.currentTransform; // Keep the transform as the identity
    swapchain_info.compositeAlpha =
        VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR; // Don't bother trying to support window managers with alpha compositing
    swapchain_info.presentMode = mode_;
    // Necessarily render to pixels even if they are offscreen as they may be read back or used in debugging
    swapchain_info.clipped = VK_FALSE;
    VkSwapchainKHR old_swapchain = swapchain_;
    swapchain_info.oldSwapchain = old_swapchain;

    DEBUG("Creating {} by {} px swapchain", width_, height_);
    ENSURE(vkCreateSwapchainKHR(context.device, &swapchain_info, nullptr, &swapchain_), "Failed to create swapchain");

    if (old_swapchain != VK_NULL_HANDLE)
    {
        vk_dispose(old_swapchain);
    }

    vkGetSwapchainImagesKHR(context.device, swapchain_, &image_count_, nullptr);
    ENSURE(image_count_ < images_.size(), "Not enough space allocated for requested swapchain size");
    ENSURE(vkGetSwapchainImagesKHR(context.device, swapchain_, &image_count_, images_.data()),
           "Failed to retrieve swapchain images");

    for (uint32_t i = 0; i != image_count_; ++i)
    {
        VkImageViewCreateInfo view_info = {};
        view_info.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
        view_info.image = images_[i];
        view_info.format = format_.format;
        view_info.flags = 0;
        view_info.viewType = VK_IMAGE_VIEW_TYPE_2D;
        view_info.components.r = VK_COMPONENT_SWIZZLE_IDENTITY;
        view_info.components.g = VK_COMPONENT_SWIZZLE_IDENTITY;
        view_info.components.b = VK_COMPONENT_SWIZZLE_IDENTITY;
        view_info.components.a = VK_COMPONENT_SWIZZLE_IDENTITY;
        view_info.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        view_info.subresourceRange.baseMipLevel = 0;
        view_info.subresourceRange.levelCount = 1;
        view_info.subresourceRange.baseArrayLayer = 0;
        view_info.subresourceRange.layerCount = 1;
        VkImageView view;
        ENSURE(vkCreateImageView(vk_device, &view_info, vk_ac, &view), "Unable to create swapchain image view");

        views_[i] = std::move(pd::image_view{static_cast<void*>(view), width_, height_});
        rts_[i]   = pd::render_target{
            &views_[i],                      // Image view
            "",                              // Debug name
            {1.0f, 0.0f, 1.0f, 1.0f},        // Clear color
            format_.format,                  // Format
            VK_IMAGE_LAYOUT_UNDEFINED,       // Initial layout
            VK_IMAGE_LAYOUT_PRESENT_SRC_KHR, // Final layout
            0b1,                             // Sample count bitfield
            4,                               // Channel count
            false,                           // Depth
            false,                           // Stencil
            true,                            // Clear
            false,                           // Load
            true                             // Store
        };
    }
}

std::optional<std::pair<pd::render_target const&, VkSemaphore>> swapchain::acquire()
{
    uint8_t frame_index = context.frame_index % 3;
    VkResult result = vkAcquireNextImageKHR(vk_device,
                                            swapchain_,
                                            std::numeric_limits<uint64_t>::max(),
                                            acquire_sems_[frame_index],
                                            VK_NULL_HANDLE,
                                            &active_index_);
    switch (result)
    {
    case VK_SUCCESS:
        return {{rts_[active_index_], acquire_sems_[frame_index]}};
    case VK_ERROR_OUT_OF_DATE_KHR:
        [[fallthrough]];
    case VK_SUBOPTIMAL_KHR:
        // Swapchain is out of date and needs to be recreated
        create();

        return {};
    default:
        break;
    }
    ENSURE(false, "Failed to acquire vk swapchain image");
    return {};
}

void swapchain::present(VkSemaphore signal)
{
    VkResult result;
    VkPresentInfoKHR present_info{};
    present_info.sType              = VK_STRUCTURE_TYPE_PRESENT_INFO_KHR;
    present_info.waitSemaphoreCount = 1;
    present_info.pWaitSemaphores    = &signal;
    present_info.swapchainCount     = 1;
    present_info.pSwapchains        = &swapchain_;
    present_info.pImageIndices      = &active_index_;
    present_info.pResults           = &result;
    vkQueuePresentKHR(context.general.q, &present_info);

    switch (result)
    {
        case VK_SUCCESS:
            return;
        case VK_ERROR_OUT_OF_DATE_KHR:
            [[fallthrough]];
        case VK_SUBOPTIMAL_KHR:
            create();
            break;
        default:
            ENSURE(result, "Failed to present swapchain image");
            break;
    }
}