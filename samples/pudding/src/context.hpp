#pragma once

#include <array>
#include <detail/uid.hpp>
#include <memory>
#include <unordered_map>
#include <vector>
#include <vk_mem_alloc.h>
#include <vulkan/vulkan.h>

#include "gc.hpp"

using pd::detail::unique_id;

// For now, we don't ever bother rendering more than one frame ahead
constexpr inline uint8_t frame_context_count = 2;

struct queue
{
    VkQueue q;
    uint32_t family;
};

struct command_pool
{
    VkCommandPool pool;
    VkFence fence;
    bool active = false;
    // At the start of a frame, if this flag is set, the fence is assumed to be submitted so it is waited on and the
    // command buffer pool is reset
    bool dirty = false;
};

struct frame
{
    command_pool general;
    command_pool dma;
    command_pool compute;
    VkSemaphore signal;
    gc garbage;
};

struct vk_context
{
    VkInstance instance;
#ifndef NDEBUG
    VkDebugUtilsMessengerEXT debug_messenger;
#endif
    VkDevice device;
    VkPhysicalDevice physical_device;
    VmaAllocator allocator;

    VkDescriptorPool descriptor_pool;

    queue general;
    queue dma;
    queue compute;
    uint32_t present_queue_family;

    VkCommandPool static_cmd_pool;
    std::array<frame, frame_context_count> frames;
    uint8_t frame_index;

    VkPhysicalDeviceLimits limits;

    struct {
        bool vsync;
    } settings;

    // Cached vulkan objects used in rendering
    struct
    {
        std::unordered_map<unique_id::hash_t, VkRenderPass> render_passes;
        std::unordered_map<unique_id::hash_t, VkFramebuffer> framebuffers;
    } objects;
};

inline vk_context context{};
inline VkInstance& vk_instance{context.instance};
inline VkDevice& vk_device{context.device};
// Placeholder allocation callbacks in case there is time to implement this later
inline VkAllocationCallbacks* vk_ac{nullptr};
inline auto& vk_objects{context.objects};

[[nodiscard]] inline frame& vk_current_frame() noexcept
{
    return context.frames[context.frame_index % frame_context_count];
}

template <typename T>
void vk_dispose(T handle) noexcept
{
    vk_current_frame().garbage.dispose(handle);
}