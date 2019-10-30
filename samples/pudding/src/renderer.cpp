#include <renderer.hpp>

#include <SDL.h>
#include <SDL_vulkan.h>
#include <stdexcept>
#include <vector>
#include <vulkan/vulkan.h>

#include "context.hpp"
#include "ensure.hpp"
#include "log.hpp"

using namespace pd;

VkBool32 VKAPI_PTR debug_utils_messenger(VkDebugUtilsMessageSeverityFlagBitsEXT severity,
                                         VkDebugUtilsMessageTypeFlagsEXT type,
                                         const VkDebugUtilsMessengerCallbackDataEXT* data,
                                         void* user_data)
{
    const char* message_id_name = data->pMessageIdName ? data->pMessageIdName : "[x]";
    switch (severity)
    {
    case VK_DEBUG_UTILS_MESSAGE_SEVERITY_VERBOSE_BIT_EXT:
        DEBUG("[VK] {} - {}", message_id_name, data->pMessage);
        break;
    case VK_DEBUG_UTILS_MESSAGE_SEVERITY_INFO_BIT_EXT:
        INFO("[VK] {} - {}", message_id_name, data->pMessage);
        break;
    case VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT:
        WARN("[VK] {} - {}", message_id_name, data->pMessage);
        break;
    case VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT:
        ERROR("[VK] {} - {}", message_id_name, data->pMessage);
        break;
    default:
        DEBUG("[VK] {}", data->pMessage);
        break;
    }
    return VK_TRUE;
}

renderer::renderer(desc d)
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Vulkan_LoadLibrary(nullptr);

    VkApplicationInfo app_info{};
    app_info.sType              = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    app_info.engineVersion      = VK_MAKE_VERSION(1, 0, 0);
    app_info.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
    app_info.apiVersion         = VK_MAKE_VERSION(1, 1, 0);
    app_info.pEngineName        = "pudding";
    app_info.pApplicationName   = d.app_name;

    VkInstanceCreateInfo create_info = {};
    create_info.sType                = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    create_info.pApplicationInfo     = &app_info;

    const char* ext[32];
    uint32_t ext_count = 0;
    ENSURE(SDL_Vulkan_GetInstanceExtensions(nullptr, &ext_count, nullptr), "Unable to query for Vulkan extensions");

    ENSURE(ext_count < 32, "Too many extensions returned");

    ENSURE(SDL_Vulkan_GetInstanceExtensions(nullptr, &ext_count, ext), "Unable to query for Vulkan extensions");

    bool use_validation  = false;
    bool use_renderdoc   = false;
    uint32_t layer_count = 0;
    std::vector<char const*> layers;
#ifndef NDEBUG
    ext[ext_count++] = VK_EXT_DEBUG_UTILS_EXTENSION_NAME;
    vkEnumerateInstanceLayerProperties(&layer_count, nullptr);
    std::vector<VkLayerProperties> layer_props;
    DEBUG("{} layers found", layer_count);
    layer_props.resize(layer_count);
    vkEnumerateInstanceLayerProperties(&layer_count, layer_props.data());

    for (auto prop : layer_props)
    {
        std::string name{prop.layerName};
        if (d.enable_validation && name == "VK_LAYER_LUNARG_standard_validation")
        {
            use_validation = true;
            layers.emplace_back("VK_LAYER_LUNARG_standard_validation");
        }
        else if (d.enable_renderdoc && name == "VK_LAYER_RENDERDOC_Capture")
        {
            use_renderdoc = true;
            layers.emplace_back("VK_LAYER_RENDERDOC_Capture");
        }
    }
#endif
    create_info.enabledLayerCount   = layers.size();
    create_info.ppEnabledLayerNames = layers.data();

    create_info.enabledExtensionCount   = ext_count;
    create_info.ppEnabledExtensionNames = ext;

    ENSURE(vkCreateInstance(&create_info, nullptr, &context.instance), "Could not create vulkan instance");
    DEBUG("Vulkan instance initialized");

    if (use_validation)
    {
        // Create debug messenger
        VkDebugUtilsMessengerCreateInfoEXT debug_info = {};
        debug_info.sType                              = VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT;
        debug_info.messageSeverity
            = VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_SEVERITY_VERBOSE_BIT_EXT
              | VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_SEVERITY_INFO_BIT_EXT;
        debug_info.messageType = VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT
                                 | VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT
                                 | VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT;
        debug_info.pfnUserCallback = debug_utils_messenger;
        auto proc_addr             = (PFN_vkGetInstanceProcAddr)SDL_Vulkan_GetVkGetInstanceProcAddr();
        auto create_debug_messenger
            = (PFN_vkCreateDebugUtilsMessengerEXT)proc_addr(context.instance, "vkCreateDebugUtilsMessengerEXT");
        if (create_debug_messenger)
        {
            ENSURE(create_debug_messenger(context.instance, &debug_info, nullptr, &context.debug_messenger),
                   "Could not create debug messenger");
            DEBUG("Vulkan debug messenger activated");
        }
        else
        {
            ERROR("Unable to create debug messenger");
        }
    }

    uint32_t device_count;
    vkEnumeratePhysicalDevices(context.instance, &device_count, nullptr);
    ENSURE(device_count > 0, "No valid gfx devices found");
    ENSURE(device_count < 16, "Too many gfx devices found! :D");

    std::array<VkPhysicalDevice, 16> physical_devices;
    vkEnumeratePhysicalDevices(context.instance, &device_count, physical_devices.data());

    for (uint32_t i = 0; i != device_count; ++i)
    {
        auto&& physical_device = physical_devices[i];
        // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkPhysicalDeviceProperties.html
        // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkPhysicalDeviceType.html
        VkPhysicalDeviceProperties props;
        vkGetPhysicalDeviceProperties(physical_device, &props);

        // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkExtensionProperties.html
        uint32_t ext_prop_count;
        std::vector<VkExtensionProperties> ext_props;
        vkEnumerateDeviceExtensionProperties(physical_device, nullptr, &ext_prop_count, nullptr);
        ext_props.resize(ext_prop_count);
        vkEnumerateDeviceExtensionProperties(physical_device, nullptr, &ext_prop_count, ext_props.data());

        bool gfx_enabled = false;
        for (auto ext_prop : ext_props)
        {
            // Check a swapchain extension is available for this device
            if (strncmp(ext_prop.extensionName, VK_KHR_SWAPCHAIN_EXTENSION_NAME, strlen(VK_KHR_SWAPCHAIN_EXTENSION_NAME)) == 0)
            {
                gfx_enabled = true;
                break;
            }
        }

        // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkPhysicalDeviceMemoryProperties.html
        VkPhysicalDeviceMemoryProperties memory_props;
        vkGetPhysicalDeviceMemoryProperties(physical_device, &memory_props);
        uint64_t vram_bytes = 0;
        for (uint32_t i = 0; i != memory_props.memoryHeapCount; ++i)
        {
            auto& heap = memory_props.memoryHeaps[i];
            if (heap.flags & VK_MEMORY_HEAP_DEVICE_LOCAL_BIT)
            {
                vram_bytes += heap.size;
            }
        }

        devices_.emplace_back(device{props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU,
                                     gfx_enabled,
                                     props.vendorID,
                                     vram_bytes,
                                     props.deviceName,
                                     static_cast<void*>(physical_device)});

        DEBUG("Found device: {}", devices_.back().to_string());
    }

    {
        // Write any relevant settings to the context
        context.settings.vsync = d.vsync;
    }
}

void renderer::init_default_device()
{
    // Not a particularly clever heuristic, but should suffice in many cases
    uint32_t candidate = 0;
    // NOTE: we're guaranteed to have at least 1 suitable device at this point
    device& candidate_device = devices_[candidate];

    for (uint32_t i = 1; i != devices_.size(); ++i)
    {
        device& dev = devices_[i];
        if ((!candidate_device.gfx_enabled && dev.gfx_enabled)    // Prioritize gfx enabled devices
            || (!candidate_device.discrete && dev.discrete) // Prioritize discrete devices
            || (candidate_device.vram_bytes < dev.vram_bytes)       // Prioritize devices with more vram
            // TODO insert other heuristics here
        )
        {
            candidate = i;
            candidate_device = dev;
        }
    }

    INFO("Device selected: {}", candidate_device.to_string());
    init_device(candidate);
}

void renderer::init_device(uint32_t index)
{
    auto& dev = devices_[index];
    auto physical_device = static_cast<VkPhysicalDevice>(dev.handle);
    context.physical_device = physical_device;

    VkPhysicalDeviceProperties props;
    vkGetPhysicalDeviceProperties(physical_device, &props);
    context.limits = props.limits;

    // Create a logical device and required queues
    VkDeviceCreateInfo device_info = {};
    device_info.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;

    // TODO consider enabling additional extensions as needed to improve performance or fidelity
    const char* swapchain_ext = VK_KHR_SWAPCHAIN_EXTENSION_NAME;
    device_info.ppEnabledExtensionNames = &swapchain_ext;
    device_info.enabledExtensionCount = 1;

    // Query physical device queue families
    // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkQueueFamilyProperties.html
    std::vector<VkQueueFamilyProperties> queue_families;
    uint32_t queue_family_count;
    vkGetPhysicalDeviceQueueFamilyProperties(physical_device, &queue_family_count, nullptr);
    queue_families.resize(queue_family_count);
    vkGetPhysicalDeviceQueueFamilyProperties(physical_device, &queue_family_count, queue_families.data());

    // Create device queues
    // The general strategy for Nvidia is to use 1 general queue and 1 DMA queue.
    // AMD cards support async compute which adds an additional compute queue (ensure all images have exclusive access).
    // Intel does not benefit from a designated DMA queue. In this case, a single general queue is used.
    // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkDeviceQueueCreateInfo.html

    float queue_priority = 1.0f; // All queues are created with the same priority
    std::array<VkDeviceQueueCreateInfo, 3> queue_infos;
    device_info.pQueueCreateInfos = queue_infos.data();

    uint32_t& general_family_index = context.general.family;
    uint32_t& dma_family_index = context.dma.family;
    uint32_t& compute_family_index = context.compute.family;

    auto& general_queue_info = queue_infos[0];
    device_info.queueCreateInfoCount = 1;

    general_queue_info.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
    general_queue_info.queueCount = 1;
    general_queue_info.flags = 0;
    general_queue_info.pNext = nullptr;
    general_queue_info.pQueuePriorities = &queue_priority;
    for (uint32_t i = 0; i != queue_families.size(); ++i)
    {
        auto& family = queue_families[i];
        if (family.queueFlags & VK_QUEUE_GRAPHICS_BIT)
        {
            general_family_index = i;
            general_queue_info.queueFamilyIndex = i;
            DEBUG("General queue family: {}", i);
            break;
        }
    }

    context.present_queue_family = general_family_index;

    if (dev.is_amd() || dev.is_nvidia())
    {
        // Try to supply a dedicated DMA queue
        for (uint32_t i = 0; i != queue_families.size(); ++i)
        {
            auto& family = queue_families[i];
            if (family.queueFlags & VK_QUEUE_TRANSFER_BIT && !(family.queueFlags & VK_QUEUE_GRAPHICS_BIT))
            {
                dma_family_index = i;
                auto& queue_info = queue_infos[device_info.queueCreateInfoCount++];
                queue_info.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
                queue_info.queueCount = 1;
                queue_info.flags = 0;
                queue_info.pNext = nullptr;
                queue_info.queueFamilyIndex = i;
                queue_info.pQueuePriorities = &queue_priority;
                DEBUG("Dedicated DMA queue family: {}", i);
                break;
            }
        }
    }

    if (dev.is_amd() && device_info.queueCreateInfoCount == 2)
    {
        // Try to supply a dedicated compute queue (provided we were also able to find a dma queue family)
        // Only AMD supports an async command queue at this time
        for (uint32_t i = 0; i != queue_families.size(); ++i)
        {
            auto& family = queue_families[i];
            if (family.queueFlags & VK_QUEUE_COMPUTE_BIT && !(family.queueFlags & VK_QUEUE_GRAPHICS_BIT))
            {
                compute_family_index = i;
                auto& queue_info = queue_infos[device_info.queueCreateInfoCount++];
                queue_info.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
                queue_info.queueCount = 1;
                queue_info.flags = 0;
                queue_info.pNext = nullptr;
                queue_info.queueFamilyIndex = i;
                queue_info.pQueuePriorities = &queue_priority;
                DEBUG("Dedicated compute queue family: {}", i);
                break;
            }
        }
    }

    uint32_t queues_used = device_info.queueCreateInfoCount;
    ENSURE(vkCreateDevice(physical_device, &device_info, nullptr, &context.device),
           "Unable to create logical vulkan device");
    DEBUG("Logical Vulkan device created");

    {
        // Fetch a single queue for each queue family requested
        // We don't bother submitting on multiple queues within a family for now
        vkGetDeviceQueue(context.device, general_family_index, 0, &context.general.q);

        if (queues_used > 1)
        {
            vkGetDeviceQueue(context.device, dma_family_index, 0, &context.dma.q);
        }

        if (queues_used > 2)
        {
            vkGetDeviceQueue(context.device, compute_family_index, 0, &context.compute.q);
        }
    }

    {
        // Initialize the memory allocator
        VmaAllocatorCreateInfo alloc_info{};
        alloc_info.physicalDevice = physical_device;
        alloc_info.device = context.device;
        ENSURE(vmaCreateAllocator(&alloc_info, &context.allocator), "Unable to create vma allocator");
    }

    {
        // Initialize the various command pools and frame contexts
        VkCommandPoolCreateInfo command_pool_info{};
        command_pool_info.sType            = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
        command_pool_info.flags            = 0;
        command_pool_info.queueFamilyIndex = general_family_index;
        vkCreateCommandPool(context.device, &command_pool_info, nullptr, &context.static_cmd_pool);

        for (auto& frame_context : context.frames)
        {
            // Every frame context needs a general command pool
            VkCommandPoolCreateInfo command_pool_info{};
            command_pool_info.sType            = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
            command_pool_info.flags            = VK_COMMAND_POOL_CREATE_TRANSIENT_BIT;
            command_pool_info.queueFamilyIndex = general_family_index;
            ENSURE(vkCreateCommandPool(context.device, &command_pool_info, nullptr, &frame_context.general.pool),
                   "Unable to create vk command pool");
            frame_context.general.active = true;

            // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkFenceCreateInfo.html
            VkFenceCreateInfo general_fence_info{};
            general_fence_info.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
            // Create fences in an initially signaled state as they are always awaited at the start of the frame
            general_fence_info.flags = VK_FENCE_CREATE_SIGNALED_BIT;
            vkCreateFence(context.device, &general_fence_info, nullptr, &frame_context.general.fence);

            if (queues_used > 1)
            {
                command_pool_info.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
                command_pool_info.flags
                    = VK_COMMAND_POOL_CREATE_TRANSIENT_BIT | VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
                command_pool_info.queueFamilyIndex = dma_family_index;
                ENSURE(vkCreateCommandPool(context.device, &command_pool_info, nullptr, &frame_context.dma.pool),
                       "Unable to create vk command pool");
                frame_context.dma.active = true;

                VkFenceCreateInfo dma_fence_info{};
                dma_fence_info.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
                dma_fence_info.flags = VK_FENCE_CREATE_SIGNALED_BIT;
                vkCreateFence(context.device, &dma_fence_info, nullptr, &frame_context.dma.fence);
            }

            if (queues_used > 2)
            {
                command_pool_info.sType            = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
                command_pool_info.flags            = VK_COMMAND_POOL_CREATE_TRANSIENT_BIT;
                command_pool_info.queueFamilyIndex = compute_family_index;
                ENSURE(vkCreateCommandPool(context.device, &command_pool_info, nullptr, &frame_context.compute.pool),
                       "Unable to create vk command pool");
                frame_context.compute.active = true;

                VkFenceCreateInfo compute_fence_info{};
                compute_fence_info.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
                compute_fence_info.flags = VK_FENCE_CREATE_SIGNALED_BIT;
                vkCreateFence(context.device, &compute_fence_info, nullptr, &frame_context.compute.fence);
            }
        }
    }

    {
        // Create primary descriptor pool heap
        VkDescriptorPoolCreateInfo dpool{};
        dpool.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        dpool.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
        // With bindless, we really should not need too many descriptor sets
        dpool.maxSets = 64;
        // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkDescriptorPoolSize.html
        VkDescriptorPoolSize pool_sizes[] = {{VK_DESCRIPTOR_TYPE_SAMPLER, 32},
                                             {VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE, 128},
                                             {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 128},
                                             {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 128},
                                             {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC, 32},
                                             {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER_DYNAMIC, 32},
                                             {VK_DESCRIPTOR_TYPE_INPUT_ATTACHMENT, 32}};
        dpool.poolSizeCount               = sizeof(pool_sizes) / sizeof(VkDescriptorPoolSize);
        dpool.pPoolSizes                  = pool_sizes;
        ENSURE(vkCreateDescriptorPool(context.device, &dpool, nullptr, &context.descriptor_pool),
               "Unable to create vulkan descriptor pool");
    }

    {
        // Create the primary render submission signal
        for (auto& frame_context : context.frames)
        {
            VkSemaphoreCreateInfo sem_info{};
            sem_info.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;
            ENSURE(vkCreateSemaphore(vk_device, &sem_info, vk_ac, &frame_context.signal),
                   "Unable to create primary render signal");
        }
    }
}