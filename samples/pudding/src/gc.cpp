#include "gc.hpp"

#include "context.hpp"

void gc::flush()
{
    for (auto sc : swapchains_)
    {
        vkDestroySwapchainKHR(vk_device, sc, vk_ac);
    }
    swapchains_.clear();

    for (auto iv : image_views_)
    {
        vkDestroyImageView(vk_device, iv, vk_ac);
    }
    image_views_.clear();
}

void gc::dispose(VkSwapchainKHR sc)
{
    swapchains_.emplace_back(sc);
}

void gc::dispose(VkImageView iv)
{
    image_views_.emplace_back(iv);
}