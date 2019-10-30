#pragma once

#include <vector>
#include <vulkan/vulkan.h>

// Garbage chute. Items and objects sent to the chute are reclaimed in a future frame. The idea
// is to have one gc per frame context such that the flush of old objects occurs at the start of that frame.
class gc
{
public:
    void flush();

    void dispose(VkSwapchainKHR sc);
    void dispose(VkImageView iv);

private:
    std::vector<VkSwapchainKHR> swapchains_;
    std::vector<VkImageView> image_views_;
};
