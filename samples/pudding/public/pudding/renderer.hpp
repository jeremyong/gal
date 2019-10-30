#pragma once

#include <vector>

#include "device.hpp"

namespace pd
{
class renderer
{
public:
    struct desc
    {
        char const* app_name;
        char const* asset_path;

        // If vsync is used, a FIFO swapchain is used and FPS will be capped to
        // the screen's refresh rate. Otherwise, the swapchain will be triple buffered
        // and rendering will occur at a possibly greater rate than the screen's
        // refresh rate.
        bool vsync = true;

        bool enable_validation = false;
        bool enable_renderdoc  = false;
    };

    renderer(desc d);
    renderer(renderer&&);
    renderer& operator=(renderer&&);

    renderer(renderer const&) = delete;
    renderer& operator=(renderer const&) = delete;

    void init_default_device();
    void init_device(uint32_t index);

private:
    std::vector<device> devices_;
};
} // namespace pd