#pragma once

#include "detail/uid.hpp"

namespace pd
{
class image_view;

// Careful. If you mutate any of these members with the exception of the image view, you are expected to invoke reset on
// the unique id so that dependent objects (e.g. render passes) can be invalidated.
struct render_target
{
    image_view* iv;
    char const* name     = "";
    union
    {
        float clear_value[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        struct
        {
            float depth;
            uint32_t stencil;
        } depth_stencil_clear;
    };
    // Opaque enum value specifying image view format
    int format = 0;
    // Opaque enum value specifying layout at the start of the render pass
    int initial_layout = 0;
    // Opaque enum value specifying the layout the image should transition to after the render pass
    int final_layout = 0;
    // Fields with a _bits suffix are bitfield representations
    uint8_t sample_count_bits = 0b1;
    uint8_t channel_count     = 4;
    // A render target that is neither depth nor stencil is a color target
    bool depth   = false;
    bool stencil = false;
    bool clear   = false;
    bool load    = false;
    bool store   = true;
    ::pd::detail::unique_id uid;
};
} // namespace pd