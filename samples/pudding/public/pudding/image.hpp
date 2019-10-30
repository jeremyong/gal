#pragma once

#include <cstdint>

namespace pd
{
enum class image_usage : uint8_t
{
    texture,
    attachment,
    depth,
    depth_stencil
};

enum class image_compression : uint8_t
{
    none,
    bc1,
    bc2,
    bc3
};

struct image_desc
{
    int width;
    int height;
    int depth = 1;

    // 4 byte block
    uint8_t dim = 2; // 1, 2, or 3
    uint8_t channel_count = 4;
    bool srgb = true;
    image_usage usage{image_usage::texture};

    // TODO: support HDR images

    int mip_count = 1;
    int slice_count = 1;
    image_compression compression{image_compression::none};
};

class image
{
public:
    image() = default;
    explicit image(image_desc const& desc);
private:
    image_desc desc_{};
    void* handle_;
};
} // namespace pd
