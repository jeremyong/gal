#include <image.hpp>

#include "ensure.hpp"
#include "context.hpp"

using namespace pd;

image::image(image_desc const& desc)
    : desc_{desc}
{
    VkImageCreateInfo create{};
    create.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
    create.flags = 0;
    switch (desc_.dim)
    {
    case 1:
        create.imageType = VK_IMAGE_TYPE_1D;
        break;
    case 2:
        create.imageType = VK_IMAGE_TYPE_2D;
        break;
    case 3:
        create.imageType = VK_IMAGE_TYPE_3D;
        break;
    default:
        ENSURE(false, "Unsupported image type {} used (use 1, 2, or 3 please)", desc_.dim);
        break;
    }

    if (desc_.srgb)
    {
        switch (desc_.channel_count)
        {
        case 1:
            create.format = VK_FORMAT_R8_SRGB;
            break;
        case 2:
            create.format = VK_FORMAT_R8G8_SRGB;
            break;
        case 3:
            create.format = VK_FORMAT_R8G8B8_SRGB;
            break;
        case 4:
            create.format = VK_FORMAT_R8G8B8A8_SRGB;
            break;
        }
    }
    else
    {
        switch (desc_.channel_count)
        {
        case 1:
            create.format = VK_FORMAT_R8_UNORM;
            break;
        case 2:
            create.format = VK_FORMAT_R8G8_UNORM;
            break;
        case 3:
            create.format = VK_FORMAT_R8G8B8_UNORM;
            break;
        case 4:
            create.format = VK_FORMAT_R8G8B8A8_UNORM;
            break;
        }
    }

    switch (desc_.usage)
    {
    case image_usage::texture:
        create.usage = VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_SAMPLED_BIT;
        break;
    case image_usage::attachment:
        create.usage = VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_INPUT_ATTACHMENT_BIT;
        break;
    case image_usage::depth:
        [[fallthrough]];
    case image_usage::depth_stencil:
        create.usage = VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT;
        break;
    default:
        ENSURE(false, "Unsupported image usage {}", static_cast<uint8_t>(desc_.usage));
        break;
    }

    create.extent.width = desc_.width;
    create.extent.height = desc_.height;
    create.extent.depth = desc_.depth;
    create.mipLevels = desc.mip_count;
    create.arrayLayers = desc.slice_count;
    create.samples = VK_SAMPLE_COUNT_1_BIT;
    create.tiling = VK_IMAGE_TILING_OPTIMAL;
    create.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    create.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
}