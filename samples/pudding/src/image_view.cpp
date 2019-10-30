#include <image_view.hpp>

#include "context.hpp"

using namespace pd;

image_view::image_view(void* handle, int width, int height)
    : handle_{handle}
    , width_{width}
    , height_{height}
{}

image_view::image_view(image_view&& other)
    : uid_{std::move(other.uid_)}
    , handle_{other.handle_}
    , width_{other.width_}
    , height_{other.height_}
{
    other.handle_ = nullptr;
}

image_view& image_view::operator=(image_view&& other)
{
    std::swap(uid_, other.uid_);
    std::swap(handle_, other.handle_);
    width_ = other.width_;
    height_ = other.height_;
    return *this;
}

image_view::~image_view()
{
    if (handle_)
    {
        vk_dispose(static_cast<VkImageView>(handle_));
    }
}