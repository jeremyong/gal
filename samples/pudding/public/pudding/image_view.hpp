#pragma once

#include "detail/uid.hpp"

namespace pd
{
class image_view
{
public:
    image_view() = default;
    image_view(void* handle, int width, int height);
    image_view(image_view&& other);
    image_view& operator=(image_view&& other);
    ~image_view();

    [[nodiscard]] auto hash() const noexcept
    {
        return uid_.hash;
    }

    [[nodiscard]] void* handle() const noexcept
    {
        return handle_;
    }

    [[nodiscard]] int width() const noexcept
    {
        return width_;
    }

    [[nodiscard]] int height() const noexcept
    {
        return height_;
    }

private:
    detail::unique_id uid_;
    void* handle_ = nullptr;
    int width_;
    int height_;
};
} // namespace pd
