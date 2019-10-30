#pragma once

#include <memory>

struct SDL_Window;
class swapchain;

namespace pd
{
struct window_desc
{
    char const* title  = "pudding!";
    int width          = 800;
    int height         = 600;
    bool fullscreen    = false;
    bool resizable     = false;
    bool grab_input    = false;
    bool high_dpi      = false;
    bool always_on_top = false;
};

// Encapsulates an OS window
class window final
{
public:
    explicit window(window_desc const& wd);
    ~window();
    window(window&&);
    window& operator=(window&&);
    window(window const&) = delete;
    window& operator=(window const&) = delete;

    [[nodiscard]] constexpr char const* title() const noexcept
    {
        return title_;
    }

    [[nodiscard]] uint32_t id() const noexcept;

    [[nodiscard]] swapchain& sc()
    {
        return *sc_;
    }

private:
    SDL_Window* window_ = nullptr;
    swapchain* sc_ = nullptr;
    char const* title_ = "";
};
} // namespace pd