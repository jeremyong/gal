#include <window.hpp>

#include <algorithm>
#include <SDL.h>
#include <SDL_vulkan.h>

#include "context.hpp"
#include "swapchain.hpp"

using namespace pd;

window::window(window_desc const& wd)
{
    uint32_t flags = SDL_WINDOW_VULKAN;

    if (wd.fullscreen)
    {
        flags |= SDL_WINDOW_FULLSCREEN;
    }
    if (wd.resizable)
    {
        flags |= SDL_WINDOW_RESIZABLE;
    }
    if (wd.grab_input)
    {
        flags |= SDL_WINDOW_INPUT_GRABBED;
    }
    if (wd.high_dpi)
    {
        flags |= SDL_WINDOW_ALLOW_HIGHDPI;
    }
    if (wd.always_on_top)
    {
        flags |= SDL_WINDOW_ALWAYS_ON_TOP;
    }

    window_ = SDL_CreateWindow(wd.title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, wd.width, wd.height, flags);
    title_  = wd.title;

    sc_ = new swapchain(*window_);
}

window::window(window&& other)
{
    std::swap(window_, other.window_);
    std::swap(sc_, other.sc_);
}

window& window::operator=(window&& other)
{
    std::swap(window_, other.window_);
    std::swap(sc_, other.sc_);
    return *this;
}

window::~window()
{
    if (window_)
    {
        SDL_DestroyWindow(window_);
    }

    if (sc_)
    {
        delete sc_;
    }
}

uint32_t window::id() const noexcept
{
    return SDL_GetWindowID(window_);
}