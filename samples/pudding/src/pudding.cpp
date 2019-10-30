#include <pudding.hpp>

#include <SDL.h>
#include <cmd_stream.hpp>
#include <cstdint>
#include <numeric>
#include <unordered_map>

#include "context.hpp"
#include "ensure.hpp"
#include "swapchain.hpp"

using namespace pd;
using window_id = uint32_t;

static std::unordered_map<window_id, pudding*> puddings;
// Timestamp in milliseconds when the key is pressed. 0 if released.
static std::unordered_map<SDL_Keycode, uint32_t> key_states;

using namespace pd;

pudding::pudding()
    : pudding{window_desc{}}
{}

pudding::pudding(window_desc const& wd)
    : window_{wd}
{
    puddings.emplace(std::make_pair(window_.id(), this));

    for (auto& frame_context : context.frames)
    {
        VkCommandBufferAllocateInfo cb_info{};
        cb_info.sType              = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
        cb_info.commandPool        = frame_context.general.pool;
        cb_info.commandBufferCount = 1;
        cb_info.level              = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        VkCommandBuffer cb;
        ENSURE(vkAllocateCommandBuffers(vk_device, &cb_info, &cb),
               "Unable to allocate vk command buffer");
        cbs_.emplace_back(static_cast<void*>(cb));

        VkSemaphoreCreateInfo sem_info{};
        sem_info.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;
        VkSemaphore sem;
        ENSURE(vkCreateSemaphore(vk_device, &sem_info, vk_ac, &sem), "Unable to create pudding signal semaphore");
        signals_.emplace_back(sem);
    }
}

pudding::~pudding()
{
    puddings.erase(window_.id());

    size_t i = 0;
    for (auto& frame_context : context.frames)
    {
        vkFreeCommandBuffers(vk_device, frame_context.general.pool, 1, reinterpret_cast<VkCommandBuffer*>(&cbs_[0]));
        vkDestroySemaphore(vk_device, static_cast<VkSemaphore>(signals_[i]), vk_ac);
        ++i;
    }
}

bool pudding::is_key_down(int key) const
{
    return key_states[key] != 0;
}

bool pudding::is_left_mouse_button_down() const
{
    auto state = SDL_GetMouseState(nullptr, nullptr);
    return state & SDL_BUTTON(SDL_BUTTON_LEFT);
}

bool pudding::is_right_mouse_button_down() const
{
    auto state = SDL_GetMouseState(nullptr, nullptr);
    return state & SDL_BUTTON(SDL_BUTTON_RIGHT);
}

bool pudding::is_middle_mouse_button_down() const
{
    auto state = SDL_GetMouseState(nullptr, nullptr);
    return state & SDL_BUTTON(SDL_BUTTON_MIDDLE);
}

gal::pga2::point<int> pudding::mouse_xy() const
{
    int x;
    int y;
    SDL_GetMouseState(&x, &y);
    return {x, y};
}

static void advance_frame()
{
    ++context.frame_index;
    auto& frame = vk_current_frame();
    frame.garbage.flush();
    bool reset_general = frame.general.active && frame.general.dirty;
    bool reset_dma = frame.dma.active && frame.dma.dirty;
    bool reset_compute = frame.compute.active && frame.compute.dirty;
    VkFence fences[3];
    uint32_t fence_count = 0;
    if (reset_general)
    {
        fences[fence_count++] = frame.general.fence;
        frame.general.dirty = false;
    }
    if (reset_dma)
    {
        fences[fence_count++] = frame.dma.fence;
        frame.dma.dirty = false;
    }
    if (reset_compute)
    {
        fences[fence_count++] = frame.compute.fence;
        frame.compute.dirty = false;
    }
    vkWaitForFences(vk_device, fence_count, fences, VK_TRUE, std::numeric_limits<uint64_t>::max());
    vkResetFences(vk_device, fence_count, fences);

    if (reset_general)
    {
        ENSURE(vkResetCommandPool(vk_device, frame.general.pool, 0), "Unable to reset general command pool");
    }
    if (reset_dma)
    {
        ENSURE(vkResetCommandPool(vk_device, frame.dma.pool, 0), "Unable to reset dma command pool");
    }
    if (reset_compute)
    {
        ENSURE(vkResetCommandPool(vk_device, frame.compute.pool, 0), "Unable to reset compute command pool");
    }
}

void pudding::start()
{
    std::vector<pudding*> pds;
    std::vector<VkSubmitInfo> submits;
    std::vector<VkCommandBuffer> cbs;
    std::vector<VkSemaphore> waits;
    std::vector<VkSemaphore> signals;
    std::vector<VkPipelineStageFlags> wait_stages;

    while (true)
    {
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            switch (event.type)
            {
            case SDL_QUIT:
                for (auto&& p : puddings)
                {
                    p.second->on_quit();
                }
                return;
            case SDL_KEYDOWN: {
                key_states[event.key.keysym.sym] = event.key.timestamp;
                auto p = puddings.find(event.key.windowID);
                if (p == puddings.end())
                {
                    break;
                }
                p->second->on_key_down(event.key.keysym.sym);
                break;
            }
            case SDL_KEYUP: {
                auto& timestamp = key_states[event.key.keysym.sym];
                auto duration_ms = event.key.timestamp - timestamp;
                timestamp = 0u;
                auto p = puddings.find(event.key.windowID);
                if (p == puddings.end())
                {
                    break;
                }
                p->second->on_key_up(event.key.keysym.sym, duration_ms);
                break;
            }
            case SDL_MOUSEMOTION: {
                auto p = puddings.find(event.motion.windowID);
                if (p == puddings.end())
                {
                    break;
                }
                p->second->on_mouse_motion(event.motion.x, event.motion.y);
            }
            case SDL_MOUSEBUTTONDOWN: {
                auto p = puddings.find(event.button.windowID);
                if (p == puddings.end())
                {
                    break;
                }
                switch (event.button.button)
                {
                    case SDL_BUTTON_LEFT:
                        p->second->on_mouse_down_left(event.button.x, event.button.y);
                        break;
                    case SDL_BUTTON_MIDDLE:
                        p->second->on_mouse_down_middle(event.button.x, event.button.y);
                        break;
                    case SDL_BUTTON_RIGHT:
                        p->second->on_mouse_down_right(event.button.x, event.button.y);
                        break;
                    default:
                        break;
                }
                break;
            }
            case SDL_MOUSEBUTTONUP: {
                auto p = puddings.find(event.button.windowID);
                if (p == puddings.end())
                {
                    break;
                }
                switch (event.button.button)
                {
                    case SDL_BUTTON_LEFT:
                        p->second->on_mouse_up_right(event.button.x, event.button.y);
                        break;
                    case SDL_BUTTON_MIDDLE:
                        p->second->on_mouse_up_middle(event.button.x, event.button.y);
                        break;
                    case SDL_BUTTON_RIGHT:
                        p->second->on_mouse_up_right(event.button.x, event.button.y);
                        break;
                    default:
                        break;
                }
                break;
            }
            case SDL_MOUSEWHEEL: {
                auto p = puddings.find(event.wheel.windowID);
                if (p == puddings.end())
                {
                    break;
                }
                // TODO: support mouse wheel scrolling on unfocused windows
                p->second->on_mouse_wheel(event.wheel.x, event.wheel.y);
            }
            case SDL_FINGERDOWN: {
                // TODO
                break;
            }
            case SDL_FINGERUP: {
                // TODO
                break;
            }
            case SDL_CLIPBOARDUPDATE: {
                // TODO
                break;
            }
            case SDL_WINDOWEVENT: {
                auto p = puddings.find(event.window.windowID);
                if (p == puddings.end())
                {
                    break;
                }
                auto& w = *p->second;
                switch (event.window.event)
                {
                case SDL_WINDOWEVENT_SHOWN:
                    w.on_shown();
                    break;
                case SDL_WINDOWEVENT_HIDDEN:
                    w.on_hidden();
                    break;
                case SDL_WINDOWEVENT_EXPOSED:
                    w.on_exposed();
                    break;
                case SDL_WINDOWEVENT_MOVED:
                    w.on_moved(event.window.data1, event.window.data2);
                    break;
                case SDL_WINDOWEVENT_RESIZED:
                    w.on_resize(event.window.data1, event.window.data2);
                    break;
                case SDL_WINDOWEVENT_SIZE_CHANGED:
                    // TODO: handle swapchain invalidation
                    // TODO: On windows, this event should be handled on a separate thread to allow repainting during a
                    // continuous resize
                    break;
                case SDL_WINDOWEVENT_MINIMIZED:
                    w.on_minimized();
                    break;
                case SDL_WINDOWEVENT_MAXIMIZED:
                    w.on_maximized();
                    break;
                case SDL_WINDOWEVENT_RESTORED:
                    w.on_restored();
                    break;
                case SDL_WINDOWEVENT_ENTER:
                    w.on_enter();
                    break;
                case SDL_WINDOWEVENT_LEAVE:
                    w.on_enter();
                    break;
                case SDL_WINDOWEVENT_FOCUS_GAINED:
                    w.on_focus_gained();
                    break;
                case SDL_WINDOWEVENT_FOCUS_LOST:
                    w.on_focus_lost();
                    break;
                case SDL_WINDOWEVENT_CLOSE:
                    w.on_close();
                    break;
                default:
                    break;
                }
            }
            default:
                break;
            }
        }

        // Events have been handled. Time to render.
        auto& frame = vk_current_frame();
        auto index = context.frame_index % frame_context_count;

        // For now, all rendering happens on the main thread
        for (auto&& [id, pd] : puddings)
        {
            VkCommandBuffer cb = static_cast<VkCommandBuffer>(pd->cbs_[index]);
            auto acq = pd->window_.sc().acquire();
            if (acq.has_value())
            {
                VkCommandBuffer cb = static_cast<VkCommandBuffer>(pd->cbs_[index]);
                VkCommandBufferBeginInfo cmd_begin_info{};
                cmd_begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
                cmd_begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
                vkBeginCommandBuffer(cb, &cmd_begin_info);

                // TODO send swapchain render target for each pudding and set appropriate semaphore on submission
                cmd_stream cs{cb};
                pd->render(cs, acq->first);

                vkEndCommandBuffer(cb);

                cbs.emplace_back(cb);
                waits.emplace_back(acq->second);
                wait_stages.emplace_back(VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT);
                signals.emplace_back(static_cast<VkSemaphore>(pd->signals_[index]));
                pds.emplace_back(pd);

                VkSubmitInfo submit{};
                submit.sType                = VK_STRUCTURE_TYPE_SUBMIT_INFO;
                submit.waitSemaphoreCount   = 1;
                submit.pWaitSemaphores      = &waits.back();
                submit.pWaitDstStageMask    = &wait_stages.back();
                submit.commandBufferCount   = 1;
                submit.pCommandBuffers      = &cbs.back();
                submit.signalSemaphoreCount = 1;
                submit.pSignalSemaphores    = &signals.back();
                submits.emplace_back(submit);
            }
        }

        if (!submits.empty())
        {
            vkQueueSubmit(context.general.q, static_cast<uint32_t>(submits.size()), submits.data(), frame.general.fence);
            frame.general.dirty = true;
        }

        for (auto* pd : pds)
        {
            // Queue a flip for each pudding that rendered to an active surface this frame
            pd->window_.sc().present(static_cast<VkSemaphore>(pd->signals_[index]));
        }

        cbs.clear();
        pds.clear();
        signals.clear();
        submits.clear();
        waits.clear();
        wait_stages.clear();
        // Advance the frame
        advance_frame();
    }
}