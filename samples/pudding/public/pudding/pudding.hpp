#pragma once

#include "window.hpp"
#include "cmd_stream.hpp"

#include <gal/pga2.hpp>
#include <tuple>
#include <vector>

namespace pd
{
class render_target;

// Inherit me! The `pudding` class encapsulates all the functionality your application needs. Each `pudding` instance
// maps to a single window, and it is up to the user to coordinate state between windows. For example, if a pudding
// window is about to close (close event received), the user can decide whether the application should shut down or not
// depending on whether the instance corresponds to the "main window" or not.
class pudding
{
public:
    // This should be invoked in your main function ONCE, at which point control is relinquished to pudding's event loop
    static void start();

    pudding();
    explicit pudding(window_desc const& wd);
    virtual ~pudding();

    // Issued once per frame after events are all removed from the event queue and handled. Applications that do not
    // animate need not request a flip. The first argument is the [command stream](cmd_stream.hpp) which is used to
    // submit render passes and commands to the gpu. The second target is the surface associated to the window attached
    // to this pudding and can be used in composing render passes within.
    virtual void render(cmd_stream& cs, render_target const& surface) {}

    // State querying
    [[nodiscard]] bool is_key_down(int key) const;
    [[nodiscard]] bool is_left_mouse_button_down() const;
    [[nodiscard]] bool is_right_mouse_button_down() const;
    [[nodiscard]] bool is_middle_mouse_button_down() const;
    [[nodiscard]] gal::pga2::point<int> mouse_xy() const;

    // Application events

    // An application quit was requested by user. The application *will* close after this handler exits. This handler
    // exists to allow time for some cleanup before exiting.
    virtual void on_quit() {}

    // Keyboard events
    virtual void on_key_down(int key) {}
    virtual void on_key_up(int key, uint32_t press_duration_ms) {}

    // Text editing events
    // TODO

    // Mouse events
    // Invoked on mouse movement when this window has focus (x and y are relative to the window)
    virtual void on_mouse_motion(int x, int y) {}
    virtual void on_mouse_down_left(int x, int y) {}
    virtual void on_mouse_down_middle(int x, int y) {}
    virtual void on_mouse_down_right(int x, int y) {}
    virtual void on_mouse_up_left(int x, int y) {}
    virtual void on_mouse_up_middle(int x, int y) {}
    virtual void on_mouse_up_right(int x, int y) {}
    virtual void on_mouse_wheel(int dx, int dy) {}

    // Window events

    // Window was shown
    virtual void on_shown() {}
    // Window was hidden
    virtual void on_hidden() {}
    // Window was exposed (you should redraw on the next opportunity)
    virtual void on_exposed() {}
    // Window was moved (arguments specify new position)
    virtual void on_moved(int x, int y) {}
    // Window was resized
    virtual void on_resize(int width, int height) {}
    // Window was minimized
    virtual void on_minimized() {}
    // Window was maximized
    virtual void on_maximized() {}
    // Window was restored
    virtual void on_restored() {}
    // User cursor has entered window
    virtual void on_enter() {}
    // User cursor has left window
    virtual void on_leave() {}
    // Window gained focus
    virtual void on_focus_gained() {}
    // Window lost focus
    virtual void on_focus_lost() {}
    // Window was closed by user
    virtual void on_close() {}

protected:

private:
    window window_;
    // Opaque handles to this pudding's primary command buffers
    std::vector<void*> cbs_;
    // Opaque handles to semaphores raised when rendering has finished and this pudding is ready for a flip
    std::vector<void*> signals_;
};
}