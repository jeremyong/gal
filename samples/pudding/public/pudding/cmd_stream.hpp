#pragma once

#include "render_target.hpp"

#include <vector>

namespace pd
{
class cmd_stream;

// A render stream is recorded to in the context of a render pass.
class render_stream final
{
public:
private:
    friend class cmd_stream;
    explicit render_stream(cmd_stream& parent);

    cmd_stream& parent_;
};

// A `cmd_stream` is passed to a pudding application at render time and encapsulates work submission to the GPU for that
// frame.
class cmd_stream final
{
public:
    // Render passes
    // Prior to starting a render pass, all render targets used as inputs and outputs of the pass are pushed onto an
    // internal render target stack via `push_render_target`. When `begin_render_pass` is finally invoked, the render
    // target stack is fully consumed.
    void push_render_target(render_target const& rt);

    template <typename L>
    void render_pass(L&& l)
    {
        begin_render_pass();
        l(render_stream{*this});
        end_render_pass();
    }

private:
    friend class pudding;
    friend class render_stream;

    cmd_stream(void* cb);

    void begin_render_pass();
    void end_render_pass();

    // Opaque handle pointing to a command buffer object
    void* cb_;

    std::vector<render_target const*> rts_;
};
} // namespace pd