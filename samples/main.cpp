#include <gal/pga.hpp>
#include <pudding/pudding.hpp>
#include <pudding/renderer.hpp>

class sample : public pd::pudding
{
    void render(pd::cmd_stream& cs, pd::render_target const& surface) override
    {
        cs.push_render_target(surface);
        cs.render_pass([](auto rs) {});
    }
};

using namespace gal;
using namespace gal::pga;

int main(int argc, char** argv)
{
    pd::renderer renderer{{"gal_samples", ""}};
    renderer.init_default_device();
    sample s{};
    pd::pudding::start();
    return 0;
}