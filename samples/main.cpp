#include <cga.hpp>

using namespace gal::cga;

template <typename T>
auto compute_self_norm(T point)
{
    return point >> point;
}
// TODO