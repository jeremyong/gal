#pragma once

namespace gal
{
// An engine is used to execute the computation expressed in a multivector expression template
struct engine
{
    template <typename... T>
    [[nodiscard]] constexpr auto compute(T mv) noexcept
    {

    }
};
}