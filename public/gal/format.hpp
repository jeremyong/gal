#pragma once

#include "entity.hpp"

#include <string>
#include <sstream>

namespace gal
{
template <typename A, typename T, elem_t... E>
[[nodiscard]] std::string to_string(entity<A, T, E...> in)
{
    if constexpr (sizeof...(E) == 0)
    {
        return "0";
    }
    else
    {
        std::stringstream stream;

        auto const& elements = in.elements;
        for (size_t i = 0; i != elements.size(); ++i)
        {
            stream << in[i];
            auto e = elements[i];
            if (e > 0)
            {
                stream << 'e';
            }
            int j = 0;
            while (e > 0)
            {
                if ((e & 1) == 1)
                {
                    stream << j;
                }
                e >>= 1;
                ++j;
            }

            if (i != elements.size() - 1)
            {
                stream << " + ";
            }
        }
        return stream.str();
    }
}
} // namespace gal