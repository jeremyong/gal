#pragma once

#include <cstdio>

template <typename T>
void print_type()
{
    std::printf("type: %s\n", __PRETTY_FUNCTION__);
}