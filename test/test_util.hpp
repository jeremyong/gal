#pragma once

#include <cstdio>
#include <formatters.hpp>
#include <doctest/doctest.h>

template <typename T>
void print_type()
{
    std::printf("type: %s\n", __PRETTY_FUNCTION__);
}

template <typename T>
void print_type(const T&)
{
    std::printf("type: %s\n", __PRETTY_FUNCTION__);
}

template <typename T>
void print(const T& in)
{
    fmt::print("{}\n", in);
}

inline auto epsilon = doctest::Approx(0.0f);