#pragma once

#include "shader.hpp"

#include <unordered_map>

class shader_cache
{
public:
private:
    std::unordered_map<std::string, shader> shaders_;
};