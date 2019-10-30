#pragma once

#include <cstdint>
#include <string>

namespace pd
{
struct device
{
    bool discrete;
    bool gfx_enabled;
    uint32_t vendor_id;
    uint64_t vram_bytes;
    std::string name;

    // Handles to opaque types such as this VkDevice are type-erased in the public interface
    void* handle;

    [[nodiscard]] std::string to_string() const noexcept;
    [[nodiscard]] bool is_amd() const noexcept;
    [[nodiscard]] bool is_nvidia() const noexcept;
    [[nodiscard]] bool is_intel() const noexcept;
};
} // namespace pd
