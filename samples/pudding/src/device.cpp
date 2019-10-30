#include <device.hpp>

#include <fmt/core.h>

using namespace pd;

#define AMD 0x1002
#define NVIDIA 0x10DE
#define INTEL 0x8086
// Unhandled vendors:
#define IMGTEC 0x1010
#define ARM 0x13B5
#define QUALCOMM 0x5143

char const* vendor_to_string(uint32_t id)
{
    switch (id)
    {
    case AMD:
        return "AMD";
    case NVIDIA:
        return "Nvidia";
    case INTEL:
        return "Intel";

    default:
        break;
    }
    return "Unknown";
}

std::string device::to_string() const noexcept
{
    return fmt::format("{} {} ({} GB)", vendor_to_string(vendor_id), name, vram_bytes >> 30);
}

bool device::is_amd() const noexcept
{
    return vendor_id == AMD;
}

bool device::is_nvidia() const noexcept
{
    return vendor_id == NVIDIA;
}

bool device::is_intel() const noexcept
{
    return vendor_id == INTEL;
}