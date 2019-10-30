#pragma once

#include <array>
#include <fmt/core.h>
#include <mutex>
#include <string>
#include <tuple>

constexpr inline size_t log_buf_size{2048};
constexpr inline bool console_log_enabled{true};

enum class log_severity
{
    debug,
    info,
    warn,
    error,
};

namespace pd
{
class log
{
public:
    void log_debug(log_severity sev, char const* file, int line, std::string const& msg);

    void log_release(log_severity sev, std::string const& msg);

private:
    std::mutex mutex_;
    std::array<std::pair<log_severity, std::string>, log_buf_size> buf_;
    size_t cursor_{0};
};
}

inline ::pd::log log_;

#ifdef NDEBUG
#    define DEBUG(...)
#    define INFO(...) log_.log_release(log_severity::info, fmt::format(__VA_ARGS__));
#    define WARN(...) log_.log_release(log_severity::warn, fmt::format(__VA_ARGS__));
#    define ERROR(...) log_.log_release(log_severity::error, fmt::format(__VA_ARGS__));
#else
#    define DEBUG(...) log_.log_debug(log_severity::debug, __FILE__, __LINE__, fmt::format(__VA_ARGS__));
#    define INFO(...) log_.log_debug(log_severity::info, __FILE__, __LINE__, fmt::format(__VA_ARGS__));
#    define WARN(...) log_.log_debug(log_severity::warn, __FILE__, __LINE__, fmt::format(__VA_ARGS__));
#    define ERROR(...) log_.log_debug(log_severity::error, __FILE__, __LINE__, fmt::format(__VA_ARGS__));
#endif