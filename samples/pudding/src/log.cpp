#include "log.hpp"

using pd::log;

constexpr char const* sev_to_string(log_severity sev) noexcept
{
    switch (sev)
    {
    case log_severity::debug:
        return "debug";
    case log_severity::info:
        return "info";
    case log_severity::warn:
        return "warn";
    case log_severity::error:
        return "error";
    }
    return "";
}

void log::log_debug(log_severity sev, char const* file, int line, std::string const& msg)
{
    std::string data = fmt::format("[{}] {} [{}@{}]", sev_to_string(sev), msg, file, line);

    if (console_log_enabled)
    {
        fmt::print("{}\n", data);
    }

    std::scoped_lock lock{mutex_};
    buf_[cursor_++ & (buf_.size() - 1)] = std::make_pair(sev, data);
}

void log::log_release(log_severity sev, std::string const& msg)
{
    std::string data = fmt::format("[{}] {}", sev_to_string(sev), msg);
    if (console_log_enabled)
    {
        fmt::print("{}\n", data);
    }

    std::scoped_lock lock{mutex_};
    buf_[cursor_++ & (buf_.size() - 1)] = std::make_pair(sev, data);
}