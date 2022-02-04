// precompiled header
#include "pch.h"

// other headers
#include "Logger.h"

#include <chrono>

void Logger::info(const std::string& msg) const
{
    using namespace std::chrono;
    const auto local = zoned_time{current_zone(), system_clock::now()};
    info_ << "[" << local << "] --- " << msg << std::endl;
}

void Logger::error(const std::string& errMsg) const
{
    using namespace std::chrono;
    const auto local = zoned_time{current_zone(), system_clock::now()};
    err_ << "[" << local << "] --- " << errMsg << std::endl;
}

void Logger::blank(const std::string& msg) const
{
    info_ << msg << std::endl;
}
