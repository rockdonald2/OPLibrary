// precompiled header
#include "pch.h"

// other headers
#include "Logger.h"

#include <chrono>

void Logger::info(const std::string_view& msg) const
{
	using namespace std::chrono;
	const auto local = zoned_time{ current_zone(), system_clock::now() };
	*info_ << "[" << local << "] [info] --- " << msg << std::endl;
}

void Logger::error(const std::string_view& errMsg) const
{
	using namespace std::chrono;
	const auto local = zoned_time{ current_zone(), system_clock::now() };
	*err_ << "[" << local << "] [error] --- " << errMsg << std::endl;
}

void Logger::blank(const std::string_view& msg) const
{
	*info_ << msg << std::endl;
}

void Logger::setInfoHandler(std::ostream* stream)
{
	info_ = stream;
}

void Logger::setErrorHandler(std::ostream* stream)
{
	err_ = stream;
}

void Logger::resetHandlers()
{
	info_ = &std::cout;
	err_ = &std::cerr;
}
