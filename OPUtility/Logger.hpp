#pragma once

#include <iostream>
#include <string>
#include <chrono>
#include <thread>
#include <syncstream>

#define LOG Logger::getInstance()

/**
 * \brief Logger singleton utility class.
 */
class Logger
{
	std::ostream* info_;
	std::ostream* err_;

public:
	/**
	 * \brief Returns the only Logger instance.
	 * \return Logger instance
	 */
	static Logger& getInstance()
	{
		static Logger instance;
		return instance;
	}

	/**
	 * \brief Writes an information message with timestamp.
	 */
	void info(const std::string_view&) const;
	/**
	 * \brief Writes an error message with timestamp.
	 */
	void error(const std::string_view&) const;
	/**
	 * \brief Writes an information message without timestamp.
	 */
	void blank(const std::string_view&) const;

	/**
	 * \brief Sets the output stream for the INFO level.
	 */
	void setInfoHandler(std::ostream* stream);

	/**
	 * \brief Sets the output stream for the ERROR level.
	 */
	void setErrorHandler(std::ostream* stream);

	/**
	 * \brief Sets the output stream to default.
	 */
	void resetHandlers();

private:
	// regular constructor
	Logger() : info_(&std::cout), err_(&std::cerr) {}

public:
	Logger(const Logger&) = delete;
	void operator=(const Logger&) = delete;

	Logger(const Logger&&) = delete;
	Logger& operator=(const Logger&&) = delete;

	~Logger() = default;
};

inline void Logger::info(const std::string_view& msg) const
{
	using namespace std::chrono;
	using namespace std::this_thread;
	const auto local = zoned_time{ current_zone(), system_clock::now() };
	std::osyncstream(*info_) << "[" << local << "] [" << get_id() << "] [info] \t " << msg << std::endl;
}

inline void Logger::error(const std::string_view& errMsg) const
{
	using namespace std::chrono;
	using namespace std::this_thread;
	const auto local = zoned_time{ current_zone(), system_clock::now() };
	std::osyncstream(*err_) << "[" << local << "] [" << get_id() << "] [error] \t " << errMsg << std::endl;
}

inline void Logger::blank(const std::string_view& msg) const
{
	std::osyncstream(*info_) << msg << std::endl;
}

inline void Logger::setInfoHandler(std::ostream* stream)
{
	info_ = stream;
}

inline void Logger::setErrorHandler(std::ostream* stream)
{
	err_ = stream;
}

inline void Logger::resetHandlers()
{
	info_ = &std::cout;
	err_ = &std::cerr;
}