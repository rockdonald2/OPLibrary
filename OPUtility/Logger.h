#pragma once

#include <iostream>
#include <string>

/**
 * \brief Logger singleton utility class.
 */
class Logger
{
	std::ostream& info_;
	std::ostream& err_;

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

private:
	// regular constructor
	Logger() : info_(std::cout), err_(std::cerr) {}

public:
	Logger(const Logger&) = delete;
	void operator=(const Logger&) = delete;

	Logger(const Logger&&) = delete;
	Logger& operator=(const Logger&&) = delete;

	~Logger() = default;
};
