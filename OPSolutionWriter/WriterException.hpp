#pragma once

#include <exception>
#include <string>

namespace OPLibrary
{
	/**
	 * \brief Representation of Writer issues, derived from std::exception.
	 */
	class WriterException : public std::exception
	{
		std::string errorMsg_;

	public:
		explicit WriterException(std::string errorMsg) : errorMsg_(std::move(errorMsg)) {}

		const char* what() const override
		{
			return errorMsg_.c_str();
		}
	};
}
