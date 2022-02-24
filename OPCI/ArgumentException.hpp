#pragma once
#include <exception>
#include <string>

namespace OPLibrary
{
	/**
	 * \brief Represents argument issues, derived from std::exception.
	 */
	class ArgumentException : public std::exception
	{
		std::string errorMsg_;

	public:
		explicit ArgumentException(std::string errorMsg) : errorMsg_(std::move(errorMsg)) {}

		const char* what() const override
		{
			return errorMsg_.c_str();
		}
	};
}
