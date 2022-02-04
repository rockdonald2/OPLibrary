#pragma once

#include <exception>
#include <string>

namespace OPLibrary
{
	/**
	 * \brief Representation of Reader issues, derived from std::exception.
	 */
	class ReaderException final : public std::exception
	{
		std::string errorMsg_;

	public:
		explicit ReaderException(std::string errorMsg) : errorMsg_(std::move(errorMsg)) {}

		const char* what() const override
		{
			return errorMsg_.c_str();
		}
	};
}
