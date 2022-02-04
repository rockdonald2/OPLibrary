#pragma once

#include <exception>
#include <string>

namespace OPLibrary
{
	/**
	 * \brief Representation of Matrix issues, derived from std::exception.
	 */
	class MatrixException final : public std::exception
	{
		std::string errorMsg_;

	public:
		explicit MatrixException(std::string errorMsg) : errorMsg_(std::move(errorMsg)) {}

		const char* what() const override
		{
			return errorMsg_.c_str();
		}
	};
}