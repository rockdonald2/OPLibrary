#pragma once

#include <exception>
#include <string>

namespace OPLibrary
{
	/**
	 * \brief Exception class for Solver issues, derived from std::exception.
	 */
	class SolverException final : public std::exception
	{
		std::string errorMsg_;

	public:
		explicit SolverException(std::string errorMsg) : errorMsg_(std::move(errorMsg)) {}

		const char* what() const override
		{
			return errorMsg_.c_str();
		}
	};
}