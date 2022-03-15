#pragma once

#include <map>
#include <string>

#include "SOCPSolver.hpp"
#include "Solver.hpp"
#include "SolverException.hpp"

namespace OPLibrary
{
	enum class SolverType
	{
		SOCP
	};

	/**
	 * \brief Creates different Solver instances.
	 */
	class SolverFactory final
	{
		inline static const std::map<std::string, SolverType> MAP_STR_TO_SOLVER{ {"SOCP", SolverType::SOCP} };

	public:
		SolverFactory() = delete;

		/**
		 * \brief Creates a Solver instance. Can throw.
		 * \param type of requested Solver as string
		 * \return Solver instance
		 */
		template <typename T>
			requires std::floating_point<T>
		static std::unique_ptr<Solver<T>> createSolver(const std::string& type)
		{
			std::string temp(type);
			std::ranges::transform(temp, temp.begin(), [](const unsigned char c) { return std::toupper(c); });

			if (!MAP_STR_TO_SOLVER.contains(temp))
				throw SolverException("Unsupported solver type.");

			switch (MAP_STR_TO_SOLVER.at(temp))
			{
			case SolverType::SOCP: return std::move(std::make_unique<SOCPSolver<T>>());
			}

			throw SolverException("Unsupported solver type.");
		}
	};
}
