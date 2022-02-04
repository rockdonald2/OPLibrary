﻿#pragma once

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
		inline static std::map<std::string, SolverType> map_Str_To_Solver_{ {"SOCP", SolverType::SOCP} };

	public:
		SolverFactory() = delete;

		/**
		 * \brief Creates a Solver instance. Can throw.
		 * \param type of requested Solver as string
		 * \return Solver instance
		 */
		template <typename T>
		static Solver<T>* createSolver(const std::string& type)
		{
			std::string temp(type);
			std::ranges::transform(temp, temp.begin(), [](const unsigned char c) { return std::tolower(c); });

			switch (map_Str_To_Solver_[temp])
			{
			case SolverType::SOCP: return new SOCPSolver<T>();
			}

			throw SolverException("Unsupported solver type.");
		}
	};
}
