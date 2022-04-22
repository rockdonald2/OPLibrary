#pragma once

#include "Problem.hpp"

namespace OPLibrary
{
	/**
	 * \brief Representation of Readers, abstract.
	 */
	template <typename T>
		requires std::floating_point<T>
	class Reader
	{
	public:
		virtual ~Reader() = default;

		/**
		 * \brief Reads in the optimization problem in order of: constraints, constraints objectives and objectives.
		 */
		virtual void readProblem(Problem<T>* problem) = 0;
		/**
		 * \brief Reads in the optimization problem in order of: constraints, constraints objectives and objectives.
		 */
		virtual void readProblem(std::shared_ptr<Problem<T>>& problem) = 0;
	};
}
