#pragma once

#include <Logger.h>
#include <string>
#include <vector>
#include <ranges>

#include "Solution.hpp"
#include "Matrix.hpp"

namespace OPLibrary
{
	/**
	 * \brief Representation of different solution statuses.
	 */
	enum class SolutionStatus
	{
		OPTIMAL,
		NONOPTIMAL
	};

	/**
	 * \brief Representation of Solver algorithm.
	 */
	template <typename T>
	class Solver
	{
	protected:
		Matrix<T>* constraints_;
		Matrix<T>* constraintObjectives_;
		Matrix<T>* objectives_;

		std::vector<std::string> INITIALIZABLE_ARGS;
		std::map<std::string, long double*> INITIALIZATOR;

	public:
		Solver() : constraints_(nullptr), constraintObjectives_(nullptr), objectives_(nullptr) {}
		virtual ~Solver() = default;

		/**
		 * \brief Returns all the arguments which can and need to be initialized before trying to solve the optimization problem.
		 * \return vector of arguments
		 */
		virtual std::vector<std::string> getInitializableArgs()
		{
			return INITIALIZABLE_ARGS;
		}
		/**
		 * \brief Sets the value of an initializable argument for the Solver.
		 */
		virtual void setInitializableArg(const std::string& arg, const long double& val)
		{
			using namespace std;

			if (ranges::find(INITIALIZABLE_ARGS.begin(), INITIALIZABLE_ARGS.end(), arg) != INITIALIZABLE_ARGS.end())
			{
				*INITIALIZATOR[arg] = val;
			}
			else
			{
				Logger::getInstance().error(format("Tried to set an invalid initializable arg {} to value {}, skipping.", arg, val));
			}
		}

		/**
		 * \brief Sets the constraints matrix.
		 */
		virtual void setConstraints(Matrix<T>*) = 0;
		/**
		 * \brief Sets the constraint objectives vector.
		 */
		virtual void setConstraintObjectives(Matrix<T>*) = 0;
		/**
		 * \brief Sets the objectives vector.
		 */
		virtual void setObjectives(Matrix<T>*) = 0;

		/**
		 * \brief Solves the optimization problem with the set parameters.
		 * \return the status of the solution
		 */
		virtual SolutionStatus solve() = 0;

		/**
		 * \brief Returns the representation of the solution for the set problem.
		 * \return Solution
		 */
		virtual Solution getSolution() = 0;
	};
}
