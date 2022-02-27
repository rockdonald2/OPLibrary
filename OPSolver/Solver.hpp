﻿#pragma once

#include <Logger.h>
#include <string>
#include <vector>
#include <ranges>
#include <map>

#include "Solution.hpp"
#include "Problem.hpp"
#include "Writer.hpp"

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
		requires std::floating_point<T>
	class Solver
	{
	protected:
		std::vector<std::string> INITIALIZABLE_ARGS;
		std::map<std::string, long double*> INITIALIZATOR;

		std::shared_ptr<Problem<T>> problem_;
		std::shared_ptr<Writer<T>> writer_;

	public:
		Solver() : problem_(nullptr), writer_(nullptr)
		{
		}

		virtual ~Solver() = default;

		/**
		 * \brief Returns all the arguments which can and need to be initialized before trying to solve the optimization problem.
		 * \return vector of arguments
		 */
		[[nodiscard]] virtual std::vector<std::string> getInitializableArgs() const;
		/**
		 * \brief Sets the value of an initializable argument for the Solver.
		 */
		virtual void setInitializableArg(const std::string& arg, const long double& val);

		/**
		 * \brief Sets the problem for this solver.
		 * \param problem instance
		 */
		virtual void setProblem(Problem<T>* problem) = 0;
		/**
		 * \brief Sets the problem for this solver.
		 * \param problem instance
		 */
		virtual void setProblem(std::shared_ptr<Problem<T>> problem) = 0;
		/**
		 * \brief Returns the set problem for this solver.
		 * \return problem instance
		 */
		[[nodiscard]] virtual std::shared_ptr<Problem<T>> getProblem() const = 0;

		/**
		 * \brief Sets the writer for this solver.
		 */
		virtual void setWriter(Writer<T>* writer) = 0;
		/**
		 * \brief Sets the writer for this solver.
		 */
		virtual void setWriter(std::shared_ptr<Writer<T>> writer) = 0;
		/**
		 * \brief Returns the writer for this solver.
		 */
		[[nodiscard]] virtual std::shared_ptr<Writer<T>> getWriter() const = 0;

		/**
		 * \brief Solves the optimization problem with the set parameters.
		 * \return the status of the solution
		 */
		virtual SolutionStatus solve() = 0;

		/**
		 * \brief Returns the representation of the solution for the set problem.
		 * \return Solution
		 */
		[[nodiscard]] virtual Solution<T> getSolution() = 0;
	};

	template <typename T>
		requires std::floating_point<T>
	std::vector<std::string> Solver<T>::getInitializableArgs() const
	{
		return INITIALIZABLE_ARGS;
	}

	template <typename T>
		requires std::floating_point<T>
	void Solver<T>::setInitializableArg(const std::string& arg, const long double& val)
	{
		using namespace std;

		if (ranges::find(INITIALIZABLE_ARGS.begin(), INITIALIZABLE_ARGS.end(), arg) != INITIALIZABLE_ARGS.end())
		{
			*INITIALIZATOR[arg] = val;
		}
		else
		{
			LOG.error(format("Tried to set an invalid initializable arg {} to value {}, skipping.",
				arg, val));
		}
	}
}
