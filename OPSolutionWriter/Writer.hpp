#pragma once

#include <concepts>

#include "Problem.hpp"
#include "Solution.hpp"

namespace OPLibrary
{
	/**
	 * \brief Base representation of Writers, abstract.
	 */
	template <typename T>
		requires std::floating_point<T>
	class Writer
	{
	public:
		virtual ~Writer() = default;

		/**
		 * \brief Writes Problem in a stylized format to the output of choice.
		 */
		virtual void writeProblem(const Problem<T>* problem) = 0;
		/**
		 * \brief Writes Problem in a stylized format to the output of choice.
		 */
		virtual void writeProblem(const std::shared_ptr<Problem<T>>& problem) = 0;

		/**
		 * \brief Writes Solution in a stylized format to the output of choice.
		 */
		virtual void writeSolution(const Solution<T>* solution) = 0;
		/**
		 * \brief Writes Solution in a stylized format to the output of choice.
		 */
		virtual void writeSolution(const std::shared_ptr<Solution<T>>& solution) = 0;

		/**
		 * \brief Writes an iteration to the output of choice.
		 */
		virtual void writeIteration(const size_t& iter, std::initializer_list<T> args) = 0;

		/**
		 * \brief Returns the current output of choice, for other writers to use, e.g., Logger.
		 */
		virtual std::ostream* getOutput() const = 0;
	};
}
