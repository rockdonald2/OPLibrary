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
	protected:
		std::ostream* output_;
		ProblemType type_;

	public:
		Writer() : output_(nullptr), type_(ProblemType::INVALID) {}
		explicit Writer(std::ostream* output, const ProblemType& type) : output_(output), type_(type) {}
		virtual ~Writer() = default;

		/**
		 * \brief Sets the iteration headers, the name of the values which will be given on each iteration. Important: iteration header don't have to be specified.
		 */
		virtual void setIterationHeaders(const std::vector<std::string>& headers) = 0;

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
		[[nodiscard]] virtual std::ostream* getOutput() const;
	};

	template <typename T>
		requires std::floating_point<T>
	std::ostream* Writer<T>::getOutput() const
	{
		return output_;
	}
}
