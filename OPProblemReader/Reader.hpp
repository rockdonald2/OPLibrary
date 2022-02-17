#pragma once

#include "Matrix.hpp"
#include "Problem.hpp"

namespace OPLibrary
{
	/**
	 * \brief Representation of Reader.
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
		 * \brief Reads in the dimensions.
		 */
		virtual void readParams() = 0;
		/**
		 * \brief Returns the dimensions.
		 * \return dimensions pair as row-column.
		 */
		[[nodiscard]] virtual std::pair<size_t, size_t> getParams() const = 0;

		/**
		 * \brief Reads in a Matrix from the input.
		 */
		virtual void readMatrix(Matrix<T>*) = 0;
		/**
		 * \brief Reads in a Vector from the input.
		 */
		virtual void readVector(Matrix<T>*) = 0;
	};
}
