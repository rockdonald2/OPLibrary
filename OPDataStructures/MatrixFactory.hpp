#pragma once

#include "Matrix.hpp"
#include "MatrixException.hpp"
#include "DenseMatrix.hpp"

namespace OPLibrary
{
	enum class MatrixType
	{
		DENSE,
	};

	/**
	 * \brief Factory for Matrix instances.
	 */
	template <typename T>
		requires std::floating_point<T>
	class MatrixFactory final
	{
		MatrixType type_;

	public:
		/**
		 * \brief Creates a MatrixFactory for the specified type of matrix.
		 * \param type of factory
		 */
		explicit MatrixFactory(const MatrixType& type = MatrixType::DENSE) : type_(type) {}

		/**
		 * \brief Creates an empty matrix.
		 * \return Matrix
		 */
		[[nodiscard]] std::unique_ptr<Matrix<T>> createMatrix() const
		{
			switch (type_)
			{
			case MatrixType::DENSE: return std::move(std::make_unique<DenseMatrix<T>>(DenseMatrix<T>()));
			}

			throw MatrixException("Invalid matrix type given.");
		}

		/**
		 * \brief Creates an empty matrix with specified row and column numbers.
		 * \param initialRows row number
		 * \param initialCols column number
		 * \return Matrix
		 */
		[[nodiscard]] std::unique_ptr<Matrix<T>> createMatrix(const size_t& initialRows, const size_t& initialCols) const
		{
			switch (type_)
			{
			case MatrixType::DENSE: return std::move(std::make_unique<DenseMatrix<T>>(DenseMatrix<T>(initialRows, initialCols)));
			}

			throw MatrixException("Invalid matrix type given.");
		}

		/**
		 * \brief Makes a deep copy of a matrix.
		 * \param matrix to make a copy of
		 * \return Matrix
		 */
		[[nodiscard]] std::unique_ptr<Matrix<T>> createMatrix(const Matrix<T>& matrix) const
		{
			switch (type_)
			{
			case MatrixType::DENSE: return std::move(std::make_unique<DenseMatrix<T>>(DenseMatrix<T>(matrix)));
			}

			throw MatrixException("Invalid matrix type given.");
		}

		/**
		 * \brief Makes a deep copy of a matrix.
		 * \param matrix to make a copy of
		 * \return Matrix
		 */
		[[nodiscard]] std::unique_ptr<Matrix<T>> createMatrix(const Matrix<T>* matrix) const
		{
			switch (type_)
			{
			case MatrixType::DENSE: return std::move(std::make_unique<DenseMatrix<T>>(DenseMatrix<T>(*matrix)));
			}

			throw MatrixException("Invalid matrix type given.");
		}

		/**
		 * \brief Makes a deep copy of a matrix.
		 * \param matrix to make a copy of
		 * \return Matrix
		 */
		[[nodiscard]] std::unique_ptr<Matrix<T>> createMatrix(const std::unique_ptr<Matrix<T>> matrix) const
		{
			switch (type_)
			{
			case MatrixType::DENSE: return std::move(std::make_unique<DenseMatrix<T>>(DenseMatrix<T>(*matrix)));
			}

			throw MatrixException("Invalid matrix type given.");
		}

		/**
		 * \brief Makes a deep copy of a matrix.
		 * \param matrix to make a copy of
		 * \return Matrix
		 */
		[[nodiscard]] std::unique_ptr<Matrix<T>> createMatrix(const std::shared_ptr<Matrix<T>> matrix) const
		{
			switch (type_)
			{
			case MatrixType::DENSE: return std::move(std::make_unique<DenseMatrix<T>>(DenseMatrix<T>(*matrix)));
			}

			throw MatrixException("Invalid matrix type given.");
		}
	};
}
