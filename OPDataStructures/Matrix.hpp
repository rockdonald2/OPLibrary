#pragma once

#include <vector>
// Eigen-hez
#include "Eigen/Dense"

namespace OPLibrary
{
	/**
	 * \brief Available decomposition types for solve method.
	 */
	enum class DecompositionType
	{
		BDCSVD,
		JACOBISVD,
		COLPIVHOUSEHOLDER
	};

	/**
	 * \brief Matrix representation, abstract users need to provide an implementation for.
	 */
	template <typename T>
		requires std::floating_point<T>
	class Matrix
	{
	protected:
		size_t rows_;
		size_t cols_;

	public:
		Matrix() : rows_(0), cols_(0)
		{
		}

		Matrix(const size_t& rows, const size_t& columns) : rows_(rows), cols_(columns)
		{
		}

		Matrix(const Matrix<T>&) = default;
		virtual Matrix<T>& operator=(const Matrix<T>&) noexcept = default;

		Matrix(Matrix<T>&&) = default;
		virtual Matrix<T>& operator=(Matrix<T>&&) noexcept = 0;

		virtual ~Matrix() = default;

		/**
		 * \brief Returns the number of rows in the Matrix.
		 * \return size_t
		 */
		[[nodiscard]] virtual size_t getRows() const = 0;
		/**
		 * \brief Sets the number of rows in the Matrix.
		 * \param rows number of rows
		 */
		virtual void setRows(const size_t& rows) = 0;
		/**
		 * \brief Returns the number of columns in the Matrix.
		 * \return size_t
		 */
		[[nodiscard]] virtual size_t getCols() const = 0;
		/**
		 * \brief Sets the number of columns in the Matrix.
		 * \param cols number of rows
		 */
		virtual void setCols(const size_t& cols) = 0;

		/**
		 * \brief Returns the size of the Matrix which equals row x cols.
		 * \return size_t
		 */
		[[nodiscard]] virtual size_t getSize() const = 0;
		/**
		 * \brief Resizes the matrix to row times cols.
		 * \param rows number of rows
		 * \param cols number of columns
		 */
		virtual void setSize(const size_t& rows, const size_t& cols) = 0;

		/**
		 * \brief Returns the value in a specified position.
		 * \param rPos position in rows
		 * \param cPos position in columns
		 * \return value
		 */
		[[nodiscard]] virtual T get(const size_t& rPos, const size_t& cPos) const = 0;
		/**
		 * \brief Sets the value in a specified position.
		 * \param rPos position in rows
		 * \param cPos position in cols
		 * \param value new value
		 */
		virtual void set(const size_t& rPos, const size_t& cPos, const T& value) = 0;

		/**
		 * \brief Returns all the values from the Matrix.
		 * \return vector of values
		 */
		[[nodiscard]] virtual std::unique_ptr<std::vector<T>> getValues() const = 0;
		/**
		 * \brief Sets all the values in the Matrix, used for reconstruction.
		 * \param values vector of values in the Matrix
		 * \param rows number of rows
		 * \param cols number of columns
		 */
		virtual void setValues(const std::vector<T>& values, const size_t& rows, const size_t& cols) = 0;

		/**
		 * \brief Returns all the values from diagonal.
		 * \return vector of values
		 */
		[[nodiscard]] virtual std::unique_ptr<std::vector<T>> getDiagonalValues() const = 0;
		/**
		 * \brief Sets all the values along the diagonal.
		 * \param values vector of values for the diagonal
		 */
		virtual void setDiagonalValues(const std::vector<T>& values) = 0;

		/**
		 * \brief Returns a column of values.
		 * \param cPos position in columns
		 * \return vector of values
		 */
		[[nodiscard]] virtual std::unique_ptr<std::vector<T>> getColumn(const size_t& cPos) const = 0;
		/**
		 * \brief Adds a new column to the right of the Matrix.
		 * \param values vector of new values
		 */
		virtual void addColumn(const std::vector<T>& values) = 0;
		/**
		 * \brief Sets a specific column to a new set of values.
		 * \param cPos position in columns
		 * \param values vector of new values
		 */
		virtual void setColumn(const size_t& cPos, const std::vector<T>& values) = 0;
		/**
		 * \brief Removes a specific column from the Matrix.
		 * \param cPos position in columns
		 * \return the removed values
		 */
		virtual std::unique_ptr<std::vector<T>> removeColumn(const size_t& cPos) = 0;

		/**
		 * \brief Returns a row of values.
		 * \param rPos position in rows
		 * \return vector of values
		 */
		[[nodiscard]] virtual std::unique_ptr<std::vector<T>> getRow(const size_t& rPos) const = 0;
		/**
		 * \brief Adds a new row to the bottom of the Matrix.
		 * \param values vector of new values
		 */
		virtual void addRow(std::vector<T> values) = 0;
		/**
		 * \brief Sets a specific row to a new set of values.
		 * \param rPos position in rows
		 * \param values vector of new values
		 */
		virtual void setRow(const size_t& rPos, const std::vector<T>& values) = 0;
		/**
		 * \brief Removes a row from the Matrix.
		 * \param rPos position in rows
		 * \return the removed values
		 */
		virtual std::unique_ptr<std::vector<T>> removeRow(const size_t& rPos) = 0;

		/**
		 * \brief Transposes the matrix, inplace or not, if not returns the new transposed Matrix.
		 * \param inplace determines if inplace or not
		 * \return transposed Matrix
		 */
		virtual std::unique_ptr<Matrix<T>> transpose(const bool& inplace = false) = 0;

		/**
		 * \brief Inverts the matrix, inplace or not, if not returns the new inverted Matrix.
		 * \param inplace determines if inplace or not
		 * \return inverted Matrix
		 */
		virtual std::unique_ptr<Matrix<T>> inverse(const bool& inplace = false) = 0;

		/**
		 * \brief Adds a value to a specific position in the Matrix.
		 * \param rPos position in rows
		 * \param cPos position in columns
		 * \param val the value to do addition with
		 */
		virtual void addToElementByPos(const size_t& rPos, const size_t& cPos, const T& val) = 0;
		/**
		 * \brief Multiplies with a value a specific position in the Matrix.
		 * \param rPos position in rows
		 * \param cPos position in columns
		 * \param val the value to do multiplication with
		 */
		virtual void multiplyElementByPos(const size_t& rPos, const size_t& cPos, const T& val) = 0;

		/**
		 * \brief Returns the Euclidean norm of the matrix.
		 */
		[[nodiscard]] virtual T norm() const = 0;

		/**
		 * \brief Multiplies with a Matrix.
		 * \param rhs second Matrix
		 * \return New Matrix as product
		 */
		virtual std::unique_ptr<Matrix<T>> operator*(const Matrix* rhs) const = 0;
		/**
		 * \brief Multiplies with a Matrix.
		 * \param rhs second Matrix
		 * \return New Matrix as product
		 */
		virtual std::unique_ptr<Matrix<T>> operator*(const Matrix& rhs) const = 0;
		/**
		 * \brief Multiplies with a Matrix.
		 * \param rhs second Matrix
		 * \return New Matrix as product
		 */
		virtual std::unique_ptr<Matrix<T>> operator*(const std::unique_ptr<Matrix<T>>& rhs) const = 0;
		/**
		 * \brief Multiplies with a Matrix.
		 * \param rhs second Matrix
		 * \return New Matrix as product
		 */
		virtual std::unique_ptr<Matrix<T>> operator*(const std::shared_ptr<Matrix<T>>& rhs) const = 0;
		/**
		 * \brief Multiplies with a vector.
		 * \param rhs second vector
		 * \return New vector as product
		 */
		virtual std::unique_ptr<std::vector<T>> operator*(const std::vector<T>* rhs) const = 0;
		/**
		 * \brief Multiplies with a vector.
		 * \param rhs second vector
		 * \return New vector as product
		 */
		virtual std::unique_ptr<std::vector<T>> operator*(const std::vector<T>& rhs) const = 0;
		/**
		 * \brief Multiplies matrix with a scalar, componentwise.
		 * \param rhs scalar value
		 * \return New Matrix as product with scalar
		 */
		virtual std::unique_ptr<Matrix<T>> operator*(const T& rhs) const = 0;
		/**
		 * \brief Matrix multiplication with assignment.
		 * \param rhs scalar value
		 * \return New Matrix as product with scalar
		 */
		virtual Matrix<T>& operator*=(const Matrix* rhs) = 0;
		/**
		 * \brief Matrix multiplication with assignment.
		 * \param rhs scalar value
		 * \return New Matrix as product with scalar
		 */
		virtual Matrix<T>& operator*=(const Matrix& rhs) = 0;
		/**
		 * \brief Matrix multiplication with assignment.
		 * \param rhs scalar value
		 * \return New Matrix as product with scalar
		 */
		virtual Matrix<T>& operator*=(const std::unique_ptr<Matrix<T>>& rhs) = 0;
		/**
		 * \brief Matrix multiplication with assignment.
		 * \param rhs scalar value
		 * \return New Matrix as product with scalar
		 */
		virtual Matrix<T>& operator*=(const std::shared_ptr<Matrix<T>>& rhs) = 0;
		/**
		 * \brief Matrix multiplication with assignment.
		 * \param rhs scalar value
		 * \return New Matrix as product with scalar
		 */
		virtual Matrix<T>& operator*=(const std::vector<T>* rhs) = 0;
		/**
		 * \brief Matrix multiplication with assignment.
		 * \param rhs scalar value
		 * \return New Matrix as product with scalar
		 */
		virtual Matrix<T>& operator*=(const std::vector<T>& rhs) = 0;
		/**
		 * \brief Matrix multiplication with assignment.
		 * \param rhs scalar value
		 * \return New Matrix as product with scalar
		 */
		virtual Matrix<T>& operator*=(const T& rhs) = 0;

		/**
		 * \brief Divides componentwise a matrix with a scalar.
		 * \param rhs scalar value
		 * \return New Matrix as result
		 */
		virtual std::unique_ptr<Matrix<T>> operator/(const T& rhs) const = 0;
		/**
		 * \brief Divides componentwise a matrix with a scalar.
		 * \param rhs scalar value
		 * \return New Matrix as result
		 */
		virtual Matrix& operator/=(const T& rhs) = 0;

		/**
		 * \brief Adds to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as addition product
		 */
		virtual std::unique_ptr<Matrix<T>> operator+(const Matrix* rhs) const = 0;
		/**
		 * \brief Adds to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as addition product
		 */
		virtual std::unique_ptr<Matrix<T>> operator+(const Matrix& rhs) const = 0;
		/**
		 * \brief Adds to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as addition product
		 */
		virtual std::unique_ptr<Matrix<T>> operator+(const std::unique_ptr<Matrix<T>>& rhs) const = 0;
		/**
		 * \brief Adds to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as addition product
		 */
		virtual std::unique_ptr<Matrix<T>> operator+(const std::shared_ptr<Matrix<T>>& rhs) const = 0;
		/**
		 * \brief Adds to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as addition product
		 */
		virtual Matrix& operator+=(const Matrix* rhs) = 0;
		/**
		 * \brief Adds to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as addition product
		 */
		virtual Matrix& operator+=(const Matrix& rhs) = 0;
		/**
		 * \brief Adds to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as addition product
		 */
		virtual Matrix& operator+=(const std::unique_ptr<Matrix<T>>& rhs) = 0;
		/**
		 * \brief Adds to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as addition product
		 */
		virtual Matrix& operator+=(const std::shared_ptr<Matrix<T>>& rhs) = 0;

		/**
		 * \brief Subtracts to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as product
		 */
		virtual std::unique_ptr<Matrix<T>> operator-(const Matrix* rhs) const = 0;
		/**
		 * \brief Subtracts to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as product
		 */
		virtual std::unique_ptr<Matrix<T>> operator-(const Matrix& rhs) const = 0;
		/**
		 * \brief Subtracts to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as product
		 */
		virtual std::unique_ptr<Matrix<T>> operator-(const std::unique_ptr<Matrix<T>>& rhs) const = 0;
		/**
		 * \brief Subtracts to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as product
		 */
		virtual std::unique_ptr<Matrix<T>> operator-(const std::shared_ptr<Matrix<T>>& rhs) const = 0;
		/**
		 * \brief Subtracts to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as product
		 */
		virtual Matrix& operator-=(const Matrix* rhs) = 0;
		/**
		 * \brief Subtracts to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as product
		 */
		virtual Matrix& operator-=(const Matrix& rhs) = 0;
		/**
		 * \brief Subtracts to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as product
		 */
		virtual Matrix& operator-=(const std::shared_ptr<Matrix<T>>& rhs) = 0;
		/**
		 * \brief Subtracts to a Matrix.
		 * \param rhs second Matrix
		 * \return New matrix as product
		 */
		virtual Matrix& operator-=(const std::unique_ptr<Matrix<T>>& rhs) = 0;

		/**
		 * \brief Prints a Matrix instance to the output stream.
		 */
		friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix);
		/**
		 * \brief Returns the string representation of the Matrix.
		 */
		[[nodiscard]] virtual std::string toString() const = 0;

		/**
		 * \brief Exchanges two rows in the Matrix.
		 * \param lRow left row
		 * \param rRow right row
		 */
		virtual void exchangeRows(const size_t& lRow, const size_t& rRow) = 0;
		/**
		 * \brief Exchanges two columns in the Matrix.
		 * \param lCol left column
		 * \param rCol right column
		 */
		virtual void exchangeCols(const size_t& lCol, const size_t& rCol) = 0;

		/**
		 * \brief Slices and returns a block of the matrix.
		 * \param sRow starting row
		 * \param sCol starting column
		 * \param eRow ending row
		 * \param eCol ending column
		 * \return (sRow, sCol) - (eRow, eCol) block of Matrix as a new Matrix
		 */
		[[nodiscard]] virtual std::unique_ptr<Matrix<T>> block(size_t sRow, size_t sCol, size_t eRow, size_t eCol) const
			= 0;
		/**
		 * \brief Slices and exchanges a block of the matrix.
		 * \param sRow starting row
		 * \param sCol starting column
		 * \param eRow ending row
		 * \param eCol ending column
		 * \param newVals the new values of the specified block
		 */
		virtual void block(size_t sRow, size_t sCol, size_t eRow, size_t eCol, Matrix<T>* newVals) = 0;
		/**
		 * \brief Slices and exchanges a block of the matrix.
		 * \param sRow starting row
		 * \param sCol starting column
		 * \param eRow ending row
		 * \param eCol ending column
		 * \param newVals the new values of the specified block
		 */
		virtual void block(size_t sRow, size_t sCol, size_t eRow, size_t eCol, const std::unique_ptr<Matrix<T>>& newVals) = 0;
		/**
		 * \brief Slices and exchanges a block of the matrix.
		 * \param sRow starting row
		 * \param sCol starting column
		 * \param eRow ending row
		 * \param eCol ending column
		 * \param newVals the new values of the specified block
		 */
		virtual void block(size_t sRow, size_t sCol, size_t eRow, size_t eCol, const std::shared_ptr<Matrix<T>>& newVals) = 0;
		/**
		 * \brief Slices and exchanges a block of the matrix.
		 * \param sRow starting row
		 * \param sCol starting column
		 * \param eRow ending row
		 * \param eCol ending column
		 * \param newVals the new values of the specified block
		 */
		virtual void block(size_t sRow, size_t sCol, size_t eRow, size_t eCol, const std::vector<T>& newVals) = 0;
		/**
		 * \brief Fills the block with a scalar value.
		 * \param sRow starting row
		 * \param sCol starting column
		 * \param eRow ending row
		 * \param eCol ending column
		 * \param scalar the new value for this block
		 */
		virtual void block(size_t sRow, size_t sCol, size_t eRow, size_t eCol, const T& scalar) = 0;

		/**
		 * \brief Solves a system of equations with lhs as this and the given rhs, clients can decide which decomposition type is used.
		 * \param rhs the right-hand side of the equation
		 * \param decomposition type
		 * \return solution matrix/vector
		 */
		[[nodiscard]] virtual std::unique_ptr<Matrix<T>> solve(Matrix<T>& rhs,
			const DecompositionType& decomposition =
			DecompositionType::BDCSVD) = 0;
		/**
		 * \brief Solves a system of equations with lhs as this and the given rhs, clients can decide which decomposition type is used.
		 * \param rhs the right-hand side of the equation
		 * \param decomposition type
		 * \return solution matrix/vector
		 */
		[[nodiscard]] virtual std::unique_ptr<Matrix<T>> solve(Matrix<T>* rhs,
			const DecompositionType& decomposition =
			DecompositionType::BDCSVD) = 0;
		/**
		 * \brief Solves a system of equations with lhs as this and the given rhs, clients can decide which decomposition type is used.
		 * \param rhs the right-hand side of the equation
		 * \param decomposition type
		 * \return solution matrix/vector
		 */
		[[nodiscard]] virtual std::unique_ptr<Matrix<T>> solve(const std::unique_ptr<Matrix<T>>& rhs,
			const DecompositionType& decomposition =
			DecompositionType::BDCSVD) = 0;
		/**
		 * \brief Solves a system of equations with lhs as this and the given rhs, clients can decide which decomposition type is used.
		 * \param rhs the right-hand side of the equation
		 * \param decomposition type
		 * \return solution matrix/vector
		 */
		[[nodiscard]] virtual std::unique_ptr<Matrix<T>> solve(const std::shared_ptr<Matrix<T>>& rhs,
			const DecompositionType& decomposition =
			DecompositionType::BDCSVD) = 0;
	};
}
