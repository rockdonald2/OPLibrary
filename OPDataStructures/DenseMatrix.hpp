#pragma once

#include "Matrix.hpp"
#include <Eigen/Dense>
#include "MatrixException.hpp"

namespace OPLibrary
{
	/**
	 * \brief Matrix with dense, rarely missing values.
	 */
	template <typename T>
		requires std::floating_point<T>
	class DenseMatrix final : public Matrix<T>
	{
		std::unique_ptr<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> matrix_;

	public:
		explicit DenseMatrix() : Matrix<T>(0, 0), matrix_(nullptr)
		{
		}

		explicit DenseMatrix(const size_t& initialRows, const size_t& initialColumns) :
			Matrix<T>(initialRows, initialColumns),
			matrix_(new Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(initialRows, initialColumns))
		{
		}

		explicit DenseMatrix(const Matrix<T>& rhs);

		~DenseMatrix() override = default;

		DenseMatrix(const DenseMatrix<T>&) noexcept;
		DenseMatrix& operator=(const DenseMatrix<T>&) noexcept;

		DenseMatrix(DenseMatrix<T>&&) noexcept;
		DenseMatrix& operator=(DenseMatrix<T>&&) noexcept;
		explicit DenseMatrix(Matrix<T>&&) noexcept;
		DenseMatrix& operator=(Matrix<T>&&) noexcept override;

		[[nodiscard]] size_t getRows() const override;
		void setRows(const size_t& rows) override;

		[[nodiscard]] size_t getCols() const override;
		void setCols(const size_t& cols) override;

		[[nodiscard]] size_t getSize() const override;
		void setSize(const size_t& rows, const size_t& cols) override;

		[[nodiscard]] T get(const size_t& rPos, const size_t& cPos) const override;
		void set(const size_t& rPos, const size_t& cPos, const T& value) override;

		[[nodiscard]] std::unique_ptr<std::vector<T>> getValues() const override;
		void setValues(const std::vector<T>& values, const size_t& rows, const size_t& cols) override;

		[[nodiscard]] std::unique_ptr<std::vector<T>> getDiagonalValues() const override;
		void setDiagonalValues(const std::vector<T>& values) override;

		[[nodiscard]] std::unique_ptr<std::vector<T>> getColumn(const size_t& cPos) const override;

		void addColumn(const std::vector<T>& values) override;
		void setColumn(const size_t& cPos, const std::vector<T>& values) override;
		std::unique_ptr<std::vector<T>> removeColumn(const size_t& cPos) override;
		[[nodiscard]] std::unique_ptr<std::vector<T>> getRow(const size_t& rPos) const override;

		void addRow(std::vector<T> values) override;
		void setRow(const size_t& rPos, const std::vector<T>& values) override;
		std::unique_ptr<std::vector<T>> removeRow(const size_t& rPos) override;

		std::unique_ptr<Matrix<T>> transpose(const bool& inplace = false) override;
		std::unique_ptr<Matrix<T>> inverse(const bool& inplace = false) override;

		void addToElementByPos(const size_t& rPos, const size_t& cPos, const T& val) override;
		void multiplyElementByPos(const size_t& rPos, const size_t& cPos, const T& val) override;

		[[nodiscard]] T norm() const override;

		std::unique_ptr<Matrix<T>> operator*(const Matrix<T>* rhs) const override;
		std::unique_ptr<Matrix<T>> operator*(const Matrix<T>& rhs) const override;

		std::unique_ptr<Matrix<T>> operator*(const std::unique_ptr<Matrix<T>>& rhs) const override;
		std::unique_ptr<Matrix<T>> operator*(const std::shared_ptr<Matrix<T>>& rhs) const override;

		std::unique_ptr<std::vector<T>> operator*(const std::vector<T>* rhs) const override;
		std::unique_ptr<std::vector<T>> operator*(const std::vector<T>& rhs) const override;

		Matrix<T>& operator*=(const Matrix<T>* rhs) override;
		Matrix<T>& operator*=(const Matrix<T>& rhs) override;

		Matrix<T>& operator*=(const std::unique_ptr<Matrix<T>>& rhs) override;
		Matrix<T>& operator*=(const std::shared_ptr<Matrix<T>>& rhs) override;

		Matrix<T>& operator*=(const std::vector<T>* rhs) override;
		Matrix<T>& operator*=(const std::vector<T>& rhs) override;

		std::unique_ptr<Matrix<T>> operator*(const T& rhs) const override;
		Matrix<T>& operator*=(const T& rhs) override;

		std::unique_ptr<Matrix<T>> operator/(const T& rhs) const override;
		Matrix<T>& operator/=(const T& rhs) override;

		std::unique_ptr<Matrix<T>> operator+(const Matrix<T>* rhs) const override;
		std::unique_ptr<Matrix<T>> operator+(const Matrix<T>& rhs) const override;
		std::unique_ptr<Matrix<T>> operator+(const std::unique_ptr<Matrix<T>>& rhs) const override;
		std::unique_ptr<Matrix<T>> operator+(const std::shared_ptr<Matrix<T>>& rhs) const override;

		Matrix<T>& operator+=(const Matrix<T>* rhs) override;
		Matrix<T>& operator+=(const Matrix<T>& rhs) override;
		Matrix<T>& operator+=(const std::unique_ptr<Matrix<T>>& rhs) override;
		Matrix<T>& operator+=(const std::shared_ptr<Matrix<T>>& rhs) override;

		std::unique_ptr<Matrix<T>> operator-(const Matrix<T>* rhs) const override;
		std::unique_ptr<Matrix<T>> operator-(const Matrix<T>& rhs) const override;
		std::unique_ptr<Matrix<T>> operator-(const std::unique_ptr<Matrix<T>>& rhs) const override;
		std::unique_ptr<Matrix<T>> operator-(const std::shared_ptr<Matrix<T>>& rhs) const override;
		Matrix<T>& operator-=(const Matrix<T>* rhs) override;
		Matrix<T>& operator-=(const Matrix<T>& rhs) override;
		Matrix<T>& operator-=(const std::shared_ptr<Matrix<T>>& rhs) override;
		Matrix<T>& operator-=(const std::unique_ptr<Matrix<T>>& rhs) override;

		friend std::ostream& operator<<(std::ostream& out, const Matrix<T>& matrix)
		{
			if (matrix.getSize() == 0) return out;

			for (auto i = 0ULL; i < matrix.getRows(); ++i)
			{
				out << "| ";

				for (auto j = 0ULL; j < matrix.getCols(); ++j)
				{
					out << matrix.get(i, j) << " ";
				}

				out << "|\n";
			}

			return out;
		}

		[[nodiscard]] std::unique_ptr<Matrix<T>>
			block(const size_t& sRow, const size_t& sCol, const size_t& eRow, const size_t& eCol) const override;
		void block(const size_t& sRow, const size_t& sCol, const size_t& eRow, const size_t& eCol, Matrix<T>* newVals) override;
		void block(const size_t& sRow, const size_t& sCol, const size_t& eRow, const size_t& eCol, const std::unique_ptr<Matrix<T>>& newVals) override;
		void block(const size_t& sRow, const size_t& sCol, const size_t& eRow, const size_t& eCol, const std::shared_ptr<Matrix<T>>& newVals) override;
		void block(const size_t& sRow, const size_t& sCol, const size_t& eRow, const size_t& eCol, const std::vector<T>& newVals) override;
		void block(const size_t& sRow, const size_t& sCol, const size_t& eRow, const size_t& eCol, const T& scalar) override;

		[[nodiscard]] std::string toString() const override;

		void exchangeRows(const size_t& lRow, const size_t& rRow) override;
		void exchangeCols(const size_t& lCol, const size_t& rCol) override;

		[[nodiscard]] std::unique_ptr<Matrix<T>> solve(Matrix<T>& rhs, const DecompositionType& decomposition) override;
		[[nodiscard]] std::unique_ptr<Matrix<T>> solve(Matrix<T>* rhs, const DecompositionType& decomposition) override;
		[[nodiscard]] std::unique_ptr<Matrix<T>>
			solve(const std::unique_ptr<Matrix<T>>& rhs, const DecompositionType& decomposition) override;
		[[nodiscard]] std::unique_ptr<Matrix<T>>
			solve(const std::shared_ptr<Matrix<T>>& rhs, const DecompositionType& decomposition) override;
	};

	template <typename T>
		requires std::floating_point<T>
	DenseMatrix<T>::DenseMatrix(const Matrix<T>& rhs) : Matrix<T>(rhs.getRows(), rhs.getCols()),
		matrix_(new Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(
			rhs.getRows(), rhs.getCols()))
	{
		this->setValues(*rhs.getValues(), rhs.getRows(), rhs.getCols());
	}

	template <typename T>
		requires std::floating_point<T>
	DenseMatrix<T>::DenseMatrix(const DenseMatrix<T>& rhs) noexcept : Matrix<T>(rhs.getRows(), rhs.getCols()),
		matrix_(
			new Eigen::Matrix<
			T, Eigen::Dynamic, Eigen::Dynamic>(
				rhs.getRows(), rhs.getCols()))
	{
		this->setValues(rhs.getValues(), rhs.getRows(), rhs.getCols());
	}

	template <typename T>
		requires std::floating_point<T>
	DenseMatrix<T>& DenseMatrix<T>::operator=(const DenseMatrix<T>& rhs) noexcept
	{
		using namespace Eigen;

		if (&rhs == this) return *this;

		this->matrix_ = std::unique_ptr<Eigen::Matrix<T, Dynamic, Dynamic>>(
			new Eigen::Matrix<T, Dynamic, Dynamic>(rhs.getRows(), rhs.getCols()));
		this->setValues(rhs.getValues(), rhs.getRows(), rhs.getCols());

		return *this;
	}

	template <typename T>
		requires std::floating_point<T>
	DenseMatrix<T>::DenseMatrix(DenseMatrix<T>&& rhs) noexcept : Matrix<T>(rhs.getRows(), rhs.getCols()),
		matrix_(
			std::unique_ptr<Eigen::Matrix<
			T, Eigen::Dynamic, Eigen::Dynamic>>(
				rhs.matrix_.release()))
	{
		rhs.setRows(0);
		rhs.setCols(0);
	}

	template <typename T>
		requires std::floating_point<T>
	DenseMatrix<T>& DenseMatrix<T>::operator=(DenseMatrix<T>&& rhs) noexcept
	{
		using namespace Eigen;

		if (&rhs == this) return *this;

		this->matrix_ = std::unique_ptr<Eigen::Matrix<T, Dynamic, Dynamic>>(rhs.matrix_.release());

		rhs.setRows(0);
		rhs.setCols(0);

		return *this;
	}

	template <typename T>
		requires std::floating_point<T>
	DenseMatrix<T>::DenseMatrix(Matrix<T>&& rhs) noexcept : Matrix<T>(rhs.getRows(), rhs.getCols()),
		matrix_(
			std::unique_ptr<Eigen::Matrix<
			T, Eigen::Dynamic, Eigen::Dynamic>>(
				dynamic_cast<DenseMatrix<T>&&>(rhs).matrix_.
				release()))
	{
		rhs.setRows(0);
		rhs.setCols(0);
	}

	template <typename T>
		requires std::floating_point<T>
	DenseMatrix<T>& DenseMatrix<T>::operator=(Matrix<T>&& rhs) noexcept
	{
		using namespace Eigen;

		if (&rhs == this) return *this;

		this->matrix_ = std::unique_ptr<Eigen::Matrix<T, Dynamic, Dynamic>>(
			dynamic_cast<DenseMatrix<T>&&>(rhs).matrix_.release());

		rhs.setRows(0);
		rhs.setCols(0);

		return *this;
	}

	template <typename T>
		requires std::floating_point<T>
	size_t DenseMatrix<T>::getRows() const
	{
		return this->rows_;
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::setRows(const size_t& rows)
	{
		this->rows_ = rows;

		if (this->matrix_ == nullptr) this->matrix_ = std::make_unique<Eigen::Matrix<
			T, Eigen::Dynamic, Eigen::Dynamic>>(1, 1);

		this->setSize(this->rows_, this->cols_);
	}

	template <typename T>
		requires std::floating_point<T>
	size_t DenseMatrix<T>::getCols() const
	{
		return this->cols_;
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::setCols(const size_t& cols)
	{
		this->cols_ = cols;

		if (this->matrix_ == nullptr) this->matrix_ = std::make_unique<Eigen::Matrix<
			T, Eigen::Dynamic, Eigen::Dynamic>>(1, 1);

		this->setSize(this->rows_, this->cols_);
	}

	template <typename T>
		requires std::floating_point<T>
	size_t DenseMatrix<T>::getSize() const
	{
		return this->rows_ * this->cols_;
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::setSize(const size_t& rows, const size_t& cols)
	{
		this->rows_ = rows;
		this->cols_ = cols;

		if (this->matrix_ == nullptr) this->matrix_ = std::make_unique<Eigen::Matrix<
			T, Eigen::Dynamic, Eigen::Dynamic>>();

		this->matrix_->conservativeResize(this->rows_, this->cols_);
	}

	template <typename T>
		requires std::floating_point<T>
	T DenseMatrix<T>::get(const size_t& rPos, const size_t& cPos) const
	{
		if (rPos >= this->rows_) throw MatrixException("Invalid row position queried: " + std::to_string(rPos));
		if (cPos >= this->cols_) throw MatrixException("Invalid column position queried: " + std::to_string(cPos));

		return this->matrix_->coeff(rPos, cPos);
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::set(const size_t& rPos, const size_t& cPos, const T& value)
	{
		if (rPos >= this->rows_) throw MatrixException(
			"Invalid row position was tried to be set: " + std::to_string(rPos));
		if (cPos >= this->cols_) throw MatrixException(
			"Invalid column position was tried to be set: " + std::to_string(cPos));

		this->matrix_->coeffRef(rPos, cPos) = value;
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<std::vector<T>> DenseMatrix<T>::getValues() const
	{
		using namespace Eigen;
		using namespace std;

		if (this->matrix_ == nullptr) throw MatrixException("Tried to retrieve a null matrix.");

		unique_ptr<vector<T>> retVector(make_unique<vector<T>>(vector<T>(this->matrix_->size())));
		Map<Eigen::Matrix<T, Dynamic, Dynamic>>(retVector->data(), this->getRows(), this->getCols()) = *this->matrix_;
		return move(retVector);
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::setValues(const std::vector<T>& values, const size_t& rows, const size_t& cols)
	{
		using namespace Eigen;

		this->rows_ = rows;
		this->cols_ = cols;

		auto tmpVals(values);
		this->matrix_ = std::unique_ptr<Eigen::Matrix<T, Dynamic, Dynamic>>(
			new Eigen::Matrix<T, Dynamic, Dynamic>(
				Map<Eigen::Matrix<T, Dynamic, Dynamic>>(tmpVals.data(), rows, cols)));
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<std::vector<T>> DenseMatrix<T>::getDiagonalValues() const
	{
		using namespace Eigen;
		using namespace std;

		if (this->matrix_ == nullptr) throw MatrixException("Tried to retrieve a null matrix.");

		const auto retrievedDiagonal = this->matrix_->diagonal();
		unique_ptr<vector<T>> retVector(make_unique<vector<T>>(vector<T>(retrievedDiagonal.size())));
		Map<Eigen::Matrix<T, Dynamic, 1>>(retVector->data(), retrievedDiagonal.size()) = retrievedDiagonal;

		return move(retVector);
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::setDiagonalValues(const std::vector<T>& values)
	{
		if (this->matrix_ == nullptr) this->matrix_ = std::make_unique<Eigen::Matrix<
			T, Eigen::Dynamic, Eigen::Dynamic>>(values.size(), 1);

		auto diagonal(this->matrix_->diagonal());

		if (static_cast<long long>(values.size()) > diagonal.size()) throw MatrixException(
			"Sent values length for diagonal is bigger than Matrix diagonal.");

		for (size_t i(0); const auto & elem : values)
		{
			diagonal[i++] = elem;
		}
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<std::vector<T>> DenseMatrix<T>::getColumn(const size_t& cPos) const
	{
		using namespace Eigen;
		using namespace std;

		if (cPos >= this->cols_) throw MatrixException("Invalid column position queried: " + std::to_string(cPos));

		const auto retrievedCol = Eigen::Matrix<T, Dynamic, 1>(this->matrix_->col(cPos).matrix());

		unique_ptr<vector<T>> retVector(make_unique<vector<T>>(vector<T>(retrievedCol.size())));

		Map<Eigen::Matrix<T, Dynamic, 1>>(retVector->data(), retrievedCol.size()) = retrievedCol;

		return move(retVector);
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::addColumn(const std::vector<T>& values)
	{
		using namespace Eigen;

		if (this->matrix_ == nullptr) this->matrix_ = std::make_unique<Eigen::Matrix<
			T, Dynamic, Dynamic>>(1, 1);
		if (values.size() > this->cols_) throw MatrixException(
			"Sent values length is bigger than Matrix columns.");

		auto tmpValues(values);

		auto& instance = this->matrix_;
		const Map<Eigen::Matrix<T, Dynamic, 1>> newVals(tmpValues.data(), tmpValues.size());

		instance->conservativeResize(instance->rows(), instance->cols() + 1);
		instance->col(instance->cols() - 1) = newVals;
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::setColumn(const size_t& cPos, const std::vector<T>& values)
	{
		using namespace Eigen;

		if (cPos >= this->cols_) throw MatrixException(
			"Invalid column position was tried to be set: " + std::to_string(cPos));
		if (values.size() > this->cols_) throw MatrixException(
			"Sent values length is bigger than Matrix columns.");

		auto tmpValues(values);

		const Map<Eigen::Matrix<T, Dynamic, 1>> newVals(tmpValues.data(), tmpValues.size());

		this->matrix_->block(0, cPos, tmpValues.size(), 1) << newVals;
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<std::vector<T>> DenseMatrix<T>::removeColumn(const size_t& cPos)
	{
		using namespace Eigen;
		using namespace std;

		if (cPos >= this->cols_) throw MatrixException(
			"Invalid column position was tried to be removed: " + std::to_string(cPos));

		auto& instance = this->matrix_;

		Eigen::Matrix<T, Dynamic, 1> removedVals;

		if (static_cast<long long>(cPos) < instance->cols())
		{
			removedVals = Eigen::Matrix<T, Dynamic, 1>(instance->col(cPos).matrix());
			instance->block(0, cPos, instance->rows(), instance->cols() - cPos - 1) = instance->rightCols(
				instance->cols() - cPos - 1);
		}

		instance->conservativeResize(instance->rows(), instance->cols() - 1);

		unique_ptr<vector<T>> retVector(make_unique<vector<T>>(vector<T>(removedVals.size())));
		Map<Eigen::Matrix<T, Dynamic, 1>>(retVector->data(), removedVals.size()) = removedVals;

		return move(retVector);
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<std::vector<T>> DenseMatrix<T>::getRow(const size_t& rPos) const
	{
		using namespace Eigen;
		using namespace std;

		if (rPos >= this->rows_) throw MatrixException("Invalid row position queried: " + std::to_string(rPos));

		const auto retrievedRow = Eigen::Matrix<T, 1, Dynamic>(this->matrix_->row(rPos).matrix());

		unique_ptr<vector<T>> retVector(make_unique<vector<T>>(vector<T>(retrievedRow.size())));

		Map<Eigen::Matrix<T, 1, Dynamic>>(retVector->data(), retrievedRow.size()) = retrievedRow;

		return move(retVector);
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::addRow(std::vector<T> values)
	{
		using namespace Eigen;

		if (this->matrix_ == nullptr) this->matrix_ = std::make_unique<Eigen::Matrix<
			T, Eigen::Dynamic, Eigen::Dynamic>>(1, 1);
		if (values.size() > this->rows_) throw MatrixException(
			"Sent values length is bigger than Matrix rows.");

		auto& instance = this->matrix_;
		const Map<Eigen::Matrix<T, 1, Dynamic>> newVals(values.data(), values.size());

		instance->conservativeResize(instance->rows() + 1, instance->cols());
		instance->row(instance->rows() - 1) = newVals;
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::setRow(const size_t& rPos, const std::vector<T>& values)
	{
		using namespace Eigen;

		if (rPos >= this->rows_) throw MatrixException(
			"Invalid row position was tried to be set: " + std::to_string(rPos));
		if (values.size() > this->rows_) throw MatrixException(
			"Sent values length is bigger than Matrix rows.");

		std::vector<T> tmpValues(values);

		const Map<Eigen::Matrix<T, Dynamic, 1>> newVals(tmpValues.data(), tmpValues.size());

		this->matrix_->block(rPos, 0, 1, tmpValues.size()) << newVals;
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<std::vector<T>> DenseMatrix<T>::removeRow(const size_t& rPos)
	{
		using namespace Eigen;
		using namespace std;

		if (rPos >= this->rows_) throw MatrixException(
			"Invalid row position was tried to be removed: " + std::to_string(rPos));

		auto& instance = this->matrix_;

		Eigen::Matrix<T, 1, Dynamic> removedVals;

		if (static_cast<long long>(rPos) < instance->rows())
		{
			removedVals = Eigen::Matrix<T, 1, Dynamic>(instance->row(rPos).matrix());
			instance->block(rPos, 0, instance->rows() - rPos - 1, instance->cols()) = instance->bottomRows(
				instance->rows() - rPos - 1);
		}

		instance->conservativeResize(instance->rows() - 1, instance->cols());

		unique_ptr<vector<T>> retVector(make_unique<vector<T>>(vector<T>(removedVals.size())));
		Map<Eigen::Matrix<T, Dynamic, 1>>(retVector->data(), removedVals.size()) = removedVals;

		return move(retVector);
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::transpose(const bool& inplace)
	{
		using namespace std;

		if (this->matrix_ == nullptr) throw MatrixException("Tried to transpose a null matrix.");

		if (inplace)
		{
			this->matrix_->transposeInPlace();
			this->rows_ = this->matrix_->rows();
			this->cols_ = this->matrix_->cols();
			return nullptr;
		}

		const Matrix<T>* upcast = this;
		unique_ptr<Matrix<T>> retVal(make_unique<DenseMatrix<T>>(DenseMatrix(*upcast)));
		retVal->transpose(true);
		return move(retVal);
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::inverse(const bool& inplace)
	{
		using namespace std;

		if (this->matrix_ == nullptr) throw MatrixException("Tried to inverse a null matrix.");

		if (inplace)
		{
			try
			{
				*this->matrix_ = this->matrix_->inverse();
			}
			catch (...)
			{
				throw MatrixException("Matrix not invertible.");
			}

			return nullptr;
		}

		const Matrix<T>* upcast = this;
		unique_ptr<Matrix<T>> retVal(make_unique<DenseMatrix<T>>(DenseMatrix(*upcast)));
		retVal->inverse(true);
		return move(retVal);
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::addToElementByPos(const size_t& rPos, const size_t& cPos, const T& val)
	{
		if (rPos >= this->rows_) throw MatrixException(
			"Invalid row position was tried to be set: " + std::to_string(rPos));
		if (cPos >= this->cols_) throw MatrixException(
			"Invalid column position was tried to be set: " + std::to_string(cPos));

		this->matrix_->coeffRef(rPos, cPos) += val;
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::multiplyElementByPos(const size_t& rPos, const size_t& cPos, const T& val)
	{
		if (rPos >= this->rows_) throw MatrixException(
			"Invalid row position was tried to be set: " + std::to_string(rPos));
		if (cPos >= this->cols_) throw MatrixException(
			"Invalid column position was tried to be set: " + std::to_string(cPos));

		this->matrix_->coeffRef(rPos, cPos) *= val;
	}

	template <typename T>
		requires std::floating_point<T>
	T DenseMatrix<T>::norm() const
	{
		if (this->matrix_ == nullptr) throw MatrixException("Tried to calculate norm for a null matrix.");

		return this->matrix_->norm();
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator*(const Matrix<T>* rhs) const
	{
		using namespace Eigen;
		using namespace std;

		if (this->cols_ != rhs->getRows()) throw MatrixException(
			"Matrices are not multipliable: l:" + std::to_string(this->cols_) + "- r:" + std::to_string(
				rhs->getRows()));

		const auto rows(this->getRows());
		const auto cols(rhs->getCols());

		auto retMatrix = make_unique<DenseMatrix<T>>(DenseMatrix(rows, cols));

		auto rhsVals(rhs->getValues());
		const auto rhsRows(rhs->getRows());
		const auto rhsCols(rhs->getCols());

		const Eigen::Matrix<T, Dynamic, Dynamic> rhsMatrix = Map<Eigen::Matrix<T, Dynamic, Dynamic>>(
			rhsVals->data(), rhsRows, rhsCols);

		const auto tmpResult = Eigen::Matrix<T, Dynamic, Dynamic>(*this->matrix_ * rhsMatrix);
		unique_ptr<vector<T>> tmpResultVec(make_unique<vector<T>>(vector<T>(tmpResult.size())));
		Eigen::Map<Eigen::Matrix<T, Dynamic, Dynamic>>(tmpResultVec->data(), tmpResult.rows(), tmpResult.cols()) =
			tmpResult;

		retMatrix->setValues(*tmpResultVec, rows, cols);

		return move(retMatrix);
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator*(const Matrix<T>& rhs) const
	{
		return std::move(*this * (&rhs));
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator*(const std::unique_ptr<Matrix<T>>& rhs) const
	{
		return std::move(*this * rhs.get());
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator*(const std::shared_ptr<Matrix<T>>& rhs) const
	{
		return std::move(*this * rhs.get());
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<std::vector<T>> DenseMatrix<T>::operator*(const std::vector<T>* rhs) const
	{
		using namespace Eigen;
		using namespace std;

		if (this->cols_ != rhs->size()) throw MatrixException(
			"Matrices are not multipliable: l:" + std::to_string(this->cols_) + "- r:" + std::to_string(rhs->size()));

		auto tmpRhs(*rhs);
		const Map<Eigen::Matrix<T, Dynamic, 1>> rhsVector(tmpRhs.data(), tmpRhs.size());

		const auto product = Eigen::Matrix<T, Dynamic, 1>(*this->matrix_ * rhsVector);

		unique_ptr<vector<T>> retVector(make_unique<vector<T>>(vector<T>(product.size())));
		Map<Eigen::Matrix<T, Dynamic, 1>>(retVector->data(), product.size()) = product;

		return move(retVector);
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<std::vector<T>> DenseMatrix<T>::operator*(const std::vector<T>& rhs) const
	{
		return std::move(*this * (&rhs));
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator*=(const Matrix<T>* rhs)
	{
		using namespace Eigen;

		if (this->cols_ != rhs->getRows()) throw MatrixException(
			"Matrices are not multipliable: l:" + std::to_string(this->cols_) + "- r:" + std::to_string(
				rhs->getRows()));

		auto rhsVals(rhs->getValues());
		const auto rhsRows(rhs->getRows());
		const auto rhsCols(rhs->getCols());

		const Eigen::Matrix<T, Dynamic, Dynamic> rhsMatrix = Map<Eigen::Matrix<T, Dynamic, Dynamic>>(
			rhsVals->data(), rhsRows, rhsCols);

		*this->matrix_ = *this->matrix_ * rhsMatrix;
		this->setSize(this->getRows(), rhs->getCols());

		return *this;
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator*=(const Matrix<T>& rhs)
	{
		return *this *= &rhs;
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator*=(const std::unique_ptr<Matrix<T>>& rhs)
	{
		return *this *= *rhs;
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator*=(const std::shared_ptr<Matrix<T>>& rhs)
	{
		return *this *= *rhs;
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator*=(const std::vector<T>* rhs)
	{
		using namespace Eigen;

		if (this->cols_ != rhs->size()) throw MatrixException(
			"Matrices are not multipliable: l:" + std::to_string(this->cols_) + "- r:" + std::to_string(rhs->size()));

		auto tmpRhs(*rhs);
		const Map<Eigen::Matrix<T, Dynamic, 1>> rhsVector(tmpRhs.data(), tmpRhs.size());

		*this->matrix_ = *this->matrix_ * rhsVector;
		this->setSize(rhs->size(), 1);

		return *this;
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator*=(const std::vector<T>& rhs)
	{
		return *this *= &rhs;
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator*(const T& rhs) const
	{
		using namespace Eigen;
		using namespace std;

		auto retMatrix(make_unique<DenseMatrix<T>>(DenseMatrix(this->rows_, this->cols_)));

		const auto tmpResult = Eigen::Matrix<T, Dynamic, Dynamic>(*this->matrix_ * rhs);
		unique_ptr<vector<T>> tmpResultVec(make_unique<vector<T>>(vector<T>(tmpResult.size())));
		Map<Eigen::Matrix<T, Dynamic, Dynamic>>(tmpResultVec->data(), tmpResult.rows(), tmpResult.cols()) = tmpResult;

		retMatrix->setValues(*tmpResultVec, tmpResult.rows(), tmpResult.cols());

		return move(retMatrix);
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator*=(const T& rhs)
	{
		*this->matrix_ = rhs * *this->matrix_;
		return *this;
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator/(const T& rhs) const
	{
		using namespace std;

		assert(this->cols_ == 1 && "Componentwise division valid for vector only.");

		auto retMatrix(make_unique<DenseMatrix<T>>(DenseMatrix(this->rows_, 1)));

		for (size_t i = 0; i < this->rows_; ++i)
		{
			retMatrix->set(i, 0, this->get(i, 0) / rhs);
		}

		return move(retMatrix);
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator/=(const T& rhs)
	{
		using namespace std;

		assert(this->cols_ == 1 && "Componentwise division valid for vector only.");

		for (size_t i = 0; i < this->rows_; ++i)
		{
			this->set(i, 0, this->get(i, 0) / rhs);
		}

		return *this;
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator+(const Matrix<T>* rhs) const
	{
		using namespace Eigen;
		using namespace std;

		if (this->cols_ != rhs->getCols() || this->rows_ != rhs->getRows())
			throw MatrixException(
				"Matrices cannot be added together: l: " + std::to_string(this->rows_) + "x" +
				std::to_string(this->cols_) + "; r: " + std::to_string(rhs->getRows()) + "x" + std::to_string(
					rhs->getCols()));

		auto retMatrix(make_unique<DenseMatrix<T>>(DenseMatrix(this->rows_, this->cols_)));

		auto& lhsMatrix(this->matrix_);

		auto rhsVals(rhs->getValues());
		const auto rhsRows(rhs->getRows());
		const auto rhsCols(rhs->getCols());

		const Eigen::Matrix<T, Dynamic, Dynamic> rhsMatrix = Map<Eigen::Matrix<T, Dynamic, Dynamic>>(
			rhsVals->data(), rhsRows, rhsCols);

		const auto tmpResult = Eigen::Matrix<T, Dynamic, Dynamic>(*lhsMatrix + rhsMatrix);
		unique_ptr<vector<T>> tmpResultVec(make_unique<vector<T>>(vector<T>(tmpResult.size())));
		Map<Eigen::Matrix<T, Dynamic, Dynamic>>(tmpResultVec->data(), tmpResult.rows(), tmpResult.cols()) = tmpResult;

		retMatrix->setValues(*tmpResultVec, rhs->getRows(), rhs->getCols());

		return move(retMatrix);
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator+(const Matrix<T>& rhs) const
	{
		return std::move(*this + (&rhs));
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator+(const std::unique_ptr<Matrix<T>>& rhs) const
	{
		return std::move(*this + rhs.get());
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator+(const std::shared_ptr<Matrix<T>>& rhs) const
	{
		return std::move(*this - rhs.get());
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator+=(const Matrix<T>* rhs)
	{
		using namespace Eigen;

		if (this->cols_ != rhs->getCols() || this->rows_ != rhs->getRows())
			throw MatrixException(
				"Matrices cannot be added together: l: " + std::to_string(this->rows_) + "x" +
				std::to_string(this->cols_) + "; r: " + std::to_string(rhs->getRows()) + "x" + std::to_string(
					rhs->getCols()));

		auto rhsVals(rhs->getValues());
		const auto rhsRows(rhs->getRows());
		const auto rhsCols(rhs->getCols());

		const Eigen::Matrix<T, Dynamic, Dynamic> rhsMatrix = Map<Eigen::Matrix<T, Dynamic, Dynamic>>(
			rhsVals->data(), rhsRows, rhsCols);

		*this->matrix_ = *this->matrix_ + rhsMatrix;

		return *this;
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator+=(const Matrix<T>& rhs)
	{
		return *this += &rhs;
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator+=(const std::unique_ptr<Matrix<T>>& rhs)
	{
		return *this += *rhs;
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator+=(const std::shared_ptr<Matrix<T>>& rhs)
	{
		return *this += *rhs;
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator-(const Matrix<T>* rhs) const
	{
		using namespace Eigen;
		using namespace std;

		if (this->cols_ != rhs->getCols() || this->rows_ != rhs->getRows())
			throw MatrixException(
				"Matrices are not substractable: l: " + std::to_string(this->rows_) + "x" + std::to_string(this->cols_)
				+ "; r: " + std::to_string(rhs->getRows()) + "x" + std::to_string(rhs->getCols()));

		auto retMatrix(make_unique<DenseMatrix<T>>(DenseMatrix(rhs->getRows(), rhs->getCols())));

		auto& lhsMatrix(this->matrix_);

		auto rhsVals(rhs->getValues());
		const auto rhsRows(rhs->getRows());
		const auto rhsCols(rhs->getCols());

		const Eigen::Matrix<T, Dynamic, Dynamic> rhsMatrix = Map<Eigen::Matrix<T, Dynamic, Dynamic>>(
			rhsVals->data(), rhsRows, rhsCols);

		const auto tmpResult = Eigen::Matrix<T, Dynamic, Dynamic>(*lhsMatrix - rhsMatrix);
		unique_ptr<vector<T>> tmpResultVec(make_unique<vector<T>>(vector<T>(tmpResult.size())));
		Map<Eigen::Matrix<T, Dynamic, Dynamic>>(tmpResultVec->data(), tmpResult.rows(), tmpResult.cols()) = tmpResult;

		retMatrix->setValues(*tmpResultVec, rhs->getRows(), rhs->getCols());

		return move(retMatrix);
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator-(const Matrix<T>& rhs) const
	{
		return std::move(*this - (&rhs));
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator-(const std::unique_ptr<Matrix<T>>& rhs) const
	{
		return std::move(*this - rhs.get());
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::operator-(const std::shared_ptr<Matrix<T>>& rhs) const
	{
		return std::move(*this - rhs.get());
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator-=(const Matrix<T>* rhs)
	{
		using namespace Eigen;

		if (this->cols_ != rhs->getCols() || this->rows_ != rhs->getRows())
			throw MatrixException(
				"Matrices are not substractable: l: " + std::to_string(this->rows_) + "x" + std::to_string(this->cols_)
				+ "; r: " + std::to_string(rhs->getRows()) + "x" + std::to_string(rhs->getCols()));

		auto rhsVals(rhs->getValues());
		const auto rhsRows(rhs->getRows());
		const auto rhsCols(rhs->getCols());

		const Eigen::Matrix<T, Dynamic, Dynamic> rhsMatrix = Map<Eigen::Matrix<T, Dynamic, Dynamic>>(
			rhsVals->data(), rhsRows, rhsCols);

		*this->matrix_ = *this->matrix_ - rhsMatrix;

		return *this;
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator-=(const Matrix<T>& rhs)
	{
		return *this -= &rhs;
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator-=(const std::unique_ptr<Matrix<T>>& rhs)
	{
		return *this -= *rhs;
	}

	template <typename T>
		requires std::floating_point<T>
	Matrix<T>& DenseMatrix<T>::operator-=(const std::shared_ptr<Matrix<T>>& rhs)
	{
		return *this -= *rhs;
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::block(const size_t& sRow, const size_t& sCol, const size_t& eRow, const size_t& eCol) const
	{
		using namespace Eigen;
		using namespace std;

		if (sRow < 0 || sRow >= this->rows_) throw MatrixException("Out of bounds exception in block for rows.");
		if (eRow < 0 || eRow >= this->rows_) throw MatrixException("Out of bounds exception in block for rows.");
		if (sRow > eRow) throw MatrixException("In block operation the start row cannot be bigger than the end row.");
		if (sCol < 0 || sCol >= this->cols_) throw MatrixException("Out of bounds exception in block for columns.");
		if (eCol < 0 || eCol >= this->cols_) throw MatrixException("Out of bounds exception in block for columns.");
		if (sCol > eCol) throw MatrixException(
			"In block operation the start column cannot be bigger than the end column.");

		const size_t nRows(eRow - sRow + 1);
		const size_t nCols(eCol - sCol + 1);

		auto retMatrix(make_unique<DenseMatrix<T>>(DenseMatrix(nRows, nCols)));

		const auto tmpMatrix(Eigen::Matrix<T, Dynamic, Dynamic>(this->matrix_->block(sRow, sCol, nRows, nCols)));
		unique_ptr<vector<T>> tmpVals(make_unique<vector<T>>(vector<T>(tmpMatrix.size())));
		Map<Eigen::Matrix<T, Dynamic, Dynamic>>(tmpVals->data(), tmpMatrix.rows(), tmpMatrix.cols()) = tmpMatrix;

		retMatrix->setValues(*tmpVals, nRows, nCols);

		return move(retMatrix);
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::block(const size_t& sRow, const size_t& sCol, const size_t& eRow, const size_t& eCol, Matrix<T>* newVals)
	{
		using namespace Eigen;

		if (sRow < 0 || sRow >= this->rows_) throw MatrixException("Out of bounds exception in block for rows.");
		if (eRow < 0 || eRow >= this->rows_) throw MatrixException("Out of bounds exception in block for rows.");
		if (sRow > eRow) throw MatrixException("In block operation the start row cannot be bigger than the end row.");
		if (sCol < 0 || sCol >= this->cols_) throw MatrixException("Out of bounds exception in block for columns.");
		if (eCol < 0 || eCol >= this->cols_) throw MatrixException("Out of bounds exception in block for columns.");
		if (sCol > eCol) throw MatrixException(
			"In block operation the start column cannot be bigger than the end column.");

		const size_t nRows(eRow - sRow + 1);
		const size_t nCols(eCol - sCol + 1);

		if (newVals->getSize() != (nRows * nCols)) throw MatrixException(
			"The supplied new value matrix's does not apply to the block.");

		auto nVals(newVals->getValues());

		const Eigen::Matrix<T, Dynamic, Dynamic> bMatrix = Map<Eigen::Matrix<T, Dynamic, Dynamic>>(
			nVals->data(), nRows, nCols);

		this->matrix_->block(sRow, sCol, nRows, nCols) = bMatrix;
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::block(const size_t& sRow, const size_t& sCol, const size_t& eRow, const size_t& eCol, const std::unique_ptr<Matrix<T>>& newVals)
	{
		this->block(sRow, sCol, eRow, eCol, newVals.get());
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::block(const size_t& sRow, const size_t& sCol, const size_t& eRow, const size_t& eCol, const std::shared_ptr<Matrix<T>>& newVals)
	{
		this->block(sRow, sCol, eRow, eCol, newVals.get());
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::block(const size_t& sRow, const size_t& sCol, const size_t& eRow, const size_t& eCol, const std::vector<T>& newVals)
	{
		const size_t nRows(eRow - sRow + 1);
		const size_t nCols(eCol - sCol + 1);

		auto tmpMatrix(std::make_unique<DenseMatrix<T>>(DenseMatrix<T>(nRows, nCols)));
		tmpMatrix->setValues(newVals, nRows, nCols);

		this->block(sRow, sCol, eRow, eCol, tmpMatrix.get());
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::block(const size_t& sRow, const size_t& sCol, const size_t& eRow, const size_t& eCol, const T& scalar)
	{
		using namespace std;

		const size_t nRows(eRow - sRow + 1);
		const size_t nCols(eCol - sCol + 1);

		auto tmpMatrix(std::make_unique<DenseMatrix<T>>(DenseMatrix<T>(nRows, nCols)));
		tmpMatrix->setValues(vector<T>(nRows * nCols, scalar), nRows, nCols);

		this->block(sRow, sCol, eRow, eCol, tmpMatrix.get());
	}

	template <typename T>
		requires std::floating_point<T>
	std::string DenseMatrix<T>::toString() const
	{
		if (this->matrix_ == nullptr) throw MatrixException("Tried to print a null matrix.");

		std::stringstream output;
		output << matrix_;
		return output.str();
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::exchangeRows(const size_t& lRow, const size_t& rRow)
	{
		if (lRow >= this->getRows() || rRow >= this->getRows()) throw MatrixException(
			"Row position out of bounds for exchanging rows: " + std::to_string(lRow) + "-" + std::to_string(rRow));

		this->matrix_->row(lRow).swap(this->matrix_->row(rRow));
	}

	template <typename T>
		requires std::floating_point<T>
	void DenseMatrix<T>::exchangeCols(const size_t& lCol, const size_t& rCol)
	{
		if (lCol >= this->getCols() || rCol >= this->getCols()) throw MatrixException(
			"Row position out of bounds for exchanging rows: " + std::to_string(lCol) + "-" + std::to_string(rCol));

		this->matrix_->col(lCol).swap(this->matrix_->col(rCol));
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::solve(Matrix<T>* rhs, const DecompositionType& decomposition)
	{
		return std::move(this->solve(*rhs, decomposition));
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::solve(const std::unique_ptr<Matrix<T>>& rhs,
		const DecompositionType& decomposition)
	{
		return std::move(this->solve(*rhs, decomposition));
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::solve(const std::shared_ptr<Matrix<T>>& rhs,
		const DecompositionType& decomposition)
	{
		return std::move(this->solve(*rhs, decomposition));
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> DenseMatrix<T>::solve(Matrix<T>& rhs, const DecompositionType& decomposition)
	{
		using namespace Eigen;
		using namespace std;

		if (this->matrix_ == nullptr) throw MatrixException("Lhs matrix of system is null.");

		auto rhsVals(rhs.getValues());
		const auto rhsRows(rhs.getRows());
		const auto rhsCols(rhs.getCols());

		const Eigen::Matrix<T, Dynamic, Dynamic> rhsE = Map<Eigen::Matrix<T, Dynamic, Dynamic>>(
			rhsVals->data(), rhsRows, rhsCols);

		auto ret(make_unique<DenseMatrix<T>>(DenseMatrix()));

		Eigen::Matrix<T, Dynamic, Dynamic> solution;
		if (decomposition == DecompositionType::BDCSVD)
		{
			solution = this->matrix_->bdcSvd(ComputeThinU | ComputeThinV).solve(rhsE).eval();
		}
		else if (decomposition == DecompositionType::JACOBISVD)
		{
			solution = this->matrix_->jacobiSvd(ComputeThinU | ComputeThinV).solve(rhsE).eval();
		}
		else if (decomposition == DecompositionType::COLPIVHOUSEHOLDER)
		{
			solution = this->matrix_->colPivHouseholderQr().solve(rhsE).eval();
		}
		else
		{
			throw MatrixException("Unsupported solver.");
		}

		unique_ptr<vector<T>> tmpResultVec(make_unique<vector<T>>(vector<T>(solution.size())));
		Map<Eigen::Matrix<T, Dynamic, Dynamic>>(tmpResultVec->data(), solution.rows(), solution.cols()) = solution;

		ret->setValues(*tmpResultVec, tmpResultVec->size(), 1);

		return move(ret);
	}
}
