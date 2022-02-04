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
	class DenseMatrix final : public Matrix<T>
	{
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_;

	public:
		explicit DenseMatrix() : Matrix<T>(0, 0), matrix_(0, 0)
		{
		}

		explicit DenseMatrix(const size_t& initialRows, const size_t& initialColumns) :
			Matrix<T>(initialRows, initialColumns), matrix_(initialRows, initialColumns)
		{
		}

		explicit DenseMatrix(const Matrix<T>& matrix)
		{
			this->setCols(matrix.getCols());
			this->setRows(matrix.getRows());

			for (auto i = 0; i < matrix.getCols(); ++i)
			{
				this->setColumn(i, matrix.getColumn(i));
			}
		}

		~DenseMatrix() override = default;

		DenseMatrix(const DenseMatrix&) = delete;
		DenseMatrix& operator=(const DenseMatrix&) = delete;

		DenseMatrix(const DenseMatrix&&) = delete;
		DenseMatrix& operator=(const DenseMatrix&&) = delete;

		size_t getRows() const override
		{
			return this->rows_;
		}

		void setRows(const size_t& rows) override
		{
			this->rows_ = rows;
			this->setSize(this->rows_, this->cols_);
		}

		size_t getCols() const override
		{
			return this->cols_;
		}

		void setCols(const size_t& cols) override
		{
			this->cols_ = cols;
			this->setSize(this->rows_, this->cols_);
		}

		size_t getSize() const override
		{
			return this->rows_ * this->cols_;
		}

		void setSize(const size_t& rows, const size_t& cols) override
		{
			this->rows_ = rows;
			this->cols_ = cols;
			this->matrix_.conservativeResize(this->rows_, this->cols_);
		}

		T get(const size_t& rPos, const size_t& cPos) const override
		{
			if (rPos >= this->rows_) throw MatrixException("Invalid row position queried: " + std::to_string(rPos));
			if (cPos >= this->cols_) throw MatrixException("Invalid column position queried: " + std::to_string(cPos));

			return this->matrix_(rPos, cPos);
		}

		void set(const size_t& rPos, const size_t& cPos, const T& value) override
		{
			if (rPos >= this->rows_) throw MatrixException("Invalid row position was tried to be set: " + std::to_string(rPos));
			if (cPos >= this->cols_) throw MatrixException("Invalid column position was tried to be set: " + std::to_string(cPos));

			this->matrix_(rPos, cPos) = value;
		}

		std::vector<T> getValues() const override
		{
			using namespace Eigen;

			std::vector<T> retVector(this->matrix_.size());
			Map<Eigen::Matrix<T, Dynamic, Dynamic>>(retVector.data(), this->getRows(), this->getCols()) = this->matrix_;
			return retVector;
		}

		void setValues(const std::vector<T>& values, const size_t& rows, const size_t& cols) override
		{
			using namespace Eigen;

			this->rows_ = rows;
			this->cols_ = cols;

			auto tmpVals(values);
			this->matrix_ = Map<Eigen::Matrix<T, Dynamic, Dynamic>>(tmpVals.data(), rows, cols);
		}

		std::vector<T> getColumn(const size_t& cPos) const override
		{
			using namespace Eigen;

			if (cPos >= this->cols_) throw MatrixException("Invalid column position queried: " + std::to_string(cPos));

			const auto retrievedCol = Eigen::Matrix<T, Dynamic, 1>(this->matrix_.col(cPos).matrix());

			std::vector<T> retVector(retrievedCol.size());

			Map<Eigen::Matrix<T, Dynamic, 1>>(retVector.data(), retrievedCol.size()) = retrievedCol;

			return retVector;
		}

		void addColumn(std::vector<T> values) override
		{
			using namespace Eigen;

			const auto instance = &this->matrix_;
			const Map<Eigen::Matrix<T, Dynamic, 1>> newVals(values.data(), values.size());

			instance->conservativeResize(instance->rows(), instance->cols() + 1);
			instance->col(instance->cols() - 1) = newVals;
		}

		void setColumn(const size_t& cPos, const std::vector<T>& values) override
		{
			using namespace Eigen;

			if (cPos >= this->cols_) throw MatrixException("Invalid column position was tried to be set: " + std::to_string(cPos));

			std::vector<T> tmpValues(values);

			const Map<Eigen::Matrix<T, Dynamic, 1>> newVals(tmpValues.data(), tmpValues.size());

			this->matrix_.block(0, cPos, tmpValues.size(), 1) << newVals;
		}

		std::vector<T> removeColumn(const size_t& cPos) override
		{
			using namespace Eigen;

			if (cPos >= this->cols_) throw MatrixException("Invalid column position was tried to be removed: " + std::to_string(cPos));

			const auto instance = &this->matrix_;

			Eigen::Matrix<T, Dynamic, 1> removedVals;

			if (static_cast<long long>(cPos) < instance->cols())
			{
				removedVals = Eigen::Matrix<T, Dynamic, 1>(instance->col(cPos).matrix());
				instance->block(0, cPos, instance->rows(), instance->cols() - cPos - 1) = instance->rightCols(instance->cols() - cPos - 1);
			}

			instance->conservativeResize(instance->rows(), instance->cols() - 1);

			std::vector<T> retVector(removedVals.size());
			Map<Eigen::Matrix<T, Dynamic, 1>>(retVector.data(), removedVals.size()) = removedVals;

			return retVector;
		}

		std::vector<T> getRow(const size_t& rPos) const override
		{
			using namespace Eigen;

			if (rPos >= this->rows_) throw MatrixException("Invalid row position queried: " + std::to_string(rPos));

			const auto retrievedRow = Eigen::Matrix<T, 1, Dynamic>(this->matrix_.row(rPos).matrix());

			std::vector<T> retVector(retrievedRow.size());

			Map<Eigen::Matrix<T, 1, Dynamic>>(retVector.data(), retrievedRow.size()) = retrievedRow;

			return retVector;
		}

		void addRow(std::vector<T> values) override
		{
			using namespace Eigen;

			const auto instance = &this->matrix_;
			const Map<Eigen::Matrix<T, 1, Dynamic>> newVals(values.data(), values.size());

			instance->conservativeResize(instance->rows() + 1, instance->cols());
			instance->row(instance->rows() - 1) = newVals;
		}

		void setRow(const size_t& rPos, const std::vector<T>& values) override
		{
			using namespace Eigen;

			if (rPos >= this->rows_) throw MatrixException("Invalid row position was tried to be set: " + std::to_string(rPos));

			std::vector<T> tmpValues(values);

			const Map<Eigen::Matrix<T, Dynamic, 1>> newVals(tmpValues.data(), tmpValues.size());

			this->matrix_.block(rPos, 0, 1, tmpValues.size()) << newVals;
		}

		std::vector<T> removeRow(const size_t& rPos) override
		{
			using namespace Eigen;

			if (rPos >= this->rows_) throw MatrixException("Invalid row position was tried to be removed: " + std::to_string(rPos));

			const auto instance = &this->matrix_;

			Eigen::Matrix<T, 1, Dynamic> removedVals;

			if (static_cast<long long>(rPos) < instance->rows())
			{
				removedVals = Eigen::Matrix<T, 1, Dynamic>(instance->row(rPos).matrix());
				instance->block(rPos, 0, instance->rows() - rPos - 1, instance->cols()) = instance->bottomRows(instance->rows() - rPos - 1);
			}

			instance->conservativeResize(instance->rows() - 1, instance->cols());

			std::vector<T> retVector(removedVals.size());
			Map<Eigen::Matrix<T, Dynamic, 1>>(retVector.data(), removedVals.size()) = removedVals;

			return retVector;
		}


		Matrix<T>* transpose(const bool& inplace = false) override
		{
			if (inplace)
			{
				this->matrix_.transposeInPlace();
				return this;
			}

			const Matrix<T>* upcast = this;
			auto* retVal = new DenseMatrix(*upcast);
			retVal->transpose(true);
			return retVal;
		}


		void addToElementByPos(const size_t& rPos, const size_t& cPos, const T& val) override
		{
			if (rPos >= this->rows_) throw MatrixException("Invalid row position was tried to be set: " + std::to_string(rPos));
			if (cPos >= this->cols_) throw MatrixException("Invalid column position was tried to be set: " + std::to_string(cPos));

			this->matrix_(rPos, cPos) += val;
		}

		void multiplyElementByPos(const size_t& rPos, const size_t& cPos, const T& val) override
		{
			if (rPos >= this->rows_) throw MatrixException("Invalid row position was tried to be set: " + std::to_string(rPos));
			if (cPos >= this->cols_) throw MatrixException("Invalid column position was tried to be set: " + std::to_string(cPos));

			this->matrix_(rPos, cPos) *= val;
		}

		Matrix<T>* operator*(const Matrix<T>* rhs) override
		{
			using namespace Eigen;

			if (this->cols_ != rhs->getRows()) throw MatrixException("Matrices are not multipliable: l:" + std::to_string(this->cols_) + "- r:" + std::to_string(rhs->getRows()));
			
			Matrix<T>* retMatrix = new DenseMatrix(rhs->getRows(), rhs->getCols());

			const auto lhsMatrix = this->matrix_;

			auto rhsVals(rhs->getValues());
			const auto rhsRows(rhs->getRows());
			const auto rhsCols(rhs->getCols());

			const Eigen::Matrix<T, Dynamic, Dynamic> rhsMatrix = Map<Eigen::Matrix<T, Dynamic, Dynamic>>(rhsVals.data(), rhsRows, rhsCols);

			const auto tmpResult = Eigen::Matrix<T, Dynamic, Dynamic>(lhsMatrix * rhsMatrix);
			std::vector<T> tmpResultVec(tmpResult.size());
			Eigen::Map<Eigen::Matrix<T, Dynamic, Dynamic>>(tmpResultVec.data(), tmpResult.rows(), tmpResult.cols()) = tmpResult;

			retMatrix->setValues(tmpResultVec, rhs->getRows(), rhs->getCols());

			return retMatrix;
		}

		Matrix<T>* operator*(const Matrix<T>& rhs) override
		{
			return *this * (&rhs);
		}

		std::vector<T> operator*(const std::vector<T>* rhs) override
		{
			using namespace Eigen;

			if (this->cols_ != rhs->size()) throw MatrixException("Matrices are not multipliable: l:" + std::to_string(this->cols_) + "- r:" + std::to_string(rhs->size()));

			auto tmpRhs(*rhs);
			const Map<Eigen::Matrix<T, Dynamic, 1>> rhsVector(tmpRhs.data(), tmpRhs.size());

			const auto product = Eigen::Matrix<T, Dynamic, 1>(this->matrix_ * rhsVector);

			std::vector<T> retVector(product.size());
			Map<Eigen::Matrix<T, Dynamic, 1>>(retVector.data(), product.size()) = product;

			return retVector;
		}

		std::vector<T> operator*(const std::vector<T>& rhs) override
		{
			return *this * (&rhs);
		}

		Matrix<T>* operator+(const Matrix<T>* rhs) const override
		{
			using namespace Eigen;

			if (this->cols_ != rhs->getCols() || this->rows_ != rhs->getRows())
				throw MatrixException("Matrices cannot be added together: l: " + std::to_string(this->rows_) + "x" + std::to_string(this->cols_) + "; r: " + std::to_string(rhs->getRows()) + "x" + std::to_string(rhs->getCols()));

			Matrix<T>* retMatrix = new DenseMatrix(rhs->getRows(), rhs->getCols());

			const auto lhsMatrix = this->matrix_;

			auto rhsVals(rhs->getValues());
			const auto rhsRows(rhs->getRows());
			const auto rhsCols(rhs->getCols());

			const Eigen::Matrix<T, Dynamic, Dynamic> rhsMatrix = Map<Eigen::Matrix<T, Dynamic, Dynamic>>(rhsVals.data(), rhsRows, rhsCols);

			const auto tmpResult = Eigen::Matrix<T, Dynamic, Dynamic>(lhsMatrix + rhsMatrix);
			std::vector<T> tmpResultVec(tmpResult.size());
			Map<Eigen::Matrix<T, Dynamic, Dynamic>>(tmpResultVec.data(), tmpResult.rows(), tmpResult.cols()) = tmpResult;

			retMatrix->setValues(tmpResultVec, rhs->getRows(), rhs->getCols());

			return retMatrix;
		}

		Matrix<T>* operator+(const Matrix<T>& rhs) const override
		{
			return *this + (&rhs);
		}

		Matrix<T>* operator-(const Matrix<T>* rhs) const override
		{
			using namespace Eigen;

			if (this->cols_ != rhs->getCols() || this->rows_ != rhs->getRows())
				throw MatrixException("Matrices are not substractable: l: " + std::to_string(this->rows_) + "x" + std::to_string(this->cols_) + "; r: " + std::to_string(rhs->getRows()) + "x" + std::to_string(rhs->getCols()));

			Matrix<T>* retMatrix = new DenseMatrix(rhs->getRows(), rhs->getCols());

			const auto lhsMatrix = this->matrix_;

			auto rhsVals(rhs->getValues());
			const auto rhsRows(rhs->getRows());
			const auto rhsCols(rhs->getCols());

			const Eigen::Matrix<T, Dynamic, Dynamic> rhsMatrix = Map<Eigen::Matrix<T, Dynamic, Dynamic>>(rhsVals.data(), rhsRows, rhsCols);

			const auto tmpResult = Eigen::Matrix<T, Dynamic, Dynamic>(lhsMatrix - rhsMatrix);
			std::vector<T> tmpResultVec(tmpResult.size());
			Map<Eigen::Matrix<T, Dynamic, Dynamic>>(tmpResultVec.data(), tmpResult.rows(), tmpResult.cols()) = tmpResult;

			retMatrix->setValues(tmpResultVec, rhs->getRows(), rhs->getCols());

			return retMatrix;
		}

		Matrix<T>* operator-(const Matrix<T>& rhs) const override
		{
			return *this - (&rhs);
		}

		friend std::ostream& operator<<(std::ostream& out, const Matrix<T>& matrix)
		{
			for (auto i = 0; i < matrix.getRows(); ++i)
			{
				for (auto j = 0; j < matrix.getCols(); ++j)
				{
					out << matrix.get(i, j) << " ";
				}

				out << "\n";
			}

			return out;
		}

		std::string toString() const override
		{
			std::stringstream output;
			output << matrix_;
			return output.str();
		}

		void exchangeRows(const size_t& lRow, const size_t& rRow) override
		{
			if (lRow >= this->getRows() || rRow >= this->getRows()) throw MatrixException("Row position out of bounds for exchanging rows: " + std::to_string(lRow) + "-" + std::to_string(rRow));

			this->matrix_.row(lRow).swap(this->matrix_.row(rRow));
		}

		void exchangeCols(const size_t& lCol, const size_t& rCol) override
		{
			if (lCol >= this->getCols() || rCol >= this->getCols()) throw MatrixException("Row position out of bounds for exchanging rows: " + std::to_string(lCol) + "-" + std::to_string(rCol));

			this->matrix_.col(lCol).swap(this->matrix_.col(rCol));
		}
	};
}
