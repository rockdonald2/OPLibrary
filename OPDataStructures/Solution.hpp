#pragma once

#include "Matrix.hpp"
#include "SolverException.hpp"

namespace OPLibrary
{
	/**
	 * \brief Representation of different solution statuses.
	 */
	enum class SolutionStatus
	{
		FEASIBLE,
		UNFEASIBLE,
		PRIMAL_UNFEASIBLE,
		DUAL_UNFEASIBLE,
		NONINITIALIZED
	};

	template <typename T>
		requires std::floating_point<T>
	class Solution final
	{
		SolutionStatus status_;
		T optimal_;
		std::shared_ptr<Matrix<T>> x_;
		std::shared_ptr<Matrix<T>> y_;
		std::shared_ptr<Matrix<T>> s_;

	public:
		Solution(const SolutionStatus& status, const T& optimal, Matrix<T>* x, Matrix<T>* y, Matrix<T>* s) :
			status_(status), optimal_(optimal), x_(x), y_(y), s_(s)
		{
		}

		Solution(const SolutionStatus& status, const T& optimal, std::unique_ptr<Matrix<T>>& x, std::unique_ptr<Matrix<T>>& y, std::unique_ptr<Matrix<T>>& s) :
			status_(status), optimal_(optimal), x_(x.release()), y_(y.release()), s_(s.release())
		{
		}

		Solution(const SolutionStatus& status, const T& optimal, const std::shared_ptr<Matrix<T>>& x, const std::shared_ptr<Matrix<T>>& y, const std::shared_ptr<Matrix<T>>& s) :
			status_(status), optimal_(optimal), x_(x), y_(y), s_(s)
		{
		}

		[[nodiscard]] T getOptimalValue() const;
		[[nodiscard]] SolutionStatus getSolutionStatus() const;
		[[nodiscard]] std::string getSolutionStatusString() const;
		[[nodiscard]] std::shared_ptr<Matrix<T>> getPrimalSolution() const;
		[[nodiscard]] std::shared_ptr<Matrix<T>> getDualSolutionY() const;
		[[nodiscard]] std::shared_ptr<Matrix<T>> getDualSolutionS() const;

		friend std::ostream& operator<<(std::ostream& out, const Solution<T>& sol)
		{
			assert(sol.x_ != nullptr && sol.y_ != nullptr && sol.s_ != nullptr && "Tried to output solution without data.");

			out << "xT:\t";
			out << *sol.x_->transpose() << "\n";
			out << "yT:\t";
			out << *sol.y_->transpose() << "\n";
			out << "sT:\t";
			out << *sol.s_->transpose() << "\n";

			return out;
		}

		[[nodiscard]] std::string toString() const;
	};

	template <typename T> requires std::floating_point<T>
	T Solution<T>::getOptimalValue() const
	{
		return optimal_;
	}

	template <typename T> requires std::floating_point<T>
	SolutionStatus Solution<T>::getSolutionStatus() const
	{
		return status_;
	}

	template <typename T> requires std::floating_point<T>
	std::string Solution<T>::getSolutionStatusString() const
	{
		switch (this->status_)
		{
		case SolutionStatus::FEASIBLE:
			return "feasible";
		case SolutionStatus::UNFEASIBLE:
			return "unfeasible";
		case SolutionStatus::PRIMAL_UNFEASIBLE:
			return "primal unfeasible";
		case SolutionStatus::DUAL_UNFEASIBLE:
			return "dual unfeasible";
		case SolutionStatus::NONINITIALIZED:
			return "uninitialized";
		}

		throw SolverException("Unknown final status.");
	}

	template <typename T>
		requires std::floating_point<T>
	std::shared_ptr<Matrix<T>> Solution<T>::getPrimalSolution() const
	{
		return x_;
	}

	template <typename T>
		requires std::floating_point<T>
	std::shared_ptr<Matrix<T>> Solution<T>::getDualSolutionY() const
	{
		return y_;
	}

	template <typename T>
		requires std::floating_point<T>
	std::shared_ptr<Matrix<T>> Solution<T>::getDualSolutionS() const
	{
		return s_;
	}

	template <typename T>
		requires std::floating_point<T>
	std::string Solution<T>::toString() const
	{
		if (this->x_ == nullptr || this->y_ == nullptr || this->s_ == nullptr) throw MatrixException("Tried to print a solution with null matrices.");

		std::stringstream output;
		output << *this;
		return output.str();
	}
}
