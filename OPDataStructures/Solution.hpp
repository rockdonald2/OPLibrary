#pragma once

#include "Matrix.hpp"

namespace OPLibrary
{
	template <typename T>
		requires std::floating_point<T>
	class Solution final
	{
		std::shared_ptr<Matrix<T>> x_;
		std::shared_ptr<Matrix<T>> y_;
		std::shared_ptr<Matrix<T>> s_;

	public:
		Solution() : x_(nullptr), y_(nullptr), s_(nullptr) {}
		Solution(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s) : x_(x), y_(y), s_(s) {}
		Solution(std::unique_ptr<Matrix<T>>& x, std::unique_ptr<Matrix<T>>& y, std::unique_ptr<Matrix<T>>& s) : x_(x.release()), y_(y.release()), s_(s.release()) {}
		Solution(const std::shared_ptr<Matrix<T>>& x, const std::shared_ptr<Matrix<T>>& y, const std::shared_ptr<Matrix<T>>& s) : x_(x), y_(y), s_(s) {}

		[[nodiscard]] std::shared_ptr<Matrix<T>> getX() const;
		[[nodiscard]] std::shared_ptr<Matrix<T>> getY() const;
		[[nodiscard]] std::shared_ptr<Matrix<T>> getS() const;

		friend std::ostream& operator<<(std::ostream& out, const Solution<T>& sol)
		{
			assert(sol.x_ != nullptr && sol.y_ != nullptr && sol.s_ != nullptr && "Tried to output solution without data.");

			out << "The x solution is for the original primal problem:\n";
			out << *sol.x_ << "\n";
			out << "The y and s solutions are for the dual problem:\n";
			out << *sol.y_ << "\n";
			out << *sol.s_ << "\n";

			return out;
		}

		[[nodiscard]] std::string toString() const;
	};

	template <typename T>
		requires std::floating_point<T>
	std::shared_ptr<Matrix<T>> Solution<T>::getX() const
	{
		return x_;
	}

	template <typename T>
		requires std::floating_point<T>
	std::shared_ptr<Matrix<T>> Solution<T>::getY() const
	{
		return y_;
	}

	template <typename T>
		requires std::floating_point<T>
	std::shared_ptr<Matrix<T>> Solution<T>::getS() const
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
