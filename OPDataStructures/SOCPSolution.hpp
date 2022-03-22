#pragma once

#include <concepts>
#include <memory>

#include "Solution.hpp"

namespace OPLibrary
{
	template <typename T>
		requires std::floating_point<T>
	class SOCPSolution : public Solution<T>
	{
		std::vector<size_t> cones_;

	public:
		SOCPSolution(const SolutionStatus& status, const T& optimal, Matrix<T>* x, Matrix<T>* y, Matrix<T>* s, const std::vector<size_t>& cones) :
			Solution<T>(status, optimal, x, y, s), cones_(cones)
		{
		}

		SOCPSolution(const SolutionStatus& status, const T& optimal, std::unique_ptr<Matrix<T>>& x, std::unique_ptr<Matrix<T>>& y, std::unique_ptr<Matrix<T>>& s, const std::vector<size_t>& cones) :
			Solution<T>(status, optimal, x.release(), y.release(), s.release()), cones_(cones)
		{
		}

		SOCPSolution(const SolutionStatus& status, const T& optimal, const std::shared_ptr<Matrix<T>>& x, const std::shared_ptr<Matrix<T>>& y, const std::shared_ptr<Matrix<T>>& s, const std::vector<size_t>& cones) :
			Solution<T>(status, optimal, x, y, s), cones_(cones)
		{
		}

		friend std::ostream& operator<<(std::ostream& out, const SOCPSolution<T>& sol)
		{
			assert(sol.x_ != nullptr && sol.y_ != nullptr && sol.s_ != nullptr && "Tried to output solution without data.");

			out << "xT:\t";
			out << *sol.x_->transpose() << "\n";
			out << "yT:\t";
			out << *sol.y_->transpose() << "\n";
			out << "sT:\t";
			out << *sol.s_->transpose() << "\n";
			out << "noc:\t";
			out << sol.getCones() << "\n";

			return out;
		}

		[[nodiscard]] std::vector<size_t> getCones() const;
	};

	template <typename T> requires std::floating_point<T>
	std::vector<size_t> SOCPSolution<T>::getCones() const
	{
		return cones_;
	}
}
