#pragma once

#include <concepts>
#include <memory>

#include "Problem.hpp"

namespace OPLibrary
{
	template <typename T>
		requires std::floating_point<T>
	class SOCPProblem final : public Problem<T>
	{
		std::vector<size_t> coneSizes_;

	public:
		SOCPProblem(const ObjectiveDirection& direction, Matrix<T>* cnstrs, Matrix<T>* cnstrsObjs, Matrix<T>* objs, const std::vector<size_t>& coneSizes)
			: Problem<T>(direction, cnstrs, cnstrsObjs, objs), coneSizes_(coneSizes) {}
		SOCPProblem(const ObjectiveDirection& direction, std::unique_ptr<Matrix<T>>& cnstrs, std::unique_ptr<Matrix<T>>& cnstrsObjs, std::unique_ptr<Matrix<T>>& objs, const std::vector<size_t>& coneSizes)
			: Problem<T>(direction, cnstrs, cnstrsObjs, objs), coneSizes_(coneSizes) {}
		SOCPProblem(const ObjectiveDirection& direction, const std::shared_ptr<Matrix<T>>& cnstrs, const std::shared_ptr<Matrix<T>>& cnstrsObjs, const std::shared_ptr<Matrix<T>>& objs, const std::vector<size_t>& coneSizes)
			: Problem<T>(direction, cnstrs, cnstrsObjs, objs), coneSizes_(coneSizes) {}

		void setConeSizes(const std::vector<size_t>& coneSizes)
		{
			this->coneSizes_ = coneSizes;
		}

		[[nodiscard]] std::vector<size_t> getConeSizes() const
		{
			return coneSizes_;
		}

		friend std::ostream& operator<<(std::ostream& out, const SOCPProblem<T>& problem)
		{
			assert(problem.constraints_ != nullptr && problem.constraintObjectives_ != nullptr && problem.objectives_ != nullptr && "Tried to output null problem.");

			out << "A:\n";
			out << *problem.constraints_ << "\n";
			out << "bT:\t";
			out << *problem.constraintObjectives_->transpose() << "\n";
			out << "cT:\t";
			out << *problem.objectives_->transpose() << "\n";
			out << "coneSizes:\t";
			out << problem.coneSizes_ << "\n";

			return out;
		}

		[[nodiscard]] std::string toString() const override;
	};

	template <typename T>
		requires std::floating_point<T>
	std::string SOCPProblem<T>::toString() const
	{
		if (this->constraintObjectives_ == nullptr || this->objectives_ == nullptr || this->constraints_ == nullptr) throw MatrixException("Tried to print a problem with null matrices.");

		std::stringstream output;
		output << *this;
		return output.str();
	}
}
