#pragma once

#include <concepts>
#include <memory>

#include "Matrix.hpp"
#include "Problem.hpp"
#include "SOCPProblem.hpp"

namespace OPLibrary
{
	template <typename T>
		requires std::floating_point<T>
	class ProblemBuilder final
	{
		inline static const std::map<std::string, ProblemType> MAP_STR_TO_PROBLEM{ {"SOCP", ProblemType::SOCP} };
		ProblemType type_;
		ObjectiveDirection direction_;
		std::shared_ptr<Matrix<T>> cnstr_;
		std::shared_ptr<Matrix<T>> cnstrObjs_;
		std::shared_ptr<Matrix<T>> objs_;

	public:
		ProblemBuilder() : type_(ProblemType::INVALID), direction_(ObjectiveDirection::INVALID), cnstr_(nullptr), cnstrObjs_(nullptr), objs_(nullptr) {}

		ProblemBuilder& setProblemType(const ProblemType& type)
		{
			type_ = type;
			return *this;
		}

		ProblemBuilder& setProblemType(const std::string& type)
		{
			std::string temp(type);
			std::ranges::transform(temp, temp.begin(), [](const unsigned char c) { return std::toupper(c); });

			if (!MAP_STR_TO_PROBLEM.contains(temp))
				throw SolverException("Unsupported problem type.");

			switch (MAP_STR_TO_PROBLEM.at(temp))
			{
			case ProblemType::SOCP: return setProblemType(ProblemType::SOCP);
			}

			throw SolverException("Unsupported problem type.");
		}

		ProblemBuilder& setObjectiveDirection(const ObjectiveDirection& direction)
		{
			direction_ = direction;
			return *this;
		}

		ProblemBuilder& setConstraints(const Matrix<T>* cnstr)
		{
			return setConstraints(std::shared_ptr<Matrix<T>>(cnstr));
		}

		ProblemBuilder& setConstraints(const std::shared_ptr<Matrix<T>>& cnstr)
		{
			this->cnstr_ = cnstr;
			return *this;
		}

		ProblemBuilder& setConstraintsObjectives(const Matrix<T>* cnstrObjs)
		{
			return setConstraintsObjectives(std::shared_ptr<Matrix<T>>(cnstrObjs));
		}

		ProblemBuilder& setConstraintsObjectives(const std::shared_ptr<Matrix<T>>& cnstrObjs)
		{
			this->cnstrObjs_ = cnstrObjs;
			return *this;
		}

		ProblemBuilder& setObjectives(const Matrix<T>* objs)
		{
			return setObjectives(std::shared_ptr<Matrix<T>>(objs));
		}

		ProblemBuilder& setObjectives(const std::shared_ptr<Matrix<T>>& objs)
		{
			this->objs_ = objs;
			return *this;
		}

		[[nodiscard]] std::shared_ptr<Problem<T>> build() const
		{
			if (this->cnstrObjs_ == nullptr || this->objs_ == nullptr || this->cnstr_ == nullptr) throw MatrixException("Tried to build problem with null values.");
			if (this->direction_ == ObjectiveDirection::INVALID) throw SolverException("Unset objective direction.");

			switch (type_)
			{
			case ProblemType::SOCP: return std::make_shared<SOCPProblem<T>>(direction_, this->cnstr_, this->cnstrObjs_, this->objs_, std::vector<size_t>());
			}

			throw SolverException("Unsupported problem type set.");
		}
	};
}
