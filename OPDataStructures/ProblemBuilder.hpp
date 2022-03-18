#pragma once

#include <concepts>
#include <memory>

#include "Matrix.hpp"
#include "Problem.hpp"

namespace OPLibrary
{
	template <typename T>
		requires std::floating_point<T>
	class ProblemBuilder final
	{
		std::shared_ptr<Matrix<T>> cnstr_;
		std::shared_ptr<Matrix<T>> cnstrObjs_;
		std::shared_ptr<Matrix<T>> objs_;

	public:
		ProblemBuilder() = default;

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
			return std::make_shared<Problem<T>>(this->cnstr_, this->cnstrObjs_, this->objs_);
		}
	};
}
