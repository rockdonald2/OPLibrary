#pragma once

#include "Matrix.hpp"
#include "MatrixException.hpp"

namespace OPLibrary
{
	enum class ProblemType
	{
		SOCP,
		INVALID
	};

	enum class ObjectiveDirection
	{
		MAXIMIZE,
		MINIMIZE,
		INVALID
	};

	/**
	 * \brief Representation of the optimization problem.
	 */
	template <typename T>
		requires std::floating_point<T>
	class Problem
	{
	protected:
		ObjectiveDirection direction_;
		std::shared_ptr<Matrix<T>> constraints_;
		std::shared_ptr<Matrix<T>> constraintObjectives_;
		std::shared_ptr<Matrix<T>> objectives_;

	public:
		Problem() : direction_(ObjectiveDirection::INVALID), constraints_(nullptr), constraintObjectives_(nullptr), objectives_(nullptr) {}
		Problem(const ObjectiveDirection& direction, Matrix<T>* cnstrs, Matrix<T>* cnstrsObjs, Matrix<T>* objs) : direction_(direction), constraints_(cnstrs), constraintObjectives_(cnstrsObjs), objectives_(objs) {}
		/**
		 * \brief Creates a new Problem instance, but releases all argument pointers, transfers ownership.
		 */
		Problem(const ObjectiveDirection& direction, std::unique_ptr<Matrix<T>>& cnstrs, std::unique_ptr<Matrix<T>>& cnstrsObjs, std::unique_ptr<Matrix<T>>& objs) : direction_(direction), constraints_(cnstrs.release()), constraintObjectives_(cnstrsObjs.release()), objectives_(objs.release()) {}
		Problem(const ObjectiveDirection& direction, const std::shared_ptr<Matrix<T>>& cnstrs, const std::shared_ptr<Matrix<T>>& cnstrsObjs, const std::shared_ptr<Matrix<T>>& objs) : direction_(direction), constraints_(cnstrs), constraintObjectives_(cnstrsObjs), objectives_(objs) {}

		virtual ~Problem() = default;

		virtual ObjectiveDirection getObjectiveDirection() const;

		virtual void setObjectiveDirection(const ObjectiveDirection& direction);

		/**
		 * \brief Sets the constraints matrix.
		 */
		virtual void setConstraints(Matrix<T>* cnstrs);
		/**
		 * \brief Sets the constraints matrix, but releases the original pointer, transfers ownership.
		 */
		virtual void setConstraints(std::unique_ptr<Matrix<T>>& cnstrs);
		/**
		 * \brief Sets the constraint objectives vector.
		 */
		virtual void setConstraintsObjectives(Matrix<T>* cnstrsObjs);
		/**
		 * \brief Sets the constraints objectives vector, but releases the original pointer, transfers ownership.
		 */
		virtual void setConstraintsObjectives(std::unique_ptr<Matrix<T>>& cnstrsObjs);
		/**
		 * \brief Sets the objectives vector.
		 */
		virtual void setObjectives(Matrix<T>* objs);
		/**
		 * \brief Sets the objectives vector, but releases the original pointer, transfers ownership.
		 */
		virtual void setObjectives(std::unique_ptr<Matrix<T>> objs);

		/**
		 * \brief Returns the constraints matrix.
		 */
		std::shared_ptr<Matrix<T>> getConstraints() const;
		/**
		 * \brief Returns the constraint objectives vector.
		 */
		std::shared_ptr<Matrix<T>> getConstraintsObjectives() const;
		/**
		 * \brief Returns the objectives vector.
		 */
		std::shared_ptr<Matrix<T>> getObjectives() const;

		friend std::ostream& operator<<(std::ostream& out, const Problem<T>& problem)
		{
			assert(problem.constraints_ != nullptr && problem.constraintObjectives_ != nullptr && problem.objectives_ != nullptr && "Tried to output null problem.");

			out << "A:\n";
			out << *problem.constraints_ << "\n";
			out << "bT:\t";
			out << *problem.constraintObjectives_->transpose() << "\n";
			out << "cT:\t";
			out << *problem.objectives_->transpose() << "\n";

			return out;
		}

		[[nodiscard]] virtual std::string toString() const;
	};

	template <typename T> requires std::floating_point<T>
	ObjectiveDirection Problem<T>::getObjectiveDirection() const
	{
		return direction_;
	}

	template <typename T> requires std::floating_point<T>
	void Problem<T>::setObjectiveDirection(const ObjectiveDirection& direction)
	{
		direction_ = direction;
	}

	template <typename T>
		requires std::floating_point<T>
	void Problem<T>::setConstraints(Matrix<T>* cnstrs)
	{
		this->constraints_ = std::shared_ptr<Matrix<T>>(cnstrs);
	}

	template <typename T>
		requires std::floating_point<T>
	void Problem<T>::setConstraints(std::unique_ptr<Matrix<T>>& cnstrs)
	{
		this->setConstraints(cnstrs.release());
	}

	template <typename T>
		requires std::floating_point<T>
	void Problem<T>::setConstraintsObjectives(Matrix<T>* cnstrsObjs)
	{
		this->constraintObjectives_ = std::shared_ptr<Matrix<T>>(cnstrsObjs);
	}

	template <typename T>
		requires std::floating_point<T>
	void Problem<T>::setConstraintsObjectives(std::unique_ptr<Matrix<T>>& cnstrsObjs)
	{
		this->setConstraintsObjectives(cnstrsObjs.release());
	}

	template <typename T>
		requires std::floating_point<T>
	void Problem<T>::setObjectives(Matrix<T>* objs)
	{
		this->objectives_ = std::shared_ptr<Matrix<T>>(objs);
	}

	template <typename T>
		requires std::floating_point<T>
	void Problem<T>::setObjectives(std::unique_ptr<Matrix<T>> objs)
	{
		this->setObjectives(objs.release());
	}

	template <typename T>
		requires std::floating_point<T>
	std::shared_ptr<Matrix<T>> Problem<T>::getConstraints() const
	{
		return this->constraints_;
	}

	template <typename T>
		requires std::floating_point<T>
	std::shared_ptr<Matrix<T>> Problem<T>::getConstraintsObjectives() const
	{
		return this->constraintObjectives_;
	}

	template <typename T>
		requires std::floating_point<T>
	std::shared_ptr<Matrix<T>> Problem<T>::getObjectives() const
	{
		return this->objectives_;
	}

	template <typename T>
		requires std::floating_point<T>
	std::string Problem<T>::toString() const
	{
		if (this->constraintObjectives_ == nullptr || this->objectives_ == nullptr || this->constraints_ == nullptr) throw MatrixException("Tried to print a problem with null matrices.");

		std::stringstream output;
		output << *this;
		return output.str();
	}
}
