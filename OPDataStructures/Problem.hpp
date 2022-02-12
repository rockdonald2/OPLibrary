#pragma once

#include "Matrix.hpp"

namespace OPLibrary
{
	/**
	 * \brief Abstract representation of the optimization problem.
	 */
	template <typename T>
	class Problem final
	{
		Matrix<T>* constraints_;
		Matrix<T>* constraintObjectives_;
		Matrix<T>* objectives_;

	public:
		Problem() : constraints_(nullptr), constraintObjectives_(nullptr), objectives_(nullptr) {}
		Problem(Matrix<T>* cnstrs, Matrix<T>* cnstrsObjs, Matrix<T>* objs) : constraints_(cnstrs), constraintObjectives_(cnstrsObjs), objectives_(objs) {}

		/**
		 * \brief Sets the constraints matrix.
		 */
		void setConstraints(Matrix<T>* cnstrs);
		/**
		 * \brief Sets the constraint objectives vector.
		 */
		void setConstraintsObjectives(Matrix<T>* cnstrsObjs);
		/**
		 * \brief Sets the objectives vector.
		 */
		void setObjectives(Matrix<T>* objs);

		/**
		 * \brief Returns the constraints matrix.
		 */
		Matrix<T>* getConstraints() const;
		/**
		 * \brief Returns the constraint objectives vector.
		 */
		Matrix<T>* getConstraintsObjectives() const;
		/**
		 * \brief Returns the objectives vector.
		 */
		Matrix<T>* getObjectives() const;
	};

	template <typename T>
	void Problem<T>::setConstraints(Matrix<T>* cnstrs)
	{
		this->constraints_ = cnstrs;
	}

	template <typename T>
	void Problem<T>::setConstraintsObjectives(Matrix<T>* cnstrsObjs)
	{
		this->constraintObjectives_ = cnstrsObjs;
	}

	template <typename T>
	void Problem<T>::setObjectives(Matrix<T>* objs)
	{
		this->objectives_ = objs;
	}

	template <typename T>
	Matrix<T>* Problem<T>::getConstraints() const
	{
		return constraints_;
	}

	template <typename T>
	Matrix<T>* Problem<T>::getConstraintsObjectives() const
	{
		return constraintObjectives_;
	}

	template <typename T>
	Matrix<T>* Problem<T>::getObjectives() const
	{
		return objectives_;
	}
}
