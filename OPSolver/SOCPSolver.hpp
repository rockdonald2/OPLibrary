﻿#pragma once

#include <format>
#include <map>
#include "Solver.hpp"
#include <numbers>
#include <vector>
#include <string>

#include "Logger.h"
#include "Solution.hpp"
#include "SolverException.hpp"

namespace OPLibrary
{
	/**
	 * \brief Solver for the SOCP algorithm.
	 */
	template <typename T>
	class SOCPSolver : public Solver<T>
	{
		long double theta_;
		long double epsilon_;
		long double tau_;
		long double alpha_;
		long double mu_;

	public:
		SOCPSolver() : theta_(std::numbers::pi_v<long double> / 4), epsilon_(1.0e-6), tau_(1.0 / 2), alpha_(0.5), mu_(1)
		{
			using namespace std;

			this->INITIALIZABLE_ARGS = { "theta", "epsilon", "tau", "alpha", "mu" };
			this->INITIALIZATOR.insert(make_pair<string, long double*>("theta", &theta_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("epsilon", &epsilon_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("tau", &tau_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("alpha", &alpha_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("mu", &mu_));
		}

		void setConstraints(Matrix<T>* cnstrs) override;
		void setConstraintObjectives(Matrix<T>* cnstrsObjs) override;
		void setObjectives(Matrix<T>* objs) override;

		SolutionStatus solve() override;

		Solution getSolution() override;
	};

	template <typename T>
	void SOCPSolver<T>::setConstraints(Matrix<T>* cnstrs)
	{
		this->constraints_ = cnstrs;
	}

	template <typename T>
	void SOCPSolver<T>::setConstraintObjectives(Matrix<T>* cnstrsObjs)
	{
		this->constraintObjectives_ = cnstrsObjs;
	}

	template <typename T>
	void SOCPSolver<T>::setObjectives(Matrix<T>* objs)
	{
		this->objectives_ = objs;
	}

	template <typename T>
	SolutionStatus SOCPSolver<T>::solve()
	{
		if (this->constraints_ == nullptr || this->objectives_ == nullptr || this->constraintObjectives_ == nullptr)
		{
			throw SolverException("Cannot solve optimization problem without correctly setting the constraints, objectives and constraint objectives.");
		}

		return SolutionStatus::NONOPTIMAL;
	}

	template <typename T>
	Solution SOCPSolver<T>::getSolution()
	{
		return {};
	}
}
