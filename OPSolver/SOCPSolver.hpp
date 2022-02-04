#pragma once

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
		const std::vector<std::string> INITIALIZABLE_ARGS = { "theta", "epsilon", "tau", "alpha", "mu" };
		std::map<std::string, long double*> INITIALIZATOR;

		long double theta_;
		long double epsilon_;
		long double tau_;
		long double alpha_;
		long double mu_;

	public:
		SOCPSolver() : theta_(std::numbers::pi_v<long double> / 4), epsilon_(1.0e-6), tau_(1.0 / 2), alpha_(0.5), mu_(1)
		{
			using namespace std;

			INITIALIZATOR.insert(make_pair<string, long double*>("theta", &theta_));
			INITIALIZATOR.insert(make_pair<string, long double*>("epsilon", &epsilon_));
			INITIALIZATOR.insert(make_pair<string, long double*>("tau", &tau_));
			INITIALIZATOR.insert(make_pair<string, long double*>("alpha", &alpha_));
			INITIALIZATOR.insert(make_pair<string, long double*>("mu", &mu_));
		}

		std::vector<std::string> getInitializableArgs() override
		{
			return INITIALIZABLE_ARGS;
		}

		void setInitializableArg(const std::string& arg, const long double& val) override
		{
			using namespace std;

			if (ranges::find(INITIALIZABLE_ARGS.begin(), INITIALIZABLE_ARGS.end(), arg) != INITIALIZABLE_ARGS.end())
			{
				*INITIALIZATOR[arg] = val;
			}
			else
			{
				Logger::getInstance().error(format("Tried to set an invalid initializable arg {} to value {}, skipping.", arg, val));
			}
		}

		void setConstraints(Matrix<T>* cnstrs) override
		{
			this->constraints_ = cnstrs;
		}

		void setConstraintObjectives(Matrix<T>* cnstrsObjs) override
		{
			this->constraintObjectives_ = cnstrsObjs;
		}
		void setObjectives(Matrix<T>* objs) override
		{
			this->objectives_ = objs;
		}

		SolutionStatus solve() override
		{
			if (this->constraints_ == nullptr || this->objectives_ == nullptr || this->constraintObjectives_ == nullptr)
			{
				throw SolverException("Cannot solve optimization problem without correctly setting the constraints, objectives and constraint objectives.");
			}

			return SolutionStatus::NONOPTIMAL;
		}

		Solution getSolution() override
		{
			return {};
		}
	};
}
