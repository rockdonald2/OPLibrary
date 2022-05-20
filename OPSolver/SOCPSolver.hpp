#pragma once

#include <format>
#include <numbers>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <execution>
#include <valarray>

#include "Solver.hpp"
#include "MatrixFactory.hpp"
#include "MatrixUtils.hpp"
#include "SOCPSolution.hpp"
#include "Solution.hpp"
#include "SolverException.hpp"
#include "VectorExtension.hpp"

namespace OPLibrary
{
	/**
	 * \brief Solver for the SOCP algorithm.
	 */
	template <typename T>
		requires std::floating_point<T>
	class SOCPSolver final : public Solver<T>
	{
		using CONE_GROUP = std::vector<std::shared_ptr<Matrix<T>>>;

		/**
		 * \brief Initializator with Strategy pattern to initialize solution vectors.
		 */
		class Initializator
		{
		public:
			Initializator() = default;
			virtual ~Initializator() = default;

			virtual void initialize(std::vector<std::shared_ptr<Matrix<T>>> x, std::shared_ptr<Matrix<T>> y, std::vector<std::shared_ptr<Matrix<T>>> s) = 0;

			Initializator(const Initializator&) = delete;
			void operator=(const Initializator&) = delete;
			Initializator(const Initializator&&) = delete;
			void operator=(const Initializator&&) = delete;
		};

		/**
		 * \brief Classic initializator, where for x0 and s0 vector their first value is 1, all other 0; for y all the values are 0.
		 */
		class ClassicInitializator final : public Initializator
		{
		public:
			void initialize(std::vector<std::shared_ptr<Matrix<T>>> x, std::shared_ptr<Matrix<T>> y, std::vector<std::shared_ptr<Matrix<T>>> s) override;
		};

		/**
		 * \brief Classic initializator, where for x0 vector first value is 1, s0 vector first value is 2, all other 0; for y all the values are 1.
		 */
		class Classic2Initializator final : public Initializator
		{
		public:
			void initialize(std::vector<std::shared_ptr<Matrix<T>>> x, std::shared_ptr<Matrix<T>> y, std::vector<std::shared_ptr<Matrix<T>>> s) override;
		};

		/**
		 * \brief Classic initializator, where for x0 and s0 vector their first value is 1, all other 0; for y all the values are 1.
		 */
		class Classic3Initializator final : public Initializator
		{
		public:
			void initialize(std::vector<std::shared_ptr<Matrix<T>>> x, std::shared_ptr<Matrix<T>> y, std::vector<std::shared_ptr<Matrix<T>>> s) override;
		};

		/**
		 * \brief Classic initializator, where for x0 and s0 vector their first value is 2, second value is 1, all other 0; for y all the values are 0.
		 */
		class Classic4Initializator final : public Initializator
		{
		public:
			void initialize(std::vector<std::shared_ptr<Matrix<T>>> x, std::shared_ptr<Matrix<T>> y, std::vector<std::shared_ptr<Matrix<T>>> s) override;
		};

		/**
		 * \brief Markowitz initializator, starting points depend on the number of variables inside a cone.
		 */
		class MarkowitzInitializator final : public Initializator
		{
		public:
			void initialize(std::vector<std::shared_ptr<Matrix<T>>> x, std::shared_ptr<Matrix<T>> y, std::vector<std::shared_ptr<Matrix<T>>> s) override;
		};

		long double epsilon_; // independent of number of cones
		std::vector<long double> alphaPrimal_; // will depend on number of cones
		std::vector<long double> alphaDual_; // will depend on number of cones
		long double mu_; // independent of noc
		long double rho_; // independent of number of cones
		long double sigma_; // independent of number of cones

		std::vector<size_t> n_; // will depend on number of cones
		size_t m_; // independent of number of cones

		size_t currIter_;
		size_t maxIters_;

		inline static const std::map<std::string, std::shared_ptr<Initializator>> INITS = {
			{"classic", std::make_shared<ClassicInitializator>()},
			{"classic2", std::make_shared<Classic2Initializator>()}, {"classic3", std::make_shared<Classic3Initializator>()}, {"classic4", std::make_shared<Classic4Initializator>()},
			{"markowitz", std::make_shared<MarkowitzInitializator>()}
		};
		std::shared_ptr<Initializator> init_;

		std::vector<std::shared_ptr<Matrix<T>>> x_; // will depend on number of cones
		std::shared_ptr<Matrix<T>> y_; // independent of number of cones
		std::vector<std::shared_ptr<Matrix<T>>> s_;

		std::vector<std::shared_ptr<Matrix<T>>> cnstrs_; // separate A matrices for individual cones
		std::vector<std::shared_ptr<Matrix<T>>> objs_; // separate b vectors for individual cones
		std::shared_ptr<Matrix<T>> cnstrsObjs_;

		static T calculateEigenMinOf(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "Eigenvalue min can only be calculated for vectors.");

			if (vec->getSize() == 1) return vec->get(0, 0);

			const auto tmpBlock(vec->block(1, 0, vec->getRows() - 1, 0));

			return vec->get(0, 0) - tmpBlock->norm();
		}

		static T calculateEigenMaxOf(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "Eigenvalue max can only be calculated for vectors.");

			if (vec->getSize() == 1) return vec->get(0, 0);

			const auto tmpBlock(vec->block(1, 0, vec->getRows() - 1, 0));

			return vec->get(0, 0) + tmpBlock->norm();
		}

		static std::shared_ptr<Matrix<T>> calculatePMatrixOf(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "P matrix can only be calculated for vectors.");

			const MatrixFactory<T> factory(MatrixType::DENSE);
			const auto n(vec->getRows());
			auto retMatrix(factory.createMatrix(n, n));

			// norm^2
			retMatrix->set(0, 0, pow(vec->norm(), 2));

			if (vec->getSize() == 1) return retMatrix;

			const auto x1(vec->get(0, 0));
			const auto x2n(vec->block(1, 0, vec->getRows() - 1, 0));
			const auto x2nT(x2n->transpose());

			// 2 * (x1 * x2:n)
			retMatrix->block(1, 0, retMatrix->getRows() - 1, 0, *(*x2n * x1) * 2);

			// 2 * (x1 * x2T:n)
			retMatrix->block(0, 1, 0, retMatrix->getCols() - 1, *(*x2nT * x1) * 2);

			// det(x) * En-1 + 2 * x2n * x2T:n
			// det(x) = x1^2 - norm(x2:n)^2
			const auto detx(calculateEigenMaxOf(vec) * calculateEigenMinOf(vec));

			const auto E(factory.createMatrix(n - 1, n - 1));
			E->setDiagonalValues(vector<T>(n - 1, 1));

			const auto helperTriangularMatrix(factory.createMatrix(n - 1, n - 1));

			for (size_t i = 0; i < (n - 1); ++i)
			{
				for (size_t j = 0; j < (n - 1); ++j)
				{
					helperTriangularMatrix->set(i, j, 2 * x2n->get(i, 0) * x2n->get(j, 0));
				}
			}

			retMatrix->block(1, 1, n - 1, n - 1, *(*E * detx) + *helperTriangularMatrix);

			return move(retMatrix);
		}

		static std::shared_ptr<Matrix<T>> calculateC1(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "C1 vector can only be calculated for vectors.");

			const auto n(vec->getRows());

			const MatrixFactory<T> factory(MatrixType::DENSE);
			auto c1(factory.createMatrix(n, 1));
			c1->set(0, 0, 0.5);

			if (vec->getSize() == 1) return move(c1);

			const auto x2n(vec->block(1, 0, n - 1, 0));
			const auto normx2n(x2n->norm());

			const auto zeroVector(vector<T>(n - 1, 0));
			if (*x2n->getValues() == zeroVector)
			{
				c1->block(1, 0, n - 1, 0, zeroVector);
			}
			else
			{
				c1->block(1, 0, n - 1, 0, *(*x2n / normx2n) / 2);
			}

			return move(c1);
		}

		static std::shared_ptr<Matrix<T>> calculateC2(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "C2 vector can only be calculated for vectors.");

			const auto n(vec->getRows());

			const MatrixFactory<T> factory(MatrixType::DENSE);
			auto c2(factory.createMatrix(n, 1));
			c2->set(0, 0, 0.5);

			if (vec->getSize() == 1) return move(c2);

			const auto x2n(vec->block(1, 0, n - 1, 0));
			const auto normx2n(x2n->norm());

			const auto zeroVector(vector<T>(n - 1, 0));
			if (*x2n->getValues() == zeroVector)
			{
				c2->block(1, 0, n - 1, 0, zeroVector);
			}
			else
			{
				c2->block(1, 0, n - 1, 0, *(*(*x2n * -1) / normx2n) / 2);
			}

			return move(c2);
		}

		static std::shared_ptr<Matrix<T>> calculatePowerOf(Matrix<T>* vec, const double power)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "Power of matrix can only be calculated for vectors.");

			const auto eigenmax(calculateEigenMaxOf(vec));
			const auto powerEigenmax(std::pow(eigenmax, power));
			const auto eigenmin(calculateEigenMinOf(vec));
			const auto powerEigenmin(std::pow(eigenmin, power));

			const auto c1(calculateC1(vec));
			const auto c2(calculateC2(vec));

			auto ret(*(*c1 * powerEigenmax) + *(*c2 * powerEigenmin));

			return move(ret);
		}

		static bool isInCone(Matrix<T>* vec)
		{
			assert(vec->getCols() == 1 && "Is In Cone method can only be applied to vectors.");
			return calculateEigenMinOf(vec) >= 0;
		}

		static bool isAllPositive(Matrix<T>* vec)
		{
			assert(vec->getCols() == 1 && "Is All Positive method can only be applied to vectors.");
			const auto vals(vec->getValues());
			return std::all_of(vals->begin(), vals->end(), [](const auto& val)
				{
					return val >= 0;
				});
		}

		static std::shared_ptr<Matrix<T>> calculateCombinedVectors(const std::vector<size_t>& individualLengths, const std::vector<std::shared_ptr<Matrix<T>>>& vectors)
		{
			const auto n(std::reduce(std::execution::par, individualLengths.begin(), individualLengths.end(), static_cast<size_t>(0)));
			auto combined(MatrixFactory<T>().createMatrix(n, 1));
			size_t lastPoz(0);

			for (const auto& v : vectors)
			{
				combined->block(lastPoz, 0, lastPoz + v->getRows() - 1, 0, *v->getValues());
				lastPoz += v->getRows();
			}

			return combined;
		}

		static std::shared_ptr<Matrix<T>> calculateHorizontallyCombinedMatrix(const std::vector<size_t>& individualRows, const std::vector<size_t>& individualCols, const std::vector<std::shared_ptr<Matrix<T>>>& matrices)
		{
			const auto m(*std::ranges::max_element(individualRows));
			const auto n(std::reduce(std::execution::par, individualCols.begin(), individualCols.end(), static_cast<size_t>(0)));

			auto combined(MatrixFactory<T>().createMatrix(m, n));
			size_t lastN(0);

			for (const auto& tm : matrices)
			{
				combined->block(0, lastN, tm->getRows() - 1, lastN + tm->getCols() - 1, tm);
				lastN += tm->getCols();
			}

			return move(combined);
		}

		[[nodiscard]] long double calculateMu() const;
		[[nodiscard]] long double calculateAlpha(const Matrix<T>* vec, const Matrix<T>* delta) const;
		[[nodiscard]] std::vector<std::shared_ptr<Matrix<T>>> calculateW() const;
		[[nodiscard]] std::vector<std::shared_ptr<Matrix<T>>> calculateV() const;
		[[nodiscard]] std::vector<std::shared_ptr<Matrix<T>>> calculatePv() const;

		[[nodiscard]] std::vector<std::shared_ptr<Matrix<T>>> calculateA_() const;

		[[nodiscard]] std::shared_ptr<Matrix<T>> calculateResidualB() const;
		[[nodiscard]] std::vector<std::shared_ptr<Matrix<T>>> calculateResidualC() const;
		[[nodiscard]] std::vector<std::shared_ptr<Matrix<T>>> calculateResidualC_() const;

		[[nodiscard]] bool checkIsTermination() const;
		[[nodiscard]] bool checkMu() const;
		[[nodiscard]] bool checkPrimalFeasibility() const;
		[[nodiscard]] bool checkDualFeasibility() const;
		[[nodiscard]] bool checkInfeasibility() const;

		SolutionStatus internalSolver();
		void internalSolverPreparations(std::shared_ptr<Matrix<T>> A, std::shared_ptr<Matrix<T>> b, std::shared_ptr<Matrix<T>> c);

	public:
		SOCPSolver() : epsilon_(1.0e-6), mu_(1.0), rho_(0.995), sigma_(1.0 / 10), m_(0),
			currIter_(0),
			maxIters_(3000), init_(std::make_shared<ClassicInitializator>()),
			y_(nullptr)
		{
			using namespace std;

			this->INITIALIZABLE_ARGS = { "epsilon", "mu", "rho", "sigma" };
			this->INITIALIZATOR.insert(make_pair<string, long double*>("epsilon", &epsilon_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("mu", &mu_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("rho", &rho_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("sigma", &sigma_));
		}

		~SOCPSolver() override = default;

		SOCPSolver(const SOCPSolver<T>&) = delete;
		SOCPSolver<T>& operator=(const SOCPSolver<T>&) = delete;

		SOCPSolver(SOCPSolver<T>&&) = delete;
		SOCPSolver<T>& operator=(SOCPSolver<T>&&) = delete;

		void setStartingPointInitializator(const std::string& init) override;

		SolutionStatus solve() override;
	};

	template <typename T>
		requires std::floating_point<T>
	bool SOCPSolver<T>::checkIsTermination() const
	{
		const auto ctx(*this->problem_->getObjectives()->transpose() * *calculateCombinedVectors(n_, x_));
		const auto bty(*this->problem_->getConstraintsObjectives()->transpose() * *y_);

		return currIter_ > maxIters_
			|| (checkMu() && checkPrimalFeasibility() && checkDualFeasibility())
			|| (containsNaN(ctx.get()) || containsNaN(bty.get()));
	}

	template <typename T> requires std::floating_point<T>
	bool SOCPSolver<T>::checkMu() const
	{
		return mu_ < epsilon_;
	}

	template <typename T>
		requires std::floating_point<T>
	long double SOCPSolver<T>::calculateMu() const
	{
		using namespace std;

		const auto n(std::reduce(std::execution::par, n_.begin(), n_.end(), static_cast<size_t>(0)));
		return (*(*calculateCombinedVectors(n_, x_)->transpose() * (2.0 * sigma_ / n)) * *calculateCombinedVectors(n_, s_))->get(0, 0);
	}

	template <typename T>
		requires std::floating_point<T>
	long double SOCPSolver<T>::calculateAlpha(const Matrix<T>* vec, const Matrix<T>* delta) const
	{
		using namespace std;

		// kiszamoljuk az A-t, B-t, C-t

		// deltax1 negyzete - deltax2n negyzetosszege
		// A = deltaX1^2 - sum(deltaX2:n^2)

		// x1 * deltax1 - x2n
		// B = 2 * (x1 * deltax1 - sum(xi * deltaxi))

		// x1 negyzete - x2n negyzetosszege
		// C = x1^2 - sum(xi^2)

		T A, B, C;
		const auto deltax1(delta->get(0, 0));
		const auto x1(vec->get(0, 0));

		if (vec->getSize() == 1)
		{
			A = pow(deltax1, 2);
			B = 2 * x1 * deltax1;
			C = pow(x1, 2);
		}
		else
		{
			const auto deltax2n(delta->block(1, 0, delta->getRows() - 1, 0)->getValues());

			{
				T sum(0);
				for_each(deltax2n->begin(), deltax2n->end(), [&sum](T n)
					{
						sum += pow(n, 2);
					});

				A = pow(deltax1, 2) - sum;
			}

			const auto x2n(vec->block(1, 0, vec->getRows() - 1, 0)->getValues());

			{
				T sum(0);
				for (size_t i = 0; i < x2n->size(); ++i)
				{
					sum += (x2n->at(i) * deltax2n->at(i));
				}

				B = x1 * deltax1 - sum;
				B *= 2;
			}
			{
				T sum(0);
				for_each(x2n->begin(), x2n->end(), [&sum](T n)
					{
						sum += pow(n, 2);
					});

				C = pow(x1, 2) - sum;
			}
		}

		vector<T> roots;
		if (A == 0.0)
		{
			// Ax2 + Bx + C = 0
			// Ha A == 0
			// akkor Bx + C = 0 => x = -C/B;

			roots.push_back(-C / B);
		}
		else
		{
			roots = solvePolynomial<T>({ A, B, C });
		}

		roots.push_back(1.0);
		roots.erase(remove_if(roots.begin(), roots.end(), [](T n)
			{
				return n <= 0.0;
			}), roots.end());

		const auto alphaP(*min_element(roots.begin(), roots.end()));

		// osszehasonlitjuk az alphaK-val, ami az x1 + alphaK * deltax1 >= 0 osszefuggesbol jon
		// ha, deltax1 >= 0, akkor alphaK = 1
		// ha, deltax1 < 0, akkor alphaK = -x1/deltax1

		T alphaK = 1;
		if (deltax1 < 0.0)
		{
			alphaK = -x1 / deltax1;
		}

		// vegso alpha = min{alphaP, alphaK}
		return min(alphaP, alphaK);
	}

	template <typename T>
		requires std::floating_point<T>
	std::vector<std::shared_ptr<Matrix<T>>> SOCPSolver<T>::calculateW() const
	{
		// for each cone
		// w = P(x^1/2) * (P(x^1/2)*s)^(-1/2)
		// NT scaling point
		// w vector
		using namespace std;

		CONE_GROUP sqrtx;
		for (const auto& x : this->x_)
		{
			sqrtx.push_back(calculatePowerOf(x.get(), 0.5));
		}
		CONE_GROUP pMatrixSqrtx;
		for (const auto& x : sqrtx)
		{
			pMatrixSqrtx.push_back(calculatePMatrixOf(x.get()));
		}

		CONE_GROUP sMultipliedBypMatrix;
		for (size_t i(0); i < pMatrixSqrtx.size(); ++i)
		{
			sMultipliedBypMatrix.push_back(*pMatrixSqrtx.at(i) * this->s_.at(i));
		}

		CONE_GROUP sqrtSMultip;
		for (const auto& x : sMultipliedBypMatrix)
		{
			sqrtSMultip.push_back(calculatePowerOf(calculatePowerOf(x.get(), 0.5).get(), -1));
		}

		CONE_GROUP w;
		for (size_t i(0); i < pMatrixSqrtx.size(); ++i)
		{
			w.push_back(*pMatrixSqrtx.at(i) * sqrtSMultip.at(i));
		}

		return move(w);
	}

	template <typename T>
		requires std::floating_point<T>
	std::vector<std::shared_ptr<Matrix<T>>> SOCPSolver<T>::calculateV() const
	{
		// (P(w^(-1/2)) * x) / sqrt(mu)
		using namespace std;

		const auto w(calculateW());

		CONE_GROUP sqrtInvW;
		for (const auto& tw : w)
		{
			sqrtInvW.push_back(calculatePowerOf(calculatePowerOf(tw.get(), 0.5).get(), -1));
		}

		CONE_GROUP pOfW;
		for (const auto& tw : sqrtInvW)
		{
			pOfW.push_back(calculatePMatrixOf(tw.get()));
		}

		CONE_GROUP xMultipOfW;
		for (size_t i(0); i < pOfW.size(); ++i)
		{
			xMultipOfW.push_back(*pOfW.at(i) * this->x_.at(i));
		}

		for (size_t i(0); i < xMultipOfW.size(); ++i)
		{
			*xMultipOfW.at(i) /= sqrt(mu_);
		}

		return move(xMultipOfW);
	}

	template <typename T>
		requires std::floating_point<T>
	std::vector<std::shared_ptr<Matrix<T>>> SOCPSolver<T>::calculatePv() const
	{
		// ahol v = P(w^-1/2) * x / sqrt(mu) = v / sqrt(mu)
		// pv = v^-1 - v, klasszikus phi(t) = t fuggveny eseten

		using namespace std;

		const auto v(calculateV());

		CONE_GROUP vInv;
		for (const auto& tv : v)
		{
			vInv.push_back(calculatePowerOf(tv.get(), -1));
		}

		CONE_GROUP ret;
		for (size_t i(0); i < v.size(); ++i)
		{
			ret.push_back(*vInv.at(i) - v.at(i));
		}

		return move(ret);
	}

	template <typename T>
		requires std::floating_point<T>
	std::vector<std::shared_ptr<Matrix<T>>> SOCPSolver<T>::calculateA_() const
	{
		// A_ = sqrt(mu) * A * P(w^1/2)
		using namespace std;

		CONE_GROUP AMultipMu;
		for (const auto& A : cnstrs_)
		{
			AMultipMu.push_back(*A * sqrt(mu_));
		}

		const auto w(calculateW());

		CONE_GROUP sqrtW;
		for (const auto& tw : w)
		{
			sqrtW.push_back(calculatePowerOf(tw.get(), 0.5));
		}

		CONE_GROUP pOfW;
		for (const auto& tw : sqrtW)
		{
			pOfW.push_back(calculatePMatrixOf(tw.get()));
		}

		for (size_t i(0); i < AMultipMu.size(); ++i)
		{
			*AMultipMu.at(i) *= *pOfW.at(i);
		}

		return move(AMultipMu);
	}

	template <typename T>
		requires std::floating_point<T>
	std::shared_ptr<Matrix<T>> SOCPSolver<T>::calculateResidualB() const
	{
		using namespace std;
		// Ax - b

		const auto A(calculateHorizontallyCombinedMatrix({ m_ }, n_, cnstrs_));
		const auto& b(cnstrsObjs_);

		// we need to concatenate each x variable groups
		const auto combinedX(calculateCombinedVectors(this->n_, this->x_));

		return std::move(*(*A * *combinedX) - *b);
	}

	template <typename T>
		requires std::floating_point<T>
	std::vector<std::shared_ptr<Matrix<T>>> SOCPSolver<T>::calculateResidualC() const
	{
		using namespace std;
		// for each cone
		// AT * y + s - c

		CONE_GROUP ret;
		for (size_t i(0); i < cnstrs_.size(); ++i)
		{
			ret.push_back(*(*(*cnstrs_.at(i)->transpose() * y_) + s_.at(i)) - *objs_.at(i));
		}

		return ret;
	}

	template <typename T>
		requires std::floating_point<T>
	std::vector<std::shared_ptr<Matrix<T>>> SOCPSolver<T>::calculateResidualC_() const
	{
		// 1/sqrt(mu) * P(w)^1/2 * rc
		const auto rc(calculateResidualC());
		const auto w(calculateW());

		CONE_GROUP sqrtW;
		for (const auto& tw : w)
		{
			sqrtW.push_back(calculatePowerOf(tw.get(), 0.5));
		}

		CONE_GROUP pOfW;
		for (const auto& tw : sqrtW)
		{
			pOfW.push_back(calculatePMatrixOf(tw.get()));
		}

		for (const auto& m : pOfW)
		{
			*m *= 1 / sqrt(mu_);
		}

		for (size_t i(0); i < pOfW.size(); ++i)
		{
			*pOfW.at(i) *= rc.at(i);
		}

		return std::move(pOfW);
	}

	template <typename T>
		requires std::floating_point<T>
	bool SOCPSolver<T>::checkPrimalFeasibility() const
	{
		// ||rb|| / (1 + ||b||) < epsilon

		const auto rb(calculateResidualB());
		const auto b(this->problem_->getConstraintsObjectives());

		return rb->norm() / (1.0 + b->norm()) < this->epsilon_;
	}

	template <typename T>
		requires std::floating_point<T>
	bool SOCPSolver<T>::checkDualFeasibility() const
	{
		// ||rc|| / (1 + ||c||)
		const auto trc(calculateResidualC());
		const auto rc(calculateCombinedVectors(n_, trc));
		const auto c(calculateCombinedVectors(n_, objs_));

		return rc->norm() / (1.0 + c->norm()) < this->epsilon_;
	}

	template <typename T>
		requires std::floating_point<T>
	bool SOCPSolver<T>::checkInfeasibility() const
	{
		auto hr(false);

		if (currIter_ > maxIters_) hr = true;
		if (!(checkMu() && checkPrimalFeasibility() && checkDualFeasibility())) hr = true;

		const auto ctx(*this->problem_->getObjectives()->transpose() * *calculateCombinedVectors(n_, x_));
		const auto bty(*this->problem_->getConstraintsObjectives()->transpose() * *y_);

		if (containsNaN(ctx.get()) || containsNaN(bty.get())) hr = true;

		return hr;
	}

	template <typename T>
		requires std::floating_point<T>
	SolutionStatus SOCPSolver<T>::internalSolver()
	{
		using namespace std;

		auto hr(SolutionStatus::FEASIBLE);
		const MatrixFactory<T> matrixFactory;
		this->currIter_ = 0;

		do
		{
			++this->currIter_;

			// calc barrier param
			this->mu_ = calculateMu();

			// calc residuals
			const auto rb(calculateResidualB());
			const auto rc_(calculateCombinedVectors(n_, calculateResidualC_()));

			// calc important metrics
			const auto pv(calculateCombinedVectors(n_, calculatePv()));
			const auto A_(calculateHorizontallyCombinedMatrix({ m_ }, n_, calculateA_()));
			const auto A_T(A_->transpose());

			// dy
			shared_ptr<Matrix<T>> dy;
			{
				// (m x nj) * (nj x m) * (m x 1) = -(m x nj) * (nj x 1) - (m x 1) - (m x nj) * (nj x 1)
				// (m x m) * (m x 1) = -(m x 1) - (m x 1) - (m x 1)
				// A_ * A_T * dy = -A_ * rc_ - rb - A_ * pv

				const auto lhs = *A_ * *A_T;
				const auto rhs = *(*(*(*A_ * -1) * *rc_) - *rb) - *(*A_ * *pv);

				dy = lhs->solve(rhs);
			}

			// ds
			shared_ptr<Matrix<T>> ds;
			{
				// A_T * dy + ds = -rc_ => ds = -rc_ - A_T * dy
				ds = *(*rc_ * -1) - *(*A_T * *dy);
			}

			//dx
			shared_ptr<Matrix<T>> dx;
			{
				// dx + ds = pv => dx = pv - ds
				dx = *pv - *ds;
			}

			// we need to separate individual cone directions
			// in case of dy we do not need to separate anything

			// deltay
			shared_ptr<Matrix<T>> deltay;
			{
				deltay = *dy * mu_;
			}

			CONE_GROUP deltax;
			CONE_GROUP deltas;

			for (size_t lastPoz(0); const auto & n : n_)
			{
				deltax.push_back(dx->block(lastPoz, 0, lastPoz + n - 1, 0));
				deltas.push_back(ds->block(lastPoz, 0, lastPoz + n - 1, 0));
				lastPoz += n;
			}

			const auto w(calculateW());
			CONE_GROUP pofsqrtw;
			CONE_GROUP pofinvsqrtw;
			for (const auto& tw : w)
			{
				const auto sqrtw(calculatePowerOf(tw.get(), 0.5));
				const auto invsqrtw(calculatePowerOf(sqrtw.get(), -1));

				pofsqrtw.push_back(calculatePMatrixOf(sqrtw.get()));
				pofinvsqrtw.push_back(calculatePMatrixOf(invsqrtw.get()));
			}

			// deltax and deltas
			for (size_t i(0); i < n_.size(); ++i)
			{
				deltax.at(i) = *(*pofsqrtw.at(i) * sqrt(mu_)) * *deltax.at(i);
				deltas.at(i) = *(*pofinvsqrtw.at(i) * sqrt(mu_)) * *deltas.at(i);

				alphaPrimal_.push_back(calculateAlpha(x_.at(i).get(), deltax.at(i).get()));
				alphaDual_.push_back(calculateAlpha(s_.at(i).get(), deltas.at(i).get()));
			}

			const auto a1(*ranges::min_element(alphaPrimal_));
			const auto a2(*ranges::min_element(alphaDual_));

			for (size_t i(0); i < n_.size(); ++i)
			{
				*deltax.at(i) *= (a1 * this->rho_);
				*x_.at(i) += *deltax.at(i);

				*deltas.at(i) *= (a2 * this->rho_);
				*s_.at(i) += *deltas.at(i);
			}

			*y_ += *(*deltay * (a2 * this->rho_));

			{
				const auto ctx(*this->problem_->getObjectives()->transpose() * *calculateCombinedVectors(n_, x_));
				const auto bty(*this->problem_->getConstraintsObjectives()->transpose() * *y_);

				this->writer_->writeIteration(this->currIter_, { ctx->get(0, 0), bty->get(0, 0), (*ctx - *bty)->get(0, 0) });
			}
		} while (!checkIsTermination());

		// check for final infeasibility
		if (checkInfeasibility()) hr = SolutionStatus::INFEASIBLE;

		return hr;
	}

	/*
	 * For some problem types there is a need for certain preparations.
	 * E.g., in case of Markowitz basic problems we need to calculate the Cholesky factorization for the covariance matrix.
	 */
	template <typename T> requires std::floating_point<T>
	void SOCPSolver<T>::internalSolverPreparations(std::shared_ptr<Matrix<T>> A, std::shared_ptr<Matrix<T>> b, std::shared_ptr<Matrix<T>> c)
	{
		using namespace std;

		if (this->problem_->getObjectiveDirection() == ObjectiveDirection::MAXIMIZE)
		{
			*c *= -1;
		}

		// complicated switches for specific initializators on which we can make assumptions what type of problem is used

		if (this->init_ == this->INITS.at("markowitz"))
		{
			const auto firstRow(A->block(0, 0, 1, A->getCols() - 1)->getValues());
			const size_t ones(count(firstRow->begin(), firstRow->end(), 1));

			// if Markowitz problem we need to replace the covariance matrix with Cholesky factorization transpose
			auto covariance(A->block(1, 0, A->getRows() - 2, ones - 1));

			// bit of a cheat
			const shared_ptr<Matrix<T>> G(calculateCholesky(covariance.get()));

			// we replace covariance with transpose of G
			A->block(1, 0, A->getRows() - 2, ones - 1, G->transpose());
		}
	}

	template <typename T> requires std::floating_point<T>
	void SOCPSolver<T>::setStartingPointInitializator(const std::string& init)
	{
		std::string temp(init);
		std::ranges::transform(temp, temp.begin(), [](const unsigned char c) { return std::tolower(c); });

		if (!INITS.contains(temp)) throw SolverException("Unsupported starting point initializator type.");

		init_ = INITS.at(temp);
	}

	template <typename T>
		requires std::floating_point<T>
	SolutionStatus SOCPSolver<T>::solve()
	{
		using namespace std;

		if (this->problem_ == nullptr)
		{
			throw SolverException("A problem was never set, try setting the problem first before trying to solve.");
		}

		if (this->problem_->getConstraints() == nullptr || this->problem_->getObjectives() == nullptr || this->problem_
			->getConstraintsObjectives() == nullptr)
		{
			throw SolverException(
				"Cannot solve optimization problem without correctly setting the constraints, objectives and constraint objectives.");
		}

		if (this->writer_ == nullptr)
		{
			throw SolverException("No writer is set.");
		}

		this->writer_->setIterationHeaders({ "cTx", "bTy", "Duality Gap" });

		const MatrixFactory<T> matrixFactory;
		vector<size_t> coneSizes;

		// we need to separate individual cones
		try {
			m_ = this->problem_->getConstraintsObjectives()->getRows();

			const auto* pr(dynamic_cast<SOCPProblem<T>*>(this->problem_.get()));

			if (!pr) throw SolverException("Failed to convert given problem to SOCP problem.");

			coneSizes = pr->getConeSizes();

			n_ = coneSizes;

			const auto& A(this->problem_->getConstraints());
			const auto& b(this->problem_->getConstraintsObjectives());
			auto c(this->problem_->getObjectives());

			internalSolverPreparations(A, b, c);

			// set A
			size_t lastPoz(0);
			for (const auto& size : coneSizes)
			{
				cnstrs_.push_back(A->block(0, lastPoz, m_ - 1, lastPoz + size - 1));
				lastPoz += size;
			}

			cnstrsObjs_ = b;

			lastPoz = static_cast<size_t>(0);
			for (const auto& size : coneSizes)
			{
				objs_.push_back(c->block(lastPoz, 0, lastPoz + size - 1, 0));
				lastPoz += size;
			}
		}
		catch (const exception& e)
		{
			throw SolverException(e.what());
		}

		// one set of x and s for each cone
		for (const auto& n : n_)
		{
			x_.push_back(matrixFactory.createMatrix(n, 1));
			s_.push_back(matrixFactory.createMatrix(n, 1));
		}
		y_ = matrixFactory.createMatrix(m_, 1); // only one

		init_->initialize(x_, y_, s_);

		const auto status(internalSolver());

		const auto cx(calculateCombinedVectors(n_, x_));
		const auto cs(calculateCombinedVectors(n_, s_));
		auto opt((*this->problem_->getObjectives()->transpose() * *cx)->get(0, 0));

		if (this->problem_->getObjectiveDirection() == ObjectiveDirection::MAXIMIZE)
		{
			opt *= -1;
		}

		this->solution_ = make_shared<SOCPSolution<T>>(status, opt, cx, y_, cs, coneSizes);

		return status;
	}

	template <typename T>
		requires std::floating_point<T>
	void SOCPSolver<T>::ClassicInitializator::initialize(std::vector<std::shared_ptr<Matrix<T>>> x, std::shared_ptr<Matrix<T>> y, std::vector<std::shared_ptr<Matrix<T>>> s)
	{
		using namespace std;

		if (x.empty() || y == nullptr || s.empty())
			throw SolverException(
				"Wrong arguments for initializer, matrices have to be allocated before initializing them.");

		for (const auto& tx : x)
		{
			tx->setValues(vector<T>(tx->getRows(), 0), tx->getRows(), 1);
			tx->set(0, 0, 1);
		}

		y->setValues(vector<T>(y->getRows(), 0), y->getRows(), 1);

		for (const auto& ts : s)
		{
			ts->setValues(vector<T>(ts->getRows(), 0), ts->getRows(), 1);
			ts->set(0, 0, 1);
		}
	}

	template <typename T>
		requires std::floating_point<T>
	void SOCPSolver<T>::Classic2Initializator::initialize(std::vector<std::shared_ptr<Matrix<T>>> x, std::shared_ptr<Matrix<T>> y, std::vector<std::shared_ptr<Matrix<T>>> s)
	{
		using namespace std;

		if (x.empty() || y == nullptr || s.empty())
			throw SolverException(
				"Wrong arguments for initializer, matrices have to be allocated before initializing them.");

		for (const auto& tx : x)
		{
			tx->setValues(vector<T>(tx->getRows(), 0), tx->getRows(), 1);
			tx->set(0, 0, 1);
		}

		y->setValues(vector<T>(y->getRows(), 1), y->getRows(), 1);

		for (const auto& ts : s)
		{
			ts->setValues(vector<T>(ts->getRows(), 0), ts->getRows(), 1);
			ts->set(0, 0, 2);
		}
	}

	template <typename T>
		requires std::floating_point<T>
	void SOCPSolver<T>::Classic3Initializator::initialize(std::vector<std::shared_ptr<Matrix<T>>> x, std::shared_ptr<Matrix<T>> y, std::vector<std::shared_ptr<Matrix<T>>> s)
	{
		using namespace std;

		if (x.empty() || y == nullptr || s.empty())
			throw SolverException(
				"Wrong arguments for initializer, matrices have to be allocated before initializing them.");

		for (const auto& tx : x)
		{
			tx->setValues(vector<T>(tx->getRows(), 0), tx->getRows(), 1);
			tx->set(0, 0, 1);
		}

		y->setValues(vector<T>(y->getRows(), 1), y->getRows(), 1);

		for (const auto& ts : s)
		{
			ts->setValues(vector<T>(ts->getRows(), 0), ts->getRows(), 1);
			ts->set(0, 0, 1);
		}
	}

	template <typename T>
		requires std::floating_point<T>
	void SOCPSolver<T>::Classic4Initializator::initialize(std::vector<std::shared_ptr<Matrix<T>>> x, std::shared_ptr<Matrix<T>> y, std::vector<std::shared_ptr<Matrix<T>>> s)
	{
		using namespace std;

		if (x.empty() || y == nullptr || s.empty())
			throw SolverException(
				"Wrong arguments for initializer, matrices have to be allocated before initializing them.");

		for (const auto& tx : x)
		{
			tx->setValues(vector<T>(tx->getRows(), 0), tx->getRows(), 1);
			tx->set(0, 0, 2);
			tx->set(1, 0, 1);
		}

		y->setValues(vector<T>(y->getRows(), 0), y->getRows(), 1);

		for (const auto& ts : s)
		{
			ts->setValues(vector<T>(ts->getRows(), 0), ts->getRows(), 1);
			ts->set(0, 0, 2);
			ts->set(1, 0, 1);
		}
	}

	template <typename T>
		requires std::floating_point<T>
	void SOCPSolver<T>::MarkowitzInitializator::initialize(std::vector<std::shared_ptr<Matrix<T>>> x, std::shared_ptr<Matrix<T>> y, std::vector<std::shared_ptr<Matrix<T>>> s)
	{
		using namespace std;

		if (x.empty() || y == nullptr || s.empty())
			throw SolverException(
				"Wrong arguments for initializer, matrices have to be allocated before initializing them.");

		for (const auto& tx : x)
		{
			tx->setValues(vector<T>(tx->getRows(), 0), tx->getRows(), 1);
			tx->set(0, 0, 1.0 / tx->getRows());
		}

		y->setValues(vector<T>(y->getRows(), 0), y->getRows(), 1);

		for (const auto& ts : s)
		{
			ts->setValues(vector<T>(ts->getRows(), 0), ts->getRows(), 1);
			ts->set(0, 0, 1.0 / ts->getRows());
		}
	}
}
