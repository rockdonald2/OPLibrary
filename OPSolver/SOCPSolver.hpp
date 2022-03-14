#pragma once

#include <format>
#include <numbers>
#include <vector>
#include <string>
#include <algorithm>

#include "Solver.hpp"
#include "MatrixFactory.hpp"
#include "MatrixUtils.hpp"
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
		/**
		 * \brief Initializator with Strategy pattern to initialize solution vectors.
		 */
		class Initializator
		{
		public:
			Initializator() = default;
			virtual ~Initializator() = default;

			virtual void initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s) = 0;

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
			void initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s) override;
		};

		/**
		 * \brief Classic initializator, where for x0 vector first value is 1, s0 vector first value is 2, all other 0; for y all the values are 1.
		 */
		class Classic2Initializator final : public Initializator
		{
		public:
			void initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s) override;
		};

		/**
		 * \brief Classic initializator, where for x0 and s0 vector their first value is 1, all other 0; for y all the values are 1.
		 */
		class Classic3Initializator final : public Initializator
		{
		public:
			void initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s) override;
		};

		/**
		 * \brief Classic initializator, where for x0 and s0 vector their first value is 2, second value is 1, all other 0; for y all the values are 0.
		 */
		class Classic4Initializator final : public Initializator
		{
		public:
			void initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s) override;
		};

		class Classic5Initializator final : public Initializator
		{
		public:
			void initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s) override;
		};

		class ManualInitializator final : public Initializator
		{
		public:
			void initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s) override;
		};

		long double epsilon_;
		long double tau_;
		long double alphaPrimal_;
		long double alphaDual_;
		long double mu_;
		long double beta_;
		long double rho_;
		long double sigma_;

		size_t n_;
		size_t nn_;

		size_t currIter_;
		size_t maxIters_;

		std::unique_ptr<Initializator> init_;

		std::unique_ptr<Matrix<T>> x_;
		std::unique_ptr<Matrix<T>> y_;
		std::unique_ptr<Matrix<T>> s_;

		SolutionStatus status_;

		[[nodiscard]] bool checkIsTermination() const;

		static T calculateEigenMinOf(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "Eigenvalue min can only be calculated for vectors.");

			const unique_ptr<Matrix<T>> tmpBlock(vec->block(1, 0, vec->getRows() - 1, 0));

			return vec->get(0, 0) - tmpBlock->norm();
		}

		static T calculateEigenMaxOf(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "Eigenvalue max can only be calculated for vectors.");

			const unique_ptr<Matrix<T>> tmpBlock(vec->block(1, 0, vec->getRows() - 1, 0));

			return vec->get(0, 0) + tmpBlock->norm();
		}

		static std::unique_ptr<Matrix<T>> calculatePMatrixOf(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "P matrix can only be calculated for vectors.");

			const MatrixFactory<T> factory(MatrixType::DENSE);
			const auto n(vec->getRows());
			auto retMatrix(factory.createMatrix(n, n));

			// norm^2
			retMatrix->set(0, 0, pow(vec->norm(), 2));

			const auto x1(vec->get(0, 0));
			const unique_ptr<Matrix<T>> x2n(vec->block(1, 0, vec->getRows() - 1, 0));
			const unique_ptr<Matrix<T>> x2nT(x2n->transpose());

			// 2 * (x1 * x2:n)
			retMatrix->block(1, 0, retMatrix->getRows() - 1, 0, *(*x2n * x1) * 2);

			// 2 * (x1 * x2T:n)
			retMatrix->block(0, 1, 0, retMatrix->getCols() - 1, *(*x2nT * x1) * 2);

			// det(x) * En-1 + 2 * x2n * x2T:n
			// det(x) = x1^2 - norm(x2:n)^2
			const auto detx(calculateEigenMaxOf(vec) * calculateEigenMinOf(vec));

			const unique_ptr<Matrix<T>> E(factory.createMatrix(n - 1, n - 1));
			E->setValues(vector<T>(static_cast<size_t>(pow((n - 1), 2)), 0), n - 1, n - 1);
			E->setDiagonalValues(vector<T>(n - 1, 1));

			const unique_ptr<Matrix<T>> helperTriangularMatrix(factory.createMatrix(n - 1, n - 1));
			helperTriangularMatrix->setValues(vector<T>(static_cast<size_t>(pow(n - 1, 2)), 0), n - 1, n - 1);

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

		static std::unique_ptr<Matrix<T>> calculateC1(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "C1 vector can only be calculated for vectors.");

			const auto n(vec->getRows());

			const MatrixFactory<T> factory(MatrixType::DENSE);
			auto c1(factory.createMatrix(n, 1));
			c1->set(0, 0, 0.5);

			const unique_ptr<Matrix<T>> x2n(vec->block(1, 0, n - 1, 0));
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

		static std::unique_ptr<Matrix<T>> calculateC2(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "C2 vector can only be calculated for vectors.");

			const auto n(vec->getRows());

			const MatrixFactory<T> factory(MatrixType::DENSE);
			auto c2(factory.createMatrix(n, 1));
			c2->set(0, 0, 0.5);

			const unique_ptr<Matrix<T>> x2n(vec->block(1, 0, n - 1, 0));
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

		static std::unique_ptr<Matrix<T>> calculatePowerOf(Matrix<T>* vec, const double power)
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

		static std::unique_ptr<Matrix<T>> oOperation(Matrix<T>* lhs, Matrix<T>* rhs)
		{
			assert(lhs->getCols() == 1 && rhs->getCols() == 1 && "Bilinear operator o can only be calculated for vectors.");
			assert(lhs->getRows() == rhs->getRows() && "Bilinear operator o can only be calculated for same dimensional vectors.");

			using namespace std;

			const auto n(lhs->getRows());

			const MatrixFactory<T> factory;

			unique_ptr<Matrix<T>> ret(factory.createMatrix(n, 1));

			ret->set(0, 0, (*lhs->transpose() * *rhs)->get(0, 0));

			const auto x1(lhs->get(0, 0));
			const auto s1(rhs->get(0, 0));

			for (size_t i = 1; i < n; ++i)
			{
				ret->set(i, 0, x1 * rhs->get(i, 0) + s1 * lhs->get(i, 0));
			}

			return ret;
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

		[[nodiscard]] T calculateMu() const;
		[[nodiscard]] T calculateAlpha(const Matrix<T>* vec, const Matrix<T>* delta) const;
		[[nodiscard]] std::unique_ptr<Matrix<T>> calculateW() const;
		[[nodiscard]] std::unique_ptr<Matrix<T>> calculateV() const;
		[[nodiscard]] std::unique_ptr<Matrix<T>> calculatePv() const;

		[[nodiscard]] std::unique_ptr<Matrix<T>> calculateA_() const;

		[[nodiscard]] T distanceFromMuCenter() const;

		[[nodiscard]] std::unique_ptr<Matrix<T>> calculateResidualB() const;
		[[nodiscard]] std::unique_ptr<Matrix<T>> calculateResidualC() const;
		[[nodiscard]] std::unique_ptr<Matrix<T>> calculateResidualC_() const;

		[[nodiscard]] bool checkPrimalFeasibility() const;
		[[nodiscard]] bool checkDualFeasibility() const;

		SolutionStatus optimizedSolver();

	public:
		SOCPSolver() : epsilon_(1.0e-6), tau_(2.0), alphaPrimal_(1.0 / 10),
			alphaDual_(1.0 / 10), mu_(1.0), beta_(1.0 / 2), rho_(0.98), sigma_(1.0 / 10), n_(0), nn_(0),
			currIter_(0),
			maxIters_(3000), init_(new Classic5Initializator()), x_(nullptr),
			y_(nullptr), s_(nullptr), status_(SolutionStatus::FEASIBLE)
		{
			using namespace std;

			this->INITIALIZABLE_ARGS = { "epsilon", "tau", "mu", "beta", "rho", "sigma" };
			this->INITIALIZATOR.insert(make_pair<string, long double*>("epsilon", &epsilon_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("tau", &tau_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("mu", &mu_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("beta", &beta_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("rho", &rho_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("sigma", &sigma_));
		}

		void setProblem(std::shared_ptr<Problem<T>> problem) override;
		void setProblem(Problem<T>* problem) override;
		[[nodiscard]] std::shared_ptr<Problem<T>> getProblem() const override;

		SolutionStatus solve() override;
		[[nodiscard]] std::shared_ptr<Solution<T>> getSolution() override;
		[[nodiscard]] std::string getStatus() override;

		void setWriter(Writer<T>* writer) override;
		void setWriter(std::shared_ptr<Writer<T>> writer) override;
		[[nodiscard]] std::shared_ptr<Writer<T>> getWriter() const override;
	};

	template <typename T>
		requires std::floating_point<T>
	bool SOCPSolver<T>::checkIsTermination() const
	{
		return
			currIter_ <= maxIters_ &&
			(mu_ >= epsilon_) &&
			isInCone(x_.get()) && isInCone(s_.get()) &&
			!checkPrimalFeasibility() && !checkDualFeasibility();
	}

	template <typename T>
		requires std::floating_point<T>
	T SOCPSolver<T>::calculateMu() const
	{
		return (*(*this->x_->transpose() * (2.0 * this->sigma_ / static_cast<long double>(this->n_))) * *this->s_)->get(0, 0);
	}

	template <typename T>
		requires std::floating_point<T>
	T SOCPSolver<T>::calculateAlpha(const Matrix<T>* vec, const Matrix<T>* delta) const
	{
		using namespace std;

		// kiszamoljuk az A-t, B-t, C-t

		// deltax1 negyzete - deltax2n negyzetosszege
		// A = deltaX1^2 - sum(deltaX2:n^2)

		// x1 * deltax1 - x2n
		// B = 2 * (x1 * deltax1 - sum(xi * deltaxi))

		// x1 negyzete - x2n negyzetosszege
		// C = x1^2 - sum(xi^2)

		const auto deltax1(delta->get(0, 0));
		const auto deltax2n(delta->block(1, 0, delta->getRows() - 1, 0)->getValues());

		T A;
		{
			T sum(0);
			for_each(deltax2n->begin(), deltax2n->end(), [&sum](T n)
				{
					sum += pow(n, 2);
				});

			A = pow(deltax1, 2) - sum;
		}

		const auto x1(vec->get(0, 0));
		const auto x2n(vec->block(1, 0, vec->getRows() - 1, 0)->getValues());

		T B;
		{
			T sum(0);
			for (size_t i = 0; i < x2n->size(); ++i)
			{
				sum += (x2n->at(i) * deltax2n->at(i));
			}

			B = x1 * deltax1 - sum;
			B *= 2;
		}

		T C;
		{
			T sum(0);
			for_each(x2n->begin(), x2n->end(), [&sum](T n)
				{
					sum += pow(n, 2);
				});

			C = pow(x1, 2) - sum;
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
			roots = solvePolynomial<T>(vector<T>({ A, B, C }));
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
	std::unique_ptr<Matrix<T>> SOCPSolver<T>::calculateW() const
	{
		// w = P(x^1/2) * (P(x^1/2)*s)^(-1/2)
		// NT skalazasi pont
		// w vektor
		using namespace std;

		const unique_ptr<Matrix<T>> sqrtx(calculatePowerOf(this->x_.get(), 0.5));
		const unique_ptr<Matrix<T>> pMatrixSqrtx(calculatePMatrixOf(sqrtx.get()));
		const unique_ptr<Matrix<T>> sMultipliedBypMatrix(*pMatrixSqrtx * *this->s_);
		const unique_ptr<Matrix<T>> sqrtSMultiplication(
			calculatePowerOf(calculatePowerOf(sMultipliedBypMatrix.get(), 0.5).get(), -1));

		return move(*pMatrixSqrtx * *sqrtSMultiplication);
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> SOCPSolver<T>::calculateV() const
	{
		// (P(w^(-1/2)) * x) / sqrt(mu)
		using namespace std;

		const auto w(calculateW());
		const auto sqrtInverseW(calculatePowerOf(calculatePowerOf(w.get(), 0.5).get(), -1));
		const auto pOfW(calculatePMatrixOf(sqrtInverseW.get()));
		// a gyokvonas feleslegesnek tunik
		const auto xMultipliedpOfW(*pOfW * *this->x_);

		return move(*xMultipliedpOfW / sqrt(this->mu_));
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> SOCPSolver<T>::calculatePv() const
	{
		// ahol v = P(w^-1/2) * x / sqrt(mu) = v / sqrt(mu)
		// pv = v^-1 - v, klasszikus phi(t) = t fuggveny eseten
		using namespace std;

		const auto v(calculateV());
		const auto vInverse(calculatePowerOf(v.get(), -1));

		return move(*vInverse - *v);
	}

	template <typename T>
		requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> SOCPSolver<T>::calculateA_() const
	{
		using namespace std;

		// A_ = sqrt(mu) * A * P(w^1/2)

		const auto AMultipliedMu(*this->problem_->getConstraints() * sqrt(this->mu_));
		const auto w(calculateW());
		const auto sqrtW(calculatePowerOf(w.get(), 0.5));
		const auto pOfW(calculatePMatrixOf(sqrtW.get()));

		return move(*AMultipliedMu * *pOfW);
	}

	template <typename T> requires std::floating_point<T>
	T SOCPSolver<T>::distanceFromMuCenter() const
	{
		// a mu centrumtol vett tavolsag egyenlo ||pv||F / 2,
		// azonban a || ||F = sqrt(2) * || ||

		const auto pv(calculatePv());
		return sqrt(2) * pv->norm() / 2;
	}

	template <typename T> requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> SOCPSolver<T>::calculateResidualB() const
	{
		// Ax - b
		const auto& A(this->problem_->getConstraints());
		const auto& b(this->problem_->getConstraintsObjectives());

		return std::move(*(*A * *x_) - *b);
	}

	template <typename T> requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> SOCPSolver<T>::calculateResidualC() const
	{
		// AT * y + s - c

		const auto AT(this->problem_->getConstraints()->transpose());
		const auto& c(this->problem_->getObjectives());

		return std::move(*(*(*AT * y_) + *s_) - *c);
	}

	template <typename T> requires std::floating_point<T>
	std::unique_ptr<Matrix<T>> SOCPSolver<T>::calculateResidualC_() const
	{
		// 1/sqrt(mu) * P(w)^1/2 * rc

		const auto mu(1.0 / sqrt(mu_));
		const auto w(calculatePowerOf(calculateW().get(), 0.5));
		const auto pOfW(calculatePMatrixOf(w.get()));
		const auto rc(calculateResidualC());

		return std::move(*(*pOfW * mu) * *rc);
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
		const auto rc(calculateResidualC());
		const auto c(this->problem_->getObjectives());

		return rc->norm() / (1.0 + c->norm()) < this->epsilon_;
	}

	template <typename T>
		requires std::floating_point<T>
	SolutionStatus SOCPSolver<T>::optimizedSolver()
	{
		using namespace std;

		auto hr(SolutionStatus::FEASIBLE);
		const MatrixFactory<T> matrixFactory;
		this->writer_->setIterationHeaders({ "cTx", "bTy", "Duality Gap" });
		this->currIter_ = 0;

		do
		{
			++this->currIter_;

			this->mu_ = calculateMu();

			const auto rb(calculateResidualB());
			const auto pv(calculatePv());
			const auto A_(calculateA_());
			const auto A_T(A_->transpose());
			const auto rc_(calculateResidualC_());

			// dy
			unique_ptr<Matrix<T>> dy;
			{
				// A_ * A_T * dy = -A_ * rc_ - rb - A_ * pv

				const auto lhs = *A_ * *A_T;
				const auto rhs = *(*(*(*A_ * -1) * *rc_) - *rb) - *(*A_ * *pv);

				dy = lhs->solve(rhs);
			}

			// ds
			unique_ptr<Matrix<T>> ds;
			{
				// A_T * dy + ds = -rc_ => ds = -rc_ - A_T * dy
				ds = *(*rc_ * -1) - *(*A_T * *dy);
			}

			//dx
			unique_ptr<Matrix<T>> dx;
			{
				// dx + ds = pv => dx = pv - ds
				dx = *pv - *ds;
			}

			// deltay
			unique_ptr<Matrix<T>> deltay;
			{
				deltay = *dy * mu_;
			}

			const auto sqrtw(calculatePowerOf(calculateW().get(), 0.5));

			// deltax
			unique_ptr<Matrix<T>> deltax;
			{
				deltax = *(*calculatePMatrixOf(sqrtw.get()) * sqrt(mu_)) * *dx;
			}

			// deltas
			unique_ptr<Matrix<T>> deltas;
			{
				const auto invsqrtw(calculatePowerOf(sqrtw.get(), -1));

				deltas = *(*calculatePMatrixOf(invsqrtw.get()) * sqrt(mu_)) * *ds;
			}

			this->alphaPrimal_ = calculateAlpha(x_.get(), deltax.get());
			this->alphaDual_ = calculateAlpha(s_.get(), deltas.get());

			*x_ += *(*deltax * (this->alphaPrimal_ * this->rho_));
			*s_ += *(*deltas * (this->alphaDual_ * this->rho_));
			*y_ += *(*deltay * (this->alphaDual_ * this->rho_));

			this->writer_->writeIteration(this->currIter_,
				{ (*this->problem_->getObjectives()->transpose() * *this->x_)->get(0, 0),
				(*this->problem_->getConstraintsObjectives()->transpose() * *this->y_)->get(0, 0),
				(*this->x_->transpose() * *this->s_)->get(0, 0) });
		} while (checkIsTermination());

		if (currIter_ > maxIters_) hr = SolutionStatus::UNFEASIBLE;
		if (status_ == SolutionStatus::FEASIBLE &&
			containsNaN(x_.get()) || containsNaN(y_.get()) || containsNaN(s_.get())) hr = SolutionStatus::UNFEASIBLE;

		return hr;
	}

	template <typename T>
		requires std::floating_point<T>
	void SOCPSolver<T>::setProblem(std::shared_ptr<Problem<T>> problem)
	{
		this->problem_ = problem;
	}

	template <typename T>
		requires std::floating_point<T>
	void SOCPSolver<T>::setProblem(Problem<T>* problem)
	{
		this->problem_ = std::shared_ptr<Problem<T>>(problem);
	}

	template <typename T>
		requires std::floating_point<T>
	std::shared_ptr<Problem<T>> SOCPSolver<T>::getProblem() const
	{
		return this->problem_;
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
			throw SolverException("No writer was set, aborting.");
		}

		const auto n(this->problem_->getObjectives()->getRows());
		this->n_ = n;

		const auto nn(this->problem_->getConstraintsObjectives()->getRows());
		this->nn_ = nn;

		const MatrixFactory<T> matrixFactory;

		this->x_ = matrixFactory.createMatrix(n, 1);
		this->y_ = matrixFactory.createMatrix(nn, 1);
		this->s_ = matrixFactory.createMatrix(n, 1);

		this->init_->initialize(x_.get(), y_.get(), s_.get());

		this->status_ = optimizedSolver();

		this->solution_ = make_shared<Solution<T>>(Solution<T>(x_, y_, s_));

		return this->status_;
	}

	template <typename T>
		requires std::floating_point<T>
	std::shared_ptr<Solution<T>> SOCPSolver<T>::getSolution()
	{
		return this->solution_;
	}

	template <typename T> requires std::floating_point<T>
	std::string SOCPSolver<T>::getStatus()
	{
		switch (status_)
		{
		case SolutionStatus::FEASIBLE:
			return "feasible";
		case SolutionStatus::UNFEASIBLE:
			return "unfeasible";
		}

		throw SolverException("Unknown final status.");
	}

	template <typename T> requires std::floating_point<T>
	void SOCPSolver<T>::setWriter(Writer<T>* writer)
	{
		this->writer_ = std::shared_ptr<Writer<T>>(writer);
	}

	template <typename T> requires std::floating_point<T>
	void SOCPSolver<T>::setWriter(std::shared_ptr<Writer<T>> writer)
	{
		this->writer_ = writer;
	}

	template <typename T> requires std::floating_point<T>
	std::shared_ptr<Writer<T>> SOCPSolver<T>::getWriter() const
	{
		return this->writer_;
	}

	template <typename T>
		requires std::floating_point<T>
	void SOCPSolver<T>::ClassicInitializator::initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s)
	{
		using namespace std;

		if (x == nullptr || y == nullptr || s == nullptr)
			throw SolverException(
				"Wrong arguments for initializer, matrices have to be allocated before initializing them.");

		const auto rows(x->getRows());

		x->setValues(vector<T>(x->getRows(), 0), x->getRows(), 1);
		x->set(0, 0, 1);

		y->setValues(vector<T>(y->getRows(), 0), y->getRows(), 1);

		s->setValues(vector<T>(s->getRows(), 0), s->getRows(), 1);
		s->set(0, 0, 1);
	}

	template <typename T>
		requires std::floating_point<T>
	void SOCPSolver<T>::Classic2Initializator::initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s)
	{
		using namespace std;

		if (x == nullptr || y == nullptr || s == nullptr)
			throw SolverException(
				"Wrong arguments for initializer, matrices have to be allocated before initializing them.");

		x->setValues(vector<T>(x->getRows(), 0), x->getRows(), 1);
		x->set(0, 0, 1);

		y->setValues(vector<T>(y->getRows(), 1), y->getRows(), 1);

		s->setValues(vector<T>(s->getRows(), 0), s->getRows(), 1);
		s->set(0, 0, 2);
	}

	template <typename T>
		requires std::floating_point<T>
	void SOCPSolver<T>::Classic3Initializator::initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s)
	{
		using namespace std;

		if (x == nullptr || y == nullptr || s == nullptr)
			throw SolverException(
				"Wrong arguments for initializer, matrices have to be allocated before initializing them.");

		x->setValues(vector<T>(x->getRows(), 0), x->getRows(), 1);
		x->set(0, 0, 1);

		y->setValues(vector<T>(y->getRows(), 1), y->getRows(), 1);

		s->setValues(vector<T>(s->getRows(), 0), s->getRows(), 1);
		s->set(0, 0, 1);
	}

	template <typename T>
		requires std::floating_point<T>
	void SOCPSolver<T>::Classic4Initializator::initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s)
	{
		using namespace std;

		if (x == nullptr || y == nullptr || s == nullptr)
			throw SolverException(
				"Wrong arguments for initializer, matrices have to be allocated before initializing them.");

		x->setValues(vector<T>(x->getRows(), 0), x->getRows(), 1);
		x->set(0, 0, 2);
		x->set(1, 0, 1);

		y->setValues(vector<T>(y->getRows(), 0), y->getRows(), 1);

		s->setValues(vector<T>(s->getRows(), 0), s->getRows(), 1);
		s->set(0, 0, 2);
		s->set(1, 0, 1);
	}

	template <typename T> requires std::floating_point<T>
	void SOCPSolver<T>::Classic5Initializator::initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s)
	{
		using namespace std;

		if (x == nullptr || y == nullptr || s == nullptr)
			throw SolverException(
				"Wrong arguments for initializer, matrices have to be allocated before initializing them.");

		x->setValues(vector<T>(x->getRows(), 1), x->getRows(), 1);

		y->setValues(vector<T>(y->getRows(), 1), y->getRows(), 1);

		s->setValues(vector<T>(s->getRows(), 1), s->getRows(), 1);
	}

	template <typename T> requires std::floating_point<T>
	void SOCPSolver<T>::ManualInitializator::initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s)
	{
		using namespace std;

		if (x == nullptr || y == nullptr || s == nullptr)
			throw SolverException(
				"Wrong arguments for initializer, matrices have to be allocated before initializing them.");

		x->setValues(vector<T>(x->getRows(), 0), x->getRows(), 1);
		y->setValues(vector<T>(y->getRows(), 0), y->getRows(), 1);
		s->setValues(vector<T>(s->getRows(), 0), s->getRows(), 1);
	}
}
