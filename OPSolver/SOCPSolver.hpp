#pragma once

#include <format>
#include <numbers>
#include <vector>
#include <string>

#include "framework.h"
#include "Solver.hpp"
#include "MatrixFactory.hpp"
#include "Solution.hpp"
#include "SolverException.hpp"

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
		 * \brief Classic initializator, where for x0 vector first value is 1, s0 vector fist value is 2, all other 0; for y all the values are 1.
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

		long double theta_;
		long double epsilon_;
		long double tau_;
		long double alpha_;
		long double mu_;
		long double beta_;

		size_t n_;
		size_t nn_;

		std::unique_ptr<Initializator> init_;

		std::unique_ptr<Matrix<T>> x_;
		std::unique_ptr<Matrix<T>> y_;
		std::unique_ptr<Matrix<T>> s_;

		std::shared_ptr<Solution<T>> solution_;

		[[nodiscard]] bool checkIsTermination() const;

		// UTILITY functions for this method

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
			//const auto detx(pow(x1, 2) - pow(x2n->norm(), 2));
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

		[[nodiscard]] T calculateMu() const;
		[[nodiscard]] std::unique_ptr<Matrix<T>> calculateW() const;
		[[nodiscard]] std::unique_ptr<Matrix<T>> calculateV() const;
		[[nodiscard]] std::unique_ptr<Matrix<T>> calculatePv() const;

		[[nodiscard]] std::unique_ptr<Matrix<T>> calculateA_() const;

		[[nodiscard]] T distanceFromMuCenter() const;

		[[nodiscard]] std::unique_ptr<Matrix<T>> calculateResidualB() const;
		[[nodiscard]] std::unique_ptr<Matrix<T>> calculateResidualC() const;
		[[nodiscard]] std::unique_ptr<Matrix<T>> calculateResidualC_() const;

		SolutionStatus classicSolver();
		SolutionStatus optimizedSolver();

	public:
		SOCPSolver() : theta_(std::numbers::pi_v<long double> / 4), epsilon_(1.0e-8), tau_(2.0), alpha_(1.0 / 10),
			mu_(1.0), beta_(1.0 / 2), n_(0), nn_(0), init_(new Classic2Initializator()), x_(nullptr),
			y_(nullptr), s_(nullptr), solution_(nullptr)
		{
			using namespace std;

			this->INITIALIZABLE_ARGS = { "theta", "epsilon", "tau", "alpha", "mu" };
			this->INITIALIZATOR.insert(make_pair<string, long double*>("theta", &theta_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("epsilon", &epsilon_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("tau", &tau_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("alpha", &alpha_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("mu", &mu_));
		}

		void setProblem(std::shared_ptr<Problem<T>> problem) override;
		void setProblem(Problem<T>* problem) override;

		[[nodiscard]] std::shared_ptr<Problem<T>> getProblem() const override;

		SolutionStatus solve() override;

		[[nodiscard]] Solution<T> getSolution() override;
	};

	template <typename T>
		requires std::floating_point<T>
	bool SOCPSolver<T>::checkIsTermination() const
	{
		// mu < epsilon; az x minden erteke pozitiv; s minden erteke pozitiv; kupfeltetel ellenorzes
		// tehat megkell nezni azt, hogy lambda_min(x) > 0
		// lambda_min(s) > 0
		// kupfeltetel -- most eltekintunk rola, ennel az atlagosabb implementacional biztos, hogy bent maradunk

		return (mu_ >= epsilon_) && (calculateEigenMinOf(x_.get()) > 0) && (calculateEigenMinOf(s_.get()) > 0);
	}

	template <typename T>
		requires std::floating_point<T>
	T SOCPSolver<T>::calculateMu() const
	{
		return (*(*this->x_->transpose() * (2.0 * this->alpha_ / (this->n_ * 1.0))) * *this->s_)->get(0, 0);
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

	template <typename T> requires std::floating_point<T>
	SolutionStatus SOCPSolver<T>::classicSolver()
	{
		/*
				  * Egyenletrendszer:
				  *	[ A_	0		0 ]		[ dx ]		[ rb ]
				  *	[ 0		A_T		I ]	 *	[ dy ]	=	[ rc_ ]
				  *	[ I		0		I ]		[ ds ]		[ pv ]
				  */

		using namespace std;

		const MatrixFactory<T> matrixFactory(MatrixType::DENSE);

		const unique_ptr<Matrix<T>> I(matrixFactory.createMatrix());
		I->setValues(vector<T>(static_cast<size_t>(pow(n_, 2)), 0), n_, n_);
		I->setDiagonalValues(vector<T>(n_, 1));

		size_t iters(0);

		while (checkIsTermination())
		{
			++iters;

			LOG.info(format("{}. iteration ----------------------", iters));

			this->mu_ = calculateMu();

			const auto lhs(matrixFactory.createMatrix());
			const auto rhs(matrixFactory.createMatrix());

			const auto A_(calculateA_());
			const auto A_T(A_->transpose());
			const auto pv(calculatePv());

			// lhs
			{
				const auto rows(A_->getRows() + A_T->getRows() + pv->getRows());
				const auto cols(A_->getCols() + A_T->getCols() + pv->getCols() * n_);

				lhs->setValues(vector<T>(rows * cols, 0), rows, cols);

				lhs->block(0, 0, A_->getRows() - 1, A_->getCols() - 1, A_);

				lhs->block(A_->getRows(), A_->getCols(), A_->getRows() + A_T->getRows() - 1,
					A_->getCols() + A_T->getCols() - 1, A_T);

				lhs->block(A_->getRows() + A_T->getRows(), 0,
					A_->getRows() + A_T->getRows() + I->getRows() - 1, I->getCols() - 1, I);

				lhs->block(A_->getRows(), A_->getCols() + A_T->getCols(),
					A_->getRows() + A_T->getRows() - 1, A_->getCols() + A_T->getCols() + I->getCols() - 1, I);

				lhs->block(A_->getRows() + A_T->getRows(), A_->getCols() + A_T->getCols(),
					A_->getRows() + A_T->getRows() + I->getRows() - 1, A_->getCols() + A_T->getCols() + I->getCols() - 1, I);
			}

			// rhs
			{
				const auto rb(calculateResidualB());
				const auto rc_(calculateResidualC_());

				rhs->setValues(vector<T>(2 * n_ + nn_, 0), 2 * n_ + nn_, 1);

				rhs->block(0, 0, nn_ - 1, 0, rb);
				rhs->block(nn_, 0, nn_ + n_ - 1, 0, rc_);
				rhs->block(n_ + nn_, 0, 2 * n_ + nn_ - 1, 0, pv);
			}

			const auto sol = lhs->solve(rhs, DecompositionType::JACOBISVD);

			cout << *rhs << endl;
			cout << *lhs << endl;
			cout << *sol << endl;

			const auto sqrtw(calculatePowerOf(calculateW().get(), 0.5));
			const auto invsqrtw(calculatePowerOf(sqrtw.get(), -1));

			const auto dx(sol->block(0, 0, n_ - 1, 0));
			const auto dy(sol->block(n_, 0, n_ + nn_ - 1, 0));
			const auto ds(sol->block(n_ + nn_, 0, 2 * n_ + nn_ - 1, 0));

			cout << *dx << endl;
			cout << *dy << endl;
			cout << *ds << endl;

			// deltay
			unique_ptr<Matrix<T>> deltay;
			{
				deltay = *dy * mu_;
			}

			cout << *deltay << endl;

			// deltax
			unique_ptr<Matrix<T>> deltax;
			{
				deltax = *(*calculatePMatrixOf(sqrtw.get()) * sqrt(mu_)) * *dx;
			}

			cout << *deltax << endl;

			// deltas
			unique_ptr<Matrix<T>> deltas;
			{
				deltas = *(*calculatePMatrixOf(invsqrtw.get()) * sqrt(mu_)) * *ds;
			}

			cout << *deltas << endl;

			const auto step(this->alpha_);

			*x_ += *(*deltax * step);
			*s_ += *(*deltas * step);
			*y_ += *(*deltay * step);

			cout << *x_ << endl;
			cout << *s_ << endl;
			cout << *y_ << endl;

			LOG.info(format("Value is: {}", (*this->problem_->getObjectives()->transpose() * this->x_)->toString()));
			LOG.info(format("END ----------------------", iters));
			LOG.info("---------------");
		}

		LOG.info(format("Solved in {} iterations.", iters));

		cout << *(*this->problem_->getConstraints() * *this->x_) << endl;

		return SolutionStatus::OPTIMAL;
	}

	template <typename T>
		requires std::floating_point<T>
	SolutionStatus SOCPSolver<T>::optimizedSolver()
	{
		using namespace std;

		const MatrixFactory<T> matrixFactory(MatrixType::DENSE);

		size_t iters(0);

		while (checkIsTermination())
		{
			++iters;

			LOG.info(format("{}. iteration ----------------------", iters));

			this->mu_ = calculateMu();
			cout << this->mu_ << endl;

			const auto rb(calculateResidualB());
			const auto rc(calculateResidualC());

			cout << *rb << endl;
			cout << *rc << endl;

			const auto w(calculateW());
			const auto v(calculateV());

			cout << *w << endl;
			cout << *v << endl;

			const auto pv(calculatePv());

			cout << *pv << endl;

			const auto A_(calculateA_());
			const auto A_T(A_->transpose());

			cout << *A_ << endl;
			cout << *A_T << endl;

			const auto rc_(calculateResidualC_());

			cout << *rc_ << endl;

			// dy
			unique_ptr<Matrix<T>> dy;
			{
				// A_ * A_T * dy = -A_ * rc_ - rb - A_ * pv

				const auto lhs(matrixFactory.createMatrix());
				const auto rhs(matrixFactory.createMatrix());

				*lhs = *(*A_ * *A_T);
				*rhs = *(*(*(*(*A_ * -1) * *rc_) - *rb) - *(*A_ * *pv));

				dy = lhs->solve(rhs, DecompositionType::JACOBISVD);
			}

			cout << *dy << endl;

			// ds
			unique_ptr<Matrix<T>> ds;
			{
				// A_T * dy + ds = -rc_ => ds = -rc_ - A_T * dy
				ds = *(*rc_ * -1) - *(*A_T * *dy);
			}

			cout << *ds << endl;

			//dx
			unique_ptr<Matrix<T>> dx;
			{
				// dx + ds = pv => dx = pv - ds
				dx = *pv - *ds;
			}

			cout << *dx << endl;

			// deltay
			unique_ptr<Matrix<T>> deltay;
			{
				deltay = *dy * mu_;
			}

			cout << *deltay << endl;

			const auto sqrtw(calculatePowerOf(calculateW().get(), 0.5));

			// deltax
			unique_ptr<Matrix<T>> deltax;
			{
				deltax = *(*calculatePMatrixOf(sqrtw.get()) * sqrt(mu_)) * *dx;
			}

			cout << *deltax << endl;

			// deltas
			unique_ptr<Matrix<T>> deltas;
			{
				const auto invsqrtw(calculatePowerOf(sqrtw.get(), -1));

				deltas = *(*calculatePMatrixOf(invsqrtw.get()) * sqrt(mu_)) * *ds;
			}

			cout << *deltas << endl;

			const auto step(this->alpha_);

			*x_ += *(*deltax * step);
			*s_ += *(*deltas * step);
			*y_ += *(*deltay * step);

			cout << *x_ << endl;
			cout << *s_ << endl;
			cout << *y_ << endl;

			LOG.info(format("Value is: {}", (*this->problem_->getObjectives()->transpose() * this->x_)->toString()));
			LOG.info(format("END ----------------------", iters));
		}

		LOG.info(format("Solved in {} iterations.", iters));

		cout << *(*this->problem_->getConstraints() * *this->x_) << endl;

		return SolutionStatus::OPTIMAL;
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

		const auto n(this->problem_->getObjectives()->getRows());
		this->n_ = n;

		const auto nn(this->problem_->getConstraintsObjectives()->getRows());
		this->nn_ = nn;

		const MatrixFactory<T> matrixFactory(MatrixType::DENSE);

		x_ = matrixFactory.createMatrix(n, 1);
		y_ = matrixFactory.createMatrix(nn, 1);
		s_ = matrixFactory.createMatrix(n, 1);

		init_->initialize(x_.get(), y_.get(), s_.get());

		const auto ret(optimizedSolver());

		solution_ = make_shared<Solution<T>>(Solution<T>(x_, y_, s_));

		return ret;
	}

	template <typename T>
		requires std::floating_point<T>
	Solution<T> SOCPSolver<T>::getSolution()
	{
		return *solution_;
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
}
