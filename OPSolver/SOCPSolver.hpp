#pragma once

#include <format>
#include "Solver.hpp"
#include <numbers>
#include <vector>
#include <string>

#include "MatrixFactory.hpp"
#include "SOCPSolver.hpp"
#include "Solution.hpp"
#include "SolverException.hpp"

namespace OPLibrary
{
	/**
	 * \brief Solver for the SOCP algorithm.
	 */
	template <typename T>
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
			void initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s) override;
		};

		long double theta_;
		long double epsilon_;
		long double tau_;
		long double alpha_;
		long double mu_;
		long double beta_;

		std::unique_ptr<Initializator> init_;

		Matrix<T>* x_;
		Matrix<T>* y_;
		Matrix<T>* s_;

		bool checkIsTermination() const;

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
		static Matrix<T>* calculatePMatrixOf(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "P matrix can only be calculated for vectors.");

			const MatrixFactory<T> factory(MatrixType::DENSE);
			const auto n(vec->getRows());
			Matrix<T>* retMatrix(factory.createMatrix(n, n));

			// norm^2
			retMatrix->set(0, 0, std::pow(vec->norm(), 2));

			const auto x1(vec->get(0, 0));
			const unique_ptr<Matrix<T>> x2n(vec->block(1, 0, vec->getRows() - 1, 0));
			const unique_ptr<Matrix<T>> x2nT(x2n->transpose(false));

			// 2 * (x1 * x2:n)
			retMatrix->block(1, 0, retMatrix->getRows() - 1, 0, *(*x2n * 2) * x1);

			// 2 * (x1 * x2T:n)
			retMatrix->block(0, 1, 0, retMatrix->getCols() - 1, *(*x2nT * 2) * x1);

			// det(x) * En-1 + 2 * x2n * x2T:n
			// det(x) = x1^2 - norm(x2:n)^2
			const auto detx(std::pow(x1, 2) - std::pow(x2n->norm(), 2));

			const unique_ptr<Matrix<T>> E(factory.createMatrix(n - 1, n - 1));
			E->setValues(vector<T>(static_cast<size_t>(std::pow((n - 1), 2)), 0), n - 1, n - 1);
			E->setDiagonalValues(vector<T>(n - 1, 1));

			const unique_ptr<Matrix<T>> helperTriangularMatrix(factory.createMatrix(n - 1, n - 1));
			helperTriangularMatrix->setValues(vector<T>(std::pow(n - 1, 2), 0), n - 1, n - 1);

			for (size_t i = 0; i < (n - 1); ++i)
			{
				for (size_t j = 0; j < (n - 1); ++j)
				{
					// TODO: kerdezd meg
					//if (i == j) continue;

					helperTriangularMatrix->set(i, j, 2 * x2n->get(i, 0) * x2n->get(j, 0));
				}
			}

			retMatrix->block(1, 1, n - 1, n - 1, *helperTriangularMatrix + *(*E * detx));

			return retMatrix;
		}
		static Matrix<T>* calculateC1(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "C1 vector can only be calculated for vectors.");

			const auto n(vec->getRows());

			const MatrixFactory<T> factory(MatrixType::DENSE);
			const auto c1(factory.createMatrix(n, 1));
			c1->set(0, 0, 0.5);

			const unique_ptr<Matrix<T>> x2n(vec->block(1, 0, n - 1, 0));
			const auto normx2n(x2n->norm());

			const auto zeroVector__(vector<T>(n - 1, 0));
			if (x2n->getValues() == zeroVector__)
			{
				c1->block(1, 0, n - 1, 0, zeroVector__);
			} else
			{
				c1->block(1, 0, n - 1, 0, *(*x2n / normx2n) * 0.5);
			}

			return c1;
		}
		static Matrix<T>* calculateC2(Matrix<T>* vec)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "C2 vector can only be calculated for vectors.");

			const auto n(vec->getRows());

			const MatrixFactory<T> factory(MatrixType::DENSE);
			const auto c2(factory.createMatrix(n, 1));
			c2->set(0, 0, 0.5);

			const unique_ptr<Matrix<T>> x2n(vec->block(1, 0, n - 1, 0));
			const auto normx2n(x2n->norm());

			const auto zeroVector__(vector<T>(n - 1, 0));
			if (x2n->getValues() == zeroVector__)
			{
				c2->block(1, 0, n - 1, 0, zeroVector__);
			} else
			{
				c2->block(1, 0, n - 1, 0, *(*(*x2n * -1) / normx2n) * 0.5);
			}

			return c2;
		}
		static Matrix<T>* calculatePowerOf(Matrix<T>* vec, const double power)
		{
			using namespace std;

			assert(vec->getCols() == 1 && "Power of matrix can only be calculated for vectors.");

			const auto eigenmax(calculateEigenMaxOf(vec));
			const auto powerEigenmax(std::pow(eigenmax, power));
			const auto eigenmin(calculateEigenMinOf(vec));
			const auto powerEigenmin(std::pow(eigenmin, power));

			const unique_ptr<Matrix<T>> c1(calculateC1(vec));
			const unique_ptr<Matrix<T>> c2(calculateC2(vec));

			const auto ret(*(*c1 * powerEigenmax) + *(*c2 * powerEigenmin));

			return ret;
		}
		Matrix<T>* calculateW() const;

	public:
		SOCPSolver() : theta_(std::numbers::pi_v<long double> / 4), epsilon_(1.0e-6), tau_(1.0 / 2), alpha_(0.5),
			mu_(1), beta_(1.0 / 2), init_(new ClassicInitializator()), x_(nullptr), y_(nullptr), s_(nullptr)
		{
			using namespace std;

			this->INITIALIZABLE_ARGS = { "theta", "epsilon", "tau", "alpha", "mu" };
			this->INITIALIZATOR.insert(make_pair<string, long double*>("theta", &theta_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("epsilon", &epsilon_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("tau", &tau_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("alpha", &alpha_));
			this->INITIALIZATOR.insert(make_pair<string, long double*>("mu", &mu_));
		}

		~SOCPSolver() override
		{
			delete x_;
			delete y_;
			delete s_;
		}

		void setProblem(Problem<T>* problem) override;
		Problem<T>* getProblem() const override;

		SolutionStatus solve() override;

		Solution getSolution() override;
	};

	template <typename T>
	bool SOCPSolver<T>::checkIsTermination() const
	{
		// mu < epsilon; az x minden erteke pozitiv; s minden erteke pozitiv; kupfeltetel ellenorzes
		// tehat megkell nezni azt, hogy lambda_min(x) > 0
		// lambda_min(s) > 0
		// kupfeltetel -- most eltekintunk rola, ennel az atlagosabb implementacional biztos, hogy bent maradunk

		return (mu_ < epsilon_) && (calculateEigenMinOf(x_) > 0) && (calculateEigenMinOf(s_) > 0);
	}

	template <typename T>
	Matrix<T>* SOCPSolver<T>::calculateW() const
	{
		// w = P(x^1/2) * (P(x^1/2)*s)^(-1/2)
		// NT skalazasi pont
		// w vektor
		using namespace std;

		const unique_ptr<Matrix<T>> sqrtx(calculatePowerOf(this->x_, 0.5));
		const unique_ptr<Matrix<T>> pMatrixSqrtx(calculatePMatrixOf(sqrtx.get()));
		const unique_ptr<Matrix<T>> sMultipliedBypMatrix(*pMatrixSqrtx * *this->s_);
		const unique_ptr<Matrix<T>> sqrtSMultiplication(calculatePowerOf(calculatePowerOf(sMultipliedBypMatrix.get(), 0.5), -1));

		return *pMatrixSqrtx * *sqrtSMultiplication;
	}

	template <typename T>
	void SOCPSolver<T>::setProblem(Problem<T>* problem)
	{
		this->problem_ = problem;
	}

	template <typename T>
	Problem<T>* SOCPSolver<T>::getProblem() const
	{
		return this->problem_;
	}

	template <typename T>
	SolutionStatus SOCPSolver<T>::solve()
	{
		using namespace std;

		if (this->problem_ == nullptr)
		{
			throw SolverException("A problem was never set, try setting the problem first before trying to solve.");
		}

		if (this->problem_->getConstraints() == nullptr || this->problem_->getObjectives() == nullptr || this->problem_->getConstraintsObjectives() == nullptr)
		{
			throw SolverException(
				"Cannot solve optimization problem without correctly setting the constraints, objectives and constraint objectives.");
		}

		[this]
		{
			const size_t rows(this->problem_->getObjectives()->getRows());
			constexpr size_t cols(1);

			const MatrixFactory<T> matrixFactory(MatrixType::DENSE);
			x_ = matrixFactory.createMatrix(rows, cols);
			y_ = matrixFactory.createMatrix(rows, cols);
			s_ = matrixFactory.createMatrix(rows, cols);
		}();

		init_->initialize(x_, y_, s_);

		size_t iters__(0);

		while (checkIsTermination())
		{
			// kiszamolni a Pv-t
			// pv = v^-1 - v, ahol a v = P(w)^(-1/2) * x

			++iters__;
		}

		Logger::getInstance().info(format("Solved in {} iterations.", iters__));

		return SolutionStatus::NONOPTIMAL;
	}

	template <typename T>
	Solution SOCPSolver<T>::getSolution()
	{
		return {};
	}

	template <typename T>
	void SOCPSolver<T>::ClassicInitializator::initialize(Matrix<T>* x, Matrix<T>* y, Matrix<T>* s)
	{
		using namespace std;

		if (x == nullptr || y == nullptr || s == nullptr)
			throw SolverException(
				"Wrong arguments for initializer, matrices have to be allocated before initializing them.");

		const size_t rows(x->getRows());
		constexpr size_t cols(1);

		const vector<T> allZero(rows * cols, 0);

		x->setValues(allZero, rows, cols);
		x->set(0, 0, 1);
		y->setValues(allZero, rows, cols);
		s->setValues(allZero, rows, cols);
		s->set(0, 0, 1);
	}
}
