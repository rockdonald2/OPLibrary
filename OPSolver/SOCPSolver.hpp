#pragma once

#include <format>
#include <map>
#include "Solver.hpp"
#include <numbers>
#include <vector>
#include <string>

#include "Logger.h"
#include "MatrixFactory.hpp"
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

    public:
        SOCPSolver() : theta_(std::numbers::pi_v<long double> / 4), epsilon_(1.0e-6), tau_(1.0 / 2), alpha_(0.5),
                       mu_(1), beta_(1.0 / 2), init_(new ClassicInitializator()), x_(nullptr), y_(nullptr), s_(nullptr)
        {
            using namespace std;

            this->INITIALIZABLE_ARGS = {"theta", "epsilon", "tau", "alpha", "mu"};
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
            throw SolverException(
                "Cannot solve optimization problem without correctly setting the constraints, objectives and constraint objectives.");
        }

        [this]
        {
            const size_t rows(this->objectives_->getRows());
            constexpr size_t cols(1);

            const MatrixFactory<T> matrixFactory(MatrixType::DENSE);
            x_ = matrixFactory.createMatrix(rows, cols);
            y_ = matrixFactory.createMatrix(rows, cols);
            s_ = matrixFactory.createMatrix(rows, cols);
        }();

        init_->initialize(x_, y_, s_);

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
