#pragma once
#include <concepts>

#include "CSVWriter.hpp"
#include "SOCPProblem.hpp"
#include "SOCPSolution.hpp"

namespace OPLibrary
{
	template <typename T>
		requires std::floating_point<T>
	class SOCPCSVWriter final : public CSVWriter<T>
	{
	public:
		SOCPCSVWriter(std::ostream* output)
			: CSVWriter<T>(output) {}

	private:
		void writeProblem(const Problem<T>* problem) override;
		void writeProblem(const std::shared_ptr<Problem<T>>& problem) override;
		void writeSolution(const Solution<T>* solution) override;
		void writeSolution(const std::shared_ptr<Solution<T>>& solution) override;
		void writeIteration(const size_t& iter, std::initializer_list<T> args) override;
	};

	template <typename T> requires std::floating_point<T>
	void SOCPCSVWriter<T>::writeProblem(const Problem<T>* problem)
	{
		const auto A(problem->getConstraints());
		const auto b(problem->getConstraintsObjectives()->getValues());
		const auto c(problem->getObjectives()->getValues());

		const auto* pr(dynamic_cast<const SOCPProblem<T>*>(problem));

		this->internalRowWrite("noc: ");
		const auto& cs(pr->getConeSizes());
		std::for_each(cs.begin(), cs.end(), [this](auto n)
			{
				this->internalRowWrite(n);
			});

		*this->output_ << this->ENDLINE;

		this->internalRowWrite("bT: ");
		std::for_each(b->begin(), b->end(), [this](auto n)
			{
				this->internalRowWrite(n);
			});
		*this->output_ << this->ENDLINE;

		this->internalRowWrite("cT: ");
		std::for_each(c->begin(), c->end(), [this](auto n)
			{
				this->internalRowWrite(n);
			});
		*this->output_ << this->ENDLINE;

		this->internalRowWrite("A: ");
		for (size_t i = 0; i < A->getRows(); ++i)
		{
			const auto currRow(A->getRow(i));
			std::for_each(currRow->begin(), currRow->end(), [this](auto n)
				{
					this->internalRowWrite(n);
				});
			*this->output_ << this->ENDLINE << this->SEPARATOR;
		}

		*this->output_ << std::endl;
	}

	template <typename T> requires std::floating_point<T>
	void SOCPCSVWriter<T>::writeProblem(const std::shared_ptr<Problem<T>>& problem)
	{
		this->writeProblem(problem.get());
	}

	template <typename T> requires std::floating_point<T>
	void SOCPCSVWriter<T>::writeSolution(const Solution<T>* solution)
	{
		if (this->wasIterationWritten_) {
			*this->output_ << this->ENDLINE;
		}

		const auto feasibility(solution->getSolutionStatusString());
		const auto optimal(solution->getOptimalValue());
		const auto x(solution->getPrimalSolution()->getValues());
		const auto y(solution->getDualSolutionY()->getValues());
		const auto s(solution->getDualSolutionS()->getValues());

		const auto* sl(dynamic_cast<const SOCPSolution<T>*>(solution));

		this->internalRowWrite("noc: ");
		const auto& c(sl->getCones());
		std::for_each(c.begin(), c.end(), [this](auto n)
			{
				this->internalRowWrite(n);
			});

		*this->output_ << this->ENDLINE;

		this->internalRowWrite("Optimal value: ");
		this->internalRowWrite(optimal);
		*this->output_ << this->ENDLINE;

		this->internalRowWrite("Feasibility: ");
		this->internalRowWrite(feasibility);
		*this->output_ << this->ENDLINE;

		this->internalRowWrite("xT: ");
		std::for_each(x->begin(), x->end(), [this](auto n)
			{
				this->internalRowWrite(n);
			});
		*this->output_ << this->ENDLINE;

		this->internalRowWrite("yT: ");
		std::for_each(y->begin(), y->end(), [this](auto n)
			{
				this->internalRowWrite(n);
			});
		*this->output_ << this->ENDLINE;

		this->internalRowWrite("sT: ");
		std::for_each(s->begin(), s->end(), [this](auto n)
			{
				this->internalRowWrite(n);
			});
		*this->output_ << this->ENDLINE;

		*this->output_ << std::endl;
	}

	template <typename T> requires std::floating_point<T>
	void SOCPCSVWriter<T>::writeSolution(const std::shared_ptr<Solution<T>>& solution)
	{
		this->writeSolution(solution.get());
	}

	template <typename T> requires std::floating_point<T>
	void SOCPCSVWriter<T>::writeIteration(const size_t& iter, std::initializer_list<T> args)
	{
		if (iter == 1)
		{
			std::for_each(this->headers_.begin(), this->headers_.end(), [this](auto n)
				{
					this->internalRowWrite(n);
				});
			*this->output_ << this->ENDLINE;

			this->wasIterationWritten_ = true;
		}

		this->internalRowWrite(iter);

		std::for_each(args.begin(), args.end(), [this](auto n)
			{
				this->internalRowWrite(n);
			});

		*this->output_ << std::endl;
	}
}
