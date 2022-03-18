#pragma once

#include <concepts>
#include <iostream>
#include <string>
#include <vector>

#include "Writer.hpp"
#include "Problem.hpp"
#include "Solution.hpp"

namespace OPLibrary
{
	template <typename T>
		requires std::floating_point<T>
	class CSVWriter final : public Writer<T>
	{
		const inline static std::string SEPARATOR = ",";
		const inline static std::string ENDLINE = "\n";

		template <typename E>
		void internalRowWrite(E elem);

		bool wasIterationWritten_;
		std::vector<std::string> headers_;
	public:
		explicit CSVWriter(std::ostream* output) : Writer<T>(output), wasIterationWritten_(false) {}

		void writeProblem(const Problem<T>* problem) override;
		void writeProblem(const std::shared_ptr<Problem<T>>& problem) override;
		void writeSolution(const Solution<T>* solution) override;
		void writeSolution(const std::shared_ptr<Solution<T>>& solution) override;
		void writeIteration(const size_t& iter, std::initializer_list<T> args) override;
		void setIterationHeaders(const std::vector<std::string>& headers) override;
	};

	template <typename T> requires std::floating_point<T>
	template <typename E>
	void CSVWriter<T>::internalRowWrite(E elem)
	{
		*this->output_ << std::format("{}" + SEPARATOR, elem);
	}

	template <typename T>
		requires std::floating_point<T>
	void CSVWriter<T>::writeProblem(const Problem<T>* problem)
	{
		const auto A(problem->getConstraints());
		const auto b(problem->getConstraintsObjectives()->getValues());
		const auto c(problem->getObjectives()->getValues());

		internalRowWrite("bT: ");
		std::for_each(b->begin(), b->end(), [this](auto n)
			{
				internalRowWrite(n);
			});
		*this->output_ << ENDLINE;

		internalRowWrite("cT: ");
		std::for_each(c->begin(), c->end(), [this](auto n)
			{
				internalRowWrite(n);
			});
		*this->output_ << ENDLINE;

		internalRowWrite("A: ");
		for (size_t i = 0; i < A->getRows(); ++i)
		{
			const auto currRow(A->getRow(i));
			std::for_each(currRow->begin(), currRow->end(), [this](auto n)
				{
					internalRowWrite(n);
				});
			*this->output_ << ENDLINE << SEPARATOR;
		}

		*this->output_ << std::endl;
	}

	template <typename T>
		requires std::floating_point<T>
	void CSVWriter<T>::writeProblem(const std::shared_ptr<Problem<T>>& problem)
	{
		return writeProblem(problem.get());
	}

	template <typename T>
		requires std::floating_point<T>
	void CSVWriter<T>::writeSolution(const Solution<T>* solution)
	{
		if (wasIterationWritten_) {
			*this->output_ << ENDLINE;
		}

		const auto optimal(solution->getOptimalValue());
		const auto x(solution->getPrimalSolution()->getValues());
		const auto y(solution->getDualSolutionY()->getValues());
		const auto s(solution->getDualSolutionS()->getValues());

		internalRowWrite("Optimal value: ");
		internalRowWrite(optimal);
		*this->output_ << ENDLINE;

		internalRowWrite("xT: ");
		std::for_each(x->begin(), x->end(), [this](auto n)
			{
				internalRowWrite(n);
			});
		*this->output_ << ENDLINE;

		internalRowWrite("yT: ");
		std::for_each(y->begin(), y->end(), [this](auto n)
			{
				internalRowWrite(n);
			});
		*this->output_ << ENDLINE;

		internalRowWrite("sT: ");
		std::for_each(s->begin(), s->end(), [this](auto n)
			{
				internalRowWrite(n);
			});
		*this->output_ << ENDLINE;

		*this->output_ << std::endl;
	}

	template <typename T>
		requires std::floating_point<T>
	void CSVWriter<T>::writeSolution(const std::shared_ptr<Solution<T>>& solution)
	{
		return writeSolution(solution.get());
	}

	template <typename T>
		requires std::floating_point<T>
	void CSVWriter<T>::writeIteration(const size_t& iter, std::initializer_list<T> args)
	{
		if (iter == 1)
		{
			std::for_each(headers_.begin(), headers_.end(), [this](auto n)
				{
					internalRowWrite(n);
				});
			*this->output_ << ENDLINE;

			wasIterationWritten_ = true;
		}

		internalRowWrite(iter);

		std::for_each(args.begin(), args.end(), [this](auto n)
			{
				internalRowWrite(n);
			});

		*this->output_ << std::endl;
	}

	template <typename T>
		requires std::floating_point<T>
	void CSVWriter<T>::setIterationHeaders(const std::vector<std::string>& headers)
	{
		this->headers_ = std::vector(headers.begin(), headers.end());
		this->headers_.insert(this->headers_.begin(), "Iteration");
	}
}
