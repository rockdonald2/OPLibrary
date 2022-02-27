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
		const std::string SEPARATOR = ",";
		const std::vector<std::string> HEADERS = { "Iteration", "cTx", "bTy", "Duality Gap" };

		std::ofstream* output_;

		template <typename E>
		void internalRowWrite(E elem);
	public:
		explicit CSVWriter(std::ofstream* output) : output_(output) {}

		void writeProblem(const Problem<T>* problem) override;
		void writeProblem(const std::shared_ptr<Problem<T>>& problem) override;
		void writeSolution(const Solution<T>* solution) override;
		void writeSolution(const std::shared_ptr<Solution<T>>& solution) override;
		void writeIteration(const size_t& iter, std::initializer_list<T> args) override;
		[[nodiscard]] std::ostream* getOutput() const override;
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

		*output_ << "bT" + SEPARATOR;
		std::for_each(b->begin(), b->end(), [this](auto n)
			{
				internalRowWrite(n);
			});
		*output_ << "\n";

		*output_ << "cT" + SEPARATOR;
		std::for_each(c->begin(), c->end(), [this](auto n)
			{
				internalRowWrite(n);
			});
		*output_ << "\n";

		*output_ << "A" + SEPARATOR;
		for (size_t i = 0; i < A->getRows(); ++i)
		{
			const auto currRow(A->getRow(i));
			std::for_each(currRow->begin(), currRow->end(), [this](auto n)
				{
					internalRowWrite(n);
				});
			*output_ << SEPARATOR;
		}

		*output_ << std::endl;
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
		const auto x(solution->getX()->getValues());
		const auto y(solution->getY()->getValues());
		const auto s(solution->getS()->getValues());

		*output_ << "xT" + SEPARATOR;
		std::for_each(x->begin(), x->end(), [this](auto n)
			{
				internalRowWrite(n);
			});
		*output_ << "\n";

		*output_ << "yT" + SEPARATOR;
		std::for_each(y->begin(), y->end(), [this](auto n)
			{
				internalRowWrite(n);
			});
		*output_ << "\n";

		*output_ << "sT" + SEPARATOR;
		std::for_each(s->begin(), s->end(), [this](auto n)
			{
				internalRowWrite(n);
			});
		*output_ << "\n";

		*output_ << std::endl;
	}

	template <typename T>
		requires std::floating_point<T>
	void CSVWriter<T>::writeSolution(const std::shared_ptr<Solution<T>>& solution)
	{
		return writeSolution(solution.get());
	}

	template <typename T> requires std::floating_point<T>
	void CSVWriter<T>::writeIteration(const size_t& iter, std::initializer_list<T> args)
	{
		if (iter == 1)
		{
			std::for_each(HEADERS.begin(), HEADERS.end(), [this](auto n)
				{
					internalRowWrite(n);
				});
			*output_ << "\n";
		}

		*output_ << std::to_string(iter) + SEPARATOR;

		std::for_each(args.begin(), args.end(), [this](auto n)
			{
				internalRowWrite(n);
			});

		*output_ << std::endl;
	}

	template <typename T>
		requires std::floating_point<T>
	std::ostream* CSVWriter<T>::getOutput() const
	{
		return output_;
	}
}
