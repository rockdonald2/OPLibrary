#pragma once

#include <concepts>
#include <iostream>

#include "Writer.hpp"
#include "Problem.hpp"
#include "Solution.hpp"

namespace OPLibrary
{
	template <typename T>
		requires std::floating_point<T>
	class FileWriter final : public Writer<T>
	{
		std::ofstream* output_;

	public:
		explicit FileWriter(std::ofstream* output) : output_(output) {}

		void writeProblem(const Problem<T>* problem) override;
		void writeProblem(const std::shared_ptr<Problem<T>>& problem) override;
		void writeSolution(const Solution<T>* solution) override;
		void writeSolution(const std::shared_ptr<Solution<T>>& solution) override;
		[[nodiscard]] std::ostream* getWriter() const override;
	};

	template <typename T>
		requires std::floating_point<T>
	void FileWriter<T>::writeProblem(const Problem<T>* problem)
	{
		*output_ << "|---------------------------Problem--------------------------|\n" << std::endl;
		*output_ << *problem;
		*output_ << "|------------------------------------------------------------|\n" << std::endl;
	}

	template <typename T>
		requires std::floating_point<T>
	void FileWriter<T>::writeProblem(const std::shared_ptr<Problem<T>>& problem)
	{
		return writeProblem(problem.get());
	}

	template <typename T>
		requires std::floating_point<T>
	void FileWriter<T>::writeSolution(const Solution<T>* solution)
	{
		*output_ << "|--------------------------Solution--------------------------|\n" << std::endl;
		*output_ << *solution;
		*output_ << "|------------------------------------------------------------|\n" << std::endl;
	}

	template <typename T>
		requires std::floating_point<T>
	void FileWriter<T>::writeSolution(const std::shared_ptr<Solution<T>>& solution)
	{
		return writeSolution(solution.get());
	}

	template <typename T>
		requires std::floating_point<T>
	std::ostream* FileWriter<T>::getWriter() const
	{
		return output_;
	}
}
