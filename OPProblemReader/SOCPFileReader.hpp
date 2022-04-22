#pragma once
#include <concepts>

#include "FileReader.hpp"

namespace OPLibrary
{
	template <typename T>
		requires std::floating_point<T>
	class SOCPFileReader final : public FileReader<T>
	{
	public:
		explicit SOCPFileReader(std::ifstream* input)
			: FileReader<T>(input) {}

		void readProblem(Problem<T>* problem) override;
		void readProblem(std::shared_ptr<Problem<T>>& problem) override;
	};

	template <typename T> requires std::floating_point<T>
	void SOCPFileReader<T>::readProblem(Problem<T>* problem)
	{
		using namespace std;

		if (problem->getObjectives() == nullptr || problem->getConstraints() == nullptr || problem->getConstraintsObjectives() == nullptr)
		{
			throw ReaderException("All problem matrices should be initialized before read.");
		}

		if (this->input_->is_open())
		{
			this->rows_ = this->readParam();

			const auto nof(this->readParam());

			auto* pr(dynamic_cast<SOCPProblem<T>*>(problem));

			vector<size_t> coneSizes;
			for (size_t i(0); i < nof; ++i)
			{
				coneSizes.push_back(this->readParam());
			}

			pr->setConeSizes(coneSizes);

			this->cols_ = std::reduce(std::execution::par, coneSizes.begin(), coneSizes.end(), static_cast<size_t>(0));

			this->readInput(this->rows_, this->cols_, problem->getConstraints().get()); // A matrix
			this->readInput(this->rows_, 1, problem->getConstraintsObjectives().get()); // b vector
			this->readInput(this->cols_, 1, problem->getObjectives().get()); // c vector
		}
		else
		{
			throw ReaderException("Input file closed or non-existing.");
		}
	}

	template <typename T> requires std::floating_point<T>
	void SOCPFileReader<T>::readProblem(std::shared_ptr<Problem<T>>& problem)
	{
		this->readProblem(problem.get());
	}
}
