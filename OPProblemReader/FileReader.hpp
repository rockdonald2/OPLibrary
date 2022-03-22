#pragma once

#include <format>
#include <fstream>
#include <string>

#include "Reader.hpp"
#include "ReaderException.hpp"
#include "Logger.hpp"
#include "SOCPProblem.hpp"
#include "SOCPSolver.hpp"

namespace OPLibrary
{
	/**
	 * \brief Reader for file input.
	 */
	template <typename T>
		requires std::floating_point<T>
	class FileReader final : public Reader<T>
	{
		size_t rows_;
		size_t cols_;

		std::ifstream* input_;

		[[nodiscard]] size_t readParam() const;
		void readInput(const size_t& maxRows, const size_t& maxCols, Matrix<T>* holder) const;

	public:
		explicit FileReader(std::ifstream* input, const ProblemType& type) : Reader<T>(type), rows_(0), cols_(0), input_(input) {}

		void readProblem(Problem<T>* problem) override;
		void readProblem(std::shared_ptr<Problem<T>>& problem) override;
	};

	template <typename T>
		requires std::floating_point<T>
	size_t FileReader<T>::readParam() const
	{
		using namespace std;
		string placeholder;

		do
		{
			*input_ >> placeholder;

			try
			{
				if (input_->is_open() && !input_->eof()) return stoll(placeholder);

				throw ReaderException("Reached end of file or suddenly closed while reading in matrices.");
			}
			catch (const exception&)
			{
				LOG.error(std::format("Failed to use {} as a number, skipping.\n", placeholder));
			}
		} while (true);
	}

	template <typename T>
		requires std::floating_point<T>
	void FileReader<T>::readInput(const size_t& maxRows, const size_t& maxCols, Matrix<T>* holder) const
	{
		using namespace std;
		string placeholder;

		if (holder->getRows() != maxRows || holder->getCols() != maxCols) holder->setSize(maxRows, maxCols);

		for (auto i = 0ULL; i < maxRows; ++i)
		{
			for (auto j = 0ULL; j < maxCols; ++j)
			{
				if (input_->is_open() && !input_->eof())
				{
					*input_ >> placeholder;

					try
					{
						holder->set(i, j, stold(placeholder));
					}
					catch (const std::exception&)
					{
						LOG.error(std::format("Failed to use {} as a number, skipping.\n", placeholder));
						--j;
					}
				}
				else
				{
					throw ReaderException("Reached end of file or suddenly closed while reading in matrices.");
				}
			}
		}
	}

	template <typename T>
		requires std::floating_point<T>
	void FileReader<T>::readProblem(Problem<T>* problem)
	{
		using namespace std;

		if (problem->getObjectives() == nullptr || problem->getConstraints() == nullptr || problem->getConstraintsObjectives() == nullptr)
		{
			throw ReaderException("All problem matrices should be initialized before read.");
		}

		if (input_->is_open())
		{
			rows_ = readParam();

			// for extra SOCP related operations
			if (this->type_ == ProblemType::SOCP)
			{
				const auto nof(readParam());

				auto* pr(dynamic_cast<SOCPProblem<T>*>(problem));

				vector<size_t> coneSizes;
				for (size_t i(0); i < nof; ++i)
				{
					coneSizes.push_back(readParam());
				}

				pr->setConeSizes(coneSizes);

				cols_ = std::reduce(std::execution::par, coneSizes.begin(), coneSizes.end(), static_cast<size_t>(0));
			}
			else
			{
				cols_ = readParam();
			}

			readInput(rows_, cols_, problem->getConstraints().get()); // A matrix
			readInput(rows_, 1, problem->getConstraintsObjectives().get()); // b vector
			readInput(cols_, 1, problem->getObjectives().get()); // c vector
		}
		else
		{
			throw ReaderException("Input file closed or non-existing.");
		}
	}

	template <typename T>
		requires std::floating_point<T>
	void FileReader<T>::readProblem(std::shared_ptr<Problem<T>>& problem)
	{
		this->readProblem(problem.get());
	}
}
