#pragma once

#include <format>
#include <fstream>
#include <string>

#include "Reader.hpp"
#include "ReaderException.hpp"
#include "Logger.h"

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
		explicit FileReader(std::ifstream* input) : rows_(0), cols_(0), input_(input) {}

		void readProblem(Problem<T>* problem) override;
		void readProblem(const std::shared_ptr<Problem<T>>& problem) override;
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
			catch (...)
			{
				Logger::getInstance().info(std::format("Failed to use {} as a number, skipping.\n", placeholder));
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
					catch (...)
					{
						Logger::getInstance().info(std::format("Failed to use {} as a number, skipping.\n", placeholder));
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
			cols_ = readParam();

			readInput(rows_, cols_, problem->getConstraints().get()); // A matrix
			readInput(rows_, 1, problem->getConstraintsObjectives().get()); // b vektor
			readInput(cols_, 1, problem->getObjectives().get()); // c vektor
		}
		else
		{
			throw ReaderException("Input file closed or non-existing.");
		}
	}

	template <typename T> requires std::floating_point<T>
	void FileReader<T>::readProblem(const std::shared_ptr<Problem<T>>& problem)
	{
		this->readProblem(problem.get());
	}
}
