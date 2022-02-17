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

		void readParam(size_t& holder) const;

		void readInput(const size_t& maxRows, const size_t& maxCols, Matrix<T>* holder) const;

	public:
		explicit FileReader(std::ifstream* input) : rows_(0), cols_(0), input_(input) {}

		void readProblem(Problem<T>* problem) override;

		void readParams() override;
		[[nodiscard]] std::pair<size_t, size_t> getParams() const override;

		void readMatrix(Matrix<T>* matrix) override;
		void readVector(Matrix<T>* vector) override;
	};

	template <typename T>
		requires std::floating_point<T>
	void FileReader<T>::readParam(size_t& holder) const
	{
		using namespace std;
		string placeholder;

		do
		{
			*input_ >> placeholder;

			try
			{
				if (input_->is_open() && !input_->eof())
				{
					holder = stoll(placeholder);
					break;
				}
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

		for (auto i = 0; i < maxRows; ++i)
		{
			for (auto j = 0; j < maxCols; ++j)
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
			readParams();
			readMatrix(problem->getConstraints().get());
			readVector(problem->getConstraintsObjectives().get());
			readVector(problem->getObjectives().get());
		}
		else
		{
			throw ReaderException("Input file closed or non-existing.");
		}
	}

	template <typename T>
		requires std::floating_point<T>
	void FileReader<T>::readParams()
	{
		readParam(rows_);
		readParam(cols_);
	}

	template <typename T>
		requires std::floating_point<T>
	std::pair<size_t, size_t> FileReader<T>::getParams() const
	{
		return std::make_pair<size_t, size_t>(static_cast<size_t>(rows_), static_cast<size_t>(cols_));
	}

	template <typename T>
		requires std::floating_point<T>
	void FileReader<T>::readMatrix(Matrix<T>* matrix)
	{
		if (rows_ == 0 || cols_ == 0) readParams();

		matrix->setSize(rows_, cols_);
		readInput(rows_, cols_, matrix);
	}

	template <typename T>
		requires std::floating_point<T>
	void FileReader<T>::readVector(Matrix<T>* vector)
	{
		if (rows_ == 0 || cols_ == 0) readParams();

		vector->setSize(rows_, 1);
		readInput(rows_, 1, vector);
	}
}
