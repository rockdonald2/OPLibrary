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
	class FileReader final : public Reader<T>
	{
		size_t rows_;
		size_t cols_;

		std::ifstream* input_;

		void readParam(size_t& holder) const
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

		void readInput(const size_t& maxRows, const size_t& maxCols, Matrix<T>* holder) const
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

	public:
		explicit FileReader(std::ifstream* input) : rows_(0), cols_(0), input_(input) {}

		void readProblem(Matrix<T>* matrix, Matrix<T>* vector1, Matrix<T>* vector2) override
		{
			using namespace std;

			if (input_->is_open())
			{
				readParams();
				readMatrix(matrix);
				readVector(vector1);
				readVector(vector2);
			}
			else
			{
				throw ReaderException("Input file closed or non-existing.");
			}
		}

		void readParams() override
		{
			readParam(rows_);
			readParam(cols_);
		}
		[[nodiscard]] std::pair<size_t, size_t> getParams() const override
		{
			return std::make_pair<size_t, size_t>(static_cast<size_t>(rows_), static_cast<size_t>(cols_));
		}

		void readMatrix(Matrix<T>* matrix) override
		{
			if (rows_ == 0 || cols_ == 0) readParams();

			matrix->setSize(rows_, cols_);
			readInput(rows_, cols_, matrix);
		}
		void readVector(Matrix<T>* vector) override
		{
			if (rows_ == 0 || cols_ == 0) readParams();

			vector->setSize(rows_, 1);
			readInput(rows_, 1, vector);
		}
	};
}