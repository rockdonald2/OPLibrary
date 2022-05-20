#pragma once

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
	class FileReader : public Reader<T>
	{
	protected:
		size_t rows_;
		size_t cols_;

		std::ifstream* input_;

		[[nodiscard]] size_t readParam() const;
		void readInput(const size_t& maxRows, const size_t& maxCols, Matrix<T>* holder) const;

	public:
		explicit FileReader(std::ifstream* input) : rows_(0), cols_(0), input_(input) {}
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
				LOG.error("Failed to use " + placeholder + " as a number, skipping.\n");
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
						LOG.error("Failed to use " + placeholder + " as a number, skipping.\n");
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
}
