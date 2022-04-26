#pragma once

#include "Reader.hpp"
#include "FileReader.hpp"
#include "ReaderException.hpp"
#include "SOCPFileReader.hpp"

namespace OPLibrary
{
	enum class ReaderType
	{
		FILE,
		INVALID
	};

	/**
	 * \brief Builder for Readers.
	 */
	template <typename T>
		requires std::floating_point<T>
	class ReaderBuilder final
	{
		inline static const std::map<std::string, ProblemType> MAP_STR_TO_PROBLEM{ {"SOCP", ProblemType::SOCP} };

		ReaderType type_;
		ProblemType problemType_;
		void* input_;

	public:
		/**
		 * \brief Creates an empty builder.
		 */
		ReaderBuilder() = default;

		/**
		 * \brief Sets type of builder.
		 * \param type of builder
		 * \return ReaderBuilder to chain
		 */
		ReaderBuilder& setType(const ReaderType& type)
		{
			type_ = type;
			return *this;
		}

		/**
		 * \brief Sets input for builder.
		 * \param input for builder
		 * \return ReaderBuilder to chain
		 */
		ReaderBuilder& setInput(void* input)
		{
			input_ = input;
			return *this;
		}

		ReaderBuilder& setProblemType(const ProblemType& type)
		{
			problemType_ = type;
			return *this;
		}

		ReaderBuilder& setProblemType(const std::string& type)
		{
			std::string temp(type);
			std::ranges::transform(temp, temp.begin(), [](const unsigned char c) { return std::toupper(c); });

			if (!MAP_STR_TO_PROBLEM.contains(temp))
				throw SolverException("Unsupported problem type.");

			switch (MAP_STR_TO_PROBLEM.at(temp))
			{
			case ProblemType::SOCP: return setProblemType(ProblemType::SOCP);
			}

			throw SolverException("Unsupported problem type.");
		}

		/**
		 * \brief Creates the Reader instance with set parameters.
		 * \return Reader instance
		 */
		[[nodiscard]] std::unique_ptr<Reader<T>> build() const
		{
			if (input_ == nullptr) throw ReaderException("Invalid input stream.");
			if (problemType_ == ProblemType::INVALID) throw ReaderException("No problem type was set for writer.");

			switch (type_)
			{
			case ReaderType::FILE:
			{
				switch (problemType_)
				{
				case ProblemType::SOCP: return std::move(std::make_unique<SOCPFileReader<T>>(static_cast<std::ifstream*>(input_)));
				}
			}
			}

			throw ReaderException("Unsupported reader type.");
		}
	};
}
