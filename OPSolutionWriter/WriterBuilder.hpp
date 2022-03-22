#pragma once

#include <concepts>
#include <memory>

#include "CSVWriter.hpp"
#include "Writer.hpp"
#include "WriterException.hpp"

namespace OPLibrary
{
	enum class WriterType
	{
		CSV,
		INVALID
	};

	template <typename T>
		requires std::floating_point<T>
	class WriterBuilder final
	{
		inline static const std::map<std::string, ProblemType> MAP_STR_TO_PROBLEM{ {"SOCP", ProblemType::SOCP} };

		WriterType type_;
		ProblemType problemType_;
		void* output_;

	public:
		WriterBuilder() = default;

		WriterBuilder& setType(const WriterType& type)
		{
			type_ = type;
			return *this;
		}

		WriterBuilder& setProblemType(const ProblemType& type)
		{
			problemType_ = type;
			return *this;
		}

		WriterBuilder& setProblemType(const std::string& type)
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

		WriterBuilder& setOutput(void* output)
		{
			output_ = output;
			return *this;
		}

		[[nodiscard]] std::shared_ptr<Writer<T>> build() const
		{
			if (output_ == nullptr) throw WriterException("Cannot build Writer with null output.");
			if (problemType_ == ProblemType::INVALID) throw WriterException("No problem type was set for writer.");

			switch (type_)
			{
			case WriterType::CSV: return std::move(std::make_shared<CSVWriter<T>>(static_cast<std::ostream*>(output_), problemType_));
			}

			throw WriterException("Unsupported writer type.");
		}
	};
}
