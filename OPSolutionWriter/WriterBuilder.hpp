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
		CSV
	};

	template <typename T>
		requires std::floating_point<T>
	class WriterBuilder final
	{
		WriterType type_;
		void* output_;

	public:
		WriterBuilder() = default;

		WriterBuilder& setType(const WriterType& type)
		{
			type_ = type;
			return *this;
		}

		WriterBuilder& setOutput(void* output)
		{
			output_ = output;
			return *this;
		}

		[[nodiscard]] std::shared_ptr<Writer<T>> build() const
		{
			if (output_ == nullptr) throw WriterException("Cannot build Writer with null output.");

			switch (type_)
			{
			case WriterType::CSV: return std::move(std::make_shared<CSVWriter<T>>(static_cast<std::ostream*>(output_)));
			}

			throw WriterException("Unsupported writer type.");
		}
	};
}
