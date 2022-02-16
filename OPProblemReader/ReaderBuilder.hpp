#pragma once

#include "Reader.hpp"
#include "FileReader.hpp"
#include "ReaderException.hpp"

namespace OPLibrary
{
	enum class ReaderType
	{
		FILE
	};

	/**
	 * \brief Builder for Readers.
	 */
	template <typename T>
	class ReaderBuilder
	{
		ReaderType type_;
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

		/**
		 * \brief Creates the Reader instance with set parameters.
		 * \return Reader instance
		 */
		[[nodiscard]] std::unique_ptr<Reader<T>> build() const
		{
			if (input_ == nullptr) throw ReaderException("Invalid input stream.");

			switch (type_)
			{
			case ReaderType::FILE: return std::move(std::make_unique<FileReader<T>>(static_cast<std::ifstream*>(input_)));
			}

			throw ReaderException("Unsupported reader type.");
		}
	};
}
