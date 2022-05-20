#pragma once

#include <concepts>
#include <iostream>
#include <string>
#include <vector>

#include "Writer.hpp"
#include "Problem.hpp"
#include "Solution.hpp"

namespace OPLibrary
{
	template <typename T>
		requires std::floating_point<T>
	class CSVWriter : public Writer<T>
	{
	protected:
		const inline static std::string SEPARATOR = ",";
		const inline static std::string ENDLINE = "\n";

		template <typename E>
		void internalRowWrite(E elem);

		bool wasIterationWritten_;
		std::vector<std::string> headers_;
	public:
		explicit CSVWriter(std::ostream* output) : Writer<T>(output), wasIterationWritten_(false) {}

		void setIterationHeaders(const std::vector<std::string>& headers) override;
	};

	template <typename T> requires std::floating_point<T>
	template <typename E>
	void CSVWriter<T>::internalRowWrite(E elem)
	{
		*this->output_ << elem << SEPARATOR;
	}

	template <typename T>
		requires std::floating_point<T>
	void CSVWriter<T>::setIterationHeaders(const std::vector<std::string>& headers)
	{
		this->headers_ = std::vector(headers.begin(), headers.end());
		this->headers_.insert(this->headers_.begin(), "Iteration");
	}
}
