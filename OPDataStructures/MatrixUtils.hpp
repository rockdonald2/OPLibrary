#pragma once

#include <concepts>

#include "Matrix.hpp"

template <typename T>
	requires std::floating_point<T>
[[nodiscard]] bool containsNaN(const OPLibrary::Matrix<T>* matrix);




template <typename T>
	requires std::floating_point<T>
[[nodiscard]] bool containsNaN(const OPLibrary::Matrix<T>* matrix)
{
	const auto vals(matrix->getValues());
	return std::any_of(vals->begin(), vals->end(), [](const T& val)
		{
			return isnan(val);
		});
}
