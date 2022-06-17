#pragma once

#include <concepts>

#include "Matrix.hpp"

template <typename T>
	requires std::floating_point<T>
[[nodiscard]] bool containsNaN(const OPLibrary::Matrix<T>* matrix);

template <typename T>
	requires std::floating_point<T>
[[nodiscard]] std::shared_ptr<OPLibrary::Matrix<T>> calculateCholesky(const OPLibrary::Matrix<T>* matrix);


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

template <typename T>
	requires std::floating_point<T>
[[nodiscard]] std::shared_ptr<OPLibrary::Matrix<T>> calculateCholesky(const OPLibrary::Matrix<T>* matrix)
{
	using namespace std;
	using namespace OPLibrary;

	auto vals(matrix->getValues());
	const auto rows(matrix->getRows());
	const auto cols(matrix->getCols());

	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigenMatrix = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(vals->data(), rows, cols);

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, 0> L(eigenMatrix.llt().matrixL());

	const MatrixFactory<T> factory;
	auto retMatrix(factory.createMatrix(L.rows(), L.cols()));

	auto tmpResultVec(make_shared<vector<T>>(L.size()));
	Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(tmpResultVec->data(), L.rows(), L.cols()) = L;

	retMatrix->setValues(*tmpResultVec, rows, cols);

	return retMatrix;
}
