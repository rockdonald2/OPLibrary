#pragma once

#include <iostream>
#include <concepts>
#include <vector>
#include <unsupported/Eigen/Polynomials>

/**
 * \brief Computes the roots of a polynom.
 * \param coeff should respect the following rule: A, B, C in row, respectively for higher polynoms. E.g., Ax^2 + Bx + C = 0.
 * \return the roots
 */
template <typename T>
	requires std::floating_point<T>
[[nodiscard]] std::vector<T> solvePolynomial(const std::vector<T>& coeff);

template <typename T>
	requires std::floating_point<T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& rhs);


/* ========================================= */

template <typename T>
	requires std::floating_point<T>
[[nodiscard]] std::vector<T> solvePolynomial(const std::vector<T>& coeff)
{
	using namespace Eigen;
	using namespace std;

	auto tmp(coeff);
	std::reverse(tmp.begin(), tmp.end());

	Matrix<T, Dynamic, 1> eCoeff = Map<Matrix<T, Dynamic, 1>, Unaligned>(tmp.data(), tmp.size());

	PolynomialSolver<T, Dynamic> solver;
	solver.compute(eCoeff);

	const auto& r = solver.roots();

	vector<T> ret;
	for (size_t i = 0; i < static_cast<size_t>(r.size()); ++i)
	{
		const auto curr(r[i]);

		// csak a valos szamokat tartsuk meg
		if (curr.imag() != 0.0) continue;

		ret.push_back(curr.real());
	}

	return ret;
}

template <typename T>
	requires std::floating_point<T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& rhs)
{
	using namespace std;

	for_each(rhs.begin(), rhs.end(), [&out](T n)
		{
			out << n << " ";
		});

	out << endl;

	return out;
}