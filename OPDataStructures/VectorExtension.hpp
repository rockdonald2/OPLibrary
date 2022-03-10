#pragma once

#include <iostream>
#include <concepts>
#include <vector>
#include <unsupported/Eigen/Polynomials>

/**
 * \brief Computes the real roots of a polynom.
 * \param coeff should respect the following rule: A, B, C in row, respectively for higher polynoms. E.g., Ax^2 + Bx + C = 0.
 * \return real roots
 */
template <typename T>
	requires std::floating_point<T>
[[nodiscard]] std::vector<T> solvePolynomial(const std::vector<T>& coeff);

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& rhs);

template <typename T>
std::string toString(const std::vector<T>& rhs, const std::string& delimiter);




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

	const auto& r(solver.roots());

	vector<T> ret;
	for (size_t i = 0; i < static_cast<size_t>(r.size()); ++i)
	{
		const complex<T>* curr = &r[i];

		// csak a valos szamokat tartsuk meg
		if (curr->imag() != 0.0) continue;

		ret.push_back(curr->real());
	}

	return ret;
}

template <typename T>
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

template <typename T>
std::string toString(const std::vector<T>& rhs, const std::string& delimiter)
{
	using namespace std;

	if (rhs.empty()) return "";

	stringstream ss;
	auto first(true);
	for (const auto& elem : rhs)
	{
		if (!first) ss << delimiter;
		ss << elem;
		if (first) first = false;
	}

	return ss.str();
}

