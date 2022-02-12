#include <format>

#include "ArgsParser.h"
#include "ReaderBuilder.hpp"
#include "SolverFactory.hpp"
#include "Logger.h"
#include "MatrixFactory.hpp"

int main(int argc, char* argv[])
{
	using namespace OPLibrary;
	using namespace std;

	int hr(EXIT_SUCCESS);

	if (ArgsParser::parseArguments(argc, argv))
	{
		try
		{
			ifstream inFile;
			inFile.open(ArgsParser::getStringArgument(ArgsParser::Args::INPUT_FILE));

			const unique_ptr<Reader<long double>> reader(ReaderBuilder<long double>().setType(ReaderType::FILE).setInput(&inFile).build());

			const MatrixFactory<long double> matrixFactory(MatrixType::DENSE);

			const unique_ptr<Matrix<long double>> matrix(matrixFactory.createMatrix());
			const unique_ptr<Matrix<long double>> vector1(matrixFactory.createMatrix());
			const unique_ptr<Matrix<long double>> vector2(matrixFactory.createMatrix());

			const unique_ptr<Problem<long double>> problem(new Problem(matrix.get(), vector1.get(), vector2.get()));

			reader->readProblem(problem.get());

			inFile.close();

			const unique_ptr<Solver<long double>> solver(SolverFactory::createSolver<long double>(ArgsParser::getStringArgument(ArgsParser::Args::SOLVER_TYPE)));

			solver->setProblem(problem.get());

			for (const auto args = solver->getInitializableArgs(); const auto & arg : args)
			{
				const auto argE = ArgsParser::getArgByString(arg);

				try
				{
					const auto argVal = ArgsParser::getLongDoubleArgument(argE);
					solver->setInitializableArg(arg, argVal);
				}
				catch (...) {}
			}

			if (solver->solve() == SolutionStatus::NONOPTIMAL)
			{
				Logger::getInstance().error("There is no optimal solution for this optimization problem.");
				hr = EXIT_FAILURE;
			}
			else
			{
				Solution solution(solver->getSolution());
			}
		}
		catch (const ReaderException& e)
		{
			Logger::getInstance().error(e.what());
			hr = EXIT_FAILURE;
		}
		catch (const MatrixException& e)
		{
			Logger::getInstance().error(e.what());
			hr = EXIT_FAILURE;
		}
		catch (const SolverException& e)
		{
			Logger::getInstance().error(e.what());
			hr = EXIT_FAILURE;
		}
		catch (...)
		{
			Logger::getInstance().error("Unknown error.");
			hr = EXIT_FAILURE;
		}
	}

	return hr;
}
