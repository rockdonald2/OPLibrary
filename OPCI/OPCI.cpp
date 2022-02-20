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

	if (!ArgsParser::parseArguments(argc, argv)) return hr;

	try
	{
		ifstream inFile;
		inFile.open(ArgsParser::getStringArgument(ArgsParser::Args::INPUT_FILE));

		const auto reader(ReaderBuilder<long double>().setType(ReaderType::FILE).setInput(&inFile).build());

		const MatrixFactory<long double> matrixFactory(MatrixType::DENSE);

		auto matrix(matrixFactory.createMatrix());
		auto vector1(matrixFactory.createMatrix());
		auto vector2(matrixFactory.createMatrix());

		const auto problem(make_shared<Problem<long double>>(Problem(matrix, vector1, vector2)));

		reader->readProblem(problem);

		cout << *problem->getConstraints() << endl;
		cout << *problem->getConstraintsObjectives() << endl;
		cout << *problem->getObjectives() << endl;

		inFile.close();

		const auto solver(
			SolverFactory::createSolver<long double>(ArgsParser::getStringArgument(ArgsParser::Args::SOLVER_TYPE)));

		solver->setProblem(problem);

		for (const auto args(solver->getInitializableArgs()); const auto & arg : args)
		{
			const auto argE(ArgsParser::getArgByString(arg));

			try
			{
				const auto argVal(ArgsParser::getLongDoubleArgument(argE));
				solver->setInitializableArg(arg, argVal);
			}
			catch (...) {}
		}

		solver->solve();

		const auto solution(solver->getSolution());

		cout << solution << endl;
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

	return hr;
}
