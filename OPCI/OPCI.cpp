#include <format>

#include "ArgsParser.h"
#include "ArgumentException.hpp"
#include "ReaderBuilder.hpp"
#include "SolverFactory.hpp"
#include "Logger.hpp"
#include "MatrixFactory.hpp"
#include "WriterBuilder.hpp"

using TYPE = double;

int optimize(int argc, char* argv[])
{
	using namespace OPLibrary;
	using namespace std;

	auto hr(EXIT_SUCCESS);

	if (!ArgsParser::parseArguments(argc, argv)) return EXIT_FAILURE;

	try
	{
		ifstream inFile;
		inFile.open(ArgsParser::getStringArgument(ArgsParser::Args::INPUT_FILE));
		ofstream outFile;
		outFile.open(ArgsParser::getStringArgument(ArgsParser::Args::OUTPUT_FILE), ios::trunc);

		if (!inFile.is_open()) throw ArgumentException("Invalid input file.");

		const auto reader(ReaderBuilder<TYPE>().setType(ReaderType::FILE).setInput(&inFile).build());
		const auto writer(WriterBuilder<TYPE>().setType(WriterType::CSV).setOutput(&outFile).build());

		const MatrixFactory<TYPE> matrixFactory;

		auto matrix(matrixFactory.createMatrix());
		auto vector1(matrixFactory.createMatrix());
		auto vector2(matrixFactory.createMatrix());

		const auto problem(make_shared<Problem<TYPE>>(Problem(matrix, vector1, vector2)));

		reader->readProblem(problem);
		writer->writeProblem(problem);

		const auto solver(
			SolverFactory::createSolver<TYPE>(ArgsParser::getStringArgument(ArgsParser::Args::SOLVER_TYPE)));

		solver->setWriter(writer);
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

		writer->writeSolution(&solution);

		LOG.info("Optimization problem successfully resolved.");
	}
	catch (const ReaderException& e)
	{
		LOG.resetHandlers();
		LOG.error(e.what());
		hr = EXIT_FAILURE;
	}
	catch (const MatrixException& e)
	{
		LOG.resetHandlers();
		LOG.error(e.what());
		hr = EXIT_FAILURE;
	}
	catch (const SolverException& e)
	{
		LOG.resetHandlers();
		LOG.error(e.what());
		hr = EXIT_FAILURE;
	}
	catch (const ArgumentException& e)
	{
		LOG.resetHandlers();
		LOG.error(e.what());
		hr = EXIT_FAILURE;
	}
	catch (const WriterException& e)
	{
		LOG.resetHandlers();
		LOG.error(e.what());
		hr = EXIT_FAILURE;
	}
	catch (...)
	{
		LOG.resetHandlers();
		LOG.error("Unknown error.");
		hr = EXIT_FAILURE;
	}

	return hr;
}

int main(int argc, char* argv[])
{
	return optimize(argc, argv);
}
