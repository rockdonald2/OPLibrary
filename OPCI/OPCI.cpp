#include <format>

#include "ArgsParser.h"
#include "ArgumentException.hpp"
#include "ReaderBuilder.hpp"
#include "SolverFactory.hpp"
#include "Logger.h"
#include "MatrixFactory.hpp"
#include "WriterBuilder.hpp"

int main(int argc, char* argv[])
{
	using TYPE = long double;

	using namespace OPLibrary;
	using namespace std;

	auto hr(EXIT_SUCCESS);

	if (!ArgsParser::parseArguments(argc, argv)) return hr;

	try
	{
		ifstream inFile;
		inFile.open(ArgsParser::getStringArgument(ArgsParser::Args::INPUT_FILE));
		ofstream outFile;
		outFile.open(ArgsParser::getStringArgument(ArgsParser::Args::OUTPUT_FILE), ios::trunc);

		if (!inFile.is_open()) throw ArgumentException("Invalid input file.");

		const auto reader(ReaderBuilder<TYPE>().setType(ReaderType::FILE).setInput(&inFile).build());
		const auto writer(WriterBuilder<TYPE>().setType(WriterType::FILE).setOutput(&outFile).build());

		//LOG.setInfoHandler(writer->getWriter());
		//LOG.setErrorHandler(writer->getWriter());

		const MatrixFactory<TYPE> matrixFactory(MatrixType::DENSE);

		auto matrix(matrixFactory.createMatrix());
		auto vector1(matrixFactory.createMatrix());
		auto vector2(matrixFactory.createMatrix());

		const auto problem(make_shared<Problem<TYPE>>(Problem(matrix, vector1, vector2)));

		reader->readProblem(problem);
		writer->writeProblem(problem);

		inFile.close();

		const auto solver(
			SolverFactory::createSolver<TYPE>(ArgsParser::getStringArgument(ArgsParser::Args::SOLVER_TYPE)));

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
