#include <format>
#include <boost/range/combine.hpp>
#include <set>
#include <vector>

#include "ArgsParser.h"
#include "ArgumentException.hpp"
#include "ReaderBuilder.hpp"
#include "SolverFactory.hpp"
#include "Logger.hpp"
#include "MatrixFactory.hpp"
#include "WriterBuilder.hpp"

using TYPE = double;

int runOptimizer(const std::string& in, const std::string& out)
{
	using namespace std;
	using namespace OPLibrary;

	auto hr(EXIT_SUCCESS);

	try
	{
		ifstream inFile;
		inFile.open(in);
		ofstream outFile;
		outFile.open(out, ios::trunc);

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
			SolverFactory::createSolver<TYPE>(*ArgsParser::getStringArgument(ArgsParser::Args::SOLVER_TYPE)));

		solver->setWriter(writer);
		solver->setProblem(problem);

		for (const auto args(solver->getInitializableArgs()); const auto & arg : args)
		{
			if (const auto argE(ArgsParser::getArgByString(arg)); argE != ArgsParser::Args::NONEXISTING)
			{
				if (const auto argVal(ArgsParser::getDoubleArgument(argE)); argVal)
				{
					solver->setInitializableArg(arg, *argVal);
				}
			}
		}

		solver->solve();

		const auto solution(solver->getSolution());

		writer->writeSolution(solution);

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

int optimize(int argc, char* argv[])
{
	using namespace OPLibrary;

	auto hr(EXIT_SUCCESS);

	if (!ArgsParser::parseArguments(argc, argv)) return EXIT_FAILURE;

	const auto inputs(ArgsParser::getListArgument(ArgsParser::Args::INPUT_FILE));
	const auto outputs(ArgsParser::getListArgument(ArgsParser::Args::OUTPUT_FILE));

	// validation
	{
		assert(inputs->size() == outputs->size() && "Individual outputs should be applied to inputs.");

		const std::set setInputs(inputs->begin(), inputs->end());
		assert(setInputs.size() == inputs->size() && "Unique inputs are asserted.");

		const std::set setOutputs(inputs->begin(), inputs->end());
		assert(setOutputs.size() == outputs->size() && "Unique outputs are asserted.");
	}


	for (const auto& [input, output] : boost::combine(*inputs, *outputs))
	{
		runOptimizer(input, output);
	}

	return hr;
}

int main(int argc, char* argv[])
{
	return optimize(argc, argv);
}
