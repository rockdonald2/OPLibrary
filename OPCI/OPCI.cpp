#include <format>
#include <boost/range/combine.hpp>
#include <set>
#include <vector>
#include <thread>

#include "ArgsParser.h"
#include "ArgumentException.hpp"
#include "ReaderBuilder.hpp"
#include "SolverFactory.hpp"
#include "Logger.hpp"
#include "MatrixFactory.hpp"
#include "ProblemBuilder.hpp"
#include "WriterBuilder.hpp"
#include "ThreadUtils.hpp"

using TYPE = double;
ExecutorContainer<NamedJThread> gexecutors;

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

		const auto problem(ProblemBuilder<TYPE>()
			.setConstraints(matrixFactory.createMatrix())
			.setConstraintsObjectives(matrixFactory.createMatrix())
			.setObjectives(matrixFactory.createMatrix())
			.build());

		reader->readProblem(problem);
		writer->writeProblem(problem);

		const auto solver(SolverFactory::createSolver<TYPE>(*ArgsParser::getStringArgument(ArgsParser::Args::SOLVER_TYPE)));

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

		const auto currThread(gexecutors.getById(this_thread::get_id()));
		LOG.info(std::format("[{}] - Optimization problem successfully resolved. It is {}.",
			currThread ? currThread.value()->getName() : "Sequential (main) Optimizer", solution->getSolutionStatusString()));
	}
	catch (const ReaderException& e)
	{
		LOG.error(e.what());
		hr = EXIT_FAILURE;
	}
	catch (const MatrixException& e)
	{
		LOG.error(e.what());
		hr = EXIT_FAILURE;
	}
	catch (const SolverException& e)
	{
		LOG.error(e.what());
		hr = EXIT_FAILURE;
	}
	catch (const ArgumentException& e)
	{
		LOG.error(e.what());
		hr = EXIT_FAILURE;
	}
	catch (const WriterException& e)
	{
		LOG.error(e.what());
		hr = EXIT_FAILURE;
	}
	catch (const exception& e)
	{
		LOG.error(e.what());
		hr = EXIT_FAILURE;
	}
	catch (...)
	{
		LOG.error("Unknown exception, possibly from Eigen.");
		hr = EXIT_FAILURE;
	}

	return hr;
}

int optimize(int argc, char* argv[])
{
	using namespace OPLibrary;
	using namespace std;

	if (!ArgsParser::parseArguments(argc, argv)) return EXIT_FAILURE;

	const auto inputs(ArgsParser::getListArgument(ArgsParser::Args::INPUT_FILE));
	const auto outputs(ArgsParser::getListArgument(ArgsParser::Args::OUTPUT_FILE));
	const auto parallel(ArgsParser::getBooleanArgument(ArgsParser::Args::PARALLEL));

	{
		if (!inputs || !outputs)
		{
			LOG.blank("Inputs/Outputs unspecified.");
			ArgsParser::printHelp();
			return EXIT_FAILURE;
		}

		assert(inputs->size() == outputs->size() && "Individual outputs should be applied to inputs.");

		if (parallel) {
			const set setInputs(inputs->begin(), inputs->end());
			assert(setInputs.size() == inputs->size() && "Unique inputs are asserted.");

			const set setOutputs(inputs->begin(), inputs->end());
			assert(setOutputs.size() == outputs->size() && "Unique outputs are asserted.");
		}
	}

	const auto combined(boost::combine(inputs.value(), outputs.value()));

	if (parallel.has_value() && parallel.value() == true)
	{
		for (const auto& [input, output] : combined)
		{
			gexecutors.createExecutor(runOptimizer, input, output);
		}

		gexecutors.waitForExecutors();
		gexecutors.clearExecutors();
	}
	else
	{
		for (const auto& [input, output] : combined)
		{
			runOptimizer(input, output);
		}
	}

	return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
	return optimize(argc, argv);
}
