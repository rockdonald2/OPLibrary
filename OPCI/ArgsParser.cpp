#include <algorithm>
#include <format>
#include <fstream>
#include <sstream>
#include <cassert>

#include "ArgsParser.h"
#include "Logger.hpp"
#include "VectorExtension.hpp"

namespace OPLibrary
{
	bool ArgsParser::doesOptionExist(const std::string& option)
	{
		return std::ranges::find(tokens_, option) != tokens_.end();
	}

	bool ArgsParser::doesOptionExist(const std::vector<std::string>& options)
	{
		return std::ranges::any_of(options.begin(), options.end(), [](const auto& option)
			{
				return doesOptionExist(option);
			});
	}

	std::vector<std::string> ArgsParser::getOptionVal(const std::string& option)
	{
		using namespace std;
		vector<string> ret;

		const auto bIt(ranges::find(tokens_, option));

		if (auto it = bIt; it != tokens_.end())
		{
			while (++it != tokens_.end() && !it->starts_with("-") && !it->starts_with("/"))
			{
				ret.push_back(*it);
			}

			const auto eIt = it;

			tokens_.erase(bIt, eIt);
		}

		return ret;
	}

	std::vector<std::string> ArgsParser::getOptionVal(const std::vector<std::string>& options)
	{
		for (const auto& option : options)
		{
			if (const auto & val(getOptionVal(option)); !val.empty()) return val;
		}

		return {};
	}

	void ArgsParser::replaceTokensFromConfig(const std::string& path)
	{
		using namespace std;
		ifstream file;
		file.open(path);

		tokens_.clear();

		if (file.is_open())
		{
			string placeholder;

			while (file.is_open())
			{
				file >> placeholder;
				if (file.eof()) break;
				tokens_.emplace_back(placeholder);
			}
		}
		else
		{
			LOG.error("Cannot find to specified config file, skipping.");
		}
	}

	void ArgsParser::printHelp()
	{
		LOG.blank("Command line interface for OPLibrary.");
		LOG.blank("Usage: OPCI.exe [OPTIONS] -f <input(s)> -o <output(s)>");
		LOG.blank("Options:");
		LOG.blank("\t -h | --help \t\t\t\t prints help message");
		LOG.blank("\t -c | --config <config file path> \t specifies config file path, replaces any CLI options");
		LOG.blank("\t -f | --file <input(s) file path> \t specifies input(s)");
		LOG.blank("\t -p | --parallel \t\t\t multiple file inputs should run parallel");
		LOG.blank("\t -o | --output <output(s) file path> \t specifies output(s)");
		LOG.blank("\t -s | --solver <solver type> \t\t specifies the solver type");
		LOG.blank("\t --init <initializator> \t\t specifies initializator for primal and dual problems");
		LOG.blank("\t --max \t\t\t\t\t specifies whether it is a maximization problem");
		LOG.blank("\t --min \t\t\t\t\t specifies whether it is a minimization problem");
		LOG.blank("\t --epsilon <value> \t\t\t specifies epsilon");
		LOG.blank("\t --rho <value> \t\t\t\t specifies rho");
		LOG.blank("\t --sigma <value> \t\t\t specifies sigma");
		LOG.blank("\t --mu <value> \t\t\t\t specifies mu");
		LOG.blank("Note:");
		LOG.blank("\t - same solver type, initializator and input parameters will be used for all input problems.");
		LOG.blank("\t - check docs for available initializators/solvers.");
	}

	bool ArgsParser::parseArguments(int argc, char** argv)
	{
		using namespace std;

		for (auto i = 1; i < argc; ++i)
		{
			// by convention options should start with / or -
			// we do not have any options with "/"
			if (argv[i][0] == '/') continue;

			tokens_.emplace_back(argv[i]);
		}

		// help option
		if (doesOptionExist(vector<string>{"-h", "--help"}) || tokens_.empty())
		{
			printHelp();
			return false;
		}

		// config option
		if (doesOptionExist(vector<string>{"-c", "--config"}))
		{
			if (const auto config = getOptionVal(vector<string>{"-c", "--config"}); config.empty())
			{
				LOG.error("Config argument was used, but no config path was specified, skipping.");
			}
			else
			{
				replaceTokensFromConfig(*config.begin());
			}
		}

		for (const auto& arg : PARAMETERIZED_ARGS_VALS)
		{
			for (const auto& option : options_[arg])
			{
				if (doesOptionExist(option))
				{
					if (const auto val(getOptionVal(option)); !val.empty())
					{
						args_.insert(make_pair(arg, val)); break;
					}

					LOG.error(format("Invalid value for argument {}.", option));
					return false;
				}
			}
		}

		for (const auto& arg : NONPARAMETERIZED_ARGS_VALS)
		{
			for (const auto& option : options_[arg])
			{
				if (doesOptionExist(option))
				{
					if (const auto val(getOptionVal(option)); val.empty())
					{
						args_.insert(make_pair(arg, vector<string>{"true"})); break;
					}

					LOG.error(format("Invalid value for argument {}, should be left empty.", option));
					return false;
				}
			}
		}

		if (args_.contains(Args::MAXIMIZE) && args_.contains(Args::MINIMIZE))
		{
			LOG.error("Invalid value set for objective direction, both max and min were specified.");
			return false;
		}

		if (!args_.contains(Args::SOLVER_TYPE))
		{
			LOG.error("No problem type was set.");
			return false;
		}

		if (!args_.contains(Args::INPUT_FILE) || !args_.contains(Args::OUTPUT_FILE))
		{
			LOG.error("No input or output was specified.");
			return false;
		}

		if (!tokens_.empty())
		{
			LOG.error(format("Unresolved parameters: {}.", toString(tokens_, " ")));
			return false;
		}

		if (args_.empty())
		{
			printHelp();
			return false;
		}

		return true;
	}

	ArgsParser::Args ArgsParser::getArgByString(const std::string& argStr)
	{
		if (STR_TO_ARG.contains(argStr)) return STR_TO_ARG.at(argStr);
		return Args::NONEXISTING;
	}

	std::optional<std::string> ArgsParser::getStringArgument(const Args& option)
	{
		if (args_.contains(option))
		{
			return { *args_.at(option).begin() };
		}

		return std::nullopt;
	}

	std::optional<int> ArgsParser::getIntegerArgument(const Args& option)
	{
		if (!args_.contains(option)) return std::nullopt;

		const auto val(args_[option]);
		const auto castedVal(std::stoi(*val.begin()));

		return { castedVal };
	}

	std::optional<double> ArgsParser::getDoubleArgument(const Args& option)
	{
		if (!args_.contains(option)) return std::nullopt;

		const auto val(args_[option]);
		const auto castedVal(std::stod(*val.begin()));

		return { castedVal };
	}

	std::optional<long> ArgsParser::getLongArgument(const Args& option)
	{
		if (!args_.contains(option)) return std::nullopt;

		const auto val(args_[option]);
		const auto castedVal(std::stol(*val.begin()));

		return { castedVal };
	}

	std::optional<long double> ArgsParser::getLongDoubleArgument(const Args& option)
	{
		if (!args_.contains(option)) return std::nullopt;

		const auto val(args_[option]);
		const auto castedVal(std::stold(*val.begin()));

		return { castedVal };
	}

	std::optional<bool> ArgsParser::getBooleanArgument(const Args& option)
	{
		if (!args_.contains(option)) return std::nullopt;

		auto val(*args_[option].begin());
		std::ranges::transform(val, val.begin(), [](const unsigned char c) { return std::tolower(c); });
		std::istringstream is(val);
		bool ret;
		is >> std::boolalpha >> ret;

		return { ret };
	}

	std::optional<std::vector<std::string>> ArgsParser::getListArgument(const Args& option)
	{
		using namespace std;

		if (!args_.contains(option)) return std::nullopt;
		if (args_.at(option).empty()) return std::nullopt;

		const auto vals(args_[option]);

		return { vals };
	}
}
