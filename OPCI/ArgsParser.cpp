#include <algorithm>
#include <format>
#include <fstream>
#include <sstream>
#include <cassert>

#include "ArgsParser.h"
#include "Logger.hpp"

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

		if (auto it = ranges::find(tokens_, option); it != tokens_.end())
		{
			while (++it != tokens_.end() && !it->starts_with("-") && !it->starts_with("/"))
			{
				ret.push_back(*it);
			}
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
		LOG.blank("Usage:");
		LOG.blank("\t -c | --config <config file path> \t -- specifies config file path, replaces any CLI arguments");
		LOG.blank("\t -f | --file <input file path> \t\t -- specifies input");
		LOG.blank("\t -o | --output <output file path> \t\t -- specifies output");
		LOG.blank("\t -s | --solver <solver type> \t\t -- specifies the solver type");
		LOG.blank("\t -e | --epsilon <value> \t\t -- specifies epsilon");
		LOG.blank("\t -t | --theta <value> \t\t\t -- specifies theta");
		LOG.blank("\t -a | --alpha <value> \t\t\t -- specifies alpha");
		LOG.blank("\t -m | --mu <value> \t\t\t -- specifies mu");
		LOG.blank("\t --tau <value> \t\t\t\t -- specifies tau");
		LOG.blank("\t -h | --help \t\t\t\t -- prints help message");
	}

	bool ArgsParser::parseArguments(int argc, char** argv)
	{
		using namespace std;

		for (auto i = 1; i < argc; ++i)
		{
			tokens_.emplace_back(argv[i]);
		}

		// help option
		if (doesOptionExist(vector<string>{"-h", "--help"}))
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

	std::shared_ptr<std::string> ArgsParser::getStringArgument(const Args& option)
	{
		if (args_.contains(option))
		{
			return std::make_shared<std::string>(*args_.at(option).begin());
		}

		return nullptr;
	}

	std::shared_ptr<int> ArgsParser::getIntegerArgument(const Args& option)
	{
		if (!args_.contains(option)) return nullptr;

		const auto val(args_[option]);
		const auto castedVal(std::stoi(*val.begin()));
		const auto ptrVal(new int);
		*ptrVal = castedVal;

		return std::shared_ptr<int>(ptrVal);
	}

	std::shared_ptr<double> ArgsParser::getDoubleArgument(const Args& option)
	{
		if (!args_.contains(option)) return nullptr;

		const auto val(args_[option]);
		const auto castedVal(std::stod(*val.begin()));
		const auto ptrVal(new double);
		*ptrVal = castedVal;

		return std::shared_ptr<double>(ptrVal);
	}

	std::shared_ptr<long> ArgsParser::getLongArgument(const Args& option)
	{
		if (!args_.contains(option)) return nullptr;

		const auto val(args_[option]);
		const auto castedVal(std::stol(*val.begin()));
		const auto ptrVal(new long);
		*ptrVal = castedVal;

		return std::shared_ptr<long>(ptrVal);
	}

	std::shared_ptr<long double> ArgsParser::getLongDoubleArgument(const Args& option)
	{
		if (!args_.contains(option)) return nullptr;

		const auto val(args_[option]);
		const auto castedVal(std::stold(*val.begin()));
		const auto ptrVal(new long double);
		*ptrVal = castedVal;

		return std::shared_ptr<long double>(ptrVal);
	}

	std::shared_ptr<bool> ArgsParser::getBooleanArgument(const Args& option)
	{
		if (!args_.contains(option)) return nullptr;

		auto val(*args_[option].begin());
		std::ranges::transform(val, val.begin(), [](const unsigned char c) { return std::tolower(c); });
		std::istringstream is(val);
		bool ret;
		is >> std::boolalpha >> ret;

		const auto ptrVal(new bool);
		*ptrVal = ret;

		return std::shared_ptr<bool>(ptrVal);
	}

	std::shared_ptr<std::vector<std::string>> ArgsParser::getListArgument(const Args& option)
	{
		using namespace std;

		if (!args_.contains(option)) return nullptr;

		const auto vals(args_[option]);
		auto heapVals = make_shared<vector<string>>();

		ranges::for_each(vals, [&heapVals](const auto& val)
			{
				heapVals->push_back(val);
			});

		return heapVals;
	}
}
