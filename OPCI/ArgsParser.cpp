#include <algorithm>
#include <format>
#include <fstream>
#include <sstream>

#include "ArgsParser.h"
#include "Logger.h"

namespace OPLibrary
{
	bool ArgsParser::doesOptionExist(const std::string& option)
	{
		return std::ranges::find(tokens_, option) != tokens_.end();
	}

	std::string ArgsParser::getOptionVal(const std::string& option)
	{
		using namespace std;
		string ret;

		if (auto it = ranges::find(tokens_, option); it != tokens_.end() && ++it != tokens_.end())
		{
			ret = *it;
		}

		return ret;
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
			Logger::getInstance().error("Cannot find to specified config file, skipping.");
		}
	}

	void ArgsParser::printHelp()
	{
		Logger::getInstance().blank("Command line interface for OPLibrary.");
		Logger::getInstance().blank("Usage:");
		Logger::getInstance().blank("\t -c | --config <config file path> \t -- specifies config file path, replaces any CLI arguments");
		Logger::getInstance().blank("\t -f | --file <input file path> \t\t -- specifies input");
		Logger::getInstance().blank("\t -s | --solver <solver type> \t\t -- specifies the solver type");
		Logger::getInstance().blank("\t -e | --epsilon <value> \t\t -- specifies epsilon");
		Logger::getInstance().blank("\t -t | --theta <value> \t\t\t -- specifies theta");
		Logger::getInstance().blank("\t -a | --alpha <value> \t\t\t -- specifies alpha");
		Logger::getInstance().blank("\t -m | --mu <value> \t\t\t -- specifies mu");
		Logger::getInstance().blank("\t --tau <value> \t\t\t\t -- specifies tau");
		Logger::getInstance().blank("\t -h | --help \t\t\t\t -- prints help message");
	}

	bool ArgsParser::parseArguments(int argc, char** argv)
	{
		using namespace std;

		for (auto i = 1; i < argc; ++i)
		{
			tokens_.emplace_back(argv[i]);
		}

		if (doesOptionExist("-h") || doesOptionExist("--help"))
		{
			printHelp();
			return false;
		}

		if (doesOptionExist("-c"))
		{
			if (const auto config = getOptionVal("-c"); config.empty())
			{
				Logger::getInstance().error("Config argument was used, but no config path was specified, skipping.");
			}
			else
			{
				replaceTokensFromConfig(config);
			}
		}
		else if (doesOptionExist("--config"))
		{
			if (const auto config = getOptionVal("--config"); config.empty())
			{
				Logger::getInstance().error("Config argument was used, but no config path was specified, skipping.");
			}
			else
			{
				replaceTokensFromConfig(config);
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
						args_.insert(make_pair(arg, val));
						continue;
					}
					else
					{
						Logger::getInstance().error(format("Invalid value {} for argument {}.", val, option));
						return false;
					}
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
		return STR_TO_ARG[argStr];
	}

	std::string ArgsParser::getStringArgument(const Args& option)
	{
		return args_[option];
	}

	int ArgsParser::getIntegerArgument(const Args& option)
	{
		const auto val = args_[option];
		return std::stoi(val);
	}

	double ArgsParser::getDoubleArgument(const Args& option)
	{
		const auto val = args_[option];
		return std::stod(val);
	}

	long ArgsParser::getLongArgument(const Args& option)
	{
		const auto val = args_[option];
		return std::stol(val);
	}

	long double ArgsParser::getLongDoubleArgument(const Args& option)
	{
		const auto val = args_[option];
		return std::stold(val);
	}

	bool ArgsParser::getBooleanArgument(const Args& option)
	{
		auto val = args_[option];
		std::ranges::transform(val, val.begin(), [](const unsigned char c) { return std::tolower(c); });

		std::istringstream is(val);

		bool ret;
		is >> std::boolalpha >> ret;

		return ret;
	}
}