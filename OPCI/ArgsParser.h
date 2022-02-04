#pragma once

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace OPLibrary
{
	class ArgsParser final
	{
	public:
		/**
		 * \brief Representation of CLI arguments.
		 */
		enum class Args
		{
			CONFIG,
			INPUT_FILE,
			SOLVER_TYPE,
			THETA,
			EPSILON,
			ALPHA,
			MU,
			TAU
		};

	private:
		const static inline std::vector<Args> PARAMETERIZED_ARGS_VALS = { Args::CONFIG, Args::INPUT_FILE, Args::SOLVER_TYPE, Args::EPSILON, Args::THETA, Args::ALPHA, Args::TAU, Args::MU };
		static inline std::unordered_map<Args, std::vector<std::string>> options_ = []
		{
			using namespace std;
			unordered_map<Args, vector<string>> options;

			options.insert(make_pair(Args::CONFIG, vector<string>({ "-c", "--config" })));
			options.insert(make_pair(Args::INPUT_FILE, vector<string>({ "-f", "--file" })));
			options.insert(make_pair(Args::SOLVER_TYPE, vector<string>({ "-s", "--solver" })));
			options.insert(make_pair(Args::EPSILON, vector<string>({ "-e", "--epsilon" })));
			options.insert(make_pair(Args::THETA, vector<string>({ "-t", "--theta" })));
			options.insert(make_pair(Args::ALPHA, vector<string>({ "-a", "--alpha" })));
			options.insert(make_pair(Args::TAU, vector<string>({ "--tau" })));
			options.insert(make_pair(Args::MU, vector<string>({ "-m", "--mu" })));

			return options;
		}();
		static inline std::unordered_map<std::string, Args> STR_TO_ARG = []
		{
			using namespace std;
			unordered_map<string, Args> strToArg;

			strToArg.insert(make_pair("config", Args::CONFIG));
			strToArg.insert(make_pair("input_file", Args::INPUT_FILE));
			strToArg.insert(make_pair("solver_type", Args::SOLVER_TYPE));
			strToArg.insert(make_pair("epsilon", Args::EPSILON));
			strToArg.insert(make_pair("theta", Args::THETA));
			strToArg.insert(make_pair("alpha", Args::ALPHA));
			strToArg.insert(make_pair("mu", Args::MU));
			strToArg.insert(make_pair("tau", Args::TAU));

			return strToArg;
		}();

		static inline std::map<Args, std::string> args_;
		static inline std::vector<std::string> tokens_;

		/**
		 * \brief Internal usage only, private member. Returns if a certain option exists inside the already parsed tokens.
		 * \param option the option to search for
		 * \return boolean
		 */
		static bool doesOptionExist(const std::string&);
		/**
		 * \brief Internal usage only, private member. Returns the value of an option, if it exists, otherwise empty string.
		 * \param option the option to search for
		 * \return string
		 */
		static std::string getOptionVal(const std::string&);
		/**
		 * \brief Internal usage only, private member. Replaces the parsed tokens to other tokens from an input config file. Erases the existing tokens.
		 * \param path on which the config file can be found
		 */
		static void replaceTokensFromConfig(const std::string&);

	public:
		ArgsParser() = delete;

		/**
		 * \brief Prints help message.
		 */
		static void printHelp();

		/**
		 * \brief Parses input CLI arguments, builds up the internal token representation.
		 * \param argc number of arguments
		 * \param argv string array of arguments
		 * \return boolean value, true means the execution can continue, false means error happened or option specified which does not permit continuation of execution.
		 */
		static bool parseArguments(int, char**);

		/**
		 * \brief Returns enum representation of string CLI argument.
		 * \param argStr string representation of CLI argument
		 * \return Args value
		 */
		static Args getArgByString(const std::string&);

		/**
		 * \brief Returns value of CLI argument as string. Can throw.
		 * \param option to search for
		 * \return string
		 */
		static std::string getStringArgument(const Args&);
		/**
		 * \brief Returns value of CLI argument as integer. Can throw.
		 * \param option to search for
		 * \return integer
		 */
		static int getIntegerArgument(const Args&);
		/**
		 * \brief Returns value of CLI argument as double. Can throw.
		 * \param option to search for
		 * \return double
		 */
		static double getDoubleArgument(const Args&);
		/**
		 * \brief Returns value of CLI argument as long. Can throw.
		 * \param option to search for
		 * \return long
		 */
		static long getLongArgument(const Args&);
		/**
		 * \brief Returns value of CLI argument as long double. Can throw.
		 * \param option to search for
		 * \return long double
		 */
		static long double getLongDoubleArgument(const Args&);
		/**
		 * \brief Returns value of CLI argument as boolean. Can throw.
		 * \param option to search for
		 * \return long double
		 */
		static bool getBooleanArgument(const Args&);
	};
}
