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
			OUTPUT_FILE,
			SOLVER_TYPE,
			THETA,
			EPSILON,
			ALPHA,
			MU,
			TAU,
			NONEXISTING
		};

	private:
		const static inline std::vector<Args> PARAMETERIZED_ARGS_VALS = { Args::CONFIG, Args::INPUT_FILE, Args::OUTPUT_FILE, Args::SOLVER_TYPE, Args::EPSILON, Args::THETA, Args::ALPHA, Args::TAU, Args::MU };
		static inline std::unordered_map<Args, std::vector<std::string>> options_ = []
		{
			using namespace std;
			unordered_map<Args, vector<string>> options;

			options.try_emplace(Args::CONFIG, vector<string>({ "-c", "--config" }));
			options.try_emplace(Args::INPUT_FILE, vector<string>({ "-f", "--file" }));
			options.try_emplace(Args::OUTPUT_FILE, vector<string>({ "-o", "--out" }));
			options.try_emplace(Args::SOLVER_TYPE, vector<string>({ "-s", "--solver" }));
			options.try_emplace(Args::EPSILON, vector<string>({ "-e", "--epsilon" }));
			options.try_emplace(Args::THETA, vector<string>({ "-t", "--theta" }));
			options.try_emplace(Args::ALPHA, vector<string>({ "-a", "--alpha" }));
			options.try_emplace(Args::TAU, vector<string>({ "--tau" }));
			options.try_emplace(Args::MU, vector<string>({ "-m", "--mu" }));


			return options;
		}();
		static inline std::unordered_map<std::string, Args> STR_TO_ARG = []
		{
			using namespace std;
			unordered_map<string, Args> strToArg;

			strToArg.try_emplace("config", Args::CONFIG);
			strToArg.try_emplace("input_file", Args::INPUT_FILE);
			strToArg.try_emplace("output_file", Args::INPUT_FILE);
			strToArg.try_emplace("solver_type", Args::SOLVER_TYPE);
			strToArg.try_emplace("epsilon", Args::EPSILON);
			strToArg.try_emplace("theta", Args::THETA);
			strToArg.try_emplace("alpha", Args::ALPHA);
			strToArg.try_emplace("mu", Args::MU);
			strToArg.try_emplace("tau", Args::TAU);

			return strToArg;
		}();

		static inline std::map<Args, std::vector<std::string>> args_;
		static inline std::vector<std::string> tokens_;

		/**
		 * \brief Internal usage only, private member. Returns if a certain option exists inside the already parsed tokens.
		 * \param option the option to search for
		 * \return boolean
		 */
		static bool doesOptionExist(const std::string&);
		/**
		 * \brief Internal usage only, private member. Returns if any given option exists inside the already parsed tokens.
		 * \return boolean
		 */
		static bool doesOptionExist(const std::vector<std::string>&);
		/**
		 * \brief Internal usage only, private member. Returns the value of an option, if it exists, otherwise empty string.
		 * \param option the option to search for
		 * \return string
		 */
		static std::vector<std::string> getOptionVal(const std::string&);
		/**
		 * \brief Internal usage only, private member. Returns the value of an option, which firstly exists in the given list, if it exists, otherwise empty string.
		 * \param option the option to search for
		 * \return string
		 */
		static std::vector<std::string> getOptionVal(const std::vector<std::string>&);
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
		static std::shared_ptr<std::string> getStringArgument(const Args&);
		/**
		 * \brief Returns value of CLI argument as integer. Can throw.
		 * \param option to search for
		 * \return integer
		 */
		static std::shared_ptr<int> getIntegerArgument(const Args&);
		/**
		 * \brief Returns value of CLI argument as double. Can throw.
		 * \param option to search for
		 * \return double
		 */
		static std::shared_ptr<double> getDoubleArgument(const Args&);
		/**
		 * \brief Returns value of CLI argument as long. Can throw.
		 * \param option to search for
		 * \return long
		 */
		static std::shared_ptr<long> getLongArgument(const Args&);
		/**
		 * \brief Returns value of CLI argument as long double. Can throw.
		 * \param option to search for
		 * \return long double
		 */
		static std::shared_ptr<long double> getLongDoubleArgument(const Args&);
		/**
		 * \brief Returns value of CLI argument as boolean. Can throw.
		 * \param option to search for
		 * \return long double
		 */
		static std::shared_ptr<bool> getBooleanArgument(const Args&);
		/**
		 * \brief Returns value of a CLI argument which contained multiple elements
		 * \return vector of strings
		 */
		static std::shared_ptr<std::vector<std::string>> getListArgument(const Args&);
	};
}
