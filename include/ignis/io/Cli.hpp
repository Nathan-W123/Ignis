// SPDX-License-Identifier: MIT
#pragma once
/// \file Cli.hpp
/// \brief Small argument parser shared by the Ignis applications.
///
/// Every application accepts the same core options so the tools feel alike:
///
///   -c, --config PATH      configuration file (required unless noted)
///   -o, --output DIR       output directory (overrides the config)
///       --prefix NAME      output file prefix (overrides the config)
///       --set KEY=VALUE    override any named parameter, repeatable
///       --quiet            suppress the human-readable report on stdout
///       --version          print version information and exit
///   -h, --help             print usage and exit
///
/// Unknown options are an error, not a warning.

#include <map>
#include <string>
#include <vector>

namespace ignis {

/// Parsed command line.
struct CommandLine {
  std::string config;
  std::string output_dir;
  std::string prefix;
  std::vector<std::pair<std::string, double>> overrides;
  bool quiet = false;
  /// Free-standing extra options, e.g. --threads 8 -> flags["threads"] = "8".
  std::map<std::string, std::string> flags;
};

/// Description of one application-specific option.
struct CliOption {
  std::string name;        ///< long name without the leading dashes
  std::string argument;    ///< placeholder, empty for a boolean flag
  std::string help;
  std::string default_value;
};

/// Parse `argc/argv`.  Returns false when --help or --version was handled, in
/// which case the caller should exit with status 0.
/// \throws ConfigError on a malformed command line.
bool parseCommandLine(int argc, char** argv, const std::string& program,
                      const std::string& description, const std::vector<CliOption>& options,
                      CommandLine& out, bool config_required = true);

/// Value of a flag as a string, or the fallback.
std::string flagString(const CommandLine& cli, const std::string& name,
                       const std::string& fallback);
/// Value of a flag as an integer.
int flagInt(const CommandLine& cli, const std::string& name, int fallback);
/// Value of a flag as a double.
double flagDouble(const CommandLine& cli, const std::string& name, double fallback);
/// True when the flag was present.
bool flagPresent(const CommandLine& cli, const std::string& name);

}  // namespace ignis
