// SPDX-License-Identifier: MIT
#include "ignis/io/Cli.hpp"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "ignis/core/Exceptions.hpp"
#include "ignis/core/Version.hpp"

namespace ignis {
namespace {

void printUsage(const std::string& program, const std::string& description,
                const std::vector<CliOption>& options, bool config_required) {
  std::ostringstream os;
  os << "Ignis " << kVersion << " -- " << program << "\n\n" << description << "\n\n"
     << "Usage: " << program << " [options]\n\n"
     << "Options:\n";
  auto line = [&](const std::string& flag, const std::string& help, const std::string& def) {
    os << "  " << std::left << std::setw(28) << flag << help;
    if (!def.empty()) os << " (default: " << def << ")";
    os << "\n";
  };
  line("-c, --config PATH",
       std::string("configuration file") + (config_required ? " [required]" : ""), "");
  line("-o, --output DIR", "output directory", "from the configuration");
  line("    --prefix NAME", "output file prefix", "from the configuration");
  line("    --set KEY=VALUE", "override a named parameter (repeatable)", "");
  for (const auto& o : options) {
    std::string f = "    --" + o.name;
    if (!o.argument.empty()) f += " " + o.argument;
    line(f, o.help, o.default_value);
  }
  line("    --quiet", "suppress the report on stdout", "");
  line("    --version", "print version information and exit", "");
  line("-h, --help", "print this message and exit", "");
  os << "\nExit status is 0 on success and non-zero on any failure.\n";
  std::cout << os.str();
}

}  // namespace

bool parseCommandLine(int argc, char** argv, const std::string& program,
                      const std::string& description, const std::vector<CliOption>& options,
                      CommandLine& out, bool config_required) {
  auto known = [&](const std::string& name) {
    for (const auto& o : options)
      if (o.name == name) return &o;
    return static_cast<const CliOption*>(nullptr);
  };

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&](const std::string& what) -> std::string {
      if (i + 1 >= argc) throw ConfigError("option " + what + " needs a value");
      return argv[++i];
    };
    if (a == "-h" || a == "--help") {
      printUsage(program, description, options, config_required);
      return false;
    }
    if (a == "--version") {
      std::cout << versionBanner() << "\n";
      return false;
    }
    if (a == "-c" || a == "--config") {
      out.config = next(a);
    } else if (a == "-o" || a == "--output") {
      out.output_dir = next(a);
    } else if (a == "--prefix") {
      out.prefix = next(a);
    } else if (a == "--quiet") {
      out.quiet = true;
    } else if (a == "--set") {
      const std::string kv = next(a);
      const auto eq = kv.find('=');
      if (eq == std::string::npos)
        throw ConfigError("--set expects KEY=VALUE, got '" + kv + "'");
      const std::string key = kv.substr(0, eq);
      char* end = nullptr;
      const double v = std::strtod(kv.c_str() + eq + 1, &end);
      if (end == kv.c_str() + eq + 1 || *end != '\0')
        throw ConfigError("--set " + key + ": '" + kv.substr(eq + 1) + "' is not a number");
      out.overrides.emplace_back(key, v);
    } else if (a.rfind("--", 0) == 0) {
      const std::string name = a.substr(2);
      const CliOption* opt = known(name);
      if (opt == nullptr) {
        std::ostringstream os;
        os << "unknown option '" << a << "'. Run " << program << " --help for the option list.";
        throw ConfigError(os.str());
      }
      out.flags[name] = opt->argument.empty() ? "1" : next(a);
    } else {
      throw ConfigError("unexpected argument '" + a + "'");
    }
  }
  if (config_required && out.config.empty())
    throw ConfigError("no configuration file given; use --config PATH (or --help)");
  return true;
}

std::string flagString(const CommandLine& cli, const std::string& name,
                       const std::string& fallback) {
  const auto it = cli.flags.find(name);
  return it == cli.flags.end() ? fallback : it->second;
}

int flagInt(const CommandLine& cli, const std::string& name, int fallback) {
  const auto it = cli.flags.find(name);
  if (it == cli.flags.end()) return fallback;
  return std::atoi(it->second.c_str());
}

double flagDouble(const CommandLine& cli, const std::string& name, double fallback) {
  const auto it = cli.flags.find(name);
  if (it == cli.flags.end()) return fallback;
  return std::atof(it->second.c_str());
}

bool flagPresent(const CommandLine& cli, const std::string& name) {
  return cli.flags.count(name) > 0;
}

}  // namespace ignis
